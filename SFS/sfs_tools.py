#!/usr/bin/env python3

"""
FOREST Thomas (thomas.forest@college-de-france.fr)

Caution : At the moment for gzipped files only.

ARGS
--------

standalone usage : vcf_to_sfs.py VCF.gz nb_indiv

TODO
_____
Externalize sfs transforms in a function
Rectify SFS comp in parsed funct.

"""

import gzip
import sys
import matplotlib.pyplot as plt
import customgraphics

def sfs_from_vcf(n, vcf_file, folded = True, diploid = True, phased = False, verbose = False,
                 strip = False, count_ext = False):

    """
    Generates a Site Frequency Spectrum from a gzipped VCF file format.

    Parameters
    ----------
    n : int
        Nb of individuals in sample.
    vcf_file : str
        SNPs in VCF file format.

        Used to generate a Site Frequency Spectrum (SFS) from a VCF.

    Returns
    -------
    dict
        Site Frequency Spectrum (SFS)


    """
    
    if diploid and not folded:
        n *= 2
    # initiate SFS_values with a zeros dict
    # if strip:
    #     # "[1" removes the 0 bin
    #     # "n-1]" crop the last bin (n or n/2 for folded)
    #     SFS_dim = [1, n-1]
    # else:
    SFS_dim = [0, n+1]
    SFS_values = dict.fromkeys(range(SFS_dim[1]),0)
    count_pluriall = 0
    with gzip.open(vcf_file, "rb") as inputgz:
        line = inputgz.readline()
        genotypes = []
        print("Parsing VCF", vcf_file, "... Please wait...")
        while line:
            # decode gzipped binary lines
            line = line.decode('utf-8').strip()
            # every snp line, not comment or header
            if not line.startswith("##") and not line.startswith("#"):
                FIELDS = line.split("\t")
                # REF is col 4 of VCF
                REF = FIELDS[3].split(",")
                # ALT is col 5 of VCF
                ALT = FIELDS[4].split(",")            
                FORMAT = line.split("\t")[8:9]
                SAMPLES = line.split("\t")[9:]
                snp_genotypes = []
                allele_counts = {}
                allele_counts_list = []
                # SKIP the SNP if :
                # 1 : missing
                # 2 : deletion among REF
                # 3 : deletion among ALT
                if "./.:." in line \
                   or len(ALT[0]) > 1 \
                   or len(REF[0]) > 1:
                    line = inputgz.readline()
                    continue
                for sample in SAMPLES:
                    if not phased:
                        # for UNPHASED data
                        smpl_genotype = [int(a) for a in sample.split(':')[0].split('/') if a != '.']
                    else:
                        # for PHASED
                        smpl_genotype = [int(a) for a in sample.split(':')[0].split('|') if a != '.']
                    nb_alleles = set(smpl_genotype)
                    snp_genotypes += smpl_genotype
                # skip if all individuals have the same genotype
                # if len(set(snp_genotypes)) == 1:
                #     if folded or (folded == False and snp_genotypes.count(1) == 0) :
                #         line = inputgz.readline()
                #         continue         
                for k in set(snp_genotypes):
                    allele_counts[snp_genotypes.count(k)] = k
                    allele_counts_list.append(snp_genotypes.count(k))
                #print(allele_counts_list)
                if  len(set(snp_genotypes)) == 1 or allele_counts_list[0] == allele_counts_list[1]:
                    # If only heterozygous sites 0/1; skip the site (equivalent to n bin or n/2 bin for folded)
                    # skip if all individuals have the same genotype
                    line = inputgz.readline()
                    continue     
                if len(ALT) >= 2:
                    #pass
                    count_pluriall +=1
                    # TODO - work in progress
                    # for al in range(len(ALT)-1):
                    #     SFS_values[min(allele_counts_list)-1] += 1/len(ALT)
                    #     allele_counts_list.remove(min(allele_counts_list))
                else:
                    if folded:
                        SFS_values[min(allele_counts_list)-SFS_dim[0]] += 1
                    else :
                        # if unfolded, count the Ones (ALT allele)
                        #print(snp_genotypes, snp_genotypes.count(1))
                        SFS_values[snp_genotypes.count(1)-SFS_dim[0]] += 1
            # all the parsing is done, change line
            line = inputgz.readline()
            if verbose:
                print("SFS=", SFS_values)
        if strip:
            del SFS_values[0]
            del SFS_values[n]
        print("Pluriallelic sites =", count_pluriall)
    return SFS_values, count_pluriall


def sfs_from_parsed_vcf(n, vcf_dict, folded = True, diploid = True, phased = False, verbose = False):

    """
    Generates a Site Frequency Spectrum from a gzipped VCF file format.

    Parameters
    ----------
    n : int
        Nb of individuals in sample.
    vcf_file : str
        SNPs in VCF file format.

        Used to generate a Site Frequency Spectrum (SFS) from a VCF.

    Returns
    -------
    dict
        Site Frequency Spectrum (SFS)


    """
    
    if diploid and not folded:
        n *= 2
    # initiate SFS_values with a zeros dict
    SFS_values = dict.fromkeys(range(n),0)
    count_pluriall = 0

    for CHROM in vcf_dict:
        for SNP in vcf_dict[CHROM]:
            snp_genotypes = []
            allele_counts = {}
            allele_counts_list = []
            print(CHROM, SNP)
            for sample in vcf_dict[CHROM][SNP]["SAMPLES"]:
                if not phased:
                    # for UNPHASED data
                    smpl_genotype = [int(a) for a in sample.split(':')[0].split('/') if a != '.']
                else:
                    # for PHASED
                    smpl_genotype = [int(a) for a in sample.split(':')[0].split('|') if a != '.']
                nb_alleles = set(smpl_genotype)
                snp_genotypes += smpl_genotype
            # skip if all individuals have the same genotype
            if len(set(snp_genotypes)) == 1:
                continue
            for k in set(snp_genotypes):
                allele_counts[snp_genotypes.count(k)] = k
                allele_counts_list.append(snp_genotypes.count(k))
            SFS_values[min(allele_counts_list)-1] += 1
        # sum pluriall counts for this CHR to the rest
        count_pluriall += vcf_dict[CHROM]['NB_PLURIALL']
                    
    if verbose:
        print("SFS=", SFS_values)
    print("Pluriallelic sites =", count_pluriall)
    
    return SFS_values, count_pluriall


def barplot_sfs(sfs,  xlab, ylab="SNP Counts", folded=True, title = "Barplot", transformed = False, normalized = False):
    sfs_val = []
    n = len(sfs.values())
    sum_sites = sum(list(sfs.values()))
    for k, ksi in sfs.items():
        #ksi = list(sfs.values())[k-1]
        # k+1 because k starts from 0
        # if folded:
        #     # ?check if 2*n or not?
        #     sfs_val.append(ksi * k * (2*n - k))
        # else:
        #     if transformed:
        #         sfs_val.append(ksi * k)
        #     else:
        #         sfs_val.append(ksi)
        if transformed:
            ylab = r'$ \phi_i $'
            if folded:
                val = ((k*(2*n - k)) / (2*n))*(ksi)
            else:
                val = ksi * k
        else:
            val = ksi
        sfs_val.append(val)

        if not transformed and not normalized:
            ylab = r'$ \eta_i $'
            
    #terminal case, same for folded or unfolded
    if transformed:
        last_bin = list(sfs.values())[n-1] * n/2
    else:
        last_bin = list(sfs.values())[n-1]
    sfs_val[-1] = last_bin
    if normalized:
        #ylab = "Fraction of SNPs "
        ylab = r'$ \phi_i $'
        sum_val = sum(sfs_val)
        for k, sfs_bin in enumerate(sfs_val):
            sfs_val[k] = sfs_bin / sum_val
        
        #print(sum(sfs_val))
    #build the plot
    title = title+" (n="+str(len(sfs_val))+") [folded="+str(folded)+"]"+" [transformed="+str(transformed)+"]"
    print("SFS =", sfs)
    if folded:
        xlab = "Minor allele frequency"
    if transformed:
        print("Transformed SFS ( n =",len(sfs_val), ") :", sfs_val)
        #plt.axhline(y=1/n, color='r', linestyle='-')
    else:
        if normalized:
            # then plot a theoritical distribution as 1/i
            expected_y = [1/(2*x+1) for x in list(sfs.keys())]
            print(sum(expected_y))
            #plt.plot([x for x in list(sfs.keys())], expected_y, color='r', linestyle='-')
            #print(expected_y)
            
    customgraphics.barplot(x = [x for x in list(sfs.keys())], y= sfs_val, xlab = xlab, ylab = ylab, title = title)
    plt.show()

if __name__ == "__main__":
            
    if len(sys.argv) != 3:
        print("Need 2 args")
        exit(0)

    # PARAM : Nb of indiv
    n = int(sys.argv[2])
    sfs = sfs_from_vcf(n, sys.argv[1], folded = True, diploid = True, phased = False, strip = True)
    print(sfs)
