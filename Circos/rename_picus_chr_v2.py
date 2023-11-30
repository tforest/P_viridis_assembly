#!/usr/bin/env python3

import csv

with open('rename_tablev2.csv', mode='r') as infile:
    reader = csv.reader(infile)
    mydict = {rows[0]:rows[1] for rows in reader}
    
# Parse fasta and split into separated files with contig name as filename

splitting_char = "_"
keep_split_parts = -2
with open("/home/thomasforest/Documents/Thèse/Picus/FINAL_ASSEMBLY/Picus_viridis_before_renamed.fa") as fasta_input:
    for line in fasta_input:
        if line.startswith(">"):
            contig_name=line[1:].strip()
            print(contig_name)
            if contig_name in mydict.keys():
                contig_name = mydict[contig_name].strip()
                if not "unplaced" in contig_name:
                    contig_name = contig_name.split('_')[0]+"_c"+contig_name.split('_')[2]
                else:
                    contig_name = "unplaced_"+contig_name.split('_')[1]+"_c"+contig_name.split('_')[3]
            #contig_name=splitting_char.join(line.strip().split(splitting_char)[keep_split_parts:])
            fasta_out = open(contig_name+".fasta", 'w')
            fasta_out.write(">"+contig_name+"\n")
        else:
            fasta_out.write(line)
                
            
