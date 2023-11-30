library(PopGenome)
library(vcfR)
library(ape)
library(pegas)
library(ggplot2)
library(svglite)

setwd("~/Documents/Thèse/Picus/Stats_picus")

vcf <- read.vcfR("./input/picus_filtered_10_380x_without_missing.vcf.gz")

x <- extract.gt(vcf, as.numeric = TRUE)
# Some stats about the number of sites in the VCF
genotyped_sites <- colSums(!is.na(x))
barplot(genotyped_sites)
hist(genotyped_sites)
# Hist is flat: same number of sites per sample, missing data were filtered out.
# Total Nb sites = 7226888

# Sequence length
L <- 1279164199
myDNA <- vcfR2DNAbin(vcf)
S <- length(seg.sites(myDNA))
s <- S/L
pi <- nuc.div(myDNA) * S / L

# Theta - Watterson
theta.s(myDNA)/L

library(pcadapt)
# to adapt text position in the PCA
library(ggrepel)

### PCA

out_bed_name = "picus_filtered"

vcf_filtered = "/home/thomasforest/Documents/Thèse/Picus/VCFs/picus_filtered_10_380x_without_missing.vcf"

# convert VCF to Plink2 BED file (not the same as UCSC BED)
#plink2_cmd <- paste("plink2 --vcf ", vcf_filtered, " --make-bed --allow-extra-chr --geno 0.1 --max-alleles 2 --snps-only --out ", out_bed_name, sep = "")

#system(command = plink2_cmd)

bed_picus = "picus_filtered.bed"
bim_picus = "picus_filtered.bim"
fam_picus = "picus_filtered.fam"

bed_file <- read.pcadapt(bed_picus, type = "bed")
x <- pcadapt(input = bed_file, ploidy = 2, K = 10)

# With names
poplist.names <- c(rep("Contemporary", 7),rep("Historical", 5))
pop.labels <- c("JF5180", "JF5191",  "JF5258",
                "JF5277",  "JF5325",
                "JF5345",  "JF5951",
                "MO19711090",      "MO19711091",
                "MO1991185",       "MO1991186",
                "MO1993190")
pop.labels <- c("ZO-2017-195", "ZO-2017-257",  "ZO-2017-316",
                "ZO-2017-334",  "ZO-2018-033",
                "ZO-2018-044",  "ZO-2020-125",
                "ZO-1971-1090",      "ZO-1971-1091",
                "ZO-1991-185",       "ZO-1991-186",
                "ZO-1993-190")

x_coord <- x$scores[,1]
y_coord <- x$scores[,2]

# Create a dataframe with the coordinates and labels
df_labels <- data.frame(x = x_coord, y = y_coord, labels = pop.labels, pop = poplist.names)
# Define the color scheme
#color_scheme <- c("Historical" = "green", "Contemporary" = "blue")

# Create the plots
ggplot(df_labels, aes(x, y)) +
  geom_point(aes(color = pop), size = 3) +
  scale_color_manual(values = c("Contemporary" = "#619CFF", "Historical" = "#F8766D")) +
  geom_text_repel(aes(label = labels), force = 0.5) +
  xlab("PC1 (11.7% Var)") +
  ylab("PC2 (11.1% Var)") +
  labs(color = "Population") 

# Explained Variance by components
##calculate EV
EV <- (x$singular.values^2)
EV2 <-EV*100
print(EV2)
scree_vals <- (x$singular.values**2/sum(x$singular.values**2))
plot(seq(1, length(scree_vals)), scree_vals, type = "b",
     lty = 4, pch = 7, xlab = "Number of Components (K)", 
     ylab = "Proportion of explained variance")
abline(v=2)

## Hierarchical clustering

library(SNPRelate)
# to run only ONCE
#snpgdsBED2GDS(bed.fn = bed_picus, bim.fn = bim_picus, fam.fn =  fam_picus, out.gdsfn = "out.gds")
gds.object <- snpgdsOpen("out.gds", allow.duplicate = TRUE)
diss.matrix <- snpgdsDiss(gdsobj = gds.object, num.thread=10, verbose=TRUE, autosome.only = FALSE)
diss.matrix$sample.id <- pop.labels
hc <- snpgdsHCluster(diss.matrix, need.mat=TRUE, hang=0.25)
clust_tree <- snpgdsCutTree(hc, label.Z = TRUE, label.H = TRUE)
#svg("dendrogram.svg", width = 591, height = 551)
snpgdsDrawTree(clust_tree, main = "Dendrogram based on dissimilarity",edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),
               y.label.kinship=T,leaflab="perpendicular",yaxis.kinship=F)
abline(h = 1, lty=2, col = "red")

### Autosomes
vcf.autosomes <- read.vcfR("~/Documents/Thèse/Picus/VCFs/Splitted_VCF/Autosomes.vcf")
myDNA_auto <- vcfR2DNAbin(vcf.autosomes)
# number of segregating sites
L_auto <- 1171329139
S_auto <- length(seg.sites(myDNA_auto))
# theta
theta.s(myDNA_auto)/L_auto
# pi
nuc.div(myDNA_auto)*S_auto/L_auto
### Z chromosome
vcf.Zchrom <- read.vcfR("~/Documents/Thèse/Picus/VCFs/Splitted_VCF/Zchrom.vcf")
myDNA_Z <- vcfR2DNAbin(vcf.Zchrom)
# number of segregating sites
L_Z <- 107835060
S_Z <- length(seg.sites(myDNA_Z))
# theta
theta.s(myDNA_Z)/L_Z
# pi 
nuc.div(myDNA_Z)*S_Z/L_Z


### HISTORICAL

calculateGeneticStats <- function(vcfFile, L) {
  # Read VCF file
  vcf <- read.vcfR(vcfFile)
  
  # Convert VCF to DNA binary format
  dna <- vcfR2DNAbin(vcf)
  
  # Calculate number of segregating sites (S)
  S <- length(seg.sites(dna))
  
  # Calculate normalized S
  S_normalized <- S / L
  
  # Calculate pi (π)
  pi <- nuc.div(dna)
  
  # Calculate theta (θ)
  theta <- theta.s(dna)
  
  # Return the results as a named list
  results <- list(S = S, S_normalized = S_normalized, pi = pi, theta = theta, vcf = vcf, dnaBin = dna)
  return(results)
}

# Total, using the custom function
picus_histo_total <- calculateGeneticStats("~/Documents/Thèse/Picus/VCFs/Splitted_VCF/Histo_Contemp/picus_histo.vcf", L=L)

# Autosomes
picus_histo_Auto <- calculateGeneticStats("~/Documents/Thèse/Picus/VCFs/Splitted_VCF/Histo_Contemp/Autosomes_histo.vcf", L=L_auto)

# Z 
picus_histo_Z <- calculateGeneticStats("~/Documents/Thèse/Picus/VCFs/Splitted_VCF/Histo_Contemp/Z_histo.vcf", L=L_Z)


# Contemp
picus_contemp_total <- calculateGeneticStats("~/Documents/Thèse/Picus/VCFs/Splitted_VCF/Histo_Contemp/picus_histo.vcf", L=L)
picus_contemp_Auto <- calculateGeneticStats("~/Documents/Thèse/Picus/VCFs/Splitted_VCF/Histo_Contemp/Autosomes_contemp.vcf", L=L_auto)
picus_contemp_Z <- calculateGeneticStats("~/Documents/Thèse/Picus/VCFs/Splitted_VCF/Histo_Contemp/Z_contemp.vcf", L=L_Z)


