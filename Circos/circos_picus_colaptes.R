library(OmicCircos)
options(stringsAsFactors=FALSE)
library(pals)
library(dplyr)
library(RColorBrewer)

setwd("/home/thomasforest/Documents/Thèse/Picus/Circos/Colaptes_auratus/")

x <- read.table("map_coords_L10000_80i.mum", sep="\t", header=T, row.names = NULL, skip = 3)

names(x) <- c("X.S1.", "X.E1.", "X.S2.", "X.E2.", "X.LEN.1.", "X.LEN.2." , 
              "X...IDY.", "X.LEN.R.", "X.LEN.Q.", "X.TAGS.", "X")

# Count the number of each combination
combination_count <- table(x$X.TAGS., x$X)

# Print the combination count
print(combination_count)

df_combinations <- as.data.frame(combination_count, stringsAsFactors = FALSE)

names(df_combinations) <- c("X.TAGS.", "X", "Occurrences")

sorted_df <- df_combinations[order(-df_combinations$Occurrences), ]

# Filter the sorted_df based on occurrences greater than 500
filtered_df <- sorted_df[sorted_df$Occurrences > 500, ]

x.1 <- filter(x, X.TAGS. %in% filtered_df$X.TAGS.)

x.1 <- filter(x.1, X %in% filtered_df$X)


manual_caur_names <- c("CM037631.1", "CM037632.1", "CM037633.1", "CM037634.1", "CM037635.1",
                       "CM037636.1", "CM037637.1", "CM037638.1", "CM037639.1", "CM037640.1",
                       "CM037641.1", "CM037642.1", "CM037643.1", "CM037644.1", "CM037645.1",
                       "CM037646.1", "CM037647.1", "CM026886.1", "CM026887.1", "CM026888.1", 
                       "CM026889.1", "CM026890.1", "CM026891.1", "CM026892.1", "CM026893.1",
                       "CM026894.1", "CM026895.1", "CM026896.1", "CM026897.1"
)

x.1 <- filter(x, X.TAGS. %in% manual_caur_names)

# rearrange the input file to the following order
# target name, start, end, query name, start, end, identity, length target, length query
x.1 <- x.1[,c(10,1,2,11,3,4,7, 9,10)]

# subset scaffold names and start sites
x.2 <- x.1[,c(1,2,4,5)]

# rename scaffolds with prefix
#x.2[,3] <- paste("Pi_v", sapply(strsplit(x.2[,3], "_"), "[[", 3), sep="")

x.2 <- x.2[order(x.2[,1],decreasing=FALSE),]

# bin the segments of the circos plot
segf_lengths <- c()
x3 <- c()
for(a in 1:length(unique(x.2[,1]))) {
  a_rep <- x.2[x.2[,1] == unique(x.2[,1])[a], ]
  test <- unique(sort(a_rep[,2]))
  a_rep[,2] <- match(a_rep[,2], test)
  x3 <- rbind(x3, a_rep)
  
  segf_lengths <- rbind(segf_lengths, c(unique(x.2[,1])[a], length(test)))
}

x4 <- c()
for(a in 1:length(unique(x.2[,3]))) {
  a_rep <- x3[x3[,3] == unique(x.2[,3])[a], ]
  test <- unique(sort(a_rep[,4]))
  a_rep[,4] <- match(a_rep[,4], test)
  x4 <- rbind(x4, a_rep)
  segf_lengths <- rbind(segf_lengths, c(unique(x.2[,3])[a], length(test)))
}
segf_lengths <- data.frame(chrom=as.character(segf_lengths[,1]), bins=as.numeric(segf_lengths[,2]))

# create segments object
seg_f <- c()
for(a in 1:nrow(segf_lengths)) {
  a_name <- segf_lengths[a,1]
  a_lengths <- 0:(segf_lengths[a,2] - 1)
  a_lengths2 <- a_lengths + 1
  a_name <- rep(a_name, length(a_lengths))
  a_v <- rep("NA", length(a_lengths))
  a_note <- rep("NA", length(a_lengths))
  a_output <- cbind(a_name, a_lengths, a_lengths2, a_v, a_note)
  seg_f <- rbind(seg_f, a_output)
}
seg_f <- data.frame(seg.name=as.character(seg_f[,1]), seg.Start=as.numeric(seg_f[,2]), seg.End=as.numeric(seg_f[,3]), 
                    the.v=as.character(seg_f[,4]), Note=as.character(seg_f[,5]))

# make the links object in the right format
link_names <- paste("n", seq(from=1, to=nrow(x4), by=1), sep="")
link_v <- data.frame(seg1=as.character(x4[,1]), po1=as.numeric(x4[,2]), name1=as.character(link_names), 
                     seg2=as.character(x4[,3]), po2=as.numeric(x4[,4]), name2=as.character(link_names))
# subset for easier plotting
link_subset <- link_v[sample(1:nrow(link_v), floor(nrow(link_v) / 100)), ]


# make the angular database
seg_names <- sort(unique(segf_lengths[,1]))

# number of gallus kept
kept_gallus_chr <- length(unique(x4$X.TAGS.))

# number of picus kept
kept_picus_chr <- length(unique(x4$X))



# reorder chr for better plotting
seg_names <- c(manual_caur_names, seg_names[length(seg_names)], sort(seg_names[1:(kept_picus_chr-1)], decreasing = TRUE))

seg_colors <- c(rep("darkorange", kept_gallus_chr), rep("navyblue", length(seg_names)-kept_gallus_chr))

db <- segAnglePo(seg_f, seg=seg_names)



# colors for plotting
#cols <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(link_subset[,4])))
cols <- polychrome(length(unique(link_subset[,4])))
cols <- c(cols, "#ff3d00", "#1a237e", "#388e3c", "#ff6f00", "#99a9ff", "#943126", "#d4ac0d")
cols_rgb <- col2rgb(cols)
new_cols <- c()
for(a in 1:ncol(cols_rgb)) {
  a_col <- rgb(cols_rgb[1,a], cols_rgb[2,a], cols_rgb[3,a], max=255, alpha=0.95*255)
  new_cols <- c(new_cols, a_col)
}
cols <- new_cols

plot_colors <- rep("un", nrow(link_subset))
for(a in 1:length(unique(link_subset[,4]))) {
  a_rep <- sort(unique(link_subset[,4]))[a]
  plot_colors[link_subset[,4] == a_rep] <- cols[a]
}



par(mar=c(2,2,2,2)) 
plot(c(1,800),  c(1,800), type="n", axes=FALSE, 
     main="")
circos(R=360, cir=db, type="chr", col=seg_colors, print.chr.lab=TRUE, W=40, scale = F, cex=0.005)
circos(R=350, cir=db, W=40, mapping=link_subset, type="link", lwd=3, col=plot_colors)

