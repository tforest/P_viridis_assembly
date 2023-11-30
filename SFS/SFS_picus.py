import sfs_tools, customgraphics

# n=11 (Folded)
sfs_picus = { 1:2739965, 2:681333, 3:370473,4:243906,5:176483, 6:136197, 7:111297, 8:95027, 9: 85945, 10:85797, 11:55880}
sfs_tools.barplot_sfs(sfs_picus, folded = True, title = "Picus viridis", xlab = "Allele frequency", transformed = True, normalized = True)
sfs_tools.barplot_sfs(sfs_picus, folded = True, title = "Picus viridis", xlab = "Allele frequency", transformed = False, normalized = False)
