#install packages
BiocManager::install("ggtree")

install.packages("remotes")
remotes::install_github("GuangchuangYu/treeio")
library(ggtree)

dna <- read.alignment("multiple_sequence_alignment/ecoli_MSA.fasta", format = "fasta")

#create distance matrix
D <- dist.alignment(dna, matrix = "similarity") #only works if all have the same length
D[is.na(D)] <- 0 #normalize d between 0 and 1
temp <- as.data.frame(as.matrix(D))

#constructs heatmap of similarity
table.paint(temp, cleg = 0, clabel.row = .5, clabel.col = .5)+
  scale_color_viridis() #order shades + colors based on similarities

#create tre with NJ algorithm
tre <- nj(D)
class(tre)
tre <- ladderize(tre)

#base R plot
plot(tre, cex = 0.6)
title("E.Coli Similarity Tree")


#cluster dendogram
h_cluster <- hclust(D, method = "average", members = NULL)
plot(h_cluster, cex = 0.6, main = "E.Coli Similarity (Dendogram)" )


#plot with ggtree - fanned plot
ggtree(tre, yscale = "NA")+
  geom_tiplab(hjust = -0.3, size = 4, align = TRUE)+
  xlim(0,0.5)
ggsave("phylogenetic_trees/results/ecoli_fanned_plot.png")



#basic tree
ggtree(tre)+
  geom_tiplab(hjust =-0.3, size = 4, align = TRUE)+
  xlim(0,0.5)
ggsave("phylogenetic_trees/results/phylogenetic_tree_1.png")







