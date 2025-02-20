dna <- read.alignment("sequences/E_Coli_Aligned.fasta", format = "fasta")


#create distance matrix
D <- dist.alignment(dna, matrix = "similarity") #only works if all have the same length
D[is.na(D)] <- 0
#normalize d between 0 and 1
temp <- as.data.frame(as.matrix(D))

table.paint(temp, cleg = 0, clabel.row = .5, clabel.col = .5)+
  scale_color_viridis()
#order shades + colors based on similarities

tre <- nj(D)
class(tre)
tre <- ladderize(tre)
#base R plot

plot(tre, cex = 0.6)
title("E.Coli Similarity based on the 16s rRNA")


#cluster dendogram
h_cluster <- hclust(D, method = "average", members = NULL)
plot(h_cluster, cex = 0.6, main = "E.Coli Similarity based on the 16s rRNA (Dendogram)" )


#plot with ggtree - fanned plot
ggtree(tre, yscale = "NA")+
  geom_tiplab(hjust = -0.3, size = 4, align = TRUE)+
  xlim(0,0.5)


#basic tree
ggtree(tre)+
  geom_tiplab(hjust =-0.3, size = 4, align = TRUE)+
  xlim(0,0.5)


#visualize alignment

alignment <- read.alignment("sequences/E_Coli_Aligned.fasta", format = "fasta")
#print
print(alignment)

#use ape
library(ape)
dna_alignment <- read.dna("sequences/E_Coli_Aligned.fasta", format = "fasta")

# Visualize alignment as an image
image(dna_alignment)


#use ggmsa
if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")
devtools::install_github("YuLab-SMU/ggmsa")

#create msa
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa")
browseVignettes("msa")
library(msa)

if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")
devtools::install_github("YuLab-SMU/ggmsa")


#put in msa format
mySequences <- readAAStringSet("sequences/DNA.fasta")
#view
mySequences

#align
myFirstAlignment <- msa(mySequences)
myFirstAlignment

myFirstAlignment <- as(myFirstAlignment, "DNAStringSet")
# Reassign the original names
names(myFirstAlignment) <- names(mySequences)

print(myFirstAlignment, show = "complete")

writeXStringSet(as(myFirstAlignment, "DNAStringSet"), "msa-alignment.fasta")



ggmsa("msa-alignment.fasta", 1080, 1340, color = "Clustal", font = "DroidSansMono", char_width = 0.5, seq_name = TRUE )



