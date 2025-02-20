## E COLI Sequence Alignment

#data from NCBI ECOLI S16 rRNA
#https://www.ncbi.nlm.nih.gov/nuccore

#install and load required packages

install.packages("remotes")
remotes::install_github("GuangchuangYu/treeio")
install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("DECIPHER")
install.packages("viridis")
install.packages("adegenet")
install.packages("seqinr")

library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)
library(ape)
library(ggmsa)

#character vector of accession numbers for the project
strains <- c('IRQBAS-238', 'Aqeel10', 'Aqeel9', 'Aqeel9', 'Aqeel7', 'Aqeel6',
             'Aqeel5', 'Aqeel4', 'Aqeel3', 'STF-8', 'VANM4', 'VANM2',
             'SSH1', 'KSE34', 'IQ1','SAS2', 'SAS1', 'ATCC:25922', 'MBG-DUTH', 'BAGh-M1',
             'JCM-16946','VL8', 'DH7', 'H7', 'JCM-20932', 'JCM-20931', 'JCM-20377', 'JCM-20376',
             'JCM-20375', 'JCM-20351', 'JCM-20350', 'JCM-20349', 'IRQBAS127', 'M30', 'E65-Zambia2018',
             'E64-Zambia2018', 'E63-Zambia2018', 'E58-Zambia2018', 'E57-Zambia2018', 'E50-Zambia2018',
             'E48-Zambia2018')

#accession numbers: LC848137.1, LC844827.1, LC844826.1, LC844825.1,LC844824.1
IDs <- c('LC848137.1', 'LC844827.1', 'LC844826.1','LC844825.1', 'LC844824.1',
         'LC844823.1','LC844822.1','LC844821.1', 'LC844820.1', 'LC796842.1', 'LC777926.1', 'LC777924.1','LC754127.1',
         'LC747145.1', 'LC738862.1', 'LC732201.1', 'LC732200.1', 'LC730904.1',
          'OU548744.1', 'LC712754.1', 'LC682250.1', 'LC682613.1', 'LC666913.1',
          'LC666912.1', 'LC654900.1', 'LC654899.1', 'LC654898.1', 'LC654897.1',
         'LC654896.1', 'LC654892.1', 'LC654891.1', 'LC654890.1', 'LC648289.1',
         'LC649234.1', 'LC599972.1', 'LC599971.1', 'LC599970.1', 'LC599967.1',
         'LC599966.1', 'LC599964.1', 'LC599963.1')

#obtain sequences using accession numbers from genbank'

s <- load_seq_from_gen_bank(IDs,
                         seq.names = IDs,
                         species.names = TRUE,
                         as.character = TRUE)
s <- DNAStringSet(unlist(s))

write.dna(s, "multiple_sequence_alignment/ecoli_strain_seqs.fasta", format = "fasta", nbcol = -1)

fasta_content <- readLines("multiple_sequence_alignment/ecoli_strain_seqs.fasta")
writeLines(fasta_content, "multiple_sequence_alignment/ecoli_strain_seqs.txt")



seqs <- readDNAStringSet("multiple_sequence_alignment/ecoli_strain_seqs.fasta", format = "fasta") #ignore warning
seqs
seqs <- OrientNucleotides(seqs)
aligned <- AlignSeqs(seqs)

#sequence analysis
seqs <- readDNAStringSet("sequences/DNA-simple.txt", format = "fasta")
seqs
seqs <- OrientNucleotides(seqs)
aligned <- AlignSeqs(seqs)
#aligns similarities and trims to same length

#browse sequence alignment
BrowseSeqs(aligned, highlight = 0)

#save as new fasta file
writeXStringSet(aligned, file = "sequences/E_Coli_Aligned.fasta")

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



