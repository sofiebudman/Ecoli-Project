#data from NCBI ECOLI S16 rRNA 
#https://www.ncbi.nlm.nih.gov/nuccore
library(ape)
#character vector of accession numbers for the project

strains <- c('IRQBAS-238', 'Aqeel10', 'Aqeel9', 'Aqeel9', 'Aqeel7', 'Aqeel6', 
             'Aqeel5', 'Aqeel4', 'Aqeel3', 'Aqeel2', 'Aqeel1', 'FaMsAh', 'KSA-89', 
             'KSA-15', 'RMA20', 'RMA19', 'RMA17', 'RMA15', 'RMA14','RMA12', 'RMA11', 
             'RMA18', 'RMA16', 'IRAQ-RAM2', 'HKA1', 'AAE', 'STF-8', 'VANM4', 'VANM2', 
             'SSH1', 'KSE34', 'IQ1','SAS2', 'SAS1', 'ATCC:25922', 'MBG-DUTH', 'BAGh-M1', 
             'JCM-16946')
#accession numbers: LC848137.1, LC844827.1, LC844826.1, LC844825.1,LC844824.1
IDs <- c('LC848137.1', 'LC844827.1', 'LC844826.1','LC844825.1', 'LC844824.1',
         'LC844823.1','LC844822.1','LC844821.1', 'LC844820.1', 'LC844819.1', 
         'LC844818.1', 'LC842012.1', 'LC834164.1', 'LC834134.1','LC817444.1', 
         'LC817443.1', 'LC817442.1','LC817441.1', 'LC817440.1', 'LC817439.1',
         'LC817438.1', 'LC817437.1', 'LC817436.1', 'PP809701.1','LC815916.1', 
         'LC764402.1', 'LC796842.1', 'LC777926.1', 'LC777924.1','LC754127.1', 
         'LC747145.1', 'LC738862.1', 'LC732201.1', 'LC732200.1', 'LC730904.1',
          'OU548744.1', 'LC712754.1', 'LC682250.1')
sequences <- read.GenBank(IDs,
                          seq.names = IDs,
                          species.names = TRUE,
                          as.character = TRUE)
write.dna(sequences, "sequences/DNA.fasta", format = "fasta")
write.dna(sequences, "sequences/DNA-simple.fasta", format = "fasta")

fasta_content <- readLines("sequences/DNA-simple.fasta")
writeLines(fasta_content, "sequences/DNA-simple.txt")
#Annotated with the strain name


#sequence analysis
install.packages("remotes")
remotes::install_github("GuangchuangYu/treeio")
install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("DECIPHER")
install.packages("viridis")
install.packages("adegenet")

library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)

seqs <- readDNAStringSet("sequences/DNA-simple.txt", format = "fasta")
seqs #view
seqs <- OrientNucleotides(seqs)
aligned <- AlignSeqs(seqs)
#aligns similarities and trims to same length

#brose sequences
BrowseSeqs(aligned, highlight = 0)

#save as new fasta
writeXStringSet(aligned, file = "sequences/E_Coli_Aligned.fasta")

dna <- read.alignment("sequences/E_Coli_Aligned.fasta", format = "fasta")


#distance matrix
D <- dist.alignment(dna, matrix = "similarity") #only works if all have the same length
D[is.na(D)] <- 0
#normalize d between 0 and 1
temp <- as.data.frame(as.matrix(D)) #data frame
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


#plot alignment with tree

  
 
  
