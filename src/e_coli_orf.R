#data from https://archive.ics.uci.edu/dataset/120/e+coli+genes
#Data giving characteristics of each ORF (potential gene) in the
#E. coli genome. Sequence, homology (similarity to other genes) and
#structural information, and function (if known) are provided.


data <- read.table(("ecoli/ecoli.data"))
data <- as.data.frame(data)
colnames(data) <- c("Sequence ID","mcg", "gvh", "lip", "chg", "aac","alm1", "alm2", "class" )
#Sequence ID - accession number from the SWISS-PROT database
#mcg McGeoch's method for singal sequence recognition
#gvh vonHeijne's method for signal sequence recognition
#lip von Heijne's Signal Peptidase II consensus sequence score
#chg Presence of charge on N-terminus of predicted lipoproteins
#aac score of discriminant analysis of the amino acid content of outer membrane and periplasmic proteins
#alm1 score of the ALOM membrane spanning region prediction program
#alm2 score of ALOM program after excluding putative cleavable signal regions from the sequence
#class
  #cp  (cytoplasm)                                    143 (number 1)
  #im  (inner membrane without signal sequence)        77 (number 2)
  #pp  (perisplasm)                                    52 (number 3)
  #imU (inner membrane, uncleavable signal sequence)   35 (number 4)
  #om  (outer membrane)                                20 (number 5)
  #omL (outer membrane lipoprotein)                     5 (number 6)
  #imL (inner membrane lipoprotein)                     2 (number 7)
  #imS (inner membrane, cleavable signal sequence)      2 (number 8)

#create hisograms
classes <- data$class
classes <- as.data.frame(classes)
classes$classes <- as.factor(classes$classes)
levels(classes$classes) <- c("cytoplasm", "inner membrane without signal sequence",
                             "perisplasm", "inner membrane, uncleavable signal sequence",
                             "outer membrane", "outer membrane lipoprotein",
                             "inner membrane lipoprotein", "inner membrane, cleavable signal sequence")




#distribution
ggplot(classes, aes(classes), color = classes)+
  geom_bar(data = classes)+
  labs(title = "Distribution of subcellular localizations")+
  ylab ("occurence")+
  xlab("location")



#create cluster graphs for each clas
#
library(dplyr)
data_process <- data[,-c(1,9)]
data$class <- as.factor(data$class)
data$class <- as.numeric(data$class)
ecoli_matrix <- as.matrix(data_process)
library(ggplot2)
library(Rtsne)
set.seed(42)#for reproducible results
tsne <- Rtsne(ecoli_matrix, pca = FALSE, perplexity = 30, theta = 0.0)
plot(tsne$Y,col=data$class, asp=1)

tsne_for_plot <- data.frame(
  tSNE1 <-tsne$Y[,1],
  tSNE2 <- tsne$Y[,2],
  class <- factor(data$class)
)
levels(tsne_for_plot$class) <- c("cytoplasm", "inner membrane without signal sequence",
                                 "perisplasm", "inner membrane, uncleavable signal sequence",
                                 "outer membrane", "outer membrane lipoprotein",
                                 "inner membrane lipoprotein", "inner membrane, cleavable signal sequence")

##plot the t-sne distribution
ggplot(tsne_for_plot, aes(x = tSNE1, y = tSNE2, color = class)) +
  geom_point(size = 2) + # Adjust point size if needed
  theme_minimal() + # Use a clean theme
  labs(
    title = "t-SNE Plot for E. Coli potential genes",
    x = "tSNE Dimension 1",
    y = "tSNE Dimension 2",
    color = "Class" # Legend title
  )
