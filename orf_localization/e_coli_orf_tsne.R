#data from https://archive.ics.uci.edu/dataset/120/e+coli+genes

#Data giving characteristics of each ORF (potential gene) in the
#E. coli genome. Sequence, homology (similarity to other genes) and
#structural information, and function (if known) are provided.

#This dataset contains various features related to protein localization in E. Coli

##load packages
library(ggplot2)
library(dplyr)
library(ggplot2)
library(Rtsne)
library(caret)


data <- read.table(("orf_localization/ecoli.data"))
head(data)
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
#levels(classes$classes) <- c("cytoplasm", "inner membrane without signal sequence",
#                            "perisplasm", "inner membrane, uncleavable signal sequence",
#                             "outer membrane", "outer membrane lipoprotein",
#                             "inner membrane lipoprotein", "inner membrane, cleavable signal sequence")




#distribution
ggplot(classes, aes(classes), color = classes)+
  geom_bar(data = classes)+
  labs(title = "Distribution of subcellular localizations")+
  ylab ("occurence")+
  xlab("location")

ggsave("orf_localization/results/distribution_bar_chart.png")



#create cluster graphs for each clas

data_process <- data[,-c(1,9)]
data$class <- as.factor(data$class)
data$class <- as.numeric(data$class)
ecoli_matrix <- as.matrix(data_process)

set.seed(42)#for reproducible results
tsne <- Rtsne(ecoli_matrix, pca = FALSE, perplexity = 30, theta = 0.0)
plot(tsne$Y,col=data$class, asp=1) #simpel graph

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
  labs(
    title = "t-SNE Plot for Localization of E. Coli ORFs",
    subtitle = "Dimensionality reduction visualization",
    x = "tSNE Dimension 1",
    y = "tSNE Dimension 2",
    color = "Class" # Legend title
  )+
  theme(
    plot.title = element_text(color = "black", face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(size = 12, hjust = 0.5, face = "italic"), # Add a subtle subtitle
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold", hjust = 0.5), # Improve legend title visibility
    legend.box = "vertical",
    legend.text = element_text(size = 5), # Improve readability of legend items
    axis.title = element_text(size = 12), # Increase axis title size for readability
    axis.text = element_text(size = 10) # Increase tick label size


  )+
  guides(
    color = guide_legend(title.position = "top",nrow = 3, byrow = TRUE) # Arrange legend items in a single row
  )

ggsave("orf_localization/results/orf_tsne.png")





