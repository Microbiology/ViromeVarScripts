library(LDheatmap)
library(genetics)

TRinput <- read.delim("/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/ForVisualization/Staph-73209-DGR1-GenotypeTable.tsv", sep="\t", header=FALSE)
VRinput <- read.delim("/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/ForVisualization/Staph-73209-DGR2-GenotypeTable.tsv", sep="\t", header=FALSE)

makeplot <- function(x){
  TR <- x
  TR[TR=="G/N"] <- NA
  TR[TR=="T/N"] <- NA
  TR[TR=="C/N"] <- NA
  TR[TR=="A/N"] <- NA
  TR <- TR[,colSums(is.na(TR))<nrow(TR)]
  gensnp <- makeGenotypes(TR, convert=c(colnames(TR)))
  LDheatmap(gensnp, color = rev(grey.colors(20)), flip=TRUE, genetic.distances=c(1:length(gensnp)))
}

write.table(output$`R^2`, "~/Desktop/R2Table.txt")

pdf("/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/ForVisualization/Staph-DGR-LD.pdf", height=10, width=10)
  makeplot(TRinput)
  makeplot(VRinput)
dev.off()
