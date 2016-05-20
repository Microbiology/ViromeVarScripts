# GenomicEvolutionOfEachGene.R
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

# Load in the libraries
library(ggplot2)
library(pgirmess)
library(RColorBrewer)

# Due to bash processing, the first row needs to be removed out of the
# data frame.
HpvInputPnps <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/EvoByGeneByTaxa/PnPsResults/HpvOverallContigPressure.tsv", sep="\t", header=TRUE)
HPVAnnotations <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/SnpHotspotEvolution/AnnotationFiles/HPVblastnResultsAnnoated.tsv", sep="\t", header=FALSE)
# Be sure to format some of the gene names for consistency
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"L1_protein","L1")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"L2_protein","L2")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"E2_protein","E2")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Replication_protein_E1","E1")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Regulatory_protein","E1")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Major_capsid_protein_L1","L1")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Early_protein_E2","E2")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"E6_protein","E6")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Capsid_protein_L1","L1")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Transforming_protein_E6","E6")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Transforming_protein","E6")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Major_capsid_protein_L1","L1")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Late_protein_L2","L2")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Major_capsid_protein_L1","L1")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Putative_E2_product","E2")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Major_capsid_protein","L1")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Protein_E7","E7")

HPVOccurList <- as.data.frame(table(HPVAnnotations$V4))
HPVOccurList <- HPVOccurList[c(which(HPVOccurList$Freq > 1)),]
HPVAnnotations <- HPVAnnotations[c(which(HPVAnnotations$V4 %in% HPVOccurList$Var1)),]


StaphInputPnps <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/EvoByGeneByTaxa/PnPsResults/StaphPhageOverallContigPressure.tsv", sep="\t", header=TRUE)
StaphPhageAnnotations <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/SnpHotspotEvolution/AnnotationFiles/StaphPhageblastnResultsAnnoated.tsv", sep="\t", header=FALSE)
# Remove all of the rows with value ORF
StaphOccurList <- as.data.frame(table(StaphPhageAnnotations$V4))
StaphOccurList <- StaphOccurList[c(which(StaphOccurList$Freq > 1)),]
StaphPhageAnnotations <- StaphPhageAnnotations[c(which(StaphPhageAnnotations$V4 %in% StaphOccurList$Var1)),]


PropInputPnps <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/EvoByGeneByTaxa/PnPsResults/PropPhageOverallContigPressure.tsv", sep="\t", header=TRUE)
PropPhageAnnotations <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/SnpHotspotEvolution/AnnotationFiles/PropPhageblastnResultsAnnoated.tsv", sep="\t", header=FALSE)


PlotGenePressure <- function (InputPnps, Annotations) {
  InputPnps$TotalID <- paste(InputPnps$ContigID, InputPnps$GeneID, sep="_.")
  
  MergeDf <- merge(InputPnps, Annotations, by.x="TotalID", by.y="V1")
  
  colourCount = length(unique(MergeDf$V4))
  getPalette = colorRampPalette(brewer.pal(8, "Set2"))
  
  
  ggplot(MergeDf, 
         aes(x=reorder(V4, PnPs, FUN=median), 
             y=log2(PnPs), 
             fill=V4)) + 
    theme_classic() + 
    geom_boxplot(notch=FALSE, outlier.size=0) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") +
    geom_jitter(position = position_jitter(width = 0.15)) +
    geom_hline(aes(yintercept=0), linetype="dashed", size=0.5) + 
    coord_flip() +
    scale_fill_manual(values = getPalette(colourCount))
}

PlotGenePressure(HpvInputPnps, HPVAnnotations)
PlotGenePressure(StaphInputPnps, StaphPhageAnnotations)
PlotGenePressure(PropInputPnps, PropPhageAnnotations)


# Save results as pdf
pdf(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/AnalysisOutput/HPV-HotspotEvolutionEachGene.pdf", width=12, height=10)
PlotGenePressure(HpvInputPnps, HPVAnnotations)
dev.off()

pdf(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/AnalysisOutput/StaphPhage-HotspotEvolutionEachGene.pdf", width=12, height=10)
PlotGenePressure(StaphInputPnps, StaphPhageAnnotations)
dev.off()

pdf(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/AnalysisOutput/PropPhage-HotspotEvolutionEachGene.pdf", width=12, height=10)
PlotGenePressure(PropInputPnps, PropPhageAnnotations)
dev.off()


