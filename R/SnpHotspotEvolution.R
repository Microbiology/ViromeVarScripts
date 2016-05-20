# SnpHotspotEvolution.R
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

# Set libraries
library(ggplot2)
library(directlabels)
library(pgirmess)
library(plyr)


# Read in the dataframe
HPVEvolution <- read.delim(file="/Users/Hannigan/Google Drive/VmTransfer/SnpHotspotEvolution/EvolutionFiles/HPVSnpHotspotPressure.gff3", header=TRUE, sep="\t")
HPVAnnotations <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/SnpHotspotEvolution/AnnotationFiles/HPVblastnResultsAnnoated.tsv", sep="\t", header=FALSE)
# Be sure to format some of the gene names
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"L1_protein","L1")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"L2_protein","L2")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"E2_protein","E2")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Replication_protein_E1","E1")
HPVAnnotations$V4 <- sub(x=HPVAnnotations$V4,"Major_capsid_protein_L1","L1")
HPVPfamAnnotations <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/PfamHmmerDomainAnnotation/PfamDomains/HPV-PfamDomainsFormat.tsv", sep="\t", header=FALSE)

StaphPhageEvolution <- read.delim(file="/Users/Hannigan/Google Drive/VmTransfer/SnpHotspotEvolution/EvolutionFiles/StaphPhageSnpHotspotPressure.gff3", header=TRUE, sep="\t")
StaphPhageAnnotations <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/SnpHotspotEvolution/AnnotationFiles/StaphPhageblastnResultsAnnoated.tsv", sep="\t", header=FALSE)
StaphPhagePfamAnnotations <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/PfamHmmerDomainAnnotation/PfamDomains/StaphPhage-PfamDomainsFormatReplace.tsv", sep="\t", header=FALSE)


PlotEvoPressure <- function(HotspotEvoInput, Annotations, PfamAnnotations) {
  HotEvolution <- HotspotEvoInput
  
  # Add column of what window the data point was taken from
  HotEvolution$Location <- sub(".*\\d_",
                               "", 
                               HotEvolution$GeneID, 
                               perl=TRUE)
  
  HotEvolution$Location <- sub("_.*",
                               "", 
                               HotEvolution$Location, 
                               perl=TRUE)
  
  HotEvolution$UniqueID <- paste(HotEvolution$ContigID, 
                                 sub("_.*",
                                     "", 
                                     HotEvolution$GeneID), 
                                 sep="_.")
  
  HotEvolution$Window <- ifelse(HotEvolution$Location %in% "HotspotWindow", 
                                "Hotspot", 
                                "Adjacent")
  
  HotEvolution$TotalSnps <- HotEvolution$SynonymousSNPs + HotEvolution$NonSynonymousSNPs
  
  HotEvolution$SnpID <- sub(".*_",
                            "", 
                            HotEvolution$GeneID, 
                            perl=TRUE)
  
  # Simply plot the pressures by Before, After, and Sliding Window
  WilcoxTest <- wilcox.test(formula=HotEvolution$PnPs~HotEvolution$Window)
  
  Boxplot <- ggplot(HotEvolution, 
                    aes(x=Window, 
                        y=PnPs, 
                        fill=Window)) + 
    theme_classic() + 
    geom_boxplot(notch=TRUE) + 
    scale_fill_brewer(palette="Set2")
  
  ViolinPlot <- ggplot(HotEvolution, 
                       aes(x=Window, 
                           y=PnPs, 
                           fill=Window
                           )
                       ) + 
    theme_classic() + 
    geom_violin(notch=TRUE) + 
    scale_fill_brewer(palette="Set2") + 
    annotate("text", 
             label = paste("p-value by wilcox:\n",
                           WilcoxTest$p.value, 
                           sep=""), 
             x = 1.5, 
             y = 0.25, 
             size = 5, 
             colour = "black") + 
    ggtitle("Evolutionary Pressure on SNP Hotspots\nAnd Adjacent Regions") + 
    geom_jitter(position = position_jitter(width = 0.15))
  
  
  
  # Now to look at only the hotspots
  HotEvolutionSub <- HotEvolution[c(which(HotEvolution$Location %in% "HotspotWindow")),]
  
  # And annotate the hotspots
  HotSpotMerge <- merge(HotEvolution,
                        Annotations,
                        by.x="UniqueID", 
                        by.y="V1")
  
  HotMergeSub <- HotSpotMerge[c(which(HotSpotMerge$Location %in% "HotspotWindow")),]
  
  HotSubOrder <- HotMergeSub[c(order(HotMergeSub$PnPs, decreasing=TRUE)),]
  
  AnnotatedPressureJitter <- ggplot(HotMergeSub, 
                                    aes(x="Evolutionary Pressure", 
                                        y=PnPs, 
                                        colour=V4)) + 
    theme_classic() + 
    geom_jitter(size=3, 
                position = position_jitter(width = 0.15)
                ) + 
    scale_fill_brewer(palette="Paired") + 
    ggtitle("Evolutionary Pressure on SNP Hotspots\nAnd Adjacent Regions")
  
  AnnotatedPressureDensity <- direct.label(ggplot(HotMergeSub, 
                                     aes(x=(PnPs), 
                                         color=V4,
                                         fill=V4)
                                     ) + 
    theme_classic() + 
    geom_density(alpha=0.2) + 
    ggtitle("Evolutionary Pressure on SNP Hotspots\nAnd Adjacent Regions")
  )
    
  BoxplotAnnotation <- ggplot(HotMergeSub, 
                        aes(x=reorder(V4, PnPs, FUN=median), 
                            y=PnPs, 
                            fill=V4)) + 
    theme_classic() + 
    geom_boxplot(notch=FALSE, outlier.size = 0) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") +
    geom_jitter(position = position_jitter(width = 0.15)) +
    coord_flip()
  
  # Check significance
  
  HotMergeSub$V4 <- factor(HotMergeSub$V4)
  
  KruskalSig <- kruskalmc(HotMergeSub$PnPs~HotMergeSub$V4)
  
  return(list(ViolinPlot,
              Boxplot))
}

# Now run through the subroutines

PlotEvoPressure(HPVEvolution,
                HPVAnnotations, 
                HPVPfamAnnotations)

PlotEvoPressure(StaphPhageEvolution, 
                StaphPhageAnnotations, 
                StaphPhagePfamAnnotations)

# Save results as pdf
pdf(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/AnalysisOutput/HPV-HotspotEvolution.pdf", width=12, height=10)
PlotEvoPressure(HPVEvolution,
                HPVAnnotations, 
                HPVPfamAnnotations)
dev.off()

pdf(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/AnalysisOutput/StaphPhage-HotspotEvolution.pdf", width=12, height=10)
PlotEvoPressure(StaphPhageEvolution, 
                StaphPhageAnnotations, 
                StaphPhagePfamAnnotations)
dev.off()


