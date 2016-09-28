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
HPVEvolution <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/SnpHotspotEvolution/EvolutionFiles/HPVSnpHotspotPressure.gff3", header=TRUE, sep="\t")

StaphPhageEvolution <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/SnpHotspotEvolution/EvolutionFiles/StaphPhageSnpHotspotPressure.gff3", header=TRUE, sep="\t")


PlotEvoPressure <- function(HotspotEvoInput) {
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
  
  Median <- ddply(HotEvolution, c("Window"), summarize, median=median(PnPs))
  
  #return(list(Boxplot, WilcoxTest, Median))
  return(Boxplot)
}

# Now run through the subroutines
# Uncomment this and return median or WilcoxTest to get those values
# PlotEvoPressure(HPVEvolution)

HPVEvolutionForPlot <- data.frame(PlotEvoPressure(HPVEvolution))
HPVEvolutionForPlot$Name <- "HPV"

StaphPhageEvolutionForPlot <- data.frame(PlotEvoPressure(StaphPhageEvolution))
StaphPhageEvolutionForPlot$Name <- "StaphPhage"

BoundForPlot <- rbind(HPVEvolutionForPlot, StaphPhageEvolutionForPlot)

# Create the plot
MergeBox <- ggplot(BoundForPlot, aes(x=Name, y=PnPs, fill=Window)) + 
  theme_classic() + 
  geom_boxplot(notch=TRUE) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))
MergeBox

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

pdf(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/AnalysisOutput/Merged-HotspotEvolution.pdf", width=12, height=6)
MergeBox
dev.off()


