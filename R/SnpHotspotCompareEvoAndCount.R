# SnpHotspotCompareEvoAndCount.R
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

# Set libraries
library(ggplot2)
library(plyr)
library(tuple)

# Read in the dataframes
StaphHotEvolution <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/SnpHotspotEvolution/EvolutionFiles/StaphPhageSnpHotspotPressure.gff3", header=TRUE, sep="\t")
HpvHotEvolution <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/SnpHotspotEvolution/EvolutionFiles/HPVSnpHotspotPressure.gff3", header=TRUE, sep="\t")

###############
# Staph Phage #
###############
# Add column of what window the data point was taken from
StaphHotEvolution$Location <- sub(".*\\d_",
                                  "", 
                                  StaphHotEvolution$GeneID, 
                                  perl=TRUE)

StaphHotEvolution$Location <- sub("_.*",
                                  "", 
                                  StaphHotEvolution$Location, 
                                  perl=TRUE)

StaphHotEvolution$UniqueID <- paste(StaphHotEvolution$ContigID, 
                                    sub("_.*",
                                        "", 
                                        StaphHotEvolution$GeneID), 
                                    sep="_.")

StaphHotEvolution$Window <- ifelse(StaphHotEvolution$Location %in% "HotspotWindow", 
                                   "Hotspot", 
                                   "Adjacent")

StaphHotEvolution$TotalSnps <- StaphHotEvolution$SynonymousSNPs + StaphHotEvolution$NonSynonymousSNPs

StaphHotEvolution$SnpID <- sub(".*_",
                               "", 
                               StaphHotEvolution$GeneID, 
                               perl=TRUE)

StaphHotEvolution$ID <- "StaphPhage"

#######
# HPV #
#######
# Add column of what window the data point was taken from
HpvHotEvolution$Location <- sub(".*\\d_",
                                "", 
                                HpvHotEvolution$GeneID, 
                                perl=TRUE)

HpvHotEvolution$Location <- sub("_.*",
                                "", 
                                HpvHotEvolution$Location, 
                                perl=TRUE)

HpvHotEvolution$UniqueID <- paste(HpvHotEvolution$ContigID, 
                                  sub("_.*",
                                      "", 
                                      HpvHotEvolution$GeneID), 
                                  sep="_.")

HpvHotEvolution$Window <- ifelse(HpvHotEvolution$Location %in% "HotspotWindow", 
                                 "Hotspot", 
                                 "Adjacent")

HpvHotEvolution$TotalSnps <- HpvHotEvolution$SynonymousSNPs + HpvHotEvolution$NonSynonymousSNPs

HpvHotEvolution$SnpID <- sub(".*_",
                             "", 
                             HpvHotEvolution$GeneID, 
                             perl=TRUE)

HpvHotEvolution$ID <- "HPV"

######################
# Merging & Analysis #
######################

# Total SNP Count for comparison
CatHotEvolution <- rbind(StaphHotEvolution,
                         HpvHotEvolution)

CatHotspotSubset <- CatHotEvolution[c(which(CatHotEvolution$Window %in% "Hotspot")),]

CatCount <- ddply(CatHotspotSubset, 
                  c("ID"), 
                  summarize, 
                  length=length(ID))

# Test significance
wilcoxResult <- wilcox.test(CatHotspotSubset$PnPs~CatHotspotSubset$ID)
wilcoxResult

median <- ddply(CatHotspotSubset, c("ID"), summarize, median=median(PnPs))

# Plot out the numbers of hotspots
SnpCountComp <- ggplot(CatHotspotSubset, aes(x=ID, y=PnPs, fill=ID)) +
  theme_classic() +
  geom_boxplot(notch=TRUE) +
  scale_colour_brewer(palette="Set1") +
  ggtitle("Hotspot Evolution Between Viruses") +
  annotate("text", 
           label = paste("p-value by wilcox:\n",
                         wilcoxResult$p.value, 
                         sep=""), 
           x = 1.5, 
           y = 1.75, 
           size = 5, 
           colour = "black") +
  scale_fill_brewer(palette="Paired")
SnpCountComp

# Save results in PDF

pdf(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/AnalysisOutput/HotspotEvolutionBetweenViruses.pdf", width=8, height=6)
SnpCountComp
dev.off()



