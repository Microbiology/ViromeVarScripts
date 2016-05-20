# HotspotSubstitutionPatterns.R
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

# Set the libraries
library(ggplot2)
library(plyr)
library(seqinr)
library(reshape2)

# Read in the files
HpvSubPat <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/SnpHotspotEvolution/EvolutionFiles/HPVSnpHotspotSubPatterns.tsv", sep="\t", header=TRUE)
StaphPhageSubPat <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/SnpHotspotEvolution/EvolutionFiles/StaphPhageSnpHotspotSubPatterns.tsv", sep="\t", header=TRUE)
# This file path is the same for all versions of the script
AminoProfile <- read.delim("/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/AminoAcidProps.tsv", sep="\t", header=FALSE)
data(pK)


PlotSubstPatterns <- function(SubPatterns, AminoProfile) {
  ###################
  # Format the Data #
  ###################
  
  # Get only the SNP hotspots
  SubPatterns$Location <- sub(".*\\d_","",
                              SubPatterns$GeneID,
                              perl=TRUE)
  
  # Format the location values
  SubPatterns$Location <- sub("_.*",
                              "",
                              SubPatterns$Location,
                              perl=TRUE)
  
  # Create value of contig ID with Gene ID
  SubPatterns$UniqueID <- paste(SubPatterns$ContigID,
                                sub("_.*",
                                    "",
                                    SubPatterns$GeneID),
                                sep="_.")
  
  # Binary value for whether it is hotspot or ajacent region
  SubPatterns$Window <- ifelse(SubPatterns$Location %in% "HotspotWindow",
                               "Hotspot",
                               "Adjacent")
  
  # Only look at the hotspots
  SubPatternsCut <- SubPatterns[
    c(which(SubPatterns$Window %in% "Hotspot"))
    ,]
  
  # Set whether the consensus amino acid values are charged or not
  SubPatternsCut$ChargeCons <- ifelse(SubPatternsCut$CodonAAConsensus %in% row.names(pK),
                                      "Yes",
                                      "No")
  
  # Set whether the variable amino acid values are charged or not
  SubPatternsCut$ChargeVar <- ifelse(SubPatternsCut$CodonAAVariable %in% row.names(pK),
                                     "Yes",
                                     "No")
  
  # Merge amion acid properties into subst pattern scripts - Consensus
  SubPatternsCut <- merge(SubPatternsCut,
                          AminoProfile,
                          by.x="CodonAAConsensus",
                          by.y="V3")
  
  # Merge amino acid properties into subst pattern scripts - Variable
  SubPatternsCut <- merge(SubPatternsCut,
                          AminoProfile,
                          by.x="CodonAAVariable",
                          by.y="V3")
  
  ##################################
  # Get Rel Abund DFs For Plotting #
  ##################################
  
  # Nucleotide pattern count
  SnpScanNuc <- ddply(SubPatternsCut,
                         c("SnpConsensus","SnpVariable"),
                         summarize,
                         count=length(SnpConsensus))
  SnpScanNuc$RelAbund <- 100 * SnpScanNuc$count / sum(SnpScanNuc$count)
  
  # Overall amino acid pattern count
  SnpScanUniqueAA <- ddply(SubPatternsCut,
                           c("CodonAAConsensus","CodonAAVariable"),
                           summarize,
                           count=length(CodonAAConsensus))
  SnpScanUniqueAA$RelAbund <- 100 * SnpScanUniqueAA$count / sum(SnpScanUniqueAA$count)
  
  # Charged amino acid pattern count
  SnpScanAACharge <- ddply(SubPatternsCut,
                           c("ChargeCons","ChargeVar"),
                           summarize,
                           count=length(ChargeCons))
  SnpScanAACharge$RelAbund <- 100 * SnpScanAACharge$count / sum(SnpScanAACharge$count)
  
  # Polarity amino acid pattern count
  SnpScanPolarAA <- ddply(SubPatternsCut,
                          c("V4.x","V4.y"),
                          summarize,
                          count=length(V4.x))
  SnpScanPolarAA$RelAbund <- 100 * SnpScanPolarAA$count / sum(SnpScanPolarAA$count)
  
  # Specific polarity amino acid pattern count
  SnpScanSpecPolarAA <- ddply(SubPatternsCut,
                              c("V5.x","V5.y"),
                              summarize,
                              count=length(V5.x))
  SnpScanSpecPolarAA$RelAbund <- 100 * SnpScanSpecPolarAA$count / sum(SnpScanSpecPolarAA$count)
  
  # Specific charged amino acid pattern count
  SnpScanUniqueAACharge <- ddply(SubPatternsCut,
                                 c("V6.x","V6.y"),
                                 summarize,
                                 count=length(V6.x))
  SnpScanUniqueAACharge$RelAbund <- 100 * SnpScanUniqueAACharge$count / sum(SnpScanUniqueAACharge$count)
  
  # Calculate the transition and transversion ratio
  SnpScanNuc$TsTvLabel <- ifelse(
    test=(SnpScanNuc$SnpConsensus %in% c("A","G") & SnpScanNuc$SnpVariable %in% c("C","T")),
    no=ifelse(
      test=(SnpScanNuc$SnpConsensus %in% c("C","T") & SnpScanNuc$SnpVariable %in% c("G","A")),
      yes="Transversion",
      no="Transition"),
    yes="Transversion")
  
  SnpTsTv <- ddply(SnpScanNuc,
                   c("TsTvLabel"),
                   summarize,
                   sum=sum(RelAbund))
  
  TsTvRatio <- 2*SnpTsTv[2,2]/SnpTsTv[1,2]
    
  #############
  # Plot Data #
  #############
  
  SnpStatsPlotNuc <- ggplot(SnpScanNuc,
                            aes(x=SnpConsensus,
                                y=SnpVariable,
                                fill=(RelAbund)
                                )
                            ) +
    theme_classic() + 
    geom_tile(colour="white") + 
    scale_fill_gradient(low = "white", 
                        high = "black",
                        limit=c(0,18)) + 
    xlab("Consensus SNP") + 
    ylab("Variable SNP") + 
    ggtitle("Patterns of Hotspot Amino Acid Mutations") + 
    geom_abline(slope=1)  
  
  # Plot rel abund as exp to highlight most abudnant
  SnpStatsPlotAA <- ggplot(SnpScanUniqueAA, 
                           aes(x=CodonAAConsensus, 
                               y=CodonAAVariable, 
                               fill=exp(RelAbund)
                               )
                           ) + 
    theme_classic() + 
    geom_tile(colour="white", 
              size=1) + 
    scale_fill_gradient(low = "white", 
                        high = "black",
                        limit=c(0,60)) + 
    xlab("Consensus AA") + 
    ylab("Variable AA") + 
    ggtitle("Patterns of Hotspot Amino Acid Mutations") + 
    geom_abline(slope=1)
   
  SnpChargePlotAA <- ggplot(SnpScanAACharge,
                            aes(x=ChargeCons, 
                                y=ChargeVar, 
                                fill=(RelAbund))) + 
    theme_classic() + 
    geom_tile(colour="white", 
              size=1) + 
    scale_fill_gradient(low = "white", 
                        high = "black") + 
    xlab("Consensus AA") + ylab("Variable AA") + 
    ggtitle("Patterns of Hotspot Amino Acid Mutations") + 
    geom_abline(slope=1)
  
  SnpChargePlotAAPolar1 <- ggplot(SnpScanPolarAA, 
                                  aes(x=V4.x, 
                                      y=V4.y, 
                                      fill=log10(RelAbund)
                                      )
                                  ) + 
    theme_classic() + 
    geom_tile(colour="white", 
              size=1) + 
    scale_fill_gradient(low = "white", 
                        high = "black",
                        limit=c(-0.5,2)) + 
    xlab("Consensus AA") + 
    ylab("Variable AA") + 
    ggtitle("Patterns of Hotspot Amino Acid Mutations") + 
    geom_abline(slope=1)

  SnpChargePlotAAPolar2 <- ggplot(SnpScanSpecPolarAA, 
                                  aes(x=V5.x, 
                                      y=V5.y, 
                                      fill=(RelAbund)
                                      )
                                  ) + 
    theme_classic() + 
    geom_tile(colour="white", 
              size=1) + 
    scale_fill_gradient(low = "white", 
                        high = "black") + 
    xlab("Consensus AA") + 
    ylab("Variable AA") + 
    ggtitle("Patterns of Hotspot Amino Acid Mutations") + 
    geom_abline(slope=1)
  
  SnpChargePlotAAChargeExp <- ggplot(SnpScanUniqueAACharge, 
                                     aes(x=V6.x, 
                                         y=V6.y, 
                                         fill=log10(RelAbund)
                                         )
                                     ) + 
    theme_classic() + 
    geom_tile(colour="white", 
              size=1) + 
    scale_fill_gradient(low = "white", 
                        high = "black",
                        limit=c(-0.5,2)) + 
    xlab("Consensus AA") + 
    ylab("Variable AA") + 
    ggtitle("Patterns of Hotspot Amino Acid Mutations") + 
    geom_abline(slope=1)  
  
  ####################
  # Export Plot Data #
  ####################
  return(list(SnpStatsPlotNuc,
              SnpStatsPlotAA,
              SnpChargePlotAA,
              SnpChargePlotAAPolar1,
              SnpChargePlotAAPolar2,
              SnpChargePlotAAChargeExp,
              TsTvRatio,
              SnpScanPolarAA))
}

# Run through the subroutine with the files of interest
PlotSubstPatterns(HpvSubPat, AminoProfile)
PlotSubstPatterns(StaphPhageSubPat, AminoProfile)



##########################
# Calculate Significance #
##########################
formatData <- function(SubPatterns, AminoProfile) {
  # Get only the SNP hotspots
  SubPatterns$Location <- sub(".*\\d_","",
                              SubPatterns$GeneID,
                              perl=TRUE)
  
  # Format the location values
  SubPatterns$Location <- sub("_.*",
                              "",
                              SubPatterns$Location,
                              perl=TRUE)
  
  # Create value of contig ID with Gene ID
  SubPatterns$UniqueID <- paste(SubPatterns$ContigID,
                                sub("_.*",
                                    "",
                                    SubPatterns$GeneID),
                                sep="_.")
  
  # Binary value for whether it is hotspot or ajacent region
  SubPatterns$Window <- ifelse(SubPatterns$Location %in% "HotspotWindow",
                               "Hotspot",
                               "Adjacent")
  
  # Only look at the hotspots
  SubPatternsCut <- SubPatterns[
    c(which(SubPatterns$Window %in% "Hotspot"))
    ,]
  
  # Set whether the consensus amino acid values are charged or not
  SubPatternsCut$ChargeCons <- ifelse(SubPatternsCut$CodonAAConsensus %in% row.names(pK),
                                      "Yes",
                                      "No")
  
  # Set whether the variable amino acid values are charged or not
  SubPatternsCut$ChargeVar <- ifelse(SubPatternsCut$CodonAAVariable %in% row.names(pK),
                                     "Yes",
                                     "No")
  
  # Merge amion acid properties into subst pattern scripts - Consensus
  SubPatternsCut <- merge(SubPatternsCut,
                          AminoProfile,
                          by.x="CodonAAConsensus",
                          by.y="V3")
  
  # Merge amion acid properties into subst pattern scripts - Variable
  SubPatternsCut <- merge(SubPatternsCut,
                          AminoProfile,
                          by.x="CodonAAVariable",
                          by.y="V3")
  return(SubPatternsCut)
}



# Create subroutine to create count data frames
CreateCountDF <- function(Input, Var1, Var2, Title1) {
  CountDf <- ddply(Input, 
                   c(Var1,Var2), 
                   .fun = function(xx){
                     length = length(xx[,Var1])
                   }
    )
  colnames(CountDf) <- c(Var1, Var2, Title1)
  # Only look at changes in value
  CountDf <- CountDf[CountDf[,Var1] != CountDf[,Var2],]
  CountDf$Pattern <- paste(CountDf[,Var1], 
                           CountDf[,Var2], 
                           sep="-")
  CountDf <- CountDf[,c("Pattern",Title1)]
  return(CountDf)
}
test <- as.data.frame(formatData(StaphPhageSubPat, AminoProfile))
CreateCountDF(test, "V6.x", "V6.y", "HPV")

SubstPatternSignf <- function(SubPatterns1, SubPatterns2, AminoProfile, Title1, Title2, Variable1, Variable2) {
  ###################
  # Format the Data #
  ###################
  SubPatternsCut1 <- as.data.frame(formatData(SubPatterns1, AminoProfile))
  SubPatternsCut2 <- as.data.frame(formatData(SubPatterns2, AminoProfile))
  
  ##################################
  # Get Rel Abund DFs For Plotting #
  ##################################
  SnpScanAAChargeMelt1 <- as.data.frame(CreateCountDF(SubPatternsCut1, 
                                                      Variable1, 
                                                      Variable2, 
                                                      Title1))
  SnpScanAAChargeMelt2 <- as.data.frame(CreateCountDF(SubPatternsCut2, 
                                                      Variable1, 
                                                      Variable2, 
                                                      Title2))
  # Merge the files  
  ChargeMeltMerge <- merge(SnpScanAAChargeMelt1, SnpScanAAChargeMelt2, by="Pattern", all=TRUE)
  ChargeMeltMerge[is.na(ChargeMeltMerge)] <- 0
  ChargeMeltMerge <- as.table(as.matrix(ChargeMeltMerge[,-1]))
  
  result <- chisq.test(ChargeMeltMerge, simulate.p.value=TRUE, B=1e5)
  
  return(list(SnpScanAAChargeMelt2, ChargeMeltMerge, result))
}

# Run through the subroutine with the files of interest
# NOTE: These p values are from Mont Carlo simulations so minor variation exists
SubstPatternSignf(HpvSubPat, StaphPhageSubPat, AminoProfile, "HPV", "Staph", "ChargeCons", "ChargeVar")
SubstPatternSignf(HpvSubPat, StaphPhageSubPat, AminoProfile, "HPV", "Staph", "V6.x", "V6.y")
SubstPatternSignf(HpvSubPat, StaphPhageSubPat, AminoProfile, "HPV", "Staph", "V5.x", "V5.y")
SubstPatternSignf(HpvSubPat, StaphPhageSubPat, AminoProfile, "HPV", "Staph", "V4.x", "V4.y")
SubstPatternSignf(HpvSubPat, StaphPhageSubPat, AminoProfile, "HPV", "Staph", "CodonAAConsensus", "CodonAAVariable")
SubstPatternSignf(HpvSubPat, StaphPhageSubPat, AminoProfile, "HPV", "Staph", "SnpConsensus", "SnpVariable")





# Save plots to PDF
pdf(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/AnalysisOutput/HPV-HotspotSubstPatterns.pdf", width=8, height=6)
PlotSubstPatterns(HpvSubPat, AminoProfile)
dev.off()

pdf(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/AnalysisOutput/StaphPhage-HotspotSubstPatterns.pdf", width=8, height=6)
PlotSubstPatterns(StaphPhageSubPat, AminoProfile)
dev.off()



