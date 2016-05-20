# CompareSnpLocations.R
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

library(plyr)
library(ggplot2)
library(reshape2)

# Load in the data frames
PrimaryDf <- read.delim("/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/Primary-UniqueSnps.vcf", sep="\t", header=FALSE)
PrimaryDf <- PrimaryDf[-nrow(PrimaryDf),]
PrimaryDf$V2 <- as.numeric(as.character(PrimaryDf$V2))
SecondaryDf <- read.delim("/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/Secondary-UniqueSnps.vcf", sep="\t", header=FALSE)
SecondaryDf <- SecondaryDf[-nrow(SecondaryDf),]
SecondaryDf$V2 <- as.numeric(as.character(SecondaryDf$V2))

# Merge together the dataframe columns
PrimaryDf$Merge <- paste(PrimaryDf$V1, PrimaryDf$V2, sep="")
PrimaryLengths <- ddply(PrimaryDf, c("V1"), summarize, length=length(V2), max=max(V2))
SecondaryDf$Merge <- paste(SecondaryDf$V1, SecondaryDf$V2, sep="")
SecondaryLengths <- ddply(SecondaryDf, c("V1"), summarize, length=length(V2), max=max(V2))

# Merge together to see how many are shared
MergedDf <- merge(PrimaryDf, SecondaryDf, by="Merge")
MergedStats <- ddply(MergedDf, c("V1.x"), summarize, length=length(V1.x))


###########################
# Generate Random Dataset #
###########################
Result <- lapply(PrimaryLengths$V1, function(i) {
  # Get length value
  LengthValue <- PrimaryLengths[c(PrimaryLengths$V1 == i),2]
  MaxValue <- PrimaryLengths[c(PrimaryLengths$V1 == i),3]
  RandomDf <- as.data.frame(sample(1:MaxValue, LengthValue, replace=FALSE))
  colnames(RandomDf) <- c("RandomPosition")
  RandomDf$Name <- i
  return(RandomDf)
})
PrimeDf <- do.call(rbind, Result)
PrimeDf$Merge <- paste(PrimeDf$Name, PrimeDf$RandomPosition, sep="")

Result <- lapply(SecondaryLengths$V1, function(i) {
  # Get length value
  LengthValue <- SecondaryLengths[c(SecondaryLengths$V1 == i),2]
  MaxValue <- SecondaryLengths[c(SecondaryLengths$V1 == i),3]
  RandomDf <- as.data.frame(sample(1:MaxValue, LengthValue, replace=FALSE))
  colnames(RandomDf) <- c("RandomPosition")
  RandomDf$Name <- i
  return(RandomDf)
})
SecondDf <- do.call(rbind, Result)
SecondDf$Merge <- paste(SecondDf$Name, SecondDf$RandomPosition, sep="")

# Merge the random data frame together
RandomMergedDf <- merge(PrimeDf, SecondDf, by="Merge")
RandomMergedStats <- ddply(RandomMergedDf, c("Name.x"), summarize, length=length(Name.x))

################################
# Comparing Random to Observed #
################################
# Merge together observed and random
MergedCompare <- merge(MergedStats, RandomMergedStats, by.x="V1.x", by.y="Name.x")
colnames(MergedCompare) <- c("ID","Obsv","Rand")
MergedCompareCalc <- merge(MergedCompare, PrimaryLengths, by.x="ID", by.y="V1")
MergedCompareCalc$RelObsv <- 100 * MergedCompareCalc$Obsv / MergedCompareCalc$length
MergedCompareCalc$RelRand <- 100 * MergedCompareCalc$Rand / MergedCompareCalc$length

# Subset for plotting
MergedSubset <- MergedCompareCalc[,c(1,6,7)]
MergedSubset <- melt(MergedSubset)

# Now plot comparing random to observed
SnpCompPlot <- ggplot(MergedSubset, aes(x=variable, y=value, fill=variable)) + 
  theme_classic() + 
  geom_boxplot(outlier.size=0, notch=TRUE) + 
  geom_jitter(position = position_jitter(width = .15)) + 
  scale_fill_brewer(palette="Paired") + 
  ylab("Percent Primary SNPs Matching Secondary") + 
  ggtitle("Comparing Percent Shared SNP Locations\nBetween Observed and Random")
SnpCompPlot

pdf(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/AnalysisOutput/SnpLocationComparison.pdf", width=12, height=10)
SnpCompPlot
dev.off()





