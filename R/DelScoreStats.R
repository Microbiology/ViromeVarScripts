# DelScoreStats.R
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

library(ggplot2)
library(RColorBrewer)
library(plyr)

HpvScores <- read.delim("/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/DeleteriousScores/HpvDeleteriousHotspots.tsv", header=TRUE, sep="\t")
HpvScores$Taxa <- "HPV"
StaphPhageScores <- read.delim("/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/DeleteriousScores/StaphPhageDeleteriousHotspots.tsv", header=TRUE, sep="\t")
StaphPhageScores$Taxa <- "StaphPhage"

summary(HpvScores)
summary(StaphPhageScores)

BoundScores <- rbind(HpvScores, StaphPhageScores)

wilcox.test(BoundScores$DeltScore ~ BoundScores$Taxa)

Density <- ggplot(BoundScores, aes(x=DeltScore, fill=Taxa)) + theme_classic() + geom_density(alpha=0.50) + scale_fill_brewer(palette="Set2") + ylab("Deleterious Score") + xlab("Density")
Density
BoxPlot <- ggplot(BoundScores, aes(x=Taxa, y=DeltScore, fill=Taxa)) + theme_classic() + geom_boxplot(notch=TRUE) + scale_fill_brewer(palette="Set2")
BoxPlot

ddply(BoundScores, c("Taxa"), summarize, median=median(DeltScore, na.rm=TRUE))

pdf("~/git/SchlossLab/notebook/Figures/2016-01/DelteriousScoresGricePaper.pdf", width=8, height=6)
Density
BoxPlot
dev.off()