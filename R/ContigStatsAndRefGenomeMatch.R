# ContigStatsAndRefGenomeMatch.R
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

# Load in libraries
library(ggplot2)

# Import Data Frame
ContigHits <- read.delim("/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/ContigStatsAndRefGenomeMatch/ContigBlastnVirusRefFormat.txt", header=FALSE, sep="\t")
ContigStats <- read.delim("/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/ContigStatsAndRefGenomeMatch/ContigStats/contig_length_with_coverage_for_graphing.tsv", header=TRUE, sep="\t")

# Provide total nucleotide length of the reference database
# Done by taking character count (wc) of reference fasta without titles
DatabaseLength <- 175373297

# Get the frequency information for the contig hits
ContigIdHits <- as.data.frame(table(ContigHits$V2))
ContigIdHitsSubsample <- ContigIdHits[c(ContigIdHits$Freq > 25),]

# Incorporate the sequence coverage information with the contigs
ContigMerge <- merge(ContigHits, ContigStats, by.x="V1", by.y="contig")
ContigMerge <- ContigMerge[which(ContigMerge$V2 %in% ContigIdHitsSubsample$Var1),]
ContigMerge$QueryCoverage <- 100 * ContigMerge$V4 / ContigMerge$contig_length
ContigMerge$SequenceCoverage <- 150 * ContigMerge$count / ContigMerge$contig_length

# Remove row duplicates
# Because ordered by e-value, taking first will provide best hit
ContigMerge$Ranked <- paste(ContigMerge$V1,ContigMerge$V2, sep="_")
ContigMerge <- ContigMerge[!duplicated(ContigMerge$Ranked), ]

# And now reclauclate the e-value because when blast does it, it stops at 0 for very low numbers
ContigMerge$RecalcEValue <- log10(ContigMerge$contig_length * DatabaseLength) - ContigMerge$V12*log10(2)

ContigQuerySeqs <- ggplot(ContigMerge, aes(x=(V4), y=SequenceCoverage, color=V2, size=(RecalcEValue*-1))) + 
  theme_classic() + 
  geom_point() + 
  scale_y_log10() + 
  scale_x_log10() +
  scale_colour_brewer(palette="Set2", name  ="Taxonomic Match") + 
  xlab("Length of Query Matching Reference Genome (log bp)") + 
  ylab("X Contig Coverage by Aligned Sequences (log)") + 
  scale_size_continuous(name  ="-1 * log10(E-Value)") +
  geom_hline(aes(yintercept=10), linetype="dashed") +
  geom_vline(aes(xintercept=750), linetype="dashed")
ContigQuerySeqs




pdf("/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/AnalysisOutput/ContigBlastnVirusRefFormat.pdf", width=8, height=6)
ContigQuerySeqs
dev.off()


