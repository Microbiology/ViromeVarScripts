# GenomicEvolutionByGene.R
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

# Load in the libraries
library(ggplot2)
library(pgirmess)

# Due to bash processing, the first row needs to be removed out of the
# data frame.
HpvInputPnps <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/EvoByGeneByTaxa/PnPsResults/HpvOverallContigPressure.tsv", sep="\t", header=TRUE)
StaphInputPnps <- read.delim(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/TransferFiles/EvoByGeneByTaxa/PnPsResults/StaphPhageOverallContigPressure.tsv", sep="\t", header=TRUE)

# Bring together the taxa data frames for a comprehensive plot
HpvInputPnps$Category <- c("HPV")
StaphInputPnps$Category <- c("StaphPhage")
TotalPressure <- rbind(HpvInputPnps, 
                       StaphInputPnps)

TotalEvoByTaxa <- ggplot(TotalPressure, 
                         aes(x=log10(as.numeric(as.character(PnPs))), 
                             fill=factor(Category)
                             )
                         ) + 
  theme_classic() + 
  geom_density(alpha = 0.2) + 
  scale_fill_brewer(palette="Set2", 
                    name="Contig ID") + 
  xlab("Evolutionary Pressure (pNpS Ratio)") + 
  geom_vline(x=log10(1), 
             linetype="dashed") + 
  geom_text(aes(log10(1),
                0,
                label="Point of Neutral Evolution", 
                hjust=0, 
                vjust=1, 
                size=3.25)
            ) + 
  ggtitle("Trends in Natural Selection Among Natural Major Communities") + 
  scale_y_sqrt()


# Use boxplot to visualize gene values by site group
ViolinPressureByVirus <- ggplot(TotalPressure, 
                                aes(x=Category, 
                                    y=PnPs, 
                                    fill=Category)) + 
  theme_classic() + 
  geom_violin() + 
  geom_jitter(width=0.5) + 
  scale_fill_brewer(palette="Set2") + 
  scale_y_log10() + 
  ylab("pNpS Evolutionary Pressure By Gene (sqrt)") + 
  geom_text(aes(1,
                10,
                label="All Sig Diff with p<0.00001")) + 
  ggtitle("Differences in Evolutionary Pressure Between Virus Taxa")

BoxplotPressureByVirus <- ggplot(TotalPressure, 
                                 aes(x=Category, 
                                     y=PnPs, 
                                     fill=Category)) + 
  theme_classic() + 
  geom_boxplot(notch=TRUE) + 
  scale_fill_brewer(palette="Set2") + 
  scale_y_log10() + 
  ylab("pNpS Evolutionary Pressure By Gene (sqrt)") + 
  geom_text(aes(1,10,label="All Sig Diff with p<0.01")) + 
  ggtitle("Differences in Evolutionary Pressure Between Virus Taxa")

# Test for significant differences
kruskalmc(data=TotalPressure, 
          resp=PnPs~Category, 
          probs=0.05)

median <- ddply(TotalPressure, c("Category"), summarize, median=median(PnPs))
median

# Print out the plots
BoxplotPressureByVirus


# Print the plots to pdf
pdf(file="/Users/Hannigan/Documents/Geoffrey_Hannigan/The_lab/Projects/Human_Skin_Virome/HumanVirome02/AnalysisOutput/TotalGenomicEvolutionComparison.pdf", width=8, height=6)
BoxplotPressureByVirus
dev.off()

