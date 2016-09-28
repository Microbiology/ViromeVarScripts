# EvoByGeneByTaxa.sh
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

######################
# Prepare the script #
######################
# Load the needed modules
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi
echo Loading module bowtie2-2.1.0...
module load bowtie2-2.1.0
echo Loading module ncbi-blast-2.2.0...
module load ncbi-blast-2.2.0
echo Loading module perl5lib...
module load perl5lib
echo Loading module samtools-1.1...
module load samtools-1.1
echo Loading module R-3.1.2...
module load R-3.1.2

# Set the path and sample variables to be used throughout the script
export WorkingDirectory=/home/ghanni/Analysis/CleanRerunHumanVirome02
export Output='EvoByGeneByTaxa'
export SnpHotspots=/home/ghanni/Analysis/CleanRerunHumanVirome02/variabilityCallGeomDist/VariableSnpRegionFiles/HypervariableLociResults.tsv
export VcfFiles=/home/ghanni/Analysis/CleanRerunHumanVirome02/AlignedReadsOfInterestFromSRA/VcfResults/VarScanOverllHits-InterestOnlySRA-NoIndels.vcf
export PerlPath=~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/
# HPV
export HpvContigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/HpvContigs.fa
export HpvOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/HpvOrfs.gff3
# Prop Phage
export PropPhageContigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/PropPhageContigs.fa
export PropPhageOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/PropPhageOrfs.gff3
# Staph Phage
export StaphPhageContigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhageContigs.fa
export StaphPhageOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhageOrfs.gff3
# Make the output directory and move to the working directory
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}

#######################
# Calculate Evolution #
#######################
mkdir ./${Output}/PnPsResults

# HPV
perl ~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/CalculatePnPsRatio.pl \
	${HpvOrfs} \
	${HpvContigs} \
	${VcfFiles} \
	./${Output}/PnPsResults/HpvOverallContigPressure.tsv \
	./${Output}/PnPsResults/HpvSnpPatterns.tsv \
	./${Output}/PnPsResults/HpvSnpPatternsSecond.tsv

# Prop Phage
perl ~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/CalculatePnPsRatio.pl \
	${PropPhageOrfs} \
	${PropPhageContigs} \
	${VcfFiles} \
	./${Output}/PnPsResults/PropPhageOverallContigPressure.tsv \
	./${Output}/PnPsResults/PropPhageSnpPatterns.tsv \
	./${Output}/PnPsResults/PropPhageSnpPatternsSecond.tsv

# Staph Phage
perl ~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/CalculatePnPsRatio.pl \
	${StaphPhageOrfs} \
	${StaphPhageContigs} \
	${VcfFiles} \
	./${Output}/PnPsResults/StaphPhageOverallContigPressure.tsv \
	./${Output}/PnPsResults/StaphPhageSnpPatterns.tsv \
	./${Output}/PnPsResults/StaphPhageSnpPatternsSecond.tsv

