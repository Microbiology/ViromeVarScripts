# IdentifyHotspotsGeomDist.sh
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

######################
# Prepare the script #
######################
# Load in the needed modules
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
export Output='variabilityCallGeomDist'
export Contigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/AlignedReadsOfInterestFromSRA/HpvRefContigsNoBlocks.fa
export PerlPath=~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/
export VcfFiles=/home/ghanni/Analysis/CleanRerunHumanVirome02/AlignedReadsOfInterestFromSRA/VcfResults/VarScanOverllHits-InterestOnlySRA-NoIndels.vcf
export Rpath=~/git/Club_Grice/scripts/ghanni/analysis-R/HumanVirome02/
# Set the working directory
echo Setting working directory ${WorkingDirectory}...
cd ${WorkingDirectory}
# Make the output directory for the results here
echo Making output directory ${Output}...
mkdir ./${Output}
#=====>

###################################################################
# Detect Hypervariable Regions Using Geometric Distribution Stats #
###################################################################
# Here I want to determine what genomic regions have the highest levels of variation across contigs
# I can use my custom perl script which calculates the number of SNPs within all windows of base pairs
# with a set length.

# For now I can run this across all of the contigs and see what kind of variation I am generally seeing.
mkdir ./${Output}/VariableSnpRegionFiles
# Use my perl script for filtering out all of the indels, leaving only the SNPs
perl ${PerlPath}RemoveIndelsFromVCF.pl \
	${VcfFiles} \
	./${Output}/VariableSnpRegionFiles/UniqueSnpsNoIndels.vcf
# First we need to get rid of the SNP duplicates in the VCF file
# These are duplicates for our purposes at least
cut -f 1,2 ./${Output}/VariableSnpRegionFiles/UniqueSnpsNoIndels.vcf \
	| sort -V \
	| uniq \
	> ./${Output}/VariableSnpRegionFiles/UniqueSnps.vcf
# Get together a single vcf file
perl ${PerlPath}variabilityCallGeomDist.pl \
	./${Output}/VariableSnpRegionFiles/UniqueSnps.vcf \
	./${Output}/VariableSnpRegionFiles/ContigVariabilityDistance.tsv \
	./${Output}/VariableSnpRegionFiles/RawSnpDistanceForPlotting.tsv
# Run R script to get quantile data from geometric distribution
Rscript ${Rpath}GeomDistForVariantRegions.R \
	./${Output}/VariableSnpRegionFiles/ContigVariabilityDistance.tsv \
	./${Output}/VariableSnpRegionFiles/GeomDistQuantileResults.tsv
# Run the second installment of the variability scripts using the R output
perl ${PerlPath}variabilityCallGeomDist2.pl \
	./${Output}/VariableSnpRegionFiles/UniqueSnps.vcf \
	./${Output}/VariableSnpRegionFiles/GeomDistQuantileResults.tsv \
	./${Output}/VariableSnpRegionFiles/HypervariableLociResults.tsv \
	./${Output}/VariableSnpRegionFiles/HypervariableLociVariants.vcf
# Pull out the distances of HPV contigs for plotting
grep --file=../HpvPhylogeny/GenomicVariability/ContigIdList.tsv  \
	> ./${Output}/VariableSnpRegionFiles/SpecificSnpDistances.tsv
