# AlignCompareHanniganSegreScripts.sh
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
export Output='AlignCompareHanniganSegreScripts'
export PerlPath=~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/
export Rpath=~/git/Club_Grice/scripts/ghanni/analysis-R/HumanVirome02/
# Contig Files
export PrimaryStaphPhageContigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhageContigs.fa
export SegreStaphPhageContigs=/home/ghanni/Analysis/SegreVirome/TaxContigs/SpecificContigsByTaxa/StaphPhageContigs.fa
# Sequence Files
export PrimarySeqs=/home/ghanni/Analysis/CleanRerunHumanVirome02/AlignedReadsOfInterestFromSRA/tmp/CategoryFasta/InterestOnlySRA.fa
export SecondarySeqs=/home/ghanni/Analysis/SegreVirome/RayAssembly/_2_for_ray_pairs.fa
# Set the working directory
echo Setting working directory ${WorkingDirectory}...
cd ${WorkingDirectory}
# Make the output directory for the results here
echo Making output directory ${Output}...
mkdir ./${Output}
#=====>



#########################################
# Align ViromeContigs to Virome Contigs #
#########################################
mkdir ./${Output}/tmp
mkdir ./${Output}/tmp/BowtieBuilds/
mkdir ./${Output}/BowtieOut/
mkdir ./${Output}/BowtieContigs/

# Build bowtie reference database		
bowtie2-build -f ${SegreStaphPhageContigs} \
	./${Output}/tmp/BowtieBuilds/BowtieReferenceGenomes

# Do the alignment	
bowtie2 -x ./${Output}/tmp/BowtieBuilds/BowtieReferenceGenomes \
	-f ${PrimaryStaphPhageContigs} \
	-S ./${Output}/BowtieOut/StaphContigAlignment.sam \
	-L 25 \
	-N 1

# Get list of contigs with matches
perl ${PerlPath}GetMappedReadIdFromSam.pl \
	./${Output}/BowtieOut/StaphContigAlignment.sam \
	./${Output}/BowtieOut/StaphContigAlignment.tsv \
	./${Output}/BowtieOut/StaphContigAlignment.log

# Use this list (with the Grice contig IDs in the second column)
# To pull out the contigs of interest
cut -f 1 ./${Output}/BowtieOut/StaphContigAlignment.tsv \
	| sort \
	| uniq \
	| sed 's/^/\>/' \
	> ./${Output}/tmp/ListOfKeepers.tsv

grep -A 1 --file=./${Output}/tmp/ListOfKeepers.tsv ${PrimaryStaphPhageContigs} \
	| grep -v '\-\-' \
	> ./${Output}/tmp/MatchingContigs.fa


#######################################
# Align Sample Seqs to Shared Contigs #
#######################################
RunAlignment () {
	if [ -f /etc/profile.d/modules.sh ]; then
		source /etc/profile.d/modules.sh
	fi
	module load bowtie2-2.1.0
	module load ncbi-blast-2.2.0
	module load perl5lib
	module load samtools-1.1
	module load R-3.1.2
	# Set the path and sample variables to be used throughout the script
	export WorkingDirectory=/home/ghanni/Analysis/CleanRerunHumanVirome02
	export Output='AlignCompareHanniganSegreScripts'

	# Variable positions
	# 1 = TaxaName
	# 2 = Sample file path

	mkdir ./${Output}/tmp/BowtieBuilds/
	mkdir ./${Output}/BowtieOut/
	mkdir ./${Output}/BowtieOutBam/

	# Now align the reads from each sample to those same contigs
	# Build bowtie reference database
	mkdir ./${Output}/tmp/BowtieMatch/
	bowtie2-build -f ./${Output}/tmp/MatchingContigs.fa \
		./${Output}/tmp/BowtieMatch/BowtieReferenceGenomes
	
	# Do the alignment	
	bowtie2 -x ./${Output}/tmp/BowtieMatch/BowtieReferenceGenomes \
		-f ${2} \
		-S ./${Output}/tmp/BowtieMatch/${1}-Alignment.sam \
		-L 25 \
		-N 1
	
	# Remove unneeded sequences and convert to bam
	samtools view -bS \
		./${Output}/tmp/BowtieMatch/${1}-Alignment.sam \
		> ./${Output}/BowtieOut/${1}-OverallContigHits-InterestOnlySRA.bam	
	
	# Sort the Bam file
	samtools sort \
		./${Output}/BowtieOut/${1}-OverallContigHits-InterestOnlySRA.bam \
		./${Output}/BowtieOutBam/${1}-OverallContigHits_Sorted-InterestOnlySRA
	
	# Index the reference genome		
	samtools faidx \
		./${Output}/tmp/MatchingContigs.fa
	
	# Use VarScan to calculate the variable nucleotides 		
	samtools mpileup \
		-f ./${Output}/tmp/MatchingContigs.fa \
		./${Output}/BowtieOutBam/${1}-OverallContigHits_Sorted-InterestOnlySRA.bam \
		> ./${Output}/BowtieOut/${1}-SortedContigs-InterestOnlySRA.mpileup
	
	# Make output directory for VCF
	mkdir ./${Output}/VcfResults/
	
	java -jar \
		/project/egricelab/ghanni_software/VarScan/VarScan.v2.3.7.jar \
		pileup2snp \
		./${Output}/BowtieOut/${1}-SortedContigs-InterestOnlySRA.mpileup \
		> ./${Output}/VcfResults/${1}-VarScanOverllHits-InterestOnlySRA.vcf
	
	# Format each file reference name		
	sed -i 's/|//g' \
		./${Output}/VcfResults/${1}-VarScanOverllHits-InterestOnlySRA.vcf \
		| sed 's/\.//g' \
		| sed 's/ .*//g'
	
	# Remove indels to prevent issues downstream		
	perl \
		~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/RemoveIndelsFromVCF.pl \
		./${Output}/VcfResults/${1}-VarScanOverllHits-InterestOnlySRA.vcf \
		./${Output}/VcfResults/${1}-VarScanOverllHits-InterestOnlySRA-NoIndels.vcf
	
	#####################
	# Identify Hotspots #
	#####################
	# For now I can run this across all of the contigs and see what kind of variation I am generally seeing.
	mkdir ./${Output}/VariableSnpRegionFiles
	# First we need to get rid of the SNP duplicates in the VCF file
	# These are duplicates for our purposes at least
	cut -f 1,2 ./${Output}/VcfResults/${1}-VarScanOverllHits-InterestOnlySRA-NoIndels.vcf \
		| sort -V \
		| uniq \
		> ./${Output}/VariableSnpRegionFiles/${1}-UniqueSnps.vcf
	# Get together a single vcf file
	perl ${PerlPath}variabilityCallGeomDist.pl \
		./${Output}/VariableSnpRegionFiles/${1}-UniqueSnps.vcf \
		./${Output}/VariableSnpRegionFiles/${1}-ContigVariabilityDistance.tsv \
		./${Output}/VariableSnpRegionFiles/${1}-RawSnpDistanceForPlotting.tsv
	# Run R script to get quantile data from geometric distribution
	Rscript ${Rpath}GeomDistForVariantRegions.R \
		./${Output}/VariableSnpRegionFiles/${1}-ContigVariabilityDistance.tsv \
		./${Output}/VariableSnpRegionFiles/${1}-GeomDistQuantileResults.tsv
	# Run the second installment of the variability scripts using the R output
	perl ${PerlPath}variabilityCallGeomDist2.pl \
		./${Output}/VariableSnpRegionFiles/${1}-UniqueSnps.vcf \
		./${Output}/VariableSnpRegionFiles/${1}-GeomDistQuantileResults.tsv \
		./${Output}/VariableSnpRegionFiles/${1}-HypervariableLociResults.tsv

	# Copy unique SNP files out of this dir for easier transfer
	cp ./${Output}/VariableSnpRegionFiles/${1}-UniqueSnps.vcf ./${1}-UniqueSnps.vcf
}

# Export as a function
export -f RunAlignment

# Now run the scripts using the viruses of interest
RunAlignment \
	"Primary" \
	${PrimarySeqs}

RunAlignment \
	"Secondary" \
	${SecondarySeqs}


