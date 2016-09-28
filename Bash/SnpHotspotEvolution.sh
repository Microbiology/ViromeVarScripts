# SnpHotspotEvolution.sh
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
export Output='SnpHotspotEvolution'
export SnpHotspots=/home/ghanni/Analysis/CleanRerunHumanVirome02/variabilityCallGeomDist/VariableSnpRegionFiles/HypervariableLociResults.tsv
export VcfFiles=/home/ghanni/Analysis/CleanRerunHumanVirome02/AlignedReadsOfInterestFromSRA/VcfResults/VarScanOverllHits-InterestOnlySRA-NoIndels.vcf
export PerlPath=~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/
# HPV
export HpvContigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/HpvContigs.fa
export HpvOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/HpvOrfs.gff3
export HpvOrfSeqs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/HpvPredOrfs.fa
export HpvRefDb=/home/ghanni/Analysis/HumanVirome02/HpvSnpHotspotEvolution/AnnotationFiles/HPVdb
# Prop Phage
export PropPhageContigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/PropPhageContigs.fa
export PropPhageOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/PropPhageOrfs.gff3
export PropPhageOrfSeqs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/PropPhagePredOrfs.fa
export PropPhageRefDb=/project/egricelab/references/UniProt-Virus-Phage/uniprot_virus_and_phage_TrEMBL_db
# Staph Phage
export StaphPhageContigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhageContigs.fa
export StaphPhageOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhageOrfs.gff3
export StaphPhageOrfSeqs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhagePredOrfs.fa
export StaphPhageRefDb=/project/egricelab/references/UniProt-Virus-Phage/uniprot_virus_and_phage_TrEMBL_db
# Make the output directory and move to the working directory
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}


#####################################
# Genomic variability and evolution #
#####################################
# Run this as a subroutine
GenomeVariability () {
	if [ -f /etc/profile.d/modules.sh ]; then
		source /etc/profile.d/modules.sh
	fi
	module load bowtie2-2.1.0
	module load ncbi-blast-2.2.0
	module load perl5lib
	module load samtools-1.1
	module load R-3.1.2
	# Set the path and sample variables to be used throughout the script
	export WorkingDirectory=/home/ghanni/Analysis/HumanVirome02
	export Output='SnpHotspotEvolution'

	# # 1 = TaxaName
	# # 2 = Contig Files
	# # 3 = Snp Hotspots
	# # 4 = Contig ORFs
	# # 5 = VCF file
	# # 6 = ORF Fasta File
	# # 7 = Genome Reference for Blast
	# # The SNP hotspots were calculated already using my custom perl scripts, so now I can pull out the hotspots using the list of contig IDs
	# mkdir ./${Output}/GenomicVariability

	# grep '>' ${2} \
	# 	> ./${Output}/GenomicVariability/${1}ContigIdList.tsv

	# # Format the SNP hot spots
	# sed 's/^/>/' ${3} \
	# 	> ./${Output}/GenomicVariability/${1}HotSpotFormat.tsv

	# # Pull out the hotspots that match the contig ID list
	# grep \
	# 	--file=./${Output}/GenomicVariability/${1}ContigIdList.tsv \
	# 	./${Output}/GenomicVariability/${1}HotSpotFormat.tsv \
	# 	> ./${Output}/GenomicVariability/${1}SpecificSNPs.tsv
	
	# # Format the contig IDs here
	# sed -i 's/>//' ./${Output}/GenomicVariability/${1}SpecificSNPs.tsv
	# sed -i 's/_//' ./${Output}/GenomicVariability/${1}SpecificSNPs.tsv

	# # Pull out the SNPs.
	# awk '$5 > 1 {print $0}' \
	# 	./${Output}/GenomicVariability/${1}SpecificSNPs.tsv \
	# 	> ./${Output}/GenomicVariability/${1}SpecificSNPsLength.tsv
	
	# # Finally add a unique ID to each of the SNPs so that they can be individually identified later
	# awk ' { print $0"\t"FNR } ' \
	# 	./${Output}/GenomicVariability/${1}SpecificSNPsLength.tsv \
	# 	> ./${Output}/GenomicVariability/${1}SpecificSNPsLengthNumbered.tsv
	
	# # Evolution
	# # I will be able to use the HPV contig gff3 and SNP hotspot position files to generate gff3 files for SNP
	# # hotspots and the surrounding regions. I will be able to use the resulting gff3 file (with postiions
	# # relative to the HPV contigs) to determine what the evolutionary pressures are at and around the hotspots.
	# # First generate the gff3 files
	# mkdir ./${Output}/EvolutionFiles

	# perl ${PerlPath}HotspotSlidingWindowGff3.pl \
	# 	${4} \
	# 	./${Output}/GenomicVariability/${1}SpecificSNPsLengthNumbered.tsv \
	# 	./${Output}/EvolutionFiles/${1}SnpHotspotWindows.gff3
	
	# Use this with the other data to get evolutionary pressure
	perl ${PerlPath}CalculatePnPsRatio.pl \
		./${Output}/EvolutionFiles/${1}SnpHotspotWindows.gff3 \
		${2} \
		${5} \
		./${Output}/EvolutionFiles/${1}SnpHotspotPressure.gff3 \
		./${Output}/EvolutionFiles/${1}SnpHotspotSubPatterns.tsv \
		./${Output}/EvolutionFiles/${1}SnpHotspotCalls.vcf
}

# Export as a function
export -f GenomeVariability

# Now run the scripts using the viruses of interest
GenomeVariability \
	"HPV" \
	${HpvContigs} \
	${SnpHotspots} \
	${HpvOrfs} \
	${VcfFiles} \
	${HpvOrfSeqs} \
	${StaphPhageRefDb}

GenomeVariability \
	"PropPhage" \
	${PropPhageContigs} \
	${SnpHotspots} \
	${PropPhageOrfs} \
	${VcfFiles} \
	${PropPhageOrfSeqs} \
	${PropPhageRefDb}

GenomeVariability \
	"StaphPhage" \
	${StaphPhageContigs} \
	${SnpHotspots} \
	${StaphPhageOrfs} \
	${VcfFiles} \
	${StaphPhageOrfSeqs} \
	${StaphPhageRefDb}


