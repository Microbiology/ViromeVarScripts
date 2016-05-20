# CalculateDeleteriousScores.sh
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
export Output='DeleteriousScores'
export PerlPath=~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/

# HPV
export HpvContigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/HpvContigs.fa
export HpvOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/HpvOrfs.gff3
export HpvVcfSubset=/home/ghanni/Analysis/CleanRerunHumanVirome02/SnpHotspotEvolution/EvolutionFiles/HPVSnpHotspotCalls.vcf
export HpvDeleteriousScoreReference=/home/ghanni/Analysis/CleanRerunHumanVirome02/DeleteriousScores/HPVDelScoreResults.tsv

# HPV
export StaphPhageContigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhageContigs.fa
export StaphPhageOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhageOrfs.gff3
export StaphPhageVcfSubset=/home/ghanni/Analysis/CleanRerunHumanVirome02/SnpHotspotEvolution/EvolutionFiles/StaphPhageSnpHotspotCalls.vcf
export StaphPhageDeleteriousScoreReference=/home/ghanni/Analysis/CleanRerunHumanVirome02/DeleteriousScores/StaphPhageDelScoreResults.tsv



# Make the output directory and move to the working directory
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}

CalcDelScores () {
	# Define variable positions
	# 1 = File Set Name
	# 2 = VCF Subset
	# 3 = ORF gff3
	# 4 = Contig Seqs
	# 5 = Del Score Reference

	# First get the VCF subset file to contain only the hotspots and not the adjacent regions
	grep "HotspotWindow" ${2} \
		| sed 's/\tUID.*//' \
		> ./${Output}/${1}HotspotCallsOnly.vcf
	# And get one with the adjacent regions only
	grep -v "HotspotWindow" ${2} \
		| sed 's/\tUID.*//' \
		> ./${Output}/${1}AdjacentCallsOnly.vcf
	
	# This is just going to be a simple run of a perl script, nothing fancy
	perl ${PerlPath}CalculateDeleteriousScores.pl \
		./${Output}/${1}HotspotCallsOnly.vcf \
		${3} \
		${4} \
		${5} \
		./${Output}/${1}DeleteriousHotspots.tsv
	
	sed -i 's/$/\tHotspot/' ./${Output}/${1}DeleteriousHotspots.tsv
}

export -f CalcDelScores

CalcDelScores \
	"Hpv" \
	${HpvVcfSubset} \
	${HpvOrfs} \
	${HpvContigs} \
	${HpvDeleteriousScoreReference}

CalcDelScores \
	"StaphPhage" \
	${StaphPhageVcfSubset} \
	${StaphPhageOrfs} \
	${StaphPhageContigs} \
	${StaphPhageDeleteriousScoreReference}
