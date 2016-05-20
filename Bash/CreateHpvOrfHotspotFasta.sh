# CreateHpvOrfHotspotFasta.sh
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
echo Loading module EMBOSS-6.6.0...
module load EMBOSS-6.6.0

# Set the path and sample variables to be used throughout the script
export WorkingDirectory=/home/ghanni/Analysis/CleanRerunHumanVirome02
export Output='DeleteriousScoreReference'
export VcfFiles=/home/ghanni/Analysis/CleanRerunHumanVirome02/AlignedReadsOfInterestFromSRA/VcfResults/VarScanOverllHits-InterestOnlySRA-NoIndels.vcf
export PerlPath=~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/
export GitBin=~/git/Club_Grice/bin/

# HPV
export HpvContigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/HpvContigs.fa
export HpvOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/HpvOrfs.gff3
export HpvOrfSeqs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/HpvPredOrfs.fa
export HpvSnpHotspots=/home/ghanni/Analysis/CleanRerunHumanVirome02/SnpHotspotEvolution/EvolutionFiles/HPVSnpHotspotWindows.gff3
# Staph Phage
export StaphPhageContigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhageContigs.fa
export StaphPhageOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhageOrfs.gff3
export StaphPhageOrfSeqs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhagePredOrfs.fa
export StaphPhageSnpHotspots=/home/ghanni/Analysis/CleanRerunHumanVirome02/SnpHotspotEvolution/EvolutionFiles/StaphPhageSnpHotspotWindows.gff3


# Make the output directory and move to the working directory
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}

GetFasta () {
	# 1 = Phage Name
	# 2 = SnpHotspots
	# 3 = ORF Seqs

	# Get a gff3 of only the SNP hotspots
	grep "HotspotWindow" ${2} \
		| sed 's/_HotspotWindow_\S\+\t/\t/' \
		> ./${Output}/${1}SnpHotspotWindowsHotspotOnly.gff3
		
	# Create a list of the unique contig and ORF IDs
	cut -f 1,6 ./${Output}/${1}SnpHotspotWindowsHotspotOnly.gff3 \
		| sed 's/\t/_\./' \
		| sed 's/^/\>/' \
		| sort \
		| uniq \
		> ./${Output}/ListOfUniqueSnpHotspotOrfs.tsv
	
	# Pull out the ORF sequences that are known to contain SNP hotspots
	grep -A 1 --file=./${Output}/ListOfUniqueSnpHotspotOrfs.tsv ${3} \
		| grep -v "\-\-" \
		> ./${Output}/${1}HotspotOrfsOnly.fa
	
	# Translate those sequences
	transeq \
		-sequence ./${Output}/${1}HotspotOrfsOnly.fa \
		-outseq ./${Output}/${1}HotspotOrfsOnlyTranslated.fa \
		-clean
	
	# Remove block formatting from that resulting file
	perl \
		${GitBin}remove_block_fasta_format.pl \
			./${Output}/${1}HotspotOrfsOnlyTranslated.fa \
			./${Output}/${1}HotspotOrfsOnlyTranslatedNoBlock.fa
}

export -f GetFasta

GetFasta \
	"HPV" \
	${HpvSnpHotspots} \
	${HpvOrfSeqs}

GetFasta \
	"StaphPhage" \
	${StaphPhageSnpHotspots} \
	${StaphPhageOrfSeqs}

