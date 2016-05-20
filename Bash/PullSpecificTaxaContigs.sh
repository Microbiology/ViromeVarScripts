# PullSpecificTaxaContigs.sh
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania


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
export Output='TaxContigs'
export Contigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/AlignedReadsOfInterestFromSRA/HpvRefContigsNoBlocks.fa
export BinPath=~/git/Club_Grice/bin/
export ToolkitPath=~/git/Microbiome_sequence_analysis_toolkit/
export PerlPath=~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/
export Rpath=~/git/Club_Grice/scripts/ghanni/analysis-R/HumanVirome02/
export TaxList=/home/ghanni/Analysis/CleanRerunHumanVirome02/ContigStatsAndRefGenomeMatch/ContigBlastnVirusRefFormat.txt
export MasterGff=/home/ghanni/Analysis/CleanRerunHumanVirome02/AlignedReadsOfInterestFromSRA/PredictedOrfs/PredictedOrfsGlimmerOutputFormat.gff3
export MasterGeneSeq=/home/ghanni/Analysis/CleanRerunHumanVirome02/AlignedReadsOfInterestFromSRA/PredictedOrfs/PredictedOrfsGlimmerOutputGenes.fa
# Set the working directory
echo Setting working directory ${WorkingDirectory}...
cd ${WorkingDirectory}
# Make the output directory for the results here
echo Making output directory ${Output}...
mkdir ./${Output}
#=====>


####################################################
# Identifying Contigs That Match Reference Genomes #
####################################################
# Get a list of those contigs that match staph
mkdir ./${Output}/SpecificContigsByTaxa

sed 's/^/\>/' ${MasterGff} \
	| sed 's/\t/_\t/' \
	> ./${Output}/SpecificContigsByTaxa/FormatMaster.gff3

# ENAHuman_papillomavirus
grep 'ENAHuman_papillomavirus' ${TaxList} \
	| cut -f 1 \
	| sort \
	| uniq \
	| sed 's/^/\>/' \
	> ./${Output}/SpecificContigsByTaxa/HpvList.tsv

grep -A 1 -f ./${Output}/SpecificContigsByTaxa/HpvList.tsv ${Contigs} \
	| grep -v '\-\-' \
	> ./${Output}/SpecificContigsByTaxa/HpvContigs.fa

grep -f ./${Output}/SpecificContigsByTaxa/HpvList.tsv ./${Output}/SpecificContigsByTaxa/FormatMaster.gff3 \
	| grep -v '\-\-' \
	| sed 's/>//' \
	| sed 's/_//' \
	> ./${Output}/SpecificContigsByTaxa/HpvOrfs.gff3

grep -A 1 -f ./${Output}/SpecificContigsByTaxa/HpvList.tsv ${MasterGeneSeq} \
	| grep -v '\-\-' \
	> ./${Output}/SpecificContigsByTaxa/HpvPredOrfs.fa


# ENAPropionibacterium_phage
grep 'ENAPropionibacterium_phage' ${TaxList} \
	| cut -f 1 \
	| sort \
	| uniq \
	| sed 's/^/\>/' \
	> ./${Output}/SpecificContigsByTaxa/PropPhageList.tsv

grep -A 1 -f ./${Output}/SpecificContigsByTaxa/PropPhageList.tsv ${Contigs} \
	| grep -v '\-\-' \
	> ./${Output}/SpecificContigsByTaxa/PropPhageContigs.fa

grep -f ./${Output}/SpecificContigsByTaxa/PropPhageList.tsv ./${Output}/SpecificContigsByTaxa/FormatMaster.gff3 \
	| grep -v '\-\-' \
	| sed 's/>//' \
	| sed 's/_//' \
	> ./${Output}/SpecificContigsByTaxa/PropPhageOrfs.gff3

grep -A 1 -f ./${Output}/SpecificContigsByTaxa/PropPhageList.tsv ${MasterGeneSeq} \
	| grep -v '\-\-' \
	> ./${Output}/SpecificContigsByTaxa/PropPhagePredOrfs.fa


# ENAStaphylococcus_phage
grep 'ENAStaphylococcus_phage' ${TaxList} \
	| cut -f 1 \
	| sort \
	| uniq \
	| sed 's/^/\>/' \
	> ./${Output}/SpecificContigsByTaxa/StaphPhageList.tsv

grep -A 1 -f ./${Output}/SpecificContigsByTaxa/StaphPhageList.tsv ${Contigs} \
	| grep -v '\-\-' \
	> ./${Output}/SpecificContigsByTaxa/StaphPhageContigs.fa

grep -f ./${Output}/SpecificContigsByTaxa/StaphPhageList.tsv ./${Output}/SpecificContigsByTaxa/FormatMaster.gff3 \
	| grep -v '\-\-' \
	| sed 's/>//' \
	| sed 's/_//' \
	> ./${Output}/SpecificContigsByTaxa/StaphPhageOrfs.gff3

grep -A 1 -f ./${Output}/SpecificContigsByTaxa/StaphPhageList.tsv ${MasterGeneSeq} \
	| grep -v '\-\-' \
	> ./${Output}/SpecificContigsByTaxa/StaphPhagePredOrfs.fa
