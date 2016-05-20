# ContigStatsAndReferenceGenomeHits.sh
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
module load ncbi-blast-2.2.0
module load bowtie2-2.1.0
module load perl5lib

#=====>
# Set the path and sample variables to be used throughout the script
export WorkingDirectory=/home/ghanni/Analysis/CleanRerunHumanVirome02
export Output='ContigStatsAndRefGenomeMatch'
export Contigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/AlignedReadsOfInterestFromSRA/HpvRefContigsNoBlocks.fa
export Reference=/project/egricelab/references/virus_and_phage_ref/virus_and_phage_db
export SampleSam=/home/ghanni/Analysis/CleanRerunHumanVirome02/AlignedReadsOfInterestFromSRA/BowtieOut/InterestOnlySRA.sam
# Set the working directory
echo Setting working directory ${WorkingDirectory}...
cd ${WorkingDirectory}
# Make the output directory for the results here
echo Making output directory ${Output}...
mkdir ./${Output}
#=====>



#########################################
# Get Abundance Data From Bowtie Output #
#########################################
# Make file directories
mkdir ./${Output}/ContigStats

#Get abundance data from the bowtie output
perl \
	/project/egricelab/bin/calculate_abundance_from_sam.pl \
	${SampleSam} \
	./${Output}/InterestOnlySRASam.txt

#Generate table with the sequence length numbers for each contig
awk 'NR % 2 {printf $0"\t"} !(NR % 2) {print length($0)}' \
	${Contigs} \
	> ./${Output}/ContigStats/contig_length.txt

sed 's/>//g' \
	./${Output}/ContigStats/contig_length.txt \
	> ./${Output}/ContigStats/contig_length_without_greater_sign.txt

#Add headers
awk 'BEGIN{print "contig\tcontig_length"}1' \
	./${Output}/ContigStats/contig_length_without_greater_sign.txt \
	> ./${Output}/ContigStats/contig_length_without_greater_sign_with_header.txt

#Merge the contig length and bowtie hit values into a single tab-delimited file
awk 'FNR==NR { a[$1]=$2; next } $1 in a { print $1"\t"$2"\t"a[$1] }' \
	./${Output}/ContigStats/contig_length_without_greater_sign_with_header.txt \
	./${Output}/InterestOnlySRASam.txt \
	> ./${Output}/ContigStats/contig_length_with_coverage_for_graphing.tsv



####################################
# Blast Contigs to Virus Reference #
####################################
# Perform blast to get the contigs that match virus reference genomes.
# This will be done using the overall contigs calculated from the previous virome study.
blastn \
	-query ${Contigs} \
	-out ./${Output}/ContigBlastnVirusRef.txt \
	-db /project/egricelab/references/virus_and_phage_ref/virus_and_phage_db \
	-evalue 1e-3 \
	-max_target_seqs 1 \
	-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen'

# Format the contig blast hits so that they are easier to interpret downstream
# WARNING: This is specific for this reference database, so be careful using a new database as a reference
sed 's/\(ENA[^_]*_[^_]*\)_[^\t]\+\t/\1\t/' \
	./${Output}/ContigBlastnVirusRef.txt \
	> ./${Output}/ContigBlastnVirusRefFormat.txt




