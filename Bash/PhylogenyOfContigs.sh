# PhylogenyOfContigs.sh
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
export Output='Phylogeny'

# HPV
export HpvOrfSeqs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/HpvPredOrfs.fa
export HpvReferenceSeqs=/home/ghanni/Analysis/HumanVirome02/HpvPhylogeny/referenceSeqs/PaveHpvL1ReferenceOf2015-06-03-Format.fa

# Prop Phage
export PropPhageOrfSeqs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/PropPhagePredOrfs.fa
export PropPhageReferenceSeqs=/home/ghanni/Analysis/HumanVirome02/HpvPhylogeny/referenceSeqs/PropPhageRefGenomesPhylogeny-LargeTerminaseSubunit-2015-09-15-Format.fa

# Staph Phage
export StaphPhageOrfSeqs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhagePredOrfs.fa
export StaphPhageReferenceSeqs=/home/ghanni/Analysis/HumanVirome02/HpvPhylogeny/referenceSeqs/StaphPhageRefGenomesPhylogeny-LargeTerminaseSubunit-2015-09-14-Format.fa

# Other helpful paths
export BinPath=~/git/Club_Grice/bin/
export ToolkitPath=~/git/Microbiome_sequence_analysis_toolkit/

# Set the working directory
echo Setting working directory ${WorkingDirectory}...
cd ${WorkingDirectory}

# Make the output directory for the results here
echo Making output directories in ${Output}...
mkdir ./${Output}


PerformPhylogeny () {
	#####################################
	# Extract Phylogenetic Marker Genes #
	#####################################
	# 1 = Taxa Name
	# 2 = Phylo Marker Gene Database
	# 3 = ORF fasta
	# 4 = minimum length of ORF for filtering

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

	# Run blast to detect what ORF sequences are most similar to the L1 capsid gene
	mkdir ./${Output}/BlastnResults
	echo ${1} - Making blast reference using phylogenetic marker genes...

	# Filter out the shorter ORFs that cannot match the entire phylogeny reference
	awk -v sizeMin=${4} '!/^>/ { next } { getline seq } length(seq) >= size { print $0 "\n" seq }' ${3} > ./${Output}/BlastnResults/${1}-FilteredOrfs.fa

	makeblastdb \
		-dbtype nucl \
		-in ${2} \
		-out ./${Output}/BlastnResults/${1}-BlastDb

	echo ${1} - Blasting to determine what ORFs match the phylogenetic marker genes...

	blastn \
		-query ./${Output}/BlastnResults/${1}-FilteredOrfs.fa \
		-out ./${Output}/BlastnResults/${1}-GeneBlastOutput.tsv \
		-db ./${Output}/BlastnResults/${1}-BlastDb \
		-outfmt 6 \
		-num_threads 2 \
		-max_target_seqs 1 \
		-evalue 1e-10

	# Now get a subset of those contig L1 genes so that they can be aligned with the rest of the genes
	mkdir ./${Output}/PhyloContigGeneSubset

	echo ${1} - Generating ORF subset for phylogenetic marker gene...

	cut \
		-f 1 \
		./${Output}/BlastnResults/${1}-GeneBlastOutput.tsv \
		> ./${Output}/PhyloContigGeneSubset/${1}-PhyloContigGeneSubsetList.tsv
	
	grep \
		-A 1 \
		-f ./${Output}/PhyloContigGeneSubset/${1}-PhyloContigGeneSubsetList.tsv \
		./${Output}/BlastnResults/${1}-FilteredOrfs.fa \
		| grep -v '\-\-' \
		> ./${Output}/PhyloContigGeneSubset/${1}-PhyloContigGenes.fa


	###################################################
	# Perform Phylogenetic Alignment and Calculations #
	###################################################
	mkdir ./${Output}/PhylogeneticResults

	# First cat together the known reference sequences with those from contigs
	echo ${1} - Cating files for MAFFT input...

	cat ${2} \
		./${Output}/PhyloContigGeneSubset/${1}-PhyloContigGenes.fa \
		> ./${Output}/PhyloContigGeneSubset/${1}-TotalPhyloSetForMafft.fa
	
	echo ${1} - Aligning phylogenetic marker genes using MAFFT...
	
	/project/egricelab/ghanni_software/bin/mafft \
		--localpair \
		--maxiterate 1000 \
		./${Output}/PhyloContigGeneSubset/${1}-TotalPhyloSetForMafft.fa \
		> ./${Output}/PhylogeneticResults/${1}-mafftAlignment.fa
	
	echo ${1} - Running MAFFT output in RAXML...
	
	/project/egricelab/ghanni_software/standard-RAxML-8.1.21/raxmlHPC-SSE3 \
		-s./${Output}/PhylogeneticResults/${1}-mafftAlignment.fa \
		-n ${1}-Phylogeny \
		-m GTRCAT \
		-p 1234 # GTRCAT is the default suggested by the manual, and is a nucleotide model, which is what I want to use.
}

# Export as a function
export -f PerformPhylogeny

# Now run through the subroutine
PerformPhylogeny \
	"HPV" \
	${HpvReferenceSeqs} \
	${HpvOrfSeqs} \
	1500

PerformPhylogeny \
	"PropPhage" \
	${PropPhageReferenceSeqs} \
	${PropPhageOrfSeqs} \
	150

PerformPhylogeny \
	"StaphPhage" \
	${StaphPhageReferenceSeqs} \
	${StaphPhageOrfSeqs} \
	150



