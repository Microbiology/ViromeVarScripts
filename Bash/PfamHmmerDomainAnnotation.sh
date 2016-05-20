# PfamHmmerDomainAnnotation.sh
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
export Output='PfamHmmerDomainAnnotation'
export PerlPath=~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/
# Don't forget to press this reference before use doing:
# /project/egricelab/bin/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmpress /project/egricelab/references/Pfam-Hmm/Pfam-A.hmm
export PfamHmmRef=/project/egricelab/references/Pfam-Hmm/Pfam-A.hmm
export PfamSeqRef=/project/egricelab/references/Pfam-fasta/Pfam-A.fasta

# HPV
export HpvOrfSeqs=/home/ghanni/Analysis/CleanRerunHumanVirome02/DeleteriousScoreReference/HPVHotspotOrfsOnly.fa
# Staph Phage
export StaphPhageOrfSeqs=/home/ghanni/Analysis/CleanRerunHumanVirome02/DeleteriousScoreReference/StaphPhageHotspotOrfsOnly.fa
# Make the output directory and move to the working directory
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}


#########################
# Annotate Pfam Domains #
#########################
PfamDomains () {
	if [ -f /etc/profile.d/modules.sh ]; then
		source /etc/profile.d/modules.sh
	fi
	module load bowtie2-2.1.0
	module load ncbi-blast-2.2.0
	module load perl5lib
	module load samtools-1.1
	module load R-3.1.2
	module load EMBOSS-6.6.0

	# 1 = Taxa Name
	# 2 = ORF Fasta (nucleotide)
	# 3 = Reference Database (pfam)
	# 4 = Reference Sequence Database (pfam)

	# Make output directory
	mkdir ./${Output}/PfamDomains

	# Translate the sequences
	transeq -sequence ${2} -outseq ./${Output}/PfamDomains/${1}-TanslatedOrfs.fa -clean

	# Perform HMM alignment against pfam HMMER database
	/project/egricelab/bin/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan \
		--cpu 8 \
		--notextw \
		--cut_ga \
		--domtblout ./${Output}/PfamDomains/${1}-PfamDomains.hmmscan \
		${3} \
		./${Output}/PfamDomains/${1}-TanslatedOrfs.fa

	# Format the data so it is easier to deal with in R analysis
	# The cut by character count works because it is space delimited
	# With the final column starting at character 181.
	grep -v '#' ./${Output}/PfamDomains/${1}-PfamDomains.hmmscan  \
		| cut -c 1-180 \
		| sed 's/\s\+/\t/g' \
		| sort -rnk22 \
		> ./${Output}/PfamDomains/${1}-PfamDomainsFormat.tsv
}

# Export as a function
export -f PfamDomains

PfamDomains \
	"HPV" \
	${HpvOrfSeqs} \
	${PfamHmmRef} \
	${PfamSeqRef}

PfamDomains \
	"StaphPhage" \
	${StaphPhageOrfSeqs} \
	${PfamHmmRef} \
	${PfamSeqRef}


