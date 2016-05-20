# GenerateDelScoreTable.sh
# Geoffrey Hannigan
# Laboratory of Elizabeth Grice
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
# The reference tables with Del scores were copied and pasted from the
# suspect website. I manually copied and pasted the sequences from the
# associated fasta to the online program interface, and copied the output
# to the tsv files referenced below.
# http://www.sbg.bio.ic.ac.uk/suspect/
export WorkingDirectory=/home/ghanni/Analysis/CleanRerunHumanVirome02
export Output='DeleteriousScores'
export PerlPath=~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/
export DelReference=/home/ghanni/Analysis/CleanRerunHumanVirome02/DeleteriousScoreReference/DelScoreReferenceTable.tsv
export StaphPhageDelReference=/home/ghanni/Analysis/CleanRerunHumanVirome02/DeleteriousScoreReference/StaphPhageDelScoreReferenceTable.tsv

# Make the output directory and move to the working directory
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}

# Temporarily move directories to make my life easier
cd ./${Output}

GenerateScoreTable () {
	# 1 = File Set Name
	# 2 = Del Score Reference File

	while read line; do
		export Contig=$(echo ${line} | awk '{print $1}')
		export ORF=$(echo ${line} | awk '{print $2}')
		export html=$(echo ${line} | awk '{print $3}')
		echo Contig is ${Contig}
		echo ORF is ${ORF}
		echo HTML is ${html}
	
		echo Loading in file...
		wget --output-document=${Contig}-${ORF}.tsv ${html}
		sed -i "s/^/${Contig}\t${ORF}\t/" ${Contig}-${ORF}.tsv
		sed -i -e '$a\' ${Contig}-${ORF}.tsv
	done < ${2}

	# Cat together the individual files
	cat ./*orf*.tsv > ./${1}DelScoreResults.tsv
	
	# Remove all of the intermediate files
	rm ./*orf*.tsv
}

export -f GenerateScoreTable

GenerateScoreTable \
	"HPV" \
	${DelReference}

GenerateScoreTable \
	"StaphPhage" \
	${StaphPhageDelReference}

