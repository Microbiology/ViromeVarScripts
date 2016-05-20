#! /bin/bash
# EvaluateDgrs.sh
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
export Output='DgrDetection'
export PerlPath=/home/ghanni/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/
# Downloaded from UniProt using their REST API
# wget "http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=%22reverse%20transcriptase%22%20(phage%20OR%20virus)&fil=&format=fasta&force=yes" -O ./virusPhageRt.fa.gz
export rtDatabase=/project/egricelab/references/UniprotRtForViromeTargetedVariability/uniprotPhageVirusRt
export HypervariableLoci=/home/ghanni/Analysis/CleanRerunHumanVirome02/variabilityCallGeomDist/VariableSnpRegionFiles/HypervariableLociResults.tsv
export ContigAlignment=/home/ghanni/Analysis/CleanRerunHumanVirome02/AlignedReadsOfInterestFromSRA/BowtieOutBam/OverallContigHits_Sorted-InterestOnlySRA.bam

# HPV
export HpvContigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/HpvContigs.fa
export HpvOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/HpvPredOrfs.fa
export HpvPfamDomains=/home/ghanni/Analysis/CleanRerunHumanVirome02/PfamHmmerDomainAnnotation/PfamDomains/HPV-PfamDomainsFormat.tsv
export HpvTransPhageOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/PfamHmmerDomainAnnotation/PfamDomains/HPV-TanslatedOrfs.fa

# HPV
export StaphPhageContigs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhageContigs.fa
export StaphPhageOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/TaxContigs/SpecificContigsByTaxa/StaphPhagePredOrfs.fa
export StaphPfamDomains=/home/ghanni/Analysis/CleanRerunHumanVirome02/PfamHmmerDomainAnnotation/PfamDomains/StaphPhage-PfamDomainsFormat.tsv
export StaphTransPhageOrfs=/home/ghanni/Analysis/CleanRerunHumanVirome02/PfamHmmerDomainAnnotation/PfamDomains/StaphPhage-TanslatedOrfs.fa


# Set the working directory
echo Setting working directory ${WorkingDirectory}...
cd ${WorkingDirectory} || exit

# Make the output directory for the results here
echo Making output directories in ${Output}...
mkdir ./${Output}

GetCandidateDgrs () {
	# 1 = Sample Name
	# 2 = Contig Fasta
	# 3 = ORF Fasta
	# 4 = RT Database

	# Detect hits to RT elements
	blastx \
		-query ${3} \
		-out ./${Output}/${1}-BlastRtHits.tsv \
		-db ${rtDatabase} \
		-outfmt 6 \
		-max_target_seqs 1 \
		-evalue 1e-5

	# Detect repeats within virus contigs
	makeblastdb \
		-dbtype nucl \
		-in ${2} \
		-out ./${Output}/${1}-db

	tblastx \
		-query ${2} \
		-out ./${Output}/${1}-tblastxRepeatHits.tsv \
		-db ./${Output}/${1}-db \
		-outfmt 6 \
		-num_threads 6 \
		-evalue 1e-50

	# Format the output and filter out regions <200 nt
	awk '$1 = $2 {print $0}' ./${Output}/${1}-tblastxRepeatHits.tsv \
		| awk '$7 != $9 {print $0}' \
		| awk 'sqrt(($8 - $7)^2) < 150 {print $0}' \
		| sed 's/ /\t/g' \
		> ./${Output}/${1}-UniqueRepeats.tsv

	# Remove duplicate candidate pairs and format
	perl ${PerlPath}ExtractDgrs.pl \
		./${Output}/${1}-BlastRtHits.tsv \
		./${Output}/${1}-UniqueRepeats.tsv \
		./${Output}/${1}-DgrCandidates.tsv

	# Get the DGRs that contain at least one HVL
	perl ${PerlPath}GetDgrHvls.pl \
		./${Output}/${1}-DgrCandidates.tsv \
		${HypervariableLoci} \
		./${Output}/${1}-DgrHvl.tsv

	# Generate gff3 file for visualizing genes with annotation
	# and annotate genes predicted to be affected by DGRs
	mkdir ./${Output}/ForVisualization

	perl ${PerlPath}GetDgrContainedGenes.pl \
		./${Output}/${1}-DgrHvl.tsv \
		${5} \
		${4} \
		./${Output}/${1}-BlastRtHits.tsv \
		${3} \
		./${Output}/${1}-DgrOutput.tsv \
		./${Output}/ForVisualization/${1}-Dgr.gff3

	# To keep things easy, lets bring in the alignment and contig
	# files into the same directory
	ContigList=$(cut -f 1 ./${Output}/ForVisualization/${1}-Dgr.gff3 | sed 's/$/_/' | sort | uniq)

	echo List is ${ContigList}

	# Get sequence alignment to those contigs by subsetting
	# existing alignment
	# Make sure it is indexed like this:
	# samtools index -b OverallContigHits_Sorted-InterestOnlySRA.bam
	for line in ${ContigList}; do
		egrep -A 1 ${line} \
			${2} \
			> ./${Output}/ForVisualization/${1}-DgrContigs.fa

		samtools \
			faidx \
			./${Output}/ForVisualization/${1}-DgrContigs.fa

		samtools view \
			-bh \
			${ContigAlignment} \
			${line} \
			> ./${Output}/ForVisualization/${1}-${line}-Contigs.bam

		# Index the bam file
		cd ./${Output}/ForVisualization
		samtools index \
			-b ./${1}-${line}-Contigs.bam
		cd ../../
	done
	# At this point I should have enough to visualize the DGRs as
	# a figure.
}

# Use pileup to create matrix to calculate LD in R
PrepareLdMatrix () {
	# 1 = Name
	# 2 = DGR GFF3 File

	# Iterate through the contig IDs by first creating a list
	cut -f 1 ${2} | sort | uniq > ./${Output}/${1}-ContigIdList.tsv

	while read line; do
		# The first nucleotide position of the first DGR element
		FIRSTone=$(awk -v name=${line} '$1 == name && $3 == "DGR-1" {print $4}' ./${Output}/ForVisualization/${1}-Dgr.gff3)
		echo ${FIRSTone}...
		# The second nucleotide position of the first DGR element
		SECONDone=$(awk -v name=${line} '$1 == name && $3 == "DGR-1" {print $5}' ./${Output}/ForVisualization/${1}-Dgr.gff3)
		echo ${SECONDone}...
		# The first nucleotide position of the second DGR element
		FIRSTtwo=$(awk -v name=${line} '$1 == name && $3 == "DGR-2" {print $4}' ./${Output}/ForVisualization/${1}-Dgr.gff3)
		echo ${FIRSTtwo}...
		# The second nucleotide position of the second DGR element
		SECONDtwo=$(awk -v name=${line} '$1 == name && $3 == "DGR-2" {print $5}' ./${Output}/ForVisualization/${1}-Dgr.gff3)
		echo ${SECONDtwo}...

		# Pull out the nucleotide matrix given the nucleotide coordinates
		if [ "${FIRSTone}" -lt "${SECONDone}" ]; then
			export COLUMNS=$(expr $FIRSTone - $SECONDone)
			samtools tview \
				-p ${line}_:${FIRSTone} \
				-d t \
				./${Output}/ForVisualization/${1}-${line}_-Contigs.bam \
				./${Output}/ForVisualization/${1}-DgrContigs.fa \
				> ./${Output}/ForVisualization/${1}-${line}-DGR1.tview
		else
			export COLUMNS=$(expr $SECONDone - $FIRSTone)
			samtools tview \
				-p ${line}_:${SECONDone} \
				-d t \
				./${Output}/ForVisualization/${1}-${line}_-Contigs.bam \
				./${Output}/ForVisualization/${1}-DgrContigs.fa \
				> ./${Output}/ForVisualization/${1}-${line}-DGR1.tview
		fi

		if [ "${FIRSTtwo}" -lt "${SECONDtwo}" ]; then
			export COLUMNS=$(expr $SECONDtwo - $FIRSTtwo)
			samtools tview \
				-p ${line}_:${FIRSTtwo} \
				-d t \
				./${Output}/ForVisualization/${1}-${line}_-Contigs.bam \
				./${Output}/ForVisualization/${1}-DgrContigs.fa \
				> ./${Output}/ForVisualization/${1}-${line}-DGR2.tview
		else
			export COLUMNS=$(expr $SECONDtwo - $FIRSTtwo)
			samtools tview \
				-p ${line}_:${SECONDtwo} \
				-d t \
				./${Output}/ForVisualization/${1}-${line}_-Contigs.bam \
				./${Output}/ForVisualization/${1}-DgrContigs.fa \
				> ./${Output}/ForVisualization/${1}-${line}-DGR2.tview
		fi

		# Convert to genotype tables
		perl ${PerlPath}mpilup2genotype.pl \
			./${Output}/ForVisualization/${1}-${line}-DGR1.tview \
			./${Output}/ForVisualization/${1}-${line}-DGR1-GenotypeTable.tsv

		perl ${PerlPath}mpilup2genotype.pl \
			./${Output}/ForVisualization/${1}-${line}-DGR2.tview \
			./${Output}/ForVisualization/${1}-${line}-DGR2-GenotypeTable.tsv

	done < ./${Output}/${1}-ContigIdList.tsv
}

export -f GetCandidateDgrs
export -f PrepareLdMatrix

GetCandidateDgrs \
	"HPV" \
	${HpvContigs} \
	${HpvOrfs} \
	${HpvPfamDomains} \
	${HpvTransPhageOrfs}

GetCandidateDgrs \
	"Staph" \
	${StaphPhageContigs} \
	${StaphPhageOrfs} \
	${StaphPfamDomains} \
	${StaphTransPhageOrfs}

PrepareLdMatrix \
	"HPV" \
	./${Output}/ForVisualization/HPV-Dgr.gff3

PrepareLdMatrix \
	"Staph" \
	./${Output}/ForVisualization/Staph-Dgr.gff3

