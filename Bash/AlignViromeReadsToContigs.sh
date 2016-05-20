# AlignViromeReadsToContigs.sh
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
export Output='AlignedReadsOfInterestFromSRA'
export Map=/home/ghanni/git/Club_Grice/mapping_files/metadata/SkinMet_and_Virome_001_metadata_for_SRA.txt
export SampleDirectory=/home/ghanni/Analysis/Human_virome_analysis/negative_clean_seqs/
export Contigs=/home/ghanni/Analysis/Human_virome_analysis/ray_contigs_from_total_cat_pairs/Contigs_no_block_with_names.fasta
export BinPath=~/git/Club_Grice/bin/
export ToolkitPath=~/git/Microbiome_sequence_analysis_toolkit/

# Set the working directory
echo Setting working directory ${WorkingDirectory}...
cd ${WorkingDirectory}

# Make the output directory for the results here
echo Making output directories in ${Output}...
mkdir ./${Output}
mkdir ./${Output}/tmp
mkdir ./${Output}/tmp/CategoryFasta



##############################################
# Concatenate Files to be used in this study #
##############################################
# Because this relies on appending to cat file, remove if it already exists
rm ./${Output}/tmp/CategoryFasta/InterestOnlySRA.fa
rm ./${Output}/tmp/CategoryFasta/InterestOnlySRA.log
# Get a file with the subject IDs		
cut -f 4 ${Map} \
	| sort \
	| uniq \
	| grep -v 1 \
	> ./${Output}/tmp/TimePoints.tsv		

cut -f 5 ${Map} \
	| sort \
	| uniq \
	| grep -v Ba \
	| grep -v Ph \
	| grep -v Vf \
	| grep -v Neg \
	| grep -v No \
	> ./${Output}/tmp/SiteCategories.tsv		

while read timepoint; do
  while read sitecat; do
    for SampleID in $(awk -v Time=${timepoint} -v Site=${sitecat} '$4 == Time && $5 == Site && $11 == "Virome" { print $1 }' ${Map}); do		
      echo Looking at ${timepoint}-${sitecat}-${SampleID}...
      # Create a log file of the samples in the cat
      printf "Added Sample ${timepoint}-${sitecat}-${SampleID}\n" >> ./${Output}/tmp/CategoryFasta/InterestOnlySRA.log
      cat ${SampleDirectory}${SampleID}* >> ./${Output}/tmp/CategoryFasta/InterestOnlySRA.fa		
    done		
  done < ./${Output}/tmp/SiteCategories.tsv		
done < ./${Output}/tmp/TimePoints.tsv		



#########################################################
# Prepare gff3 from phylogenetically assignable contigs #
#########################################################
perl \
	~/git/Club_Grice/bin/remove_block_fasta_format.pl \
	${Contigs} \
	./${Output}/HpvRefContigsNoBlocks.fa

mkdir ./${Output}/PredictedOrfs

echo 'Use glimmer wrapper to detect ORFs...'
echo 'Building icm...'

# Use the cat fasta file of the HPV reference genomes with the contigs 
# containing phylogenetic information
/project/egricelab/meiselj_software/glimmer3.02/bin/build-icm \
	./${Output}/PredictedOrfs/PredictedOrfs.icm \
	< ./${Output}/HpvRefContigsNoBlocks.fa

echo 'Running glimmer...'

/project/egricelab/meiselj_software/glimmer3.02/bin/glimmer3 \
	--linear \
	-g 50 \
	./${Output}/HpvRefContigsNoBlocks.fa \
	./${Output}/PredictedOrfs/PredictedOrfs.icm \
	./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutput

# Format the output file
cat ./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutput.predict | while read line;
    do if [ "${line:0:1}" == ">" ]
        then seqname=${line#'>'}
        else
        orf="$seqname.${line%%' '*}"
        coords="${line#*' '}"
        echo -e "$orf\t$seqname\t$coords"
        fi
    done > ./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutput.predict.formatted

# Extract the gene sequences

echo 'Performing multi-extract...'

/project/egricelab/meiselj_software/glimmer3.02/bin/multi-extract \
	-l 50 \
	--nostop \
	./${Output}/HpvRefContigsNoBlocks.fa \
	./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutput.predict.formatted \
	> ./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutput.genes

# Remove the block formatting
perl \
	${BinPath}remove_block_fasta_format.pl \
	./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutput.genes \
	./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutputGenes.fa

# Convert the predict file to gff3 format for downstream analysis
# This needs to be run through a loop
mkdir ./${Output}/PredictedOrfs/Gff3ConvertDocs
for contigID \
	in $(grep '>' ./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutput.predict | sed 's/>//' | sed 's/_//g')
	do
  		echo Loop contig ID is ${contigID}...
  		perl \
  			${ToolkitPath}GlimmerPredict2Gff3/GlimmerPredict2Gff3.pl \
  				./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutput.predict \
  				${contigID} \
  				./${Output}/PredictedOrfs/Gff3ConvertDocs/PredictedOrfsGlimmerOutput${contigID}.gff3
done

# Cat together the gff3 files to one big file, which will be used downstream in pnps ratio calculations
cat \
	./${Output}/PredictedOrfs/Gff3ConvertDocs/PredictedOrfsGlimmerOutput*.gff3 \
	> ./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutput.gff3

sed 's/|//g' \
	./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutput.gff3 \
	| sed 's/\.//g' \
	| sed 's/ .*//' \
	> ./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutputFormat.gff3



#################################
# Align Fasta to Virome Contigs #
#################################
mkdir ./${Output}/tmp/BowtieBuilds/
mkdir ./${Output}/BowtieOut/
mkdir ./${Output}/BowtieOutBam/
# Build bowtie reference database		
bowtie2-build -f ./${Output}/HpvRefContigsNoBlocks.fa \
	./${Output}/tmp/BowtieBuilds/BowtieReferenceGenomes
	
# Do the alignment	
bowtie2 -x ./${Output}/tmp/BowtieBuilds/BowtieReferenceGenomes \
	-f ./${Output}/tmp/CategoryFasta/InterestOnlySRA.fa \
	-S ./${Output}/BowtieOut/InterestOnlySRA.sam \
	-L 25 \
	-N 1 \
	-p 16

# Remove unneeded sequences and convert to bam
samtools view -bS \
	./${Output}/BowtieOut/InterestOnlySRA.sam \
	> ./${Output}/BowtieOut/OverallContigHits-InterestOnlySRA.bam		

# Sort the Bam file
samtools sort \
	./${Output}/BowtieOut/OverallContigHits-InterestOnlySRA.bam \
	./${Output}/BowtieOutBam/OverallContigHits_Sorted-InterestOnlySRA		

# Index the reference genome		
samtools faidx \
	./${Output}/HpvRefContigsNoBlocks.fa	

# Use VarScan to calculate the variable nucleotides 		
samtools mpileup \
	-f ./${Output}/HpvRefContigsNoBlocks.fa \
	./${Output}/BowtieOutBam/OverallContigHits_Sorted-InterestOnlySRA.bam \
	> ./${Output}/BowtieOut/SortedContigs-InterestOnlySRA.mpileup

# Make output directory for VCF
mkdir ./${Output}/VcfResults/

java -jar \
	/project/egricelab/ghanni_software/VarScan/VarScan.v2.3.7.jar \
	pileup2snp \
	./${Output}/BowtieOut/SortedContigs-InterestOnlySRA.mpileup \
	> ./${Output}/VcfResults/VarScanOverllHits-InterestOnlySRA.vcf

# Format each file reference name		
sed -i 's/|//g' \
	./${Output}/VcfResults/VarScanOverllHits-InterestOnlySRA.vcf \
	| sed 's/\.//g' \
	| sed 's/ .*//g'

# Remove indels to prevent issues downstream		
perl \
	~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/RemoveIndelsFromVCF.pl \
	./${Output}/VcfResults/VarScanOverllHits-InterestOnlySRA.vcf \
	./${Output}/VcfResults/VarScanOverllHits-InterestOnlySRA-NoIndels.vcf

# Make directory for pNpS calculation output
mkdir ./${Output}/PnPsResults

# Run through pNpS calculator		
perl \
	~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/CalculatePnPsRatio.pl \
	./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutputFormat.gff3 \
	./${Output}/HpvRefContigsNoBlocks.fa \
	./${Output}/VcfResults/VarScanOverllHits-InterestOnlySRA-NoIndels.vcf \
	./${Output}/PnPsResults/OverallHpvContigs-InterestOnlySRA.tsv \
	./${Output}/PnPsResults/SnpPatterns-InterestOnlySRA.tsv

grep \
	-v "ContigID" \
	./${Output}/PnPsResults/OverallHpvContigs-InterestOnlySRA.tsv \
	> ./${Output}/PnPsResults/OverallHpvContigsFormat-InterestOnlySRA.tsv		

#Merge the gff3 information with the pnps ratios		
awk 'NR==FNR {a[$1$6] = $0; next} {print $0"\t"a[$1$2]}' \
	./${Output}/PredictedOrfs/PredictedOrfsGlimmerOutputFormat.gff3 \
	./${Output}/PnPsResults/OverallHpvContigsFormat-InterestOnlySRA.tsv \
	| awk -v Name=InterestOnlySRA ' { print $0"\t"Name } ' \
	| sed 's/-//1' \
	| sed 's/-/\t/g' \
	| sed 's/\t\t/\tNA\t/' \
	| sed 's/\t\t/\t/g' \
	> ./${Output}/PnPsResults/OverallHpvContigsWithGff3-InterestOnlySRA.tsv		

