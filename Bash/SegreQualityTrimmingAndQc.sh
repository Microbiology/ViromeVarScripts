# SegreQualityTrimmingAndQc.sh
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
echo Loading module fastx_toolkit-0.0.14...
module load fastx_toolkit-0.0.14

# Set the path and sample variables to be used throughout the script
export WorkingDirectory=/home/ghanni/Analysis/SegreVirome
export Output='QualityTrimmingAndQc'
export PerlPath=~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/
export Rpath=~/git/Club_Grice/scripts/ghanni/analysis-R/HumanVirome02/
export SampleDirectory=/home/ghanni/Analysis/HumanVirome02/DownloadSegreSequences/FastqFilesFromSra

# For potential debugging purposes, echo back the parameters selected
echo Working directory is ${WorkingDirectory}...
echo Output directory name is ${Output}...
echo The name of the perl path is ${PerlPath}...
echo The name of the R path is ${Rpath}...
echo ***+++***
# =====================================================

# First we need to gunzip the files so that they are easier to work with.
# Do this after making an output directory within the working directory.
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}
echo Making direcotry for error files...
mkdir ./${Output}/StdErr

# Copy over the Segre virome files into this directory
mkdir ./${Output}/AdapterTrimmedFastq
echo Moving files to new working directory from processing directory...
cp ${SampleDirectory}/* ./${Output}/AdapterTrimmedFastq/

#==============================================
#FastaX quality trimming of the raw fastq files
#==============================================
#Make directory for the quality trimmed fastq files
mkdir ./${Output}/TrimmedFastq/
echo Quality trimming fastq sequences with FASTX...
ls ./${Output}/AdapterTrimmedFastq/* \
    | sed -e 's/^.*\///g' \
    | xargs \
        -I {} \
        --max-procs=48 \
        fastq_quality_trimmer \
            -t 33 \
            -Q 33 \
            -i ./${Output}/AdapterTrimmedFastq/{} \
            -o ./${Output}/TrimmedFastq/{}

#==================================
#DeconSeq clean trimmed fastq files
#==================================
mkdir ./${Output}/deconseqFastq
mkdir ./${Output}/StdErr/deconseqFastq
echo Removing human sequences from the dataset using DeconSeq...
RunDeconseqForHuman () {
    #Modules need to be loaded in each submitted function to work properly
    if [ -f /etc/profile.d/modules.sh ]; then
        source /etc/profile.d/modules.sh
    fi
    module load perl5lib
    perl /project/egricelab/meiselj_software/deconseq-standalone-0.4.3/deconseq.pl \
        -f ./${Output}/TrimmedFastq/${1} \
        -dbs hsref \
        -out_dir ./${Output}/deconseqFastq/${1}
}
export -f RunDeconseqForHuman

ls ./${Output}/TrimmedFastq/* \
    | sed -e 's/^.*\///g' \
    | xargs \
        -I {} \
        --max-procs=48 \
        sh -c 'bsub -K -M 64000 -eo ./${Output}/StdErr/deconseqFastq/{} "RunDeconseqForHuman {}"'
wait

# Automatically rename the output files so that they are specific for each sample.
# Here I am only renaming the 'clean' files because those are the ones I will continue to work with downstream.
echo Renaming clean and contaminated files and moving to different directory...
mkdir ./${Output}/DeconSeqCleanFastq
mkdir ./${Output}/DeconSeqContFastq

ls -d ./${Output}/deconseqFastq/* \
    | sed 's/^.*\///g' \
    | xargs \
        -I {} \
        --max-procs=1 \
        sh -c 'cp ./${Output}/deconseqFastq/{}/*_clean.fq ./${Output}/DeconSeqCleanFastq/{}'

ls -d ./${Output}/deconseqFastq/* \
    | sed 's/^.*\///g' \
    | xargs \
        -I {} \
        --max-procs=1 \
        sh -c 'cp ./${Output}/deconseqFastq/{}/*_cont.fq ./${Output}/DeconSeqContFastq/{}'

#Quantify the numbers of reads in the clean and contaminated files, and calculate the percent human contamination.
echo Quantifying the levels of human contamination...
mkdir ./${Output}/HumanContaminationStats
# Be sure to remove the files we will be creating because we are appending and preexisting files can mess it up
rm ./${Output}/HumanContaminationStats/CleanRawSeqCounts.txt
rm ./${Output}/HumanContaminationStats/ContRawSeqCounts.txt

for file in $(ls ./${Output}/DeconSeqCleanFastq); do
    LINES=`cat ./${Output}/DeconSeqCleanFastq/$file | wc -l`
    READS=`expr $LINES / 4`
    echo ${file} $READS >> ./${Output}/HumanContaminationStats/cleanRawSeqCounts.txt
done
for file in $(ls ./${Output}/DeconSeqContFastq); do
    LINES=`cat ./${Output}/DeconSeqContFastq/$file | wc -l`
    READS=`expr $LINES / 4`
    echo ${file} $READS >> ./${Output}/HumanContaminationStats/ContRawSeqCounts.txt
done

# These files can be used to plot the percent human contamination in R
#Get information on decontamination efforts, and generate graphs using R
echo Using human_decontamination_plot.R to generate decontamination stats graphs...
Rscript ${Rpath}human_decontamination_plot.R \
    ./${Output}/HumanContaminationStats/cleanRawSeqCounts.txt \
    ./${Output}/HumanContaminationStats/ContRawSeqCounts.txt \
    ./${Output}/HumanContaminationStats

# Convert the output fastq files to fasta for downstream processes
mkdir ./${Output}/clean_R2_fasta
for file in $(ls ${Output}/DeconSeqCleanFastq/); do
    /project/egricelab/meiselj_software/idba_ud-1.0.9/bin/fq2fa ./${Output}/DeconSeqCleanFastq/$file ./${Output}/clean_R2_fasta/${file}
done
for file in $(ls ./${Output}/clean_R2_fasta/); do
        mv ./${Output}/clean_R2_fasta/"${file}" ./${Output}/clean_R2_fasta/"${file/%.fastq/.fa}"
done

# This dataset does not have negative control data, so I will be skipping that part of the analysis
# Move the R2 (forward) files to the next directory so that I can easily move on with the rest of my scripts
mkdir ./${Output}/NegCleanedFasta
cp ./${Output}/clean_R2_fasta/*_2.fa ./${Output}/NegCleanedFasta/








