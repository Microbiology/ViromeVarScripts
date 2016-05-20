# SegreAssembleContigsRay.sh
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
echo Loading module openmpi-1.5.4-x86_64...
module load openmpi-1.5.4-x86_64
echo Loading module python-2.7.5...
module load python-2.7.5

# Set the path and sample variables to be used throughout the script
export WorkingDirectory=/home/ghanni/Analysis/SegreVirome
export Output='RayAssembly'
export PerlPath=~/git/Club_Grice/scripts/ghanni/analysis-perl/HumanVirome02/
export Rpath=~/git/Club_Grice/scripts/ghanni/analysis-R/HumanVirome02/
export SampleDirectory=/home/ghanni/Analysis/SegreVirome/QualityTrimmingAndQc/DeconSeqCleanFastq/

# For potential debugging purposes, echo back the parameters selected
echo Working directory is ${WorkingDirectory}...
echo Output directory name is ${Output}...
echo The name of the perl path is ${PerlPath}...
echo The name of the R path is ${Rpath}...
echo ***+++***
# =====================================================

# Do this after making an output directory within the working directory.
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}

# Assemble contigs by sample library
echo Assembling individual libraries...
### Place libraries in different directories bases on _1 vs _2 (forward vs reverse)
cat ${SampleDirectory}*_1* > ./${Output}/_1_for_ray.fq
cat ${SampleDirectory}*_2* > ./${Output}/_2_for_ray.fq
### Get sequence pairs
echo Getting sequence pairs...
python /project/egricelab/bin/get_trimmed_pairs.py \
	-f ./${Output}/_1_for_ray.fq -s ./${Output}/_2_for_ray.fq \
	-o ./${Output}/_1_for_ray_pairs.fq \
	-t ./${Output}/_2_for_ray_pairs.fq
### Convert the sequence files to fasta from fastq
echo Converting sequences to fasta from fastq...
/project/egricelab/meiselj_software/idba_ud-1.0.9/bin/fq2fa \
	./${Output}/_1_for_ray_pairs.fq \
	./${Output}/_1_for_ray_pairs.fa

/project/egricelab/meiselj_software/idba_ud-1.0.9/bin/fq2fa \
	./${Output}/_2_for_ray_pairs.fq \
	./${Output}/_2_for_ray_pairs.fa
### Perform assembly for each set of sequences
echo Assembling sequences per individual sample...

mpiexec \
	-n 9 \
	/project/egricelab/ghanni_software/Ray-2.3.1/Ray \
	-minimum-contig-length 500 \
	-p ./${Output}/_1_for_ray_pairs.fa \
	./${Output}/_2_for_ray_pairs.fa \
	-o ./${Output}/ray_contigs

#Add name of sample to the end of each contig ID to ensure they are all unique IDs when they are cat together
#Remove block format in contig fasta file | Next three part of same thing | Replace the spaces in the contig names with underscores | Add sample ID to the end of each name
sed -r 's/\s/_/g' ./${Output}/ray_contigs/Contigs.fasta \
	| sed 's/^\([A,T,G,C,n]*\)$/\1\@/g' \
	| sed ':a;N;$!ba;s/\@\n\([A,C,G,T,n]\)/\1/g' \
	| sed 's/\@//g' \
	| sed '/\>/s/ /_/g' \
	| sed "/>/s/$/\_$file/" \
	> ./${Output}/ContigsFormat.fasta
# Also number the contigs with easy unique identifiers
perl /home/ghanni/git/Club_Grice/scripts/ghanni/ViromeAnalysis/scripts/RenameContigs.pl \
	./${Output}/ContigsFormat.fasta \
	./${Output}/ContigsFormatNumbered.fasta

# Add underscore at the end for easy termination recognition of the IDs later on in various regex's
sed -i 's/\(^>.*$\)/\1_/' ./${Output}/ContigsFormatNumbered.fasta

