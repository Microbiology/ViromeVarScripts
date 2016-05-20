# PrepareReferencesForPhylogeny.sh
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania


# The HPV L1 reference genes are easily downloaded from the
# PAVE website. All HPV L1 genes were downloaded on 2015-06-03.

# The Staphylococcus phage and Propionibacterium phage large terminase
# subunit genes were taken from the NCBI gene website. The search terms
# used were as follows:
# Staph:
# ((phage terminase large subunit staphylococcus)) AND 
# "viruses"[porgn:__txid10239] NOT "ORF" NOT "hypothetical" NOT "putative" 
# Prop:
# ((phage terminase large subunit propionibacterium)) AND 
# "viruses"[porgn:__txid10239] 
# The gene sequences (fasta nucleotides) were manually copied and pasted
# into the reference file to be used.

# Further formatting was performed on the reference databases after
# they were downloaded.
# This is done in the directory where the databases where downloaded
cd /home/ghanni/Analysis/HumanVirome02/HpvPhylogeny/referenceSeqs

sed 's/_L1.*//g' PaveHpvL1ReferenceOf2015-06-03.fa \
	| sed 's/.*HPV/>HPV/g' \
	> PaveHpvL1ReferenceOf2015-06-03-Format.fa

sed 's/\,//g' StaphPhageRefGenomesPhylogeny-LargeTerminaseSubunit-2015-09-14.fa 
	| sed 's/\://g' 
	| sed 's/ //g' 
	| sed 's/\;//g' 
	| perl -pe 's/.*?-//' 
	| sed 's/complete.*//' 
	| sed 's/phage/-phage-/' 
	| sed 's/\([0-9]\)/>\1/' 
	| sed 's/Staph/-Staph/' 
	> StaphPhageRefGenomesPhylogeny-LargeTerminaseSubunit-2015-09-14-Format.fa

sed 's/\,//g' PropPhageRefGenomesPhylogeny-LargeTerminaseSubunit-2015-09-15.fa \
	| sed 's/\://g' \
	| sed 's/ //g' \
	| sed 's/\;//g' \
	> PropPhageRefGenomesPhylogeny-LargeTerminaseSubunit-2015-09-15-Format.fa

# At this point the files can be used in the phylogeny analysis scripts.
