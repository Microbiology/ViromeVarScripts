#!/usr/bin/perl
# CalculatePnPsRatio.pl
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

# Set use
use strict;
use warnings;
use Data::Dumper;
# Also measure time to run script
my $start_run = time();

# Set files to scalar variables
my $usage = "Usage: perl $0 <GFF3> <FASTA> <VCF> <OUT> <SEQOUT> <VCFOUT>";
# The gff3 file should include all genes of interest
# The fasta file should include all contigs or genomes
#	described in the gff3 file.
# VCF file should include SNPs in contig of interest
# OUT should be the desired output file 
my $gff3 = shift or die $usage;
my $fasta = shift or die $usage;
my $vcf = shift or die $usage;
my $outfile = shift or die $usage;
my $seqoutfile = shift or die $usage;
my $vcfout = shift or die $usage;
open(GFF3, "<$gff3") || die "Unable to open $gff3: $!";
open(FASTA, "<$fasta") || die "Unable to open $fasta: $!";
open(VCF, "<$vcf") || die "Unable to open $vcf: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
open(SEQOUT, ">$seqoutfile") || die "Unable to write to $seqoutfile: $!";
# This file will contain the VCF information from the hotspot SNP loci
open(VCFOUT, ">$vcfout") || die "Unable to write to $vcfout: $!";

# Set variables for the script
my %FastaHash;
my %AminoAcidHash;
my %VcfHash;
my @gffArray;
my @VcfArray;
my @gff3Array;
my @vcfArray;
my @CoverageArray;
my $FastaID = 0;
my $FastaErrorCounter = 0;
my $ContigID = 0;
my $geneID = 0;
my $fastaSequence = 0;
my $GeneStart = 0;
my $GeneEnd = 0;
my $PositionPositive = 2;
my $GeneLength = 0;
my $SnpPosition = 0;
my $SnpConsensus = 0;
my $SnpVariable = 0;
my $vcfLine = 0;
my $GeneSequence = 0;
my $SnpCodonPosition = 0;
my $SnpCodonConsensus = 0;
my $SnpCodonVariable = 0;
my $CodonAAConsensus = 0;
my $SnpCodon = 0;
my $CodonAAVariable = 0;
my $SynMutation = 0;
my $NonSynMutation = 0;
my $PnPsRatio = 0;
my $VcfContigID = 0;
my $vcfline = 0;
my $VcfPosition = 0;
my $nucleotideConsensus = 0;
my $FirstAlter = 0;
my $SecondAlter = 0;
my $ThirdAlter = 0;
my $NucleotideCodonPosition = 0;
my $NucleotideCodonConsensus = 0;
my $NucleotideConsensus = 0;
my $NucleotideCodonFirst = 0;
my $NucleotideCodonSecond = 0;
my $NucleotideCodonThird = 0;
my $CodonAAFirst = 0;
my $CodonAASecond = 0;
my $CodonAAThird = 0;
my $SynNucleotideMut = 0;
my $NonSynNucleotideMut = 0;
my $SynMutationValue = 0;
my $NonSynMutationValue = 0;
my $SynMutationCounter = 0;
my $NonSynMutationCounter = 0;
my $NormFactor = 0;
my $SnpCoverage = 0;
my $GeneSequenceInt = 0;
my $SnpPositionCorrNuc = 0;
my $skipFlag = 0;

# Set subroutine to create hash sequences.
# High praise to the monks -> http://www.perlmonks.org/?node_id=904666
sub SetCodonHashInMemory {
	print "PROGRESS: Loading in amino acid codon hash.\n";
	my %prot = (
		'TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S',
		'TTC'=>'F','TTT'=>'F',
		'TTA'=>'L','TTG'=>'L',
		'TAC'=>'Y','TAT'=>'Y',
		'TAA'=>'_','TAG'=>'_','TGA'=>'_',
		'TGC'=>'C','TGT'=>'C',
		'TGG'=>'W',
		'CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L',
		'CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P',
		'CAC'=>'H','CAT'=>'H',
		'CAA'=>'Q','CAG'=>'Q',
		'CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R',
		'ATA'=>'I','ATC'=>'I','ATT'=>'I',
		'ATG'=>'M',
		'ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T',
		'AAC'=>'N','AAT'=>'N',
		'AAA'=>'K','AAG'=>'K',
		'AGC'=>'S','AGT'=>'S',
		'AGA'=>'R','AGG'=>'R',
		'GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V',
		'GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A',
		'GAC'=>'D','GAT'=>'D',
		'GAA'=>'E','GAG'=>'E',
		'GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G'
	);
	return %prot;
}

# Read the fasta file in with the ID as the key, and the sequence
# as the value.

# Get ready to calculate the median
sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub ReadInFasta {
	print "PROGRESS: Loading in fasta file.\n";
	# Set the variable for the fasta input file
	my $fastaInput = shift;
	# Setup fasta hash to return at the end
	while (my $line = <$fastaInput>) {
		# print $line;
		# print "$FastaErrorCounter\n";
		if ($line =~ /\>/ && $FastaErrorCounter == 0) {
			# print "Made it to ID!\n";
			chomp $line;
			$FastaID = $line;
			# Get rid of the arrow from the ID
			$FastaID =~ s/\>//;
			$FastaID =~ s/_//g;
			$FastaErrorCounter = 1;
		} elsif ($line =~ /\>/ && $FastaErrorCounter == 1) {
			die "KILLED BY BAD FASTA! There was no sequence before ID $line: $!";
		} elsif ($line !~ /\>/ && $FastaErrorCounter == 0) {
			print STDERR "Yikes, is this in block format? That is totally not allowed!\n";
			die "KILLED BY BAD FASTA! There was a missing ID: $!";
		} elsif ($line !~ /\>/ && $FastaErrorCounter == 1) {
			chomp $line;
			# print "Made it to sequence!\n";
			# Change out the lower case letters so they match the codon hash
			$line =~ s/g/G/g;
			$line =~ s/a/A/g;
			$line =~ s/c/C/g;
			$line =~ s/t/T/g;
			$FastaHash{$FastaID} = $line;
			$FastaErrorCounter = 0;
		}
	}
	return %FastaHash;
}

sub ReadInGff3 {
	print "PROGRESS: Loading in gene annotation file.\n";
	# Save the gff3 file into memory for faster access
	my $gff3Input = shift;
	while (my $line = <$gff3Input>) {
		chomp $line;
		push @gffArray, $line;
	}
	return @gffArray;
}

sub ReadInVCF {
	print "PROGRESS: Loading in VCF.\n";
	# Save the VCF file into memory for faster access
	my $VcfInput = shift;
	while (my $line = <$VcfInput>) {
		chomp $line;
		$line =~ s/_//g;
		$VcfContigID = (split /\t/, $line)[0];
		$VcfPosition = (split /\t/, $line)[1];
		# print "$VcfContigID\n";
		next if ($VcfContigID eq "Chrom");
		# Save this in as a nested hash within a hash
		$VcfHash{$VcfContigID}{$VcfPosition} = $line;
	}
	return %VcfHash;
}

sub CalculatePnPsRatio {
	print "PROGRESS: Calculating pNpS ratio.\n";
	# Print out a header for the output file
	print OUT "ContigID\tGeneID\tSynonymousSNPs\tNonSynonymousSNPs\tPotentialSynonymousMutations\tPotentialNonSynonymousMutations\tPnPs\tPositiveDirection\n";
	# Also print out the header for the alternative output file
	print SEQOUT "ContigID\tGeneID\tSnpPosition\tSnpConsensus\tSnpVariable\tCodonAAConsensus\tCodonAAVariable\tPositiveDirection\n";
	# Get the variables from the subroutine call
	my ($FastaHash, $gff3Array, $vcfArray, $AminoAcidHash) = @_;
	foreach my $gff3Line (@$gff3Array) {
		# Remember to reset the counters for the next ratio calculation
		$SynMutation = 0;
		$NonSynMutation = 0;
		$SynMutationCounter = 0;
		$NonSynMutationCounter = 0;
		$SnpCoverage = 0;
		$NormFactor = 0;
		$GeneSequence = 0;
		undef @CoverageArray;
		# print "$gff3Line\n";
		# Set the contig ID
		my $contigID = (split /\t/, $gff3Line)[0];
		# Skip if the sequence is not found in the fasta file
		print STDERR "No ContigID match for $contigID in fasta.\n" unless exists($FastaHash -> {$contigID});
		next unless exists($FastaHash -> {$contigID});
		$geneID = (split /\t/, $gff3Line)[5];
		# print "$geneID\n";
		# Pull out the fasta sequence from hash into variable
		print STDERR "Looking at $contigID $geneID.\n";
		$fastaSequence = $FastaHash -> {$contigID};
		# print "$fastaSequence\n";
		# Set the gene coordinates
		# Be sure to add one because of zero offset
		if ((split /\t/, $gff3Line)[3] < (split /\t/, $gff3Line)[4]) {
			$GeneStart = (split /\t/, $gff3Line)[3] - 1;
			$GeneEnd = (split /\t/, $gff3Line)[4] - 1;
			$PositionPositive = 1;
		} elsif ((split /\t/, $gff3Line)[3] > (split /\t/, $gff3Line)[4]) {
			# Here I also want the start to be in reference to the positive strand
			# This was start is always less than the end
			$GeneStart = (split /\t/, $gff3Line)[4] - 1;
			$GeneEnd = (split /\t/, $gff3Line)[3] - 1;
			$PositionPositive = 0;
		} else {
			# Kill if these conditions are not met because there is an error
			# in the gff3 file.
			die "HOLD ON! The gene start and end positions are the same: $!";
		}
		# Get the length of the gene
		# Add one because length will be offset by one at start
		$GeneLength = $GeneEnd - $GeneStart + 1;
		# Having a gene length that is not divisible by three can result from hotspot IDs (not whole genes)
		# and can cause downstream issues. Therefore make sure it is divisible by three.
		if ($GeneLength % 3 == 0) {
			print STDERR "Gene length is $GeneLength and divisible by three.\n";
		} elsif ($GeneLength % 3 == 2) {
			print STDERR "Gene length is $GeneLength and is not divisible by three.\n";
			print STDERR "This is a problem I can fix, but I want you to know about it.\n";
			$GeneLength = $GeneLength + 1;
			print STDERR "Fixed! Now gene length is $GeneLength and divisible by three.\n";
		} elsif ($GeneLength % 3 == 1) {
			print STDERR "Gene length is $GeneLength and is not divisible by three.\n";
			print STDERR "This is a problem I can fix, but I want you to know about it.\n";
			$GeneLength = $GeneLength + 2;
			print STDERR "Fixed! Now gene length is $GeneLength and divisible by three.\n";
		} else {
			die "LAME! Something bad happened when validating gene length.\n"
		}
		while (my ($VcfKey, $vcfline) = each(%{ $vcfArray -> {$contigID} })) {
			# Set SNP position to variable
			# WARNING: Changing this can affect the modulo operation below for codon pos
			$SnpPosition = (split /\t/, $vcfline)[1];
			# Skip if it is the VCF header
			next if ($SnpPosition eq "Position");
			# print "Gene starts at $GeneStart\n";
			# print "Gene ends at $GeneEnd\n";
			# print "SNP $SnpPosition is outside of gene.\n" if ($SnpPosition < $GeneStart || $SnpPosition > $GeneEnd);
			$SnpConsensus = (split /\t/, $vcfline)[2];
			$SnpVariable = (split /\t/, $vcfline)[18];
			# Print out the fasta to confirm
			# print STDERR "$fastaSequence\n";
			# Extract the gene sequence from the gene's parent fasta
			$GeneSequence = substr($fastaSequence, $GeneStart, $GeneLength);
			# Determine whether SNP position is in current ORF
			# If not, move onto the next SNP
			next if ($SnpPosition < $GeneStart || $SnpPosition > $GeneEnd);
			# And for downstream analysis, set SNP position with respect to gene
			# This is the SNP position count within the gene, so need to add 1
			$SnpPositionCorrNuc = $SnpPosition - $GeneStart + 1;
			# Use a modulo operation to determine SNP position within codon
			# 0 = third position
			# 2 = second position
			# 1 = first position
			$SnpCodonPosition = $SnpPositionCorrNuc % 3;
			# Extract the codon containing the SNP
			if ($SnpCodonPosition == 0) {
				# Gene is in third codon position
				# Extract two sequences before the SNP
				# Add one to offset zero count
				# Subtract two to start at fisrt nucleotide of codon
				# I left the math there as a reminder
				# Cat the first two nucleotides with the third as SNP consensus
				$SnpCodonConsensus = substr($GeneSequence, $SnpPositionCorrNuc-1-2, 2) . $SnpConsensus;
				$SnpCodonVariable = substr($GeneSequence, $SnpPositionCorrNuc-1-2, 2) . $SnpVariable;
			} elsif ($SnpCodonPosition == 2) {
				# Gene is in second codon position
				$SnpCodonConsensus = substr($GeneSequence, $SnpPositionCorrNuc-1-1, 1) . $SnpConsensus . substr($GeneSequence, $SnpPositionCorrNuc-1+1, 1);
				$SnpCodonVariable = substr($GeneSequence, $SnpPositionCorrNuc-1-1, 1) . $SnpVariable . substr($GeneSequence, $SnpPositionCorrNuc-1+1, 1);
			} elsif ($SnpCodonPosition == 1) {
				$SnpCodonConsensus = $SnpConsensus . substr($GeneSequence, $SnpPositionCorrNuc-1+1, 2);
				$SnpCodonVariable = $SnpVariable . substr($GeneSequence, $SnpPositionCorrNuc-1+1, 2);
			} else {
				die "ARG! Killed because of error in codon modulo operation: $!";
			}
			# Reverse compliment the codon if it is on the neg strand
			if ($PositionPositive == 0) {
				# Get the compliment sequences of the gene
				$SnpCodonConsensus =~ tr/ACGT/TGCA/;
				$SnpCodonVariable =~ tr/ACGT/TGCA/;
				# Reverse to finish getting reverse compliment
				$SnpCodonConsensus = reverse $SnpCodonConsensus;
				$SnpCodonVariable = reverse $SnpCodonVariable;
			} elsif ($PositionPositive == 2) {
				die "Error processing ORF directionality: $!";
			}
			$CodonAAConsensus = $AminoAcidHash -> {$SnpCodonConsensus};
			$CodonAAVariable = $AminoAcidHash -> {$SnpCodonVariable};
			# print "Codon consensus is $SnpCodonConsensus for position $SnpPosition.\n";
			# print "Codon variable is $SnpCodonVariable for position $SnpPosition.\n";
			# print "AA consensus is $CodonAAConsensus\n";
			# print "AA variable is $CodonAAVariable\n";
			# Print out the additional output with the substitution patterns
			print SEQOUT "$contigID\t$geneID\t$SnpPosition\t$SnpConsensus\t$SnpVariable\t$CodonAAConsensus\t$CodonAAVariable\t$PositionPositive\n";
			# Print out the VCF line to seperate file for alternate analysis
			# Add a unigue ID to the end for parsing that can be removed later
			print VCFOUT "$vcfline\tUID=$geneID\n";
			# Add a count depending on the mutation type
			++$SynMutation if ($CodonAAConsensus eq $CodonAAVariable);
			++$NonSynMutation if ($CodonAAConsensus ne $CodonAAVariable);
			# Also start adding the number of reads (coverage) to an array to average
			# for normalization.
			$SnpCoverage = (split /\t/, $vcfline)[5];
			push @CoverageArray, $SnpCoverage;
		}
	# If there are no SNPs, then make a note of it (with neutral evolution) and move on
	unless (@CoverageArray) {
		print STDERR "There were no SNPs in $contigID $geneID.\n";
		$PnPsRatio = 1;
		print OUT "$contigID\t$geneID\t$SynMutation\t$NonSynMutation\t$SynMutationCounter\t$NonSynMutationCounter\t$PnPsRatio\t$PositionPositive\n";
		# Reset position so it will throw error if not defined as 1 or 0
		$PositionPositive = 2;
		next;
	}
	# Now to normalize the data, I have to calculate the expected ratios of
	# synonymous and non-synonymous mutations and use this in the final
	# formula.
	# I will start by iterating through the gene sequence one nucleotide at a time
	# Here I want to move through the gene subset from the start to the end
	foreach my $NucleotideIteration (1 .. $GeneLength) {
		# print "Gene start is $GeneStart.\n";
		# print "Gene end is $GeneEnd.\n";
		# print "Position is $NucleotideIteration.\n";
		$nucleotideConsensus = substr($GeneSequence, $NucleotideIteration-1, 1);
		print STDERR "Gene sequence is: $GeneSequence\n";
		print STDERR "Nuc consensus is: $nucleotideConsensus\n";
		# Assign the other three nucleotide letters depending on the initial letter
		if ($nucleotideConsensus eq 'A') {
			$FirstAlter = 'G';
			$SecondAlter = 'C';
			$ThirdAlter = 'T';
		} elsif ($nucleotideConsensus eq 'T') {
			$FirstAlter = 'G';
			$SecondAlter = 'C';
			$ThirdAlter = 'A';
		} elsif ($nucleotideConsensus eq 'G') {
			$FirstAlter = 'T';
			$SecondAlter = 'C';
			$ThirdAlter = 'A';
		}elsif ($nucleotideConsensus eq 'C') {
			$FirstAlter = 'G';
			$SecondAlter = 'T';
			$ThirdAlter = 'A';
		} else {
			print STDERR "This does not look like a nucleotide to me: $nucleotideConsensus\n";
			die "WOAH! There are strange symbols in this fasta\n";
		}
		# print "Nuc consensus is $nucleotideConsensus.\n";
		# Again use a modulo operation to determine nucleotide of interest's position within codon
		# 0 = third position
		# 2 = second position
		# 1 = first position
		$NucleotideCodonPosition = $NucleotideIteration % 3;
		# Extract the codon containing the nucleotide of interest
		if ($NucleotideCodonPosition == 0) {
			# Gene is in third codon position
			# Extract two sequences before the SNP
			# Add one to offset zero count
			# Subtract two to start at fisrt nucleotide of codon
			# I left the math there as a reminder
			# Cat the first two nucleotides with the third as SNP consensus
			$NucleotideCodonConsensus = substr($GeneSequence, $NucleotideIteration-1-2, 2) . $nucleotideConsensus;
			$NucleotideCodonFirst = substr($GeneSequence, $NucleotideIteration-1-2, 2) . $FirstAlter;
			$NucleotideCodonSecond = substr($GeneSequence, $NucleotideIteration-1-2, 2) . $SecondAlter;
			$NucleotideCodonThird = substr($GeneSequence, $NucleotideIteration-1-2, 2) . $ThirdAlter;
		} elsif ($NucleotideCodonPosition == 2) {
			# Gene is in second codon position
			$NucleotideCodonConsensus = substr($GeneSequence, $NucleotideIteration-1-1, 1) . $nucleotideConsensus . substr($GeneSequence, $NucleotideIteration-1+1, 1);
			$NucleotideCodonFirst = substr($GeneSequence, $NucleotideIteration-1-1, 1) . $FirstAlter . substr($GeneSequence, $NucleotideIteration-1+1, 1);
			$NucleotideCodonSecond = substr($GeneSequence, $NucleotideIteration-1-1, 1) . $SecondAlter . substr($GeneSequence, $NucleotideIteration-1+1, 1);
			$NucleotideCodonThird = substr($GeneSequence, $NucleotideIteration-1-1, 1) . $ThirdAlter . substr($GeneSequence, $NucleotideIteration-1+1, 1);
		} elsif ($NucleotideCodonPosition == 1) {
			$NucleotideCodonConsensus = $nucleotideConsensus . substr($GeneSequence, $NucleotideIteration-1+1, 2);
			$NucleotideCodonFirst = $FirstAlter . substr($GeneSequence, $NucleotideIteration-1+1, 2);
			$NucleotideCodonSecond = $SecondAlter . substr($GeneSequence, $NucleotideIteration-1+1, 2);
			$NucleotideCodonThird = $ThirdAlter . substr($GeneSequence, $NucleotideIteration-1+1, 2);
		} else {
			die "ARG! Killed because of error in codon modulo operation 2: $!";
		}
		# print "Consensus is $NucleotideCodonConsensus.\n";
		# If this was a reverse compliment, we need to rev comp the codons before matching
		if ($PositionPositive == 0) {
				# Get the compliment sequences of the gene
				$NucleotideCodonConsensus =~ tr/ACGT/TGCA/;
				$NucleotideCodonFirst =~ tr/ACGT/TGCA/;
				$NucleotideCodonSecond =~ tr/ACGT/TGCA/;
				$NucleotideCodonThird =~ tr/ACGT/TGCA/;
			# Reverse to finish getting reverse compliment
				$NucleotideCodonConsensus = reverse $SnpCodonConsensus;
				$NucleotideCodonFirst = reverse $NucleotideCodonFirst;
				$NucleotideCodonSecond = reverse $NucleotideCodonSecond;
				$NucleotideCodonThird = reverse $NucleotideCodonThird;
			}
		# Get the amino acids that match the codons
		$CodonAAConsensus = $AminoAcidHash -> {$NucleotideCodonConsensus};
		$CodonAAFirst = $AminoAcidHash -> {$NucleotideCodonFirst};
		$CodonAASecond = $AminoAcidHash -> {$NucleotideCodonSecond};
		$CodonAAThird = $AminoAcidHash -> {$NucleotideCodonThird};
		# Count the value for potential non-synonymous and synonymous mutations
		++$SynNucleotideMut if ($CodonAAConsensus eq $CodonAAFirst);
		++$NonSynNucleotideMut if ($CodonAAConsensus ne $CodonAAFirst);
		++$SynNucleotideMut if ($CodonAAConsensus eq $CodonAASecond);
		++$NonSynNucleotideMut if ($CodonAAConsensus ne $CodonAASecond);
		++$SynNucleotideMut if ($CodonAAConsensus eq $CodonAAThird);
		++$NonSynNucleotideMut if ($CodonAAConsensus ne $CodonAAThird);
		# Calculate the value of syn and non-syn mutation rate for this SNP
		$SynMutationValue = $SynNucleotideMut / 3;
		$NonSynMutationValue = $NonSynNucleotideMut / 3;
		# Append the value of each to the counter
		# print STDOUT "Counter is $SynMutationCounter.\n";
		$SynMutationCounter = $SynMutationCounter + $SynMutationValue;
		$NonSynMutationCounter = $NonSynMutationCounter + $NonSynMutationValue;
	}
	# # Get the ratio between the two mutations
	# print "SynMutation is $SynMutation.\n";
	# print "NonSynMutation is $NonSynMutation.\n";
	# print "SynMutationCounter is $SynMutationCounter.\n";
	# print "NonSynMutationCounter is $NonSynMutationCounter.\n";
	# Set the normalization factor as the median of the SNP coverage
	# print join("\n", @CoverageArray);
	# print "Median is ".median(@CoverageArray)."\n";
	# Set the evolutionary pressure values
	$NormFactor = sqrt(median(@CoverageArray))/2.000;
	# print "Norm factor is ".$NormFactor."\n";
	$PnPsRatio = (($NonSynMutation + $NormFactor) / $NonSynMutationCounter) / (($SynMutation + $NormFactor) / $SynMutationCounter);
	print OUT "$contigID\t$geneID\t$SynMutation\t$NonSynMutation\t$SynMutationCounter\t$NonSynMutationCounter\t$PnPsRatio\t$PositionPositive\n";
	# Reset position so it will throw error if not defined as 1 or 0
	$PositionPositive = 2;
	}
}

# Call the subroutines to process the files
my %CodonHash = SetCodonHashInMemory();
# print Dumper \%CodonHash;
my %Fasta = ReadInFasta(\*FASTA);
# print Dumper \%Fasta;
my @Gff = ReadInGff3(\*GFF3);
# print join("\n", @Gff);
my %Vcf = ReadInVCF(\*VCF);
# print Dumper \%Vcf;
CalculatePnPsRatio(\%Fasta, \@Gff, \%Vcf, \%CodonHash);


close(GFF3);
close(FASTA);
close(VCF);
close(OUT);
close(SEQOUT);
close(VCFOUT);
# Get the toal time to run the script
my $end_run = time();
my $run_time = $end_run - $start_run;
print STDERR "\nThe ratios were calculated in $run_time seconds.\n";
