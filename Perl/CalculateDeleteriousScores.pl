#!/usr/bin/perl
# CalculateDeleteriousScores.pl
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

use strict;
use warnings;
use Data::Dumper;

# Also measure time to run script
my $start_run = time();

# Set the unput variables
my $usage = "Usage: perl $0 <VCF> <GFF3> <FASTA> <SCORES> <OUTPUT>";
# And set the variables to open and write
my $vcf = shift or die $usage;
my $gff3 = shift or die $usage;
my $fasta = shift or die $usage;
my $scores = shift or die $usage;
my $output = shift or die $usage;
# Open the files for reading or writing
open(VCF, "<$vcf") || die "I was not able to open the file $vcf: $!";
open(GFF3, "<$gff3") || die "I was not able to open the file $gff3: $!";
open(FASTA, "<$fasta") || die "I was not able to open the file $fasta: $!";
open(SCORES, "<$scores") || die "I was not able to open the file $scores: $!";
open(OUTPUT, ">$output") || die "I am not able to write to the file $output: $!";

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

sub ReadInFasta {
	my $FastaID = 0;
	my %FastaHash;
	my $FastaErrorCounter = 0;
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
	my @gffArray;
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
	my %VcfHash;
	print "PROGRESS: Loading in VCF.\n";
	# Save the VCF file into memory for faster access
	my $VcfInput = shift;
	while (my $line = <$VcfInput>) {
		chomp $line;
		$line =~ s/_//g;
		my $VcfContigID = (split /\t/, $line)[0];
		my $VcfPosition = (split /\t/, $line)[1];
		# print "$VcfContigID\n";
		next if ($VcfContigID eq "Chrom");
		# Save this in as a nested hash within a hash
		$VcfHash{$VcfContigID}{$VcfPosition} = $line;
	}
	return %VcfHash;
}

sub ReadInMatrix {
	my %ScoreReference;
	my %MatrixScores;
	my $MatrixScore;
	my $InterpScore;
	print "PROGRESS: Loading in score matrix.\n";
	my $scoreInput = shift;
	while (my $line = <$scoreInput>) {
		chomp $line;
		my $MatrixContigID = (split /\t/, $line)[0];
		my $MatrixOrfID = (split /\t/, $line)[1];
		my $position = (split /\t/, $line)[2];
		if ($position =~ /pos/) {
			foreach my $iterator (3 .. 22) {
				my $AminoAcidID = (split /\t/, $line)[$iterator];
				# Set this into a hash
				$ScoreReference{$iterator} = $AminoAcidID;
			}
		} else {
			# Iterate through each of the potential substitutions
			foreach my $iterator (3 .. 22) {
				$MatrixScore = (split /\t/, $line)[$iterator];
				$InterpScore = $ScoreReference{$iterator};
				# Set the appropriate amino acid ID
				$MatrixScores{$MatrixContigID}{$MatrixOrfID}{$position}{$InterpScore} = $MatrixScore;
				# This results in a hash with:
				# contig id containing orf ->
				# orf id containing amino acid variant ->
				# amino acid position ->
				# amino acid identifier ->
				# associated deleterious score
			}
		}
	}
	return %MatrixScores;
}

sub GetGeneFasta {
	# Print header for output
	print OUTPUT "ContigID\tGeneID\tContigPosition\tAminoAcidPosition\tDeltScore\n";
	# Input the VCF, GFF3, and Fasta to get amino acid associated with SNP
	# Also don't forget the amino acid hash
	my %GeneFastaHash;
	my ($FastaHash, $gff3Array, $SubVcfHash, $AminoAcidHash, $SavScoreHash) = @_;
	foreach my $gff3Line (@$gff3Array) {
		my $GeneStart = 0;
		my $GeneEnd = 0;
		my $PositionPositive;
		# print "$gff3Line\n";
		# Set the contig ID
		my $contigID = (split /\t/, $gff3Line)[0];
		# Skip if the sequence is not found in the fasta file
		print STDERR "No ContigID match for $contigID in fasta.\n" unless exists($FastaHash -> {$contigID});
		next unless exists($FastaHash -> {$contigID});
		my $geneID = (split /\t/, $gff3Line)[5];
		# Skip this iteration if the gene is not in the Score Matrix hash
		print "$contigID $geneID did not have del scores.\n" unless exists($SavScoreHash -> {$contigID}{$geneID});
		next unless exists($SavScoreHash -> {$contigID}{$geneID});
		# Pull out the fasta sequence from hash into variable
		print STDERR "Looking at $contigID $geneID.\n";
		my $fastaSequence = $FastaHash -> {$contigID};
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
			die "HOLD ON! The gene start and end positions are the same: $!\n";
		}
		# Get the length of the gene
		# Add one because length will be offset by one at start
		my $GeneLength = $GeneEnd - $GeneStart + 1;
		# Get the gene sequence
		my $GeneSequence = substr($fastaSequence, $GeneStart, $GeneLength);
		# Now incorporte the variant information
		while (my ($VcfKey, $vcfline) = each(%{ $SubVcfHash -> {$contigID} })) {
			my $SnpPosition = (split /\t/, $vcfline)[1];
			# Skip if it is the VCF header
			next if ($SnpPosition eq "Position");
			# Determine whether SNP position is in current ORF
			# If not, move onto the next SNP
			next if ($SnpPosition < $GeneStart || $SnpPosition > $GeneEnd);
			# Set SNP position with respect to gene
			# This is the SNP position count within the gene, so need to add 1
			my $SnpPositionCorrNuc = $SnpPosition - $GeneStart + 1;
			# Get the SNP nucleotides
			my $SnpConsensus = (split /\t/, $vcfline)[2];
			my $SnpVariable = (split /\t/, $vcfline)[18];
			# Use a modulo operation to determine SNP position within codon
			# 0 = third position
			# 2 = second position
			# 1 = first position
			my $SnpCodonPosition = $SnpPositionCorrNuc % 3;
			# Extract the codon containing the SNP
			my $SnpCodonConsensus;
			my $SnpCodonVariable;
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
			}
			my $CodonAAConsensus = $AminoAcidHash -> {$SnpCodonConsensus};
			my $CodonAAVariable = $AminoAcidHash -> {$SnpCodonVariable};
			# Get the gene position as amino acid
			# Remainers are removed here
			my $GeneAminoAcidPosition = int( $SnpPositionCorrNuc / 3);
			print STDERR "Amino acid position is $GeneAminoAcidPosition\n";
			print STDERR "Amino acid value is $CodonAAVariable\n";
			print STDERR "Gene ID is $geneID\n";
			# Use the amino acid to get the score associated with the variant
			my $VariantScore = $SavScoreHash -> {$contigID}{$geneID}{$GeneAminoAcidPosition}{$CodonAAVariable};
			print OUTPUT "$contigID\t$geneID\t$SnpPosition\t$GeneAminoAcidPosition\t$VariantScore\n";
		}
	}
}

# Note that when this runs, it will be able to deal with stop codons since those are not included in the
# delterious score prediction.


# Call the subroutines to process the files
my %CodonHash = SetCodonHashInMemory();
# print Dumper \%CodonHash;
my %Fasta = ReadInFasta(\*FASTA);
# print Dumper \%Fasta;
my @Gff = ReadInGff3(\*GFF3);
# print join("\n", @Gff);
my %Vcf = ReadInVCF(\*VCF);
# print Dumper \%Vcf;
my %ScoreMatrixResult = ReadInMatrix(\*SCORES);
# print Dumper \%ScoreMatrixResult;
GetGeneFasta(\%Fasta, \@Gff, \%Vcf, \%CodonHash, \%ScoreMatrixResult);

# Get the toal time to run the script
my $end_run = time();
my $run_time = $end_run - $start_run;
print STDERR "\nThe scores were calculated in $run_time seconds.\n";

