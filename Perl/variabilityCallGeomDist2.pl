# variabilityCallGeomDist2.pl
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

# Set use
use strict;
use warnings;
use Data::Dumper;
use List::Util qw(sum);

# WARNING: This script required the VCF input be ordered

# Set files to scalar variables
my $usage = "Usage: perl $0 <INFILE> <GEOMDIST> <OUTFILE> <SUBVCF>";
my $infile = shift or die $usage;
my $geomDist = shift or die $usage;
my $outfile = shift or die $usage;
my $subvcf = shift or die $usage;
open(IN, "<$infile") || die "Unable to open $infile: $!";
# This assumes VarScan VCF output
open(GD, "<$geomDist") || die "Unable to open $infile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
# This last one is unfortunately not all that helpful
open(SUBVCF, ">$subvcf") || die "Unable to write to $subvcf: $!";

# Set the variables to be used in the script
my %MasterInput;
my %ReferenceInput;
my %ContigIdHash;
my %ReferenceHash;
my $SnpPosition = 0;
my $SnpDistance = 0;
my $SnpDistMean = 0;
my $flag = 0;
my $SnpCounter = 0;
my $SnpHolder = 1;
my $VariableStart = 0;
my $VariableEnd = 0;
my $VariableLength = 0;

# Set subroutine to read input file into memory as array
# Add the array to a hash with the contig ID as the key
sub ReadIntoMemory {
	while (my $line = <IN>) {
		next if ((split /\t/, $line)[0] eq "Chrom");
		chomp $line;
		my $HashContigID = (split /\t/, $line)[0];
		push ( @{$MasterInput{$HashContigID}}, $line);
	}
	return %MasterInput;
}

# Set subroutine to read input file into memory as array
# Add the array to a hash with the contig ID as the key
sub ReadReferenceIntoMem {
	while (my $line = <GD>) {
		next if ((split /\t/, $line)[0] eq "ContigID");
		chomp $line;
		my $HashContigID = (split /\t/, $line)[0];
		my $QuantileValue = (split /\t/, $line)[1];
		$ReferenceInput{$HashContigID} = $QuantileValue;
	}
	return %ReferenceInput;
}

sub GetSnpHotspots {
	# Print header for the output file
	print OUT "ContigID\t".
		"VariableLociStart\t".
		"VariableLociEnd\t".
		"NumberOfSnps\t".
		"VariableLociLength\n";
	my ($ContigIdHash, $ReferenceHash) = @_;
	# print Dumper \%ContigIdHash;
	# print Dumper \%ReferenceHash;
	while (my ($ContigKey, $ContigValue) = each %$ContigIdHash) {
		foreach my $line (@$ContigValue) {
			my $ContigGD = %$ReferenceHash -> {$ContigKey};			
			$SnpPosition = (split /\t/, $line)[1];
			# Add the distance between SNPs to an array
			$SnpDistance = $SnpPosition - $SnpHolder;
			# Test whether this distance is less than GD
			if ($SnpDistance <= $ContigGD && $flag == 0) {
				$VariableStart = $SnpHolder;
				$SnpHolder = $SnpPosition;
				$flag = 1;
				++$SnpCounter;
			} elsif ($SnpDistance <= $ContigGD && $flag == 1) {
				$SnpHolder = $SnpPosition;
				++$SnpCounter;
				$flag = 1;
			} elsif ($SnpDistance > $ContigGD && $flag == 1) {
				# This ends the string of hypervariable SNP loci
				++$SnpCounter;
				$VariableEnd = $SnpHolder;
				$VariableLength = $VariableEnd - $VariableStart + 1;
				print OUT "$ContigKey\t".
					"$VariableStart\t".
					"$VariableEnd\t".
					"$SnpCounter\t".
					"$VariableLength\n";
				print SUBVCF "$line\n";
				$flag = 0;
				$SnpHolder = $SnpPosition;
				$SnpCounter = 0;
			} elsif ($SnpDistance > $ContigGD && $flag == 0) {
				$SnpHolder = $SnpPosition;
				$SnpCounter = 0;
				# Need to reset the flag back to zero in case the SNP
				# region ended with the end of the contig, to prevent
				# it from carying over.
				$flag = 0;
			}
		}
		$SnpHolder = 1;
		$SnpPosition = 0;
		$flag = 0;
		print STDERR "Completed SNP analysis for contig $ContigKey.\n"
	}
}

my %MasterInputOutside = ReadIntoMemory();
# print Dumper \%MasterInputOutside;
my %ReferenceInputOutside = ReadReferenceIntoMem();
# print Dumper \%ReferenceInputOutside;
GetSnpHotspots(\%MasterInputOutside,\%ReferenceInputOutside);

close(IN);
close(GD);
close(OUT);
close(SUBVCF);
print STDERR "Hypervariable SNP loci have been calculated.\n";
