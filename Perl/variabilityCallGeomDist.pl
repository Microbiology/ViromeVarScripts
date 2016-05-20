# variabilityCallGeomDist.pl
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
my $usage = "Usage: perl $0 <INFILE> <OUTFILE> <DISTOUT>";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $distout = shift or die $usage;
open(IN, "<$infile") || die "Unable to open $infile: $!";
# This assumes VarScan VCF output
open(OUTFILE, ">$outfile") || die "Unable to write to $outfile: $!";
open(DISTOUT, ">$distout") || die "Unable to write to $outfile: $!";


# Set the variables to be used in the script
my %MasterInput;
my %DistanceHolder;
my @MeanCalulator;
my $SnpHolder = 1; # Genome position is 1 based counting
my $SnpPosition = 0;
my $SnpDistance = 0;
my $SnpDistMean = 0;


# Set subroutine for getting mean of SNP distance array
sub mean {
    return sum(@_)/@_;
}

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

# Iterate through the master hash by contig ID, and get the
# average distance between SNPs.
sub GetSnpDistances {
	my %ContigIdHash = @_;
	# Print out the distout header
	print DISTOUT "ContigID\tSnpPosition\tDistanceFromPreviousSnp\n";
	# print Dumper \%ContigIdHash;
	while (my ($ContigKey, $ContigValue) = each %ContigIdHash) {
		# print join("\n", @$ContigValue);
		foreach my $line (@$ContigValue) {
			$SnpPosition = (split /\t/, $line)[1];
			# Add the distance between SNPs to an array
			$SnpDistance = $SnpPosition - $SnpHolder;
			# print "Position is $SnpPosition\n";
			# print "Distance is $SnpDistance\n";
			push @MeanCalulator, $SnpDistance;
			$SnpHolder = $SnpPosition;
			# print distances out to different folder
			print DISTOUT "$ContigKey\t$SnpPosition\t$SnpDistance\n";
		}
		my $SnpDistMean = mean(@MeanCalulator);
		# print join("\n", @MeanCalulator);
		# print "Mean is $SnpDistMean\n";
		# Save distance into new hash
		$DistanceHolder{$ContigKey} = $SnpDistMean;
		# Reset back to one for next iteration
		$SnpHolder = 1;
		@MeanCalulator = ();
		print STDERR "Completed SNP analysis for contig $ContigKey.\n"
	}
	return %DistanceHolder;
}

# Set sub to print final hash to output file
sub PrintHashToOutFile {
	my %InputHash = @_;
	print OUTFILE "ContigID\tAverageSnpDistance\n";
	while (my ($Key, $Value) = each %InputHash) {
		print OUTFILE "$Key\t$Value\n";
	}
}

# Run the subroutines called above
# Print hashes for debugging script
my %MasterInputOutside = ReadIntoMemory();
# print Dumper \%MasterInputOutside;
my %DistanceHashOutside = GetSnpDistances(%MasterInputOutside);
# print Dumper \%DistanceHashOutside;
PrintHashToOutFile(%DistanceHashOutside);

close(IN);
close(OUTFILE);
close(DISTOUT);
print STDERR "SNP distance averages have been calculated.\n";


