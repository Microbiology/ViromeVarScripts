#!/usr/bin/perl
# GetMappedReadIdFromSam.pl
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

# Set use
use strict;
use warnings;
# Also measure time to run script
my $start_run = time();

# Set files to scalar variables
my $usage = "Usage: perl $0 <SAM> <OUT> <LOG>";
my $sam = shift or die $usage;
my $out = shift or die $usage;
my $log = shift or die $usage;
open(SAM, "<$sam") || die "Unable to open $sam: $!";
open(OUT, ">$out") || die "Unable to write to $out: $!";
open(LOG, ">$log") || die "Unable to write to $log: $!";

# Set variables to use in script
my $lineID = 0;

# Print a header for the output
print LOG "ContigID\tMatchStatus\tSamCode\n";
print OUT "ContigID\tMatchID\n";

# Iterate through the Sam file (just read off disk once, and probably
# a small file here anyways)
while (my $line = <SAM>) {
	chomp $line;
	# Convert file to tab delimited if that is not the case
	$lineID = (split /\t/, $line)[0];
	print STDERR "LineID is $lineID.\n";
	# Skip SAM header lines
	next if ($lineID eq '@HD');
	next if ($lineID eq '@SQ');
	next if ($lineID eq '@PG');
	# Log reasult of whether the read had a match
	if ((split /\t/, $line)[2] eq '*') {
		print LOG "$lineID\tNoMatch\t".(split /\t/, $line)[1]."\n";
	} else {
		print LOG "$lineID\tFoundMatch\t".(split /\t/, $line)[1]."\n";
		print OUT "$lineID\t".(split /\t/, $line)[2]."\n";
	}
}

close(SAM);
close(OUT);
close(LOG);
# Get the toal time to run the script
my $end_run = time();
my $run_time = $end_run - $start_run;
print STDERR "Contig IDs were pulled in $run_time seconds.\n";

