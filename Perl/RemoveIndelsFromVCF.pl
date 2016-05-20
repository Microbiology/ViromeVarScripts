#!/usr/bin/perl
# RemoveIndelsFromVCF.pl
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

# Set use
use strict;
use warnings;

# Set the input files to scalar variables
my $usage = "Usage: perl $0 <INPUT> <OUTPUT>";
my $input = shift or die $usage;
my $output = shift or die $usage;
open (INPUT, "<$input") || die "Unable to open file $input: $!";
open(OUT, ">$output") || die "Unable to write to $output: $!";

# Set variables to be used in the script

# Filter these out in a pretty straight-forward subroutine
sub FilterForIndels {
	while (my $line = <INPUT>) {
		chomp $line;
		my $NucOfInterest = (split /\t/, $line)[18];
		if ($NucOfInterest =~ /\+/ || $NucOfInterest =~ /\-/) {
			next;
		} else {
			print OUT "$line\n";
		}
	}
}

# Run the subroutine
FilterForIndels();

# Close out
close(INPUT);
close(OUT);
print STDERR "The indels have been filtered out of your VCF!\n";
