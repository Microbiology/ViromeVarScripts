#!/usr/bin/perl -w
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

# Use this script to get the HVLs present in DGRs

# Set use for script
use strict;
use warnings;
use Data::Dumper;
# Also measure time to run script
my $start_run = time();

# Set the input from the command line
my $usage = "Usage: perl <PU> <OUT>";
my $mpileup = shift or die $usage;
my $out = shift or die $usage;

open (my $PU, "<", "$mpileup") || die "Unable to open $mpileup: $!";
open (my $OUT, ">", "$out") || die "Unable to open $out: $!";

# Set variables
my $counter = 0;
my @consensus;
my $charactercount = 0;

# Get Array of the contig IDs
while (my $line = <$PU>) {
	chomp $line;
	++$counter;
	next if ($counter == 1);
	my $position = 0;
	print $OUT "\n" if ($counter > 3);
	foreach my $character (split //, $line) {
		my $totallength = scalar @consensus;
		if ($counter == 2) {
			($character = "N") if ($character =~ /\*/);
			push @consensus, $character;
			$charactercount = 0;
			next;
		} elsif ($character =~ /\./) {
			print $OUT "$consensus[$position]\/$consensus[$position]\t";
			$charactercount = 0;
		} elsif ($character =~ /\,/) {
			print $OUT "$consensus[$position]\/$consensus[$position]\t";
			$charactercount = 0;
		} elsif ($character =~ /\s/ && $position == 0) {
			print $OUT "NA\t";
			$charactercount = 2;
		}  elsif ($character =~ /\s/ && $charactercount == 0) {
			print $OUT "NA\t"x($totallength - $position)."\n"."NA\t"x($position + 1);
			$charactercount = 1;
		} elsif ($character =~ /\s/ && $charactercount == 1) {
			print $OUT "NA\t";
		} elsif ($character =~ /\s/ && $charactercount == 2) {
			print $OUT "NA\t";
			$charactercount = 2;
		}  elsif ($character =~ /[ATGCatgc]/) {
			$character =~ tr/atgc/ATGC/;
			print $OUT "$character\/$consensus[$position]\t";
			$charactercount = 0;
		} elsif ($character =~ /\*/) {
			print $OUT "NA\t";
		}
		++$position;
	}
}

close($OUT);
