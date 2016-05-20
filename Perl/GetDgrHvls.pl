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
my $usage = "Usage: perl <DGR> <HVL> <OUT>";
my $dgr = shift or die $usage;
my $hvl = shift or die $usage;
my $out = shift or die $usage;

open (my $DGR, "<", "$dgr") || die "Unable to open $dgr: $!";
open (my $HVL, "<", "$hvl") || die "Unable to open $hvl: $!";
open (my $OUT, ">", "$out") || die "Unable to open $out: $!";

# Set other variables
my %dgrhash;
my %hvlhash;
my $itrid = 0;
my $contigid;

# Read in DGR file
while (my $line = <$DGR>) {
	chomp $line;
	my $contigid = (split "\t", $line)[0];
	my $dgrid = (split "\t", $line)[1];
	my @positions = (split /\t/, $line)[2,3,4,5];
	$dgrhash{$contigid}{$dgrid} = [@positions];
}

# Reset contig variable
$contigid = '';

# Read in the HVL file
while (my $line = <$HVL>) {
	chomp $line;
	$itrid = 0 unless ((split "\t", $line)[0] eq $contigid."_");
	my $dgrid = $itrid;
	++$itrid;
	$contigid = $1 if $line =~ /^(\d+)_.+/;
	next if ($contigid eq "ContigID");
	my @positions = (split /\t/, $line)[1,2];
	$hvlhash{$contigid}{$dgrid} = [@positions];
}

# Get DGRs that contain HVLs
foreach my $key (keys %dgrhash) {
	foreach my $keyid (keys %{ $dgrhash{$key} }) {
		my @dgrposition = @{ $dgrhash{$key}{$keyid} };
		my $flagtwo = 0;
		my $flagone = 0;
		my $finalflag = 0;
		my $firstDGR = 0;
		my $secondDGR = 0;
		foreach my $hvlkey (keys %{ $hvlhash{$key} }) {
			my @hvlposition = @{ $hvlhash{$key}{$hvlkey} };
			# Determine whether one or both pairs contain HVL
			if ($dgrposition[0] < $dgrposition[1]) {
				# Running forward
				$flagone = 1 if ($dgrposition[0] <= $hvlposition[0] && $hvlposition[0] <= $dgrposition[1] && $dgrposition[0] <= $hvlposition[1] && $hvlposition[1] <= $dgrposition[1]);
				++$firstDGR if ($dgrposition[0] <= $hvlposition[0] && $hvlposition[0] <= $dgrposition[1] && $dgrposition[0] <= $hvlposition[1] && $hvlposition[1] <= $dgrposition[1]);
			} elsif ($dgrposition[0] > $dgrposition[1]) {
				# Running backwards
				$flagone = 1 if ($dgrposition[1] <= $hvlposition[0] && $hvlposition[0] <= $dgrposition[0] && $dgrposition[1] <= $hvlposition[1] && $hvlposition[1] <= $dgrposition[0]);
				++$firstDGR if ($dgrposition[1] <= $hvlposition[0] && $hvlposition[0] <= $dgrposition[0] && $dgrposition[1] <= $hvlposition[1] && $hvlposition[1] <= $dgrposition[0]);
			}
			if ($dgrposition[2] < $dgrposition[3]) {
				# Running forward
				$flagtwo = 2 if ($dgrposition[2] <= $hvlposition[0] && $hvlposition[0] <= $dgrposition[3] && $dgrposition[2] <= $hvlposition[1] && $hvlposition[1] <= $dgrposition[3]);
				++$secondDGR if ($dgrposition[2] <= $hvlposition[0] && $hvlposition[0] <= $dgrposition[3] && $dgrposition[2] <= $hvlposition[1] && $hvlposition[1] <= $dgrposition[3]);
			} elsif ($dgrposition[2] > $dgrposition[3]) {
				# Running backwards
				$flagtwo = 2 if ($dgrposition[3] <= $hvlposition[0] && $hvlposition[0] <= $dgrposition[2] && $dgrposition[3] <= $hvlposition[1] && $hvlposition[1] <= $dgrposition[2]);
				++$secondDGR if ($dgrposition[3] <= $hvlposition[0] && $hvlposition[0] <= $dgrposition[2] && $dgrposition[3] <= $hvlposition[1] && $hvlposition[1] <= $dgrposition[2]);
			}
		}
		$finalflag = $flagone + $flagtwo;
		next if ($finalflag == 0);
		my $arraystring = join("\t", @dgrposition);
		print $OUT "$key\tFlag=$finalflag\tHVL1=$firstDGR\tHVL2=$secondDGR\t$arraystring\n";
	}
}

close($DGR);
close($HVL);
close($OUT);
