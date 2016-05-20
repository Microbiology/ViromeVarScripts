#!/usr/bin/perl -w
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

# Set use for script
use strict;
use warnings;
use Data::Dumper;
# Also measure time to run script
my $start_run = time();

# Set the input from the command line
my $usage = "Usage: perl <REVT> <REPEAT> <OUT>";
my $revt = shift or die $usage;
my $repeat = shift or die $usage;
my $out = shift or die $usage;
# And open the files
open (my $REVT, "<", "$revt") || die "Unable to open $revt: $!";
open (my $REPEAT, "<", "$repeat") || die "Unable to open $repeat: $!";
open (my $OUT, ">", "$out") || die "Unable to open $out: $!";

# Assign variables
my %rtHash;
my %candidatepositions;
my $candidatepositions;
my %exclusion;
my @positions;
my $match;
my $contigid;
my $idcounter = 0;
my $flag = 0;

# Get the contig IDs that contain reverse transcriptase elements
print "PROGRESS: Parsing RT.\n";
while (my $line = <$REVT>) {
	chomp $line;
	$match = $1 if $line =~ /^(\d+)_\..+/;
	$rtHash{$match} = 1;
}

# Get the positions of the potential DGRs
while (my $line = <$REPEAT>) {
	chomp $line;
	$contigid = $1 if $line =~ /^(\d+)_.+/;
	next unless (exists $rtHash{$contigid});
	$idcounter = 0 unless (exists $candidatepositions{$contigid});
	@positions = (split /\t/, $line)[6,7,8,9];
	# Don't save if either end is within it's pair
	$flag = 0;
	if ($positions[0] < $positions[1]) {
		# Running forward
		$flag = 1 if ($positions[0] <= $positions[2] && $positions[2] <= $positions[1]);
		$flag = 1 if ($positions[0] <= $positions[3] && $positions[3] <= $positions[1]);
	} elsif ($positions[1] < $positions[0]) {
		# Running backwards
		$flag = 1 if ($positions[1] <= $positions[2] && $positions[2] <= $positions[0]);
		$flag = 1 if ($positions[1] <= $positions[3] && $positions[3] <= $positions[0]);
	}
	# Do the same but switch compared pairs
	if ($positions[2] < $positions[3]) {
		# Running forwards
		$flag = 1 if ($positions[2] <= $positions[0] && $positions[0] <= $positions[3]);
		$flag = 1 if ($positions[2] <= $positions[1] && $positions[1] <= $positions[3]);
	} elsif ($positions[1] < $positions[0]) {
		# Running backwards
		$flag = 1 if ($positions[3] <= $positions[0] && $positions[0] <= $positions[2]);
		$flag = 1 if ($positions[3] <= $positions[1] && $positions[1] <= $positions[2]);
	}
	next if ($flag == 1);
	$candidatepositions{$contigid}{$idcounter} = [@positions];
	++$idcounter;
}

# Make sure the pair is not contained in any of the other pairs from
# the same contig
foreach my $contigidhash (keys %candidatepositions) {
	# while (my ($key, @value) = each %{ $candidatepositions{$contigidhash} }) {
	undef %exclusion;
	foreach my $key (keys %{ $candidatepositions{$contigidhash} }) {
		my @value = @{ $candidatepositions{$contigidhash}{$key} };
		my $counter = 0;
		foreach my $keyitr (keys %{ $candidatepositions{$contigidhash} }) {
			$flag = 0;
			my @valueitr = @{ $candidatepositions{$contigidhash}{$keyitr} };
			# Pass if comparing same key values
			next if ($key == $keyitr);
			if ($valueitr[0] < $valueitr[1]) {
				# Running forward
				$flag = 1 if ($valueitr[0] <= $value[0] && $value[0] <= $valueitr[1]);
				$flag = 1 if ($valueitr[0] <= $value[1] && $value[1] <= $valueitr[1]);
				$flag = 1 if ($valueitr[0] <= $value[2] && $value[2] <= $valueitr[1]);
				$flag = 1 if ($valueitr[0] <= $value[3] && $value[3] <= $valueitr[1]);
			} elsif ($valueitr[1] < $valueitr[0]) {
				# Running backwards
				$flag = 1 if ($valueitr[1] <= $value[0] && $value[0] <= $valueitr[0]);
				$flag = 1 if ($valueitr[1] <= $value[1] && $value[1] <= $valueitr[0]);
				$flag = 1 if ($valueitr[1] <= $value[2] && $value[2] <= $valueitr[0]);
				$flag = 1 if ($valueitr[1] <= $value[3] && $value[3] <= $valueitr[0]);
			}
			# Do the same with second of reference pair
			if ($valueitr[2] < $valueitr[3]) {
				# Running forward
				$flag = 1 if ($valueitr[2] <= $value[0] && $value[0] <= $valueitr[3]);
				$flag = 1 if ($valueitr[2] <= $value[1] && $value[1] <= $valueitr[3]);
				$flag = 1 if ($valueitr[2] <= $value[2] && $value[2] <= $valueitr[3]);
				$flag = 1 if ($valueitr[2] <= $value[3] && $value[3] <= $valueitr[3]);
			} elsif ($valueitr[3] < $valueitr[2]) {
				# Running backwards
				$flag = 1 if ($valueitr[3] <= $value[0] && $value[0] <= $valueitr[2]);
				$flag = 1 if ($valueitr[3] <= $value[1] && $value[1] <= $valueitr[2]);
				$flag = 1 if ($valueitr[3] <= $value[2] && $value[2] <= $valueitr[2]);
				$flag = 1 if ($valueitr[3] <= $value[3] && $value[3] <= $valueitr[2]);
			}
			$exclusion{$keyitr} = 1 if ($flag == 1);
			$counter = 1 if ($flag == 1);
		}
		my $arraystring = join("\t", @value);
		print $OUT "$contigidhash\t$key\t$arraystring\n" if ($counter == 0);
		next if ($counter == 0);
		print $OUT "$contigidhash\t$key\t$arraystring\n" unless ($exclusion{$key});
	}
}

