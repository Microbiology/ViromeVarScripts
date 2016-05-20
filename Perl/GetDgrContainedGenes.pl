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
my $genes = shift or die $usage;
my $ant = shift or die $usage;
my $rt = shift or die $usage;
my $generef = shift or die $usage;
my $out = shift or die $usage;
my $gff3 = shift or die $usage;

open (my $DGR, "<", "$dgr") || die "Unable to open $dgr: $!";
open (my $GENES, "<", "$genes") || die "Unable to open $genes: $!";
open (my $ANT, "<", "$ant") || die "Unable to open $ant: $!";
open (my $RT, "<", "$rt") || die "Unable to open $rt: $!";
open (my $GENEREF, "<", "$generef") || die "Unable to open $generef: $!";
open (my $OUT, ">", "$out") || die "Unable to open $out: $!";
open (my $GFF3, ">", "$gff3") || die "Unable to open $gff3: $!";

# Set other variables
my %dgrhash;
my %hvlhash;
my %annotationHash;
my %RtHash;
my %GeneReferenceHash;
my $itrid = 0;
my $contigid;
my $dgrid = 0;
my $annotation;

# Read in DGR file
while (my $line = <$DGR>) {
	chomp $line;
	my $contigid = (split "\t", $line)[0];
	$dgrid = $dgrid++;
	my @positions = (split /\t/, $line)[4,5,6,7];
	$dgrhash{$contigid}{$dgrid} = [@positions];
}

# Read in RT blast file
while (my $line = <$RT>) {
	chomp $line;
	$line =~ s/\s+/\t/g;
	my $SplitId = (split "\t", $line)[0];
	(my $contigid = $SplitId) =~ s/_\..+//;
	(my $orfId = $SplitId) =~ s/^\d+_\.//;
	$RtHash{$contigid} = $orfId;
}

# Read in annotation reference
while (my $line = <$ANT>) {
	chomp $line;
	$line =~ s/\s+/\t/g;
	my $SplitId = (split "\t", $line)[3];
	(my $contigid = $SplitId) =~ s/_\..+//;
	(my $orfId = $SplitId) =~ s/^\d+_\.//;
	my $annotation = (split "\t", $line)[0];
	$annotationHash{$contigid}{$orfId} = $annotation ;
}

# Reset contig variable
$contigid = '';

# Neet to make a gene position reference
while (my $line = <$GENEREF>) {
	next unless ($line =~ /^\>/);
	chomp $line;
	$line =~ s/\s+/\t/g;
	my $contigid = (split "\t", $line)[1];
	$contigid =~ s/_//g;
	my $orfId = (split "\t", $line)[0];
	$orfId =~ s/^\S+_\.//;
	$orfId =~ s/_\d$//;
	my @positions = (split "\t", $line)[2,3];
	$GeneReferenceHash{$contigid}{$orfId} = [@positions];
}

# Get the DGRs inside a predicted gene
while (my $line = <$GENES>) {
	next unless ($line =~ /^\>/);
	$line =~ s/\s+/\t/g;
	my $contigid = (split "\t", $line)[1];
	$contigid =~ s/_//g;
	my $geneId = (split "\t", $line)[0];
	$geneId =~ s/^\S+_\.//;
	my @positions = (split "\t", $line)[2,3];
	foreach my $keyid (keys %{ $dgrhash{$contigid} }) {
		my @dgrposition = @{ $dgrhash{$contigid}{$keyid} };
		my $flagtwo = 0;
		my $flagone = 0;
		my $finalflag = 0;
		if ($positions[0] < $positions[1]) {
			if ($dgrposition[0] < $dgrposition[1]) {
				$flagone = 1 if ($positions[0] <= $dgrposition[0] && $dgrposition[0] <= $positions[1] && $positions[0] <= $dgrposition[1] && $dgrposition[1] <= $positions[1]);
			} elsif ($dgrposition[0] > $dgrposition[1]) {
				$flagone = 1 if ($positions[1] <= $dgrposition[0] && $dgrposition[0] <= $positions[0] && $positions[1] <= $dgrposition[1] && $dgrposition[1] <= $positions[0]);
			}
			if ($dgrposition[2] < $dgrposition[3]) {
				$flagone = 1 if ($positions[0] <= $dgrposition[2] && $dgrposition[2] <= $positions[1] && $positions[0] <= $dgrposition[3] && $dgrposition[3] <= $positions[1]);
			} elsif ($dgrposition[2] > $dgrposition[3]) {
				$flagone = 1 if ($positions[1] <= $dgrposition[2] && $dgrposition[2] <= $positions[0] && $positions[1] <= $dgrposition[3] && $dgrposition[3] <= $positions[0]);
			}
		} elsif ($positions[1] < $positions[0]) {
			if ($dgrposition[0] < $dgrposition[1]) {
				$flagtwo = 1 if ($positions[0] <= $dgrposition[0] && $dgrposition[0] <= $positions[1] && $positions[0] <= $dgrposition[1] && $dgrposition[1] <= $positions[1]);
			} elsif ($dgrposition[0] > $dgrposition[1]) {
				$flagtwo = 1 if ($positions[1] <= $dgrposition[0] && $dgrposition[0] <= $positions[0] && $positions[1] <= $dgrposition[1] && $dgrposition[1] <= $positions[0]);
			}
			if ($dgrposition[2] < $dgrposition[3]) {
				$flagtwo = 1 if ($positions[0] <= $dgrposition[2] && $dgrposition[2] <= $positions[1] && $positions[0] <= $dgrposition[3] && $dgrposition[3] <= $positions[1]);
			} elsif ($dgrposition[2] > $dgrposition[3]) {
				$flagtwo = 1 if ($positions[1] <= $dgrposition[2] && $dgrposition[2] <= $positions[0] && $positions[1] <= $dgrposition[3] && $dgrposition[3] <= $positions[0]);
			}
		}
		$finalflag = $flagone + $flagtwo;
		next if ($finalflag == 0);
		# Annotate the genes
		if (exists $annotationHash{$contigid}{$geneId}) {
			$annotation = $annotationHash{$contigid}{$geneId};
		} else {
			$annotation = "Unknown";
		}
		my $arraystring = join("\t", @dgrposition);
		print $OUT "$contigid\t$keyid\tFlag=$finalflag\t$annotation\t$geneId\t$arraystring\n";
		# Print out the gene with the DGR
		print $GFF3 "$contigid\tGLIMMER\t$annotation\t$positions[0]\t$positions[1]\t$annotation\t.\t1\tID=$annotation; NOTE: DGR Prediction;\n";
		# Print out the first and second DGR array
		print $GFF3 "$contigid\tGLIMMER\tDGR-1\t$dgrposition[0]\t$dgrposition[1]\tDGR-1\t.\t1\tID=DGR-1; NOTE: DGR Prediction;\n";
		print $GFF3 "$contigid\tGLIMMER\tDGR-2\t$dgrposition[2]\t$dgrposition[3]\tDGR-2\t.\t1\tID=DGR-2; NOTE: DGR Prediction;\n";
		# Get the RT gene ORF ID
		my $rtOrf = $RtHash{$contigid};
		print "Contig ID is $contigid\n";
		print "Gene is $rtOrf\n";
		print "$GeneReferenceHash{$contigid}{$rtOrf}\n";
		# Use that ID to get coordinates
		my @RTposition = @{ $GeneReferenceHash{$contigid}{$rtOrf} };
		my $arraystringRT = join("\t", @RTposition);
		# Print this into the gff3 file
		print $GFF3 "$contigid\tGLIMMER\tReverseTranscriptase\t$arraystringRT\tReverseTranscriptase\t.\t1\tID=ReverseTranscriptase; NOTE: DGR prediction;\n";
	}
}

close($DGR);
close($GENES);
close($ANT);
close($RT);
close($OUT);
close($GFF3);
