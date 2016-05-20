#!/usr/bin/perl
# HotspotSlidingWindowGff3.pl
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
my $usage = "Usage: perl <GFF3> <HOTSPOTS> <OUT>";
my $gff3 = shift or die $usage;
my $hotspots = shift or die $usage;
my $out = shift or die $usage;
# And open the files
open (GFF3, "<$gff3") || die "Unable to open $gff3: $!";
open (HOTSPOTS, "<$hotspots") || die "Unable to open $hotspots: $!";
open (OUT, ">$out") || die "Unable to open $out: $!";
# NOTE: The resulting gff3 will have the sliding window coordinated within the entire contig,
# given the corresponding gene it is found in.

# Set how long the adjacent region should be over the HVL
my $windowFactor = 5;

# Set the variables to be used in the script
my @gffArray;
my @hotspotArray;
my @SubHotspotArray;
my @SubGffArray;
my $SubGffArray = 0;
my $SubHotspotArray = 0;
my $HotID = 0;
my $HotStart = 0;
my $HotEnd = 0;
my $HotLength = 0;
my $GffContigID = 0;
my $GffOrfID = 0;
my $GffStart = 0;
my $GffEnd = 0;
my $WindowStart = 0;
my $WindowEnd = 0;
my $WindowLength = 0;
my $SurroundWindowLength = 0;
my $WindowSlideBeforeEnd = 0;
my $WindowSlideBeforeStart = 0;
my $WindowSlideAfterStart = 0;
my $WindowSlideAfterEnd = 0;
my $HotUnique = 0;


# Read the gff3 file into memory as an array
sub ReadInGff3 {
	print "PROGRESS: Loading in gene annotation file.\n";
	while (my $line = <GFF3>) {
		chomp $line;
		push @gffArray, $line;
	}
	return @gffArray;
}

# Read the hotspot coordinates into memory as a array
sub ReadInHotspots {
	print "PROGRESS: Loading in hotspot coordinate file.\n";
	while (my $line = <HOTSPOTS>) {
		chomp $line;
		push @hotspotArray, $line;
	}
	return @hotspotArray;
}

# Generate the gff3 file with the hotspot and surrounding regions
sub GenerateHotspotGff3 {
	print "PROGRESS: Generating gff3 with hotspot coordinates.\n";
	# Get the variables for the subroutine
	my ($SubGffArray, $SubHotspotArray) = @_;
	# print join("\n", @$SubGffArray);
	# Start by iterating through the Hotspot array to get hotspot coordinates
	foreach my $HotspotLine (@$SubHotspotArray) {
		# Skip the header row
		next if ($HotspotLine =~ /ContigID/);
		# Pull out the contig ID that we are dealing with here
		$HotID = (split /\t/, $HotspotLine)[0];
		# Pull out the start and end coordinates
		$HotStart = (split /\t/, $HotspotLine)[1];
		$HotEnd = (split /\t/, $HotspotLine)[2];
		# Also get the unique ID for that SNP
		$HotUnique = (split /\t/, $HotspotLine)[5];
		# Set the hotspot size using the data from the file
		$HotLength = $HotEnd - $HotStart + 1;
		# Use this information to now iterate through the gff3 file
		foreach my $GffLine (@$SubGffArray) {
			# Set the Gff3 line contig ID
			$GffContigID = (split /\t/, $GffLine)[0];
			# Also get the ID of the gene containing the hotspot
			$GffOrfID = (split /\t/, $GffLine)[5];
			# Skip if the ID does not match (I know, I know, a hash would be better here)
			next if ($HotID != $GffContigID);
			# Get the coordinates for the Gff3 gene line
			$GffStart = (split /\t/, $GffLine)[3];
			$GffEnd = (split /\t/, $GffLine)[4];
			# Set coordinates if forward gene
			if ($GffStart < $GffEnd) {
				# Skip if the hotspot is not found within the gene
				next unless ($GffStart <= $HotStart && $GffEnd >= $HotEnd);
				# Get the hotspot window as in frame codons
				# In frame codon position in start
				$WindowStart = $HotStart - (($HotStart - $GffStart) % 3);
				# In frame codon position in end
				$WindowEnd = $HotEnd + (2 - (($HotEnd - $GffStart) % 3));
				# Get the length of the corrected window
				$WindowLength = $WindowEnd - $WindowStart;
				$SurroundWindowLength = $WindowLength * $windowFactor;
				# Now that we have the hotspot window coords, we need to build the
				# adjacent sliding windows
				# Start by adding the window before the hot window
				$WindowSlideBeforeEnd = $WindowStart - 1;
				$WindowSlideBeforeStart =  $WindowSlideBeforeEnd - $SurroundWindowLength + 1;
				# Set sliding window following hot window
				$WindowSlideAfterStart = $WindowEnd + 1;
				$WindowSlideAfterEnd = $WindowSlideAfterStart + $SurroundWindowLength - 1;
			} elsif ($GffStart > $GffEnd) {
				# Skip if the hotspot is not found within the gene
				next unless ($GffEnd <= $HotStart && $GffStart >= $HotEnd);
				# Set coords if the gene is reversed
				# Get the hotspot window as in frame codons
				# In frame codon position in start
				$WindowStart = $HotEnd + (($GffStart - $HotEnd) % 3);
				# In frame codon position in end
				$WindowEnd = $HotStart - (3 - (($HotStart - $GffStart) % 3));
				# Get the length of the corrected window
				$WindowLength = $WindowStart - $WindowEnd;
				$SurroundWindowLength = $WindowLength * $windowFactor;
				# Now that we have the hotspot window coords, we need to build the
				# adjacent sliding windows
				# Start by adding the window before the hot window
				$WindowSlideBeforeEnd = $WindowStart + 1;
				$WindowSlideBeforeStart =  $WindowSlideBeforeEnd + $SurroundWindowLength -1;
				# Set sliding window following hot window
				$WindowSlideAfterStart = $WindowEnd - 1;
				$WindowSlideAfterEnd = $WindowSlideAfterStart - $SurroundWindowLength + 1;
			} elsif ($GffStart == $GffEnd) {
				die "KILLED BECAUSE EQUAL GFF START AND END: $!\n";
			}
			# Now that we have the coordinates we can plug them into the gff3 file
			print OUT "$GffContigID\tSlidingWindowPosition\tGene\t$WindowSlideBeforeStart\t$WindowSlideBeforeEnd\t$GffOrfID"."_BeforeHotspot_$HotUnique\t.\t1\tID=BeforeHotspot;\n";
			print OUT "$GffContigID\tSlidingWindowPosition\tGene\t$WindowStart\t$WindowEnd\t$GffOrfID"."_HotspotWindow_$HotUnique\t.\t1\tID=HotspotWindow;\n";
			print OUT "$GffContigID\tSlidingWindowPosition\tGene\t$WindowSlideAfterStart\t$WindowSlideAfterEnd\t$GffOrfID"."_AfterHotspot_$HotUnique\t.\t1\tID=AfterHotspot;\n";
		}
	}
}

# Now run the subroutines
ReadInGff3();
# print join("\n", @gffArray);
ReadInHotspots();
# print join("\n", @hotspotArray);
GenerateHotspotGff3(\@gffArray, \@hotspotArray);

close(GFF3);
close(HOTSPOTS);
close(OUT);

# Get the toal time to run the script
my $end_run = time();
my $run_time = $end_run - $start_run;
print STDERR "\nThe gff3 sliding window files were calculated in $run_time seconds.\n";




