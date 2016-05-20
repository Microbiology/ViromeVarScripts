Evolutionary and Functional Implications of Hypervariable Loci Within 
the Skin Virome

Hannigan, et al.



README file for the virome variability analysis workflow.
This file explains the analysis workflow.

NOTE: In this analysis, we did a primary and validation analysis.
This describes the workflow for the primary analysis. The same
scripts were used for the validation dataset from the Segre
dataset. These scripts are designated as starting with 'Segre'.

=====================================
=           SHELL SCRIPTS           =
=====================================

##########################
# Align Reads to Contigs #
##########################
-                              -
- AlignViromeReadsToContigs.sh -
-                              -

This script will:

+ Combine all of the negative control cleaned files
of the sites and time of interest, and place them into a single
fasta file for bowtie alignment

+ Predict the ORFs present in each contig

+ Call SNPs within each contig

+ Predict evolutionary pressure on each contig


|
V

##########################
# Calculate Contig Stats #
##########################
-                                      -
- ContigStatsAndReferenceGenomeHits.sh -
-                                      -

Use this script to calculate contig length, the number of
sequences mapping to each contig, and the taxonomic ID
of the contig using simple blast similarity to know virus
reference genomes.


|
V

##############################################
# Identify SNP Hotspots Using Geometric Dist #
##############################################
-                             -
- IdentifyHotspotsGeomDist.sh -
-                             -

Use this script to implement my perl and R scripts for
identifying SNP hotspots using the geometric distribution
of the nucleotide distances between SNPs on each contig.


|
V

######################################################
# Pull Contigs With Taxa IDs For Downstream Analysis #
######################################################
-                            -
- PullSpecificTaxaContigs.sh -
-                            -

Use this script to pull out the skin virome contigs
with specific taxonomic IDs of interest. The IDs of
interest can be determined using contig stats to see
what taxonomic IDs had highest coverage and greatest
identity to virus/phage references.


|
V

####################################################
# Determine Hotspot Evolution Pressure and Gene ID #
####################################################
-                        -
- SnpHotspotEvolution.sh -
-                        -

Use this script to calculate the evolutionary pressures
on the virome contigs by type (i.e. HPV, Staph phage).
Also use this script for annotating the ORFs and
determining what functionalities containt SNP hotspots,
and what their selective pressures are.


|
V

#################################################
# Determine Genomic Evo Pressure Using All SNPs #
#################################################
-                    -
- EvoByGeneByTaxa.sh -
-                    -

Use this script to calculate the evolutionary pressures
on the virome contigs by type (i.e. HPV, Staph phage)
across all of the genomic contigs using all of the SNPs.
This provides a "hotspot-free" view of overall pressure.


|
V

##########################################################
# Prepare Phylogeny References for Phylogenetic Analysis #
##########################################################
-                                  -
- PrepareReferencesForPhylogeny.sh -
-                                  -

This script is used to prepare the phylogenetic marker gene
reference dataset to be used in the next step (the actual
phylogentic analysis). The formatting is for removing illegal
characters and for making the names easier to read (important
for visualization step).


|
V

######################################################
# Calculate Phylogenetic Relationships Among Viruses #
######################################################
-                       -
- PhylogenyOfContigs.sh -
-                       -

Use this script to evaluate the phylogenetic
relationships among the virus contigs and the reference
genomes. The end result here is a phylogenetic tree
based on alignments of phylogenetic marker genes from
both the assembled contigs and the virus reference
genomes.


|
V

#####################################################
# Predict Phenotypic Impacat of SAVs on Virus Genes #
#####################################################
-                               -
- CalculateDeleteriousScores.sh -
-                               -

This allows us to evaluate the deleterious scores associated
with the hypervariable loci SAVs.


|
V

##############################################
# Predict Diversity Generating Retroelements #
##############################################
-                 -
- EvaluateDgrs.sh -
-                 -

Use this script to identify candidate DGR cassettes.


|
V

#################################
# Comparing Dataset to Oh et al #
#################################
-                                     -
- AlignCompareHanniganSegreScripts.sh -
-                                     -

Use this script to compare the contigs between the
Segre and Oh dataset and our dataset. This includes
aligning the assembled contigs between the datasets.


|
V

###################
# Quality Control #
###################
-                              -
- SegreQualityTrimmingAndQc.sh -
-                              -

Use this script to perform our standard quality
control processed on the validation dataset. This
needs to be done before further processing.


|
V


########################################
# Assemble Contigs Using Ray Assembler #
########################################
-                            -
- SegreAssembleContigsRay.sh -
-                            -

Use this script to compare the contigs between the
Segre and Oh dataset and our dataset. This includes
aligning the assembled contigs between the datasets.


|
V

##########################
# HVL ORF Fasta Creation #
##########################
-                             -
- CreateHpvOrfHotspotFasta.sh -
-                             -

Use this script to create a fasta file for only the
ORFs that contain HVLs. This will create a seperate
fasta for each virus type, as specified in the script.


|
V

##########################
# Pfam Domain Annotation #
##########################
-                              -
- PfamHmmerDomainAnnotation.sh -
-                              -

Use this script to annotate the hypervariable loci
containing ORFS using Pfam Hmmer. The result is a
table of annotations.


|
V

#############################################
# Create Deleterious Score Reference Matrix #
#############################################
-                          -
- GenerateDelScoreTable.sh -
-                          -

Use this script to generate a reference matrix of the
delterious scores associated with the HVL ORFS. This will
download the individual ORF references from the SuSPect web
browser interface, annotate the tables, and then combine
them.


|
V

########################################
# Calculate HVL ORF Deleterious Scores #
########################################
-                               -
- CalculateDeleteriousScores.sh -
-                               -

This script will calculate the delterious scores of the HVL
ORFs using the matrix generated in the previous step. The data
can be moved to analysis in R.



==========================================
=           R ANALYSIS SCRIPTS           =
==========================================

################
# Contig Stats #
################
--                                          --
--      ContigStatsAndRefGenomeMatch.R      --
--                                          --

-    ContigStatsAndReferenceGenomeHits.sh    - 

Use this script to plot the contig stats including
coverage, taxonomy, and similarity to reference
genomes.


|
V

######################################
# SNP Hotspot Evolution & Annotation #
######################################
--                                         --
--          SnpHotspotEvolution.R          --
--                                         --

Use this script to plot the evolutionary pressures
present on the SNP hotspots. This will compare the
pressures of the hotspots to the immediately
adjacent regions. It is also used to show the
hotspot pressures grouped by gene annotation (each
dot is a hotspot within that gene, showing
evolutionary pressure).


|
V

#########################################
# SNP Hotspot Evolution Between Viruses #
#########################################
--                                                  --
--          SnpHotspotCompareEvoAndCount.R          --
--                                                  --

Use this script to plot the evolutionary pressures
present on the SNP hotspots between viruses. In this
specific case, we are looking at the differences in
hotspot evo pressure between Staph phage and HPV.


|
V

#########################
# Substitution Patterns #
#########################
--                                         --
--      HotspotSubstitutionPatterns.R      --
--                                         --

Use this script to plot the substitution patterns
of the SNP hotspots within HPV and Staph phage
reference genomes.


|
V

##########################################
# Compare SNP Locations Between Datasets #
##########################################
--                                 --
--      CompareSnpLocations.R      --
--                                 --

Use this script to compare the locations of SNPs
between the primary and secondary datasets.


|
V

#############################
# Overall Genomic Evolution #
#############################
--                                    --
--      GenomicEvolutionByGene.R      --
--                                    --

Use this script to calculate overall genomic
evolution pressure of each virus, using a gy-gene approach.


|
V

###################################
# Evolution Pressure on Each Gene #
###################################
--                                        --
--      GenomicEvolutionOfEachGene.R      --
--                                        --

Use this script to calculate selective pressure on the genes
present in the virus genomes.


|
V

##########################################################
# Use Geom Distribution for Hypervariable Loci Detection #
##########################################################
--                                       --
--      GeomDistForVariantRegions.R      --
--                                       --

Use this script to calculate the geometric distributions of each
contig based on the distances between SNPs. This is used in
combination with the related perl scripts.


|
V

#####################################################
# Calculate Stats Associated with Delterious Scores #
#####################################################
--                           --
--      DelScoreStats.R      --
--                           --

Use this script to calculate the stats of the delterious scores
between viruses and to visualize their distribution.


|
V

####################################
# Calculate Linkage Disequilibrium #
####################################
--                    --
--      DGR-LD.R      --
--                    --

Use this script to calculate the linkage disequilibrium of
DGRs using the pileup matrix generated above with EvaluateDgrs.sh.




=============================================
=           PERL ANALYSIS SCRIPTS           =
=============================================

#########################
# Calculate pNpS Ratios #
#########################
--                                 --
--      CalculatePnPsRatio.pl      --
--                                 --

Use this script to calculate the pNpS ratios of the genes
within a contig.


####################################
# Extract Contig IDs From SAM File #
####################################
--                                     --
--      GetMappedReadIdFromSam.pl      --
--                                     --

Use this script to get the contig IDs from alignments
in Bowtie using the SAM output format.


##############################
# Remove Indels From Dataset #
##############################
--                                  --
--      RemoveIndelsFromVCF.pl      --
--                                  --

Use this script to remove indels from the dataset because we will not
be focusing on them in this analysis.


####################################
# First Step of Geom Dist For HVLs #
####################################
--                                      --
--      variabilityCallGeomDist.pl      --
--                                      --

Use this script for the first step in calculating hypervariable loci
using a geometric distribution.


#####################################
# Second Step of Geom Dist For HVLs #
#####################################
--                                       --
--      variabilityCallGeomDist2.pl      --
--                                       --

Use this script for the second step in calculating hypervariable loci
using a geometric distribution.


##################################
# Get Hypervariable Loci as Gff3 #
##################################
--                                       --
--      HotspotSlidingWindowGff3.pl      --
--                                       --

Use this to get the hypervariable loci as gff3 files.


###############################
# Calculate Delterious Scores #
###############################
--                                         --
--      CalculateDeleteriousScores.pl      --
--                                         --

This will let you calculate the delterious scores assocaited with your
genes of interest.

##################################
# Convert mpileup to genome file #
##################################
--                              --
--      mpilup2genotype.pl      --
--                              --

This script allows for conversion from .mpileup format to genotype, which
will be used in the R LD analysis.

#################################
# Get DGRs that are Within ORFs #
#################################
--                                   --
--      GetDgrContainedGenes.pl      --
--                                   --

This script allows us to filter out DGRs that are not found within
predicted ORFs.

##############################
# Get DGRs that Contain HVLs #
##############################
--                         --
--      GetDgrHvls.pl      --
--                         --

This script filters out DGRs that do not have at least one hypervariable
loci in the template/variable region pair.

###############################
# Get DGR positions in Contig #
###############################
--                          --
--      ExtractDgrs.pl      --
--                          --

This script gets the positions of the extracted DGRs.
