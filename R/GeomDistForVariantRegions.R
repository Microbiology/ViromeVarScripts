# GeomDistForVariantRegions.R
# Geoffrey Hannigan
# Elizabeth Grice Lab
# University of Pennsylvania

# Read in the command line args
args <- commandArgs(trailingOnly=TRUE)
input = read.delim(args[1],sep='\t',header=TRUE)
output = args[2]

GetGeomDistStats <- function(inputDF) {
  # Set empty data frame to append data into
  dataFrame = c()
  # Run loop across each line of the input file
  for (counter in 1:nrow(inputDF)) {
    line <- 1/input[counter,2]
    ContigID <- as.character(input[counter,1])
    # The p-value can be altered to customize predictions
    geomResult <- qgeom(prob=line, p=c(0.05))
    dataFrameLine <- c(ContigID, geomResult)
    dataFrame <- rbind(dataFrame,dataFrameLine)
  }
  return(dataFrame)
}

SubOutput <- as.data.frame(GetGeomDistStats(input))
colnames(SubOutput) <- c("ContigID","QuantileFunction")

# Send data frame to specified output file
write.table(x=SubOutput, file=output, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)




