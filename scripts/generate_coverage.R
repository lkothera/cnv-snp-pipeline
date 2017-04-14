library(Sushi)

# get command line arguments
args <- commandArgs(trailingOnly=TRUE)
gene_coords_file <- args[1]
sample_names_file <- args[2]
wgx_path <- args[3]
coverage_path <- args[4]

# read gene information
chroms <- read.table(file=gene_coords_file, header=FALSE, sep="\t")
factors <- sapply(chroms, is.factor)
chroms[factors] <- lapply(chroms[factors], as.character)

# read all sample names
sample_names <- scan(file=sample_names_file, character());

# read all coverage files
d <- vector(mode='list', length=length(sample_names))
for (i in 1:length(sample_names)) {
  print(sample_names[i])
  data <- read.table(paste(wgx_path,sample_names[i],'.wgx',sep=''), header=FALSE, sep="\t")
  d[[i]] <- data
}

# for each gene
for (j in 1:nrow(chroms)) {
  # create png
  png(file=paste(coverage_path,chroms[j,"V1"],'.png',sep=''), width=11, height=8.5, units="in", res=150)
  par(mfrow=c(11,4), mar=c(1.75,3,0.82,0.42))

  # for each sample
  for (i in 1:length(sample_names)) {
    plotBedgraph(d[[i]], 'chr1', chroms[j,"V2"], chroms[j,"V3"], color=SushiColors(2)(2)[1])
    axis(side=2, las=2, tcl=0.2)
    labelgenome(sample_names[i], chroms[j,"V2"], chroms[j,"V3"], side=1, n=3, chromcex=0.75, scale="Kb", scalecex=0.75)
  }

  dev.off()
}
