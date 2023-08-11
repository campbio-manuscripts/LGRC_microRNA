## Process targetscan results into a matrix

## Read in genes
UTRs = read.table(gzfile("UTR_Sequences_human_for_ts.txt.gz"), sep="\t", header=FALSE, stringsAsFactors=FALSE)
my.gene = unique(UTRs[UTRs[,2] == 9606,1])

## Read in microRNAs
seeds = read.table("IsomiR_Seeds_Filtered_For_TS.txt", sep="\t", stringsAsFactors=FALSE)
my.mir = unique(seeds[,1])

## Set up full matrix
target.matrix = matrix(0, nrow=length(my.gene), ncol=length(my.mir))
colnames(target.matrix) = my.mir
rownames(target.matrix) = my.gene

## Parse targetscan results file, one line at a time, and count the number of seeds of each miRNA in each UTR

cc = rep(c("character", "NULL"), c(2, 11))
ts.results = read.table("IsomiR_seed_hits.txt", header=TRUE, sep="\t", colClasses=cc, as.is=TRUE, check.names=FALSE, fill=TRUE)
for(i in 1:nrow(ts.results)) {
  target.matrix[ts.results[i,1], ts.results[i,2]] = target.matrix[ts.results[i,1], ts.results[i,2]] + 1
}


## Intersect genes miRs with those in mRNA and miRNA data
library(affy)
mRNA.data = readRDS("../../Data/Gene/LGRC_Gene_Overlap_MicroRNA.rds")
mRNA.annot = fData(mRNA.data)

miRNA.data = readRDS("../../Data/MicroRNA/LGRC_MicroRNA_Isoform_RPM_Log2_Combat.rds")
miRNA.annot = fData(miRNA.data)


refseq.match = match(mRNA.annot$REFSEQ, rownames(target.matrix))
mir.match = match(rownames(miRNA.annot), colnames(target.matrix))
target.matrix.matched = target.matrix[refseq.match, mir.match]

target.matrix.matched[is.na(target.matrix.matched)] = 0
rownames(target.matrix.matched) = rownames(mRNA.annot)

saveRDS(target.matrix, "130826_IsomiR_Targets_All_Refseq.rds")
saveRDS(target.matrix.matched, "130826_IsomiR_Targets_Matched.rds")


