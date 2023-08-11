
# Get list of files and IDs
fs = list.files(path = "temp_fastq/", pattern=".fq", full.names=TRUE)
id = gsub(".fq", "", basename(fs))

# Read in sequences
data.seq <- list()
data.id <- list()
for(i in 1:length(fs)) {
	d = scan(fs[i], what="")
	s = d[seq(2,length(d), by=4)]

    data.id[[i]] = rep(id[i], length(s))
    data.seq[[i]] = s
}
data.id.all <- do.call(c, data.id)
data.seq.all <- do.call(c, data.seq)

# Create table and save to file
ta = table(data.seq.all, data.id.all)
class(ta) = "matrix"
write.table(data.frame(Sequence = rownames(ta), ta), "mir34c_sequence_counts.txt", sep="\t", quote=FALSE, row.names = FALSE)
