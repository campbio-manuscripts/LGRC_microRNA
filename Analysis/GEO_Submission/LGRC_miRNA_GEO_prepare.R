
# Set up directory to transfer
dir.create("LGRC_microRNA", showWarnings = FALSE)
dir.create("LGRC_microRNA/fastq", showWarnings = FALSE)
dir.create("LGRC_microRNA/matrices", showWarnings = FALSE)

# Read in raw counts and normalized data
mir.counts.full = read.table("/protected/projects/lgrc/small_RNA/Final_Analysis/Data/MicroRNA/130720_LGRC_miRNA_Counts.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
mir.annotation.full = mir.counts.full[,1:7]
mir.counts.full = as.matrix(mir.counts.full[,-(1:7)])
rownames(mir.counts.full) = mir.annotation.full[,5]
rownames(mir.annotation.full) = mir.annotation.full[,5]
colnames(mir.annotation.full) = c("Chrom", "Start", "Stop", "Name", "Mature_ID", "Strand", "miRBase_Version")

library(affy)
data <- readRDS("../../Data/MicroRNA/LGRC_MicroRNA_Counts.rds")
pdata <- pData(data)
mir.counts.qc <- mir.counts.full[,colnames(data)]
mir.counts.with.annot <- data.frame(mir.annotation.full, mir.counts.qc)

mir.rpm = sweep(mir.counts.qc, 2, pdata$Mirna_Reads / 1e6, "/")
mir.rpm.log2 = log2(mir.rpm+1)
mir.rpm.with.annot <- data.frame(mir.annotation.full, mir.rpm.log2)

# Write tables/matrices
write.table(mir.counts.with.annot, "LGRC_microRNA/matrices/microRNA_raw_counts.txt", quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
write.table(mir.rpm.with.annot, "LGRC_microRNA/matrices/microRNA_normalized.txt", quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

# Read table and copy/rename fastq files
fastq.input <- read.table("../../../Unprotected_lgrc_365/Alignment_Pipeline/130720_param_input.txt", header=FALSE, sep="\t")
fastq.input <- fastq.input[fastq.input[,2] %in% rownames(pdata),]
new.name <- paste0("LGRC_microRNA/fastq/", fastq.input[,2], ".fq.gz")
file.copy(fastq.input[,1], new.name, overwrite = TRUE)

# Get md5 checksum
library(tools)
md5 <- md5sum(new.name)
write.table(cbind(basename(new.name), md5), "md5sum.txt", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

# Create metadata table
demographics = read.table("../../Annotation/Sample_Demographics/111228_LGRC_MicroRNA_Demographics.txt", sep="\t", header=TRUE, quote="", row.names=1, check.names=FALSE)
demo <- demographics[colnames(data),]

colnames(pdata) <- tolower(colnames(pdata))
library_name <- rownames(pdata)
status <- c("Control", "ILD", "COPD")[pdata$disease + 1]
pdata$disease <- status
pdata$copd <- NULL
pdata$control <- NULL
pdata$ild <- NULL
pdata$flow_cell <- as.numeric(as.factor(pdata$flow_cell))

pdata$cancer_status <- demo$f010qI.1K..Lung.Cancer
pdata$smoking_status <- demo$Cigarette.smoking
pdata$smoking_status[grep("^2-Ever", pdata$smoking_status)] <- "2-Ever"
pdata$gender <- demo$Gender
pdata$race_ethnicity <- demo$Race.or.Ethnicity
pdata$cancer <- NULL
pdata$smoke <- NULL
pdata$race <- NULL
pdata$percent_emphysema <- pdata$emphysema
pdata$emphysema <- NULL

title <- paste0(library_name, "_", status)
org <- "Homo sapiens"
tissue <- "lung"
ltrc_id <- rownames(pdata)
raw_file <- basename(new.name)
processed_data_file1 <- "microRNA_raw_counts.txt"
processed_data_file2 <- "microRNA_normalized.txt"
instrument <- ifelse(pdata$protocol == "Singleplex", "Illumina Genome Analyzer IIx", "Illumina HiSeq 2000")

metadata <- data.frame("library name" = library_name, title = title, organism = org, tissue = tissue, pdata, LTRC_ID = ltrc_id, "single or paired-end" = "single", molecule = "Total RNA", "instrument model" = instrument, "processed data file" = processed_data_file1, "processed data file" = processed_data_file2, "raw file" = raw_file, check.names = FALSE)
write.table(metadata, "metadata_table.txt", col.names = TRUE, sep="\t", quote=FALSE, row.names=FALSE)

sessionInfo()
