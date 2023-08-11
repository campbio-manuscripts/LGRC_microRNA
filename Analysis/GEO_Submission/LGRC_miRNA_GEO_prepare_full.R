
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

# Read demo and QC metrics (need to align LTRC ids by removing duplicate indicators)
ltrc.id <- colnames(mir.counts.full)
ltrc.id <- gsub("-2$", "", ltrc.id)
ltrc.id <- gsub("_HiSeq$", "", ltrc.id)

demographics = read.table("../../Annotation/Sample_Demographics/111228_LGRC_MicroRNA_Demographics.txt", sep="\t", header=TRUE, quote="", row.names=1, check.names=FALSE)
demo <- demographics[ltrc.id,]
qc <- read.table("../..//Data/MicroRNA/QC_table.txt", sep="\t", header=TRUE, row.names=1)
qc <- qc[colnames(mir.counts.full),]

# Create raw count and normalized matrix
mir.counts.with.annot <- data.frame(mir.annotation.full, mir.counts.full)
mir.rpm = sweep(mir.counts.full, 2, qc$Mirna_Reads / 1e6, "/")
mir.rpm.log2 = log2(mir.rpm+1)
mir.rpm.with.annot <- data.frame(mir.annotation.full, mir.rpm.log2)

# Write tables/matrices 
library(R.utils)
write.table(mir.counts.with.annot, "LGRC_microRNA/matrices/microRNA_raw_counts.txt", quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
write.table(mir.rpm.with.annot, "LGRC_microRNA/matrices/microRNA_normalized.txt", quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
gzip("LGRC_microRNA/matrices/microRNA_raw_counts.txt", overwrite = TRUE)
gzip("LGRC_microRNA/matrices/microRNA_normalized.txt", overwrite = TRUE)


# Read table and copy/rename fastq files
fastq.input <- read.table("../../../Unprotected_lgrc_365/Alignment_Pipeline/130720_param_input.txt", header=FALSE, sep="\t")
rownames(fastq.input) <- fastq.input[,2]
fastq.input <- fastq.input[colnames(mir.counts.full),]
new.name <- paste0("LGRC_microRNA/fastq/", fastq.input[,2], ".fq.gz")
file.copy(fastq.input[,1], new.name, overwrite = TRUE)

# Get md5 checksum
library(tools)
md5 <- md5sum(new.name)
write.table(cbind(basename(new.name), md5), "md5sum.txt", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

# Process demo and qc data to produce metadata file
disease.status = demo[,"Major.Diagnosis..Final.Clinical."]
status <- substring(disease.status, 3)
status[status == "COPD/Emphysema"] <- "COPD"
disease2 <- demo[,"diagmin"]
smoke.status = demo[,"Cigarette.smoking"]
smoke.status[grep("^2-Ever", smoke.status)] <- "2-Ever"
race = demo[,"Race.or.Ethnicity"]
cancer.status = demo[,"f010qI.1K..Lung.Cancer"]
gender = demo[,"Gender"]
pemphy = as.numeric(demo[,90])
pack.years = as.numeric(demo[,92])
age = as.numeric(demo[,94])
fev1.pp = as.numeric(demo[,98])
fev1.fvc = as.numeric(demo[,99])
bode = as.numeric(demo[,101])
dlco = as.numeric(demo[,53])
tissue.type <- demo[,"Type.of.tissue"]
tissue.type <- gsub('"', "", tissue.type)

demo.simple = data.frame(disease=disease.status, minor_diagnosis=disease2,
                         smoking_status=smoke.status, race=race, cancer_status=cancer.status, gender=gender, age=age,
						 percent_emphysema=pemphy, pack_years=pack.years, FEV1_percent_predicted=fev1.pp, FEV1_FVC_ratio=fev1.fvc,
						 BODE_score=bode, DLCO=dlco, tissue_type = tissue.type,
						 check.names=FALSE)
rownames(demo.simple) = rownames(demo)


# Add QC metrics
replicate <- ifelse(grepl("-2$|_HiSeq$", colnames(mir.counts.full)), 2, 1)

qc.fail <- qc$Align_Read_Len_JS_Cluster != 1 & qc$Align_Read_Len_JS_Cluster != 3
included <- ifelse(!qc.fail & replicate == 1, "yes", "no")

flow.cell <- as.integer(as.factor(qc$Flow_Cell))
lane <- qc$Lane
index <- as.numeric(as.factor(qc$Index))

mir.reads <- qc$Mirna_Reads
tot.reads <- qc$Total_Reads
align.reads <- qc$Reads_aligned



qc_final <- data.frame(flow_cell = flow.cell, lane = lane, index = index,
            total_reads = tot.reads, total_reads_aligned = align.reads, total_reads_aligned_mirna = mir.reads,
            replicate = replicate, qc_pass = !qc.fail, included = included)

# Other metadata fields
title <- paste0(ltrc.id, "_", status)
title[replicate == 2] <- paste0(title[replicate == 2], "_Replicate")
org <- "Homo sapiens"
tissue <- "lung"

raw_file <- basename(new.name)
processed_data_file1 <- "microRNA_raw_counts.txt.gz"
processed_data_file2 <- "microRNA_normalized.txt.gz"

instrument <- rep(c("Illumina Genome Analyzer IIx", "Illumina HiSeq 2000"), c(45, 326))
multiplex <- rep(c("Singleplex", "Multiplex"), c(51, 320))
protocol <- rep(c("Illumina_1.5", "Illumina_Truseq"), c(51, 320))

metadata <- data.frame("library name" = rownames(qc), title = title, organism = org, tissue = tissue, LTRC_ID = ltrc.id, demo.simple, qc_final, "single or paired-end" = "single", molecule = "Total RNA", "instrument model" = instrument, protocol = protocol, multiplex = multiplex, "processed data file" = processed_data_file1, "processed data file" = processed_data_file2, "raw file" = raw_file, check.names = FALSE)
write.table(metadata, "metadata_table.txt", col.names = TRUE, sep="\t", quote=FALSE, row.names=FALSE)

sessionInfo()
