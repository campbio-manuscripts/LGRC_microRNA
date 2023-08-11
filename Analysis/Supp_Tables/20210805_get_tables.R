# Gets the data and formats some tables for the supplement

qc <- read.table("../../Data/MicroRNA/QC_table.txt", header = TRUE)
qc$Exclude_From_Analysis <- ifelse(qc$Align_Read_Len_JS_Cluster %in% c(2,4,5),"Exclude", "Keep")
write.table(qc, "Supp_Table_QC_Stats.txt", sep="\t", quote=FALSE)

library(SummarizedExperiment)
d <- readRDS("../Differential_Expression//LGRC_MicroRNA_Isoform_RPM_Log2_Batch_Corrected_With_Cluster.rds")
fd <- fData(d)
fd <- data.frame(ID = rownames(fd), fd)
write.table(fd, "Supp_Table_MicroRNA_Information_with_DE_and_Cluster.txt", sep="\t", quote=FALSE)

pd <- pData(d)
pd <- data.frame(ID = rownames(pd), pd)
write.table(pd, "Supp_Table_Sample_Information_with_Pheno_and_Cluster.txt", sep="\t", quote=FALSE)
