mirs <- c("MI0003673_MIMAT0003327", "MI0001648_MIMAT0001541", "MI0003823_MIMAT0010251", "MI0000743_MIMAT0004677", "MI0000742_MIMAT0004676", "MI0000743_MIMAT0000686", "MI0000742_MIMAT0000685")

triplets <- read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/triplets_cit_ord_unique_COPD_new.txt', header=TRUE, sep="\t")
triplets2 <- read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_ILD/triplets_cit_ord_unique_ILD_new.txt', header=TRUE, sep="\t")

mir449b_5p_genes_pos <- unique(triplets[which(triplets[,"Mir"]=="MI0003673_MIMAT0003327" & triplets$"Mir_Gene_Cor">0), "Gene"])
as.character(mir449b_5p_genes_pos)  

mir449b_5p_genes_neg <- unique(triplets[which(triplets[,"Mir"]=="MI0003673_MIMAT0003327" & triplets$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir449b_5p_genes_neg)  

mir449b_5p_genes_pos2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0003673_MIMAT0003327" & triplets2$"Mir_Gene_Cor">0), "Gene"])
as.character(mir449b_5p_genes_pos2)  

mir449b_5p_genes_neg2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0003673_MIMAT0003327" & triplets2$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir449b_5p_genes_neg2)  

m1 <- union(as.character(mir449b_5p_genes_pos), as.character(mir449b_5p_genes_pos2))
m1n <- union(as.character(mir449b_5p_genes_neg), as.character(mir449b_5p_genes_neg2))


mir449a_genes_pos <- unique(triplets[which(triplets[,"Mir"]=="MI0001648_MIMAT0001541" & triplets$"Mir_Gene_Cor">0), "Gene"])
as.character(mir449a_genes_pos)  

mir449a_genes_neg <- unique(triplets[which(triplets[,"Mir"]=="MI0001648_MIMAT0001541" & triplets$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir449a_genes_neg)  

mir449a_genes_pos2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0001648_MIMAT0001541" & triplets2$"Mir_Gene_Cor">0), "Gene"])
as.character(mir449a_genes_pos2)  

mir449a_genes_neg2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0001648_MIMAT0001541" & triplets2$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir449a_genes_neg2)  

m2 <- union(as.character(mir449a_genes_pos), as.character(mir449a_genes_pos2))
m2n <- union(as.character(mir449a_genes_neg), as.character(mir449a_genes_neg2))



mir449c_5p_genes_pos <- unique(triplets[which(triplets[,"Mir"]=="MI0003823_MIMAT0010251" & triplets$"Mir_Gene_Cor">0), "Gene"])
as.character(mir449c_5p_genes_pos)  

mir449c_5p_genes_neg <- unique(triplets[which(triplets[,"Mir"]=="MI0003823_MIMAT0010251" & triplets$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir449c_5p_genes_neg)

mir449c_5p_genes_pos2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0003823_MIMAT0010251" & triplets2$"Mir_Gene_Cor">0), "Gene"])
as.character(mir449c_5p_genes_pos2)  

mir449c_5p_genes_neg2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0003823_MIMAT0010251" & triplets2$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir449c_5p_genes_neg2)  

m3 <- union(as.character(mir449c_5p_genes_pos), as.character(mir449c_5p_genes_pos2))
m3n <- union(as.character(mir449c_5p_genes_neg), as.character(mir449c_5p_genes_neg2))


mir34c_3p_genes_pos <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0004677" & triplets$"Mir_Gene_Cor">0), "Gene"])
as.character(mir34c_3p_genes_pos)  

mir34c_3p_genes_neg <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0004677" & triplets$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir34c_3p_genes_neg)

mir34c_3p_genes_pos2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0000743_MIMAT0004677" & triplets2$"Mir_Gene_Cor">0), "Gene"])
as.character(mir34c_3p_genes_pos2)  

mir34c_3p_genes_neg2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0000743_MIMAT0004677" & triplets2$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir34c_3p_genes_neg2)  

m4 <- union(as.character(mir34c_3p_genes_pos), as.character(mir34c_3p_genes_pos2))
m4n <- union(as.character(mir34c_3p_genes_neg), as.character(mir34c_3p_genes_neg2))

mir34b_3p_genes_pos <- unique(triplets[which(triplets[,"Mir"]=="MI0000742_MIMAT0004676" & triplets$"Mir_Gene_Cor">0), "Gene"])
as.character(mir34b_3p_genes_pos)  

mir34b_3p_genes_neg <- unique(triplets[which(triplets[,"Mir"]=="MI0000742_MIMAT0004676" & triplets$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir34b_3p_genes_neg)

mir34b_3p_genes_pos2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0000742_MIMAT0004676" & triplets2$"Mir_Gene_Cor">0), "Gene"])
as.character(mir34b_3p_genes_pos2)  

mir34b_3p_genes_neg2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0000742_MIMAT0004676" & triplets2$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir34b_3p_genes_neg2)  

m5 <- union(as.character(mir34b_3p_genes_pos), as.character(mir34b_3p_genes_pos2))
m5n <- union(as.character(mir34b_3p_genes_neg), as.character(mir34b_3p_genes_neg2))


mir34c_5p_genes_pos <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0000686" & triplets$"Mir_Gene_Cor">0), "Gene"])
as.character(mir34c_5p_genes_pos)  

mir34c_5p_genes_neg <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0000686" & triplets$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir34c_5p_genes_neg)

mir34c_5p_genes_pos2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0000743_MIMAT0000686" & triplets2$"Mir_Gene_Cor">0), "Gene"])
as.character(mir34c_5p_genes_pos2)  

mir34c_5p_genes_neg2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0000743_MIMAT0000686" & triplets2$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir34c_5p_genes_neg2)  

m6 <- union(as.character(mir34c_5p_genes_pos), as.character(mir34c_5p_genes_pos2))
m6n <- union(as.character(mir34c_5p_genes_neg), as.character(mir34c_5p_genes_neg2))

mir34b_5p_genes_pos <- unique(triplets[which(triplets[,"Mir"]=="MI0000742_MIMAT0000685" & triplets$"Mir_Gene_Cor">0), "Gene"])
as.character(mir34b_5p_genes_pos)  

mir34b_5p_genes_neg <- unique(triplets[which(triplets[,"Mir"]=="MI0000742_MIMAT0000685" & triplets$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir34b_5p_genes_neg)

mir34b_5p_genes_pos2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0000742_MIMAT0000685" & triplets2$"Mir_Gene_Cor">0), "Gene"])
as.character(mir34b_5p_genes_pos2)  

mir34b_5p_genes_neg2 <- unique(triplets2[which(triplets2[,"Mir"]=="MI0000742_MIMAT0000685" & triplets2$"Mir_Gene_Cor"<0), "Gene"])
as.character(mir34b_5p_genes_neg2)  

m7 <- union(as.character(mir34b_5p_genes_pos), as.character(mir34b_5p_genes_pos2))
m7n <- union(as.character(mir34b_5p_genes_neg), as.character(mir34b_5p_genes_neg2))

genelist_global <- union(m7, union(m6, union(m5, union(m4,union(m1, union(m2, m3))))))
genelist_global_n <- union(m7n, union(m6n, union(m5n, union(m4n,union(m1n, union(m2n, m3n))))))


#gene_list = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/COPD_Gene_miR449c_5p.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#gene_list = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/COPD_Gene_miR34c_5p.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#gene_list = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/COPD_Gene_miR34c_3p.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#gene_list = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/COPD_Gene_miR34b_3p.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#gene_list = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/COPD_Gene_miR449a.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#gene_list = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/COPD_Gene_miR449b_5p.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#gene_list = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_ILD/ILD_Gene_miR34b_5p.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

#gene_list2 = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_ILD/ILD_Gene_miR34b_3p.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#gene_list2 = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_ILD/ILD_Gene_miR34c_5p.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#gene_list2 = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_ILD/ILD_Gene_miR449a.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#gene_list2 = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_ILD/ILD_Gene_miR449b_5p.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#gene_list2 = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_ILD/ILD_Gene_miR449c_5p.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

ann_gene = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/121212_LGRC_Gene_Array_Annotation2.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(ann_gene) <- ann_gene$NAME
#ann_gene[gene_list$Gene, "GENE"]

exp_matrix <- as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/Ross_GSE5264_entrezgcdf_v16_noat.txt", sep="\t", stringsAsFactors=FALSE))
#final_list <- intersect(rownames(exp_matrix), ann_gene[gene_list$Gene, "GENE"])
#final_list <- intersect(rownames(exp_matrix), union(ann_gene[gene_list$Gene, "GENE"], ann_gene[gene_list2$Gene, "GENE"]))
final_list <- intersect(rownames(exp_matrix), ann_gene[genelist_global, "GENE"])
#final_list <- intersect(rownames(exp_matrix), ann_gene[genelist_global_n, "GENE"])

#Run GSVA
library(GSVA)
final_list_gsva <- list(set1=final_list)
mygsva <- gsva(exp_matrix, final_list_gsva, min.sz=10, max.sz=700, verbose=TRUE, rnaseq=FALSE)
gsva.val <- mygsva$es.obs
#write.table(gsva.val, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/gsva.val.ILD.COPD.global.neg.txt", sep="\t")

colnames(gsva.val) <- c( "P1_d0", "P1_d1", "P1_d4", "P1_d8", "P1_d10", "P1_d14", "P1_d21", "P1_d28", 
"P2_d0", "P2_d1", "P2_d2", "P2_d4", "P2_d8", "P2_d10", "P2_d12", "P2_d14", "P2_d17", "P2_d21", "P2_d28",
"P3_d0", "P3_d1", "P3_d2", "P3_d4", "P3_d8", "P3_d10", "P3_d12", "P3_d14", "P3_d17", "P3_d21", "P3_d28")
pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/TargetsILD.COPD.Ross.global.pos.dots.patient1.pdf", height=6, width=6)
#barplot(gsva.val[order(gsva.val)], cex.names=0.5)
#barplot(gsva.val, names.arg=substr(colnames(gsva.val),4,9), cex.names=0.5)
plot(gsva.val[1:8], xlab="Days", ylab="GSVA scores", cex.names=0.5, col="orange", pch=19, main="Patient 1")
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/TargetsILD.COPD.Ross.global.pos.dots.patient2.pdf", height=6, width=6)
#barplot(gsva.val[order(gsva.val)], cex.names=0.5)
#barplot(gsva.val, names.arg=substr(colnames(gsva.val),4,9), cex.names=0.5)
plot(gsva.val[9:19], xlab="Days", ylab="GSVA scores", cex.names=0.5, col="blue", pch=19, main="Patient 2")
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/TargetsILD.COPD.Ross.global.pos.dots.patient3.pdf", height=6, width=6)
#barplot(gsva.val[order(gsva.val)], cex.names=0.5)
#barplot(gsva.val, names.arg=substr(colnames(gsva.val),4,9), cex.names=0.5)
plot(gsva.val[20:30], xlab="Days", ylab="GSVA scores", cex.names=0.5, col="green", pch=19, main="Patient 3")
dev.off()
