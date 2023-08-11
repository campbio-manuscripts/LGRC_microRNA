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

y <- c()
yn <- c()
w <- c()
wn <- c()

y[[1]] <- as.character(mir449b_5p_genes_pos)
yn[[1]] <- as.character(mir449b_5p_genes_neg)
w[[1]] <- as.character(mir449b_5p_genes_pos2)
wn[[1]] <- as.character(mir449b_5p_genes_neg2)

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

y[[2]] <- as.character(mir449a_genes_pos)
yn[[2]] <- as.character(mir449a_genes_neg)
w[[2]] <- as.character(mir449a_genes_pos2)
wn[[2]] <- as.character(mir449a_genes_neg2)

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

y[[3]] <- as.character(mir449c_5p_genes_pos)
yn[[3]] <- as.character(mir449c_5p_genes_neg)
w[[3]] <- as.character(mir449c_5p_genes_pos2)
wn[[3]] <- as.character(mir449c_5p_genes_neg2)
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

y[[4]] <- as.character(mir34c_3p_genes_pos)
yn[[4]] <- as.character(mir34c_3p_genes_neg)
w[[4]] <- as.character(mir34c_3p_genes_pos2)
wn[[4]] <- as.character(mir34c_3p_genes_neg2)
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

y[[5]] <- as.character(mir34b_3p_genes_pos)
yn[[5]] <- as.character(mir34b_3p_genes_neg)
w[[5]] <- as.character(mir34b_3p_genes_pos2)
wn[[5]] <- as.character(mir34b_3p_genes_neg2)
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

y[[6]] <- as.character(mir34c_5p_genes_pos)
yn[[6]] <- as.character(mir34c_5p_genes_neg)
w[[6]] <- as.character(mir34c_5p_genes_pos2)
wn[[6]] <- as.character(mir34c_5p_genes_neg2)
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

y[[7]] <- as.character(mir34b_5p_genes_pos)
yn[[7]] <- as.character(mir34b_5p_genes_neg)
w[[7]] <- as.character(mir34b_5p_genes_pos2)
wn[[7]] <- as.character(mir34b_5p_genes_neg2)
m7 <- union(as.character(mir34b_5p_genes_pos), as.character(mir34b_5p_genes_pos2))
m7n <- union(as.character(mir34b_5p_genes_neg), as.character(mir34b_5p_genes_neg2))

#genelist_global <- union(m7, union(m6, union(m5, union(m4,union(m1, union(m2, m3))))))
#genelist_global_n <- union(m7n, union(m6n, union(m5n, union(m4n,union(m1n, union(m2n, m3n))))))

genelist_global_copd <- union(y[[7]], union(y[[6]], union(y[[5]], union(y[[4]],union(y[[3]], union(y[[2]], y[[1]]))))))
genelist_global_copd_n <- union(yn[[7]], union(yn[[6]], union(yn[[5]], union(yn[[4]],union(yn[[3]], union(yn[[2]], yn[[1]]))))))

genelist_global_ild <- union(w[[7]], union(w[[6]], union(w[[5]], union(w[[4]],union(w[[3]], union(w[[2]], w[[1]]))))))
genelist_global_ild_n <- union(wn[[7]], union(wn[[6]], union(wn[[5]], union(wn[[4]],union(wn[[3]], union(wn[[2]], wn[[1]]))))))

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
#final_list <- intersect(rownames(exp_matrix), ann_gene[genelist_global, "GENE"])
#final_list <- intersect(rownames(exp_matrix), ann_gene[genelist_global_n, "GENE"])

final_list_copd <- intersect(rownames(exp_matrix), ann_gene[genelist_global_copd, "GENE"])
final_list_ild <- intersect(rownames(exp_matrix), ann_gene[genelist_global_ild, "GENE"])

#Run GSVA
library(GSVA)
final_list_gsva <- list(set1=final_list_copd)
mygsva <- gsva(exp_matrix, final_list_gsva, min.sz=10, max.sz=700, verbose=TRUE, rnaseq=FALSE)
gsva.val <- mygsva$es.obs
#write.table(gsva.val, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/gsva.val.COPD.pos.txt", sep="\t")

#write.table(gsva.val, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/gsva.val.ILD.pos.txt", sep="\t")

colnames(gsva.val) <- c( "P1_d0", "P1_d1", "P1_d4", "P1_d8", "P1_d10", "P1_d14", "P1_d21", "P1_d28", 
"P2_d0", "P2_d1", "P2_d2", "P2_d4", "P2_d8", "P2_d10", "P2_d12", "P2_d14", "P2_d17", "P2_d21", "P2_d28",
"P3_d0", "P3_d1", "P3_d2", "P3_d4", "P3_d8", "P3_d10", "P3_d12", "P3_d14", "P3_d17", "P3_d21", "P3_d28")
write.table(gsva.val, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/gsva.val.COPD.pos.txt", sep="\t")

days <- c(0,1,4,8,10,14,21,28,0,1,2,4,8,10,12,14,17,21,28,0,1,2,4,8,10,12,14,17,21,28)
#days <- c(c(1:8), c(1:11), c(1:11))

days <- as.numeric(days)
color <- c(rep("magenta", 8), rep("blue", 11), rep("green", 11))
pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/gsva_results_ross_abl.copd.pdf", height=6, width=6)
#jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/gsva_results_ross_abl.copd.jpg")
plot(days, as.numeric(gsva.val), xlab="Time of airway cells differentiation (days)", ylab="miR-449/34 gene set family", col=color, pch=19, cex.lab=1.3, cex.axis=1.3, main=("r=0.74, p<1e-05"), las=1, font.lab=2, font.axis=2)
#abline(lm(as.numeric(gsva.val)~days))
abline(lm(as.numeric(gsva.val[1:8])~days[1:8]), col="magenta", lwd=2)
abline(lm(as.numeric(gsva.val[9:19])~days[9:19]), col="blue", lwd=2)
abline(lm(as.numeric(gsva.val[20:30])~days[20:30]), col="green", lwd=2)
par(font=2)
legend("bottomright", legend=c("Patient 1: r=0.78, p=0.02", "Patient 2: r=0.88, p=0.003", "Patient 3: r=0.89, p=0.002"), col=c("magenta", "blue", "green"), cex=1.2, pch=19)
dev.off()


###Hogg airway
airway <- as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/rehoggdata/airway.exp.txt", sep="\t", stringsAsFactors=FALSE))
ann_gene_all = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/121212_LGRC_Gene_Array_Annotation2.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(ann_gene_all) <- ann_gene$NAME
rownames(airway) <- as.character(sapply(rownames(airway), function(x){(gsub("_at", "", x))}))
final_list_ild <- ann_gene_all[genelist_global_ild, "GENE"]
final_list_copd <- ann_gene_all[genelist_global_copd, "GENE"]

library(GSVA)
final_list_gsva <- list(set1=final_list_copd)
mygsva <- gsva(airway, final_list_gsva, min.sz=10, max.sz=700, verbose=TRUE, rnaseq=FALSE)
gsva.val <- mygsva$es.obs

#write.table(gsva.val, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/gsva.val.ILD.pos.hogg_airway.txt", sep="\t")
#write.table(gsva.val,"/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/gsva.val.COPD.pos.hogg_airway.txt", sep="\t")


airway.pheno <- as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/rehoggdata/airway.pheno.txt", sep="\t", stringsAsFactors=FALSE))

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/gsva_results_hogg_airway_abl.copd.pdf", height=6, width=6)
#jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/gsva_results_hogg_airway_abl.copd.jpg")
#plot(as.numeric(airway.pheno[, "airway.thickness"]), as.numeric(gsva.val), xlab="Airway thickness", ylab="GSVA scores", cex.names=0.5, col="orange", pch=19, cex.lab=1.3, cex.axis=1.3, main=("r=0.51, p<1e-04"))
plot(as.numeric(airway.pheno[, "airway.thickness"]), as.numeric(gsva.val), xlab="Airway thickness", ylab="miR-449/34 gene set family", cex.names=0.5, col="purple", pch=19, cex.lab=1.3, cex.axis=1.3, main=("r=0.50, p<1e-04"), las=1, font.lab=2, font.axis=2)
abline(lm(as.numeric(gsva.val)~as.numeric(airway.pheno[, "airway.thickness"])), col="purple", lwd=2)
par(font=2)
legend("bottomright", legend=c("r=0.5, p=5e-5"), cex=1.2)
dev.off()

#Hogg parenchyma
parenchyma <- as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/rehoggdata/parenchima.exp.txt", sep="\t", stringsAsFactors=FALSE))
ann_gene_all = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/121212_LGRC_Gene_Array_Annotation2.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(ann_gene_all) <- ann_gene$NAME
rownames(parenchyma) <- as.character(sapply(rownames(parenchyma), function(x){(gsub("_at", "", x))}))
final_list_ild <- ann_gene_all[genelist_global_ild, "GENE"]
final_list_copd <- ann_gene_all[genelist_global_copd, "GENE"]

final_list <- union (final_list_copd, final_list_ild)
library(GSVA)
final_list_gsva <- list(set1=final_list)
mygsva <- gsva(parenchyma, final_list_gsva, min.sz=10, max.sz=700, verbose=TRUE, rnaseq=FALSE)
gsva.val <- mygsva$es.obs

write.table(gsva.val, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/gsva.val.ILD.pos.hogg_parenchima.txt", sep="\t")
#write.table(gsva.val, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/gsva.val.COPD.pos.hogg_parenchima.txt", sep="\t")


parenchima.pheno <- as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/rehoggdata/parenchima.pheno.txt", sep="\t", stringsAsFactors=FALSE))

#pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/gsva_results_hogg_parenchima_abl.ild.pdf", height=6, width=6)
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/gsva_results_hogg_parenchima_abl.ild.jpg")
#plot(as.numeric(parenchima.pheno[, "lm"]), as.numeric(gsva.val), xlab="LM", ylab="GSVA scores", cex.names=0.5, col="magenta", pch=19, cex.lab=1.3, cex.axis=1.3, main=("r=0.15, p=0.23"))
plot(as.numeric(parenchima.pheno[, "lm"]), as.numeric(gsva.val), xlab="LM", ylab="GSVA scores", cex.names=0.5, col="orange", pch=19, cex.lab=1.3, cex.axis=1.3, main=("r=0.17, p=0.18"))

abline(lm(as.numeric(gsva.val)~as.numeric(parenchima.pheno[, "lm"])))
dev.off()


# pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/TargetsILD.COPD.Ross.global.pos.dots.patient1.pdf", height=6, width=6)
# #barplot(gsva.val[order(gsva.val)], cex.names=0.5)
# #barplot(gsva.val, names.arg=substr(colnames(gsva.val),4,9), cex.names=0.5)
# plot(gsva.val[1:8], xlab="Days", ylab="GSVA scores", cex.names=0.5, col="orange", pch=19, main="Patient 1")
# dev.off()

# pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/TargetsILD.COPD.Ross.global.pos.dots.patient2.pdf", height=6, width=6)
# #barplot(gsva.val[order(gsva.val)], cex.names=0.5)
# #barplot(gsva.val, names.arg=substr(colnames(gsva.val),4,9), cex.names=0.5)
# plot(gsva.val[9:19], xlab="Days", ylab="GSVA scores", cex.names=0.5, col="blue", pch=19, main="Patient 2")
# dev.off()

# pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Ross/TargetsILD.COPD.Ross.global.pos.dots.patient3.pdf", height=6, width=6)
# #barplot(gsva.val[order(gsva.val)], cex.names=0.5)
# #barplot(gsva.val, names.arg=substr(colnames(gsva.val),4,9), cex.names=0.5)
# plot(gsva.val[20:30], xlab="Days", ylab="GSVA scores", cex.names=0.5, col="green", pch=19, main="Patient 3")
# dev.off()
