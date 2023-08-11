triplets <- read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Preprocess/eQTL_262_output_PC_COPD/Triplets_eQTL_with_CIT_COPD.txt', header=TRUE, sep="\t")
nrow(triplets)

primary.key=paste(as.character(triplets[, "SNP"], as.character(triplets[, "Mir"]), as.character(triplets[, "Gene"])))
filtered.triplets=unique(triplets)
nrow(filtered.triplets)

filtered.triplets.cit <- filtered.triplets[which(filtered.triplets$Causal_Call==1),]
nrow(filtered.triplets.cit)
filtered.triplets.cit.ord <- order(filtered.triplets.cit$P_Causal)

primary.key=paste(as.character(filtered.triplets.cit[, "Mir"]), as.character(filtered.triplets.cit[, "Gene"]))
mir_gene = subset(filtered.triplets.cit, !duplicated(primary.key))
nrow(mir_gene)

mirs_unique <- unique(as.character(mir_gene$Mir))
mir_t <- matrix(nrow=length(mirs_unique), ncol=1)
colnames(mir_t) <- c("Freq")
rownames(mir_t) <- mirs_unique

for (i in mirs_unique){
    nr <- 0
    for (j in 1:nrow(mir_gene)){
        if (as.character(mir_gene[j,"Mir"])==i){
            nr <- nr + 1
        }
     #mir_t[i, "Freq"] <- nr/nrow(mir_gene)   
     mir_t[i, "Freq"] <- nr        
    }
}
mir_t_ord <- order(-mir_t[,1])
mir_t_ord_desc <- mir_t[mir_t_ord,1]


###PLOT DC mirs COPD
freq_COPD <- read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Preprocess/eQTL_262_output_PC_COPD/No_miRNA_COPD_new.txt", header=FALSE, sep="\t")

colnames(freq_COPD) <- c("Mir_COPD", "Freq_COPD")
rownames(freq_COPD) <- freq_COPD$Mir_COPD

freq_Control <- read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Preprocess/eQTL_262_output_PC_COPD/No_miRNA_Control_new.txt", header=FALSE, sep="\t")
#freq_Control <- read.table("Freq_miRNA_Control_new.txt", header=FALSE, sep="\t")

colnames(freq_Control) <- c("Mir_Control", "Freq_Control")
rownames(freq_Control) <- freq_Control$Mir_Control

length(intersect(rownames(freq_Control), rownames(freq_COPD)))
length(union(rownames(freq_Control), rownames(freq_COPD)))

comb <- matrix(0, nrow=length(union(rownames(freq_Control), rownames(freq_COPD))), ncol=2)
colnames(comb) <- c("Freq_Control", "Freq_COPD")
rownames(comb) <- union(rownames(freq_Control), rownames(freq_COPD))

comb[rownames(freq_Control),"Freq_Control"] <- freq_Control$Freq_Control
comb[rownames(freq_COPD),"Freq_COPD"] <- freq_COPD$Freq_COPD


mirnames <- matrix("", ncol = 1, nrow = nrow(comb))
rownames(mirnames) <- rownames(comb)
colnames(mirnames) <- "name"
 # mirnames["MI0003673_MIMAT0003327", "name"] <- "hsa-miR-449b-5p" 
 # mirnames["MI0001648_MIMAT0001541", "name"] <- "hsa-miR-449a" 
 # mirnames["MI0003823_MIMAT0010251", "name"] <- "hsa-miR-449c-5p" 
 # mirnames["MI0000743_MIMAT0004677", "name"] <- "hsa-miR-34c-3p" 
 # mirnames["MI0000742_MIMAT0004676", "name"] <- "hsa-miR-34b-3p" 
 # #mirnames["MI0000742_MIMAT0000685", "name"] <- "hsa-miR-34b-5p" 
# mirnames["MI0000743_MIMAT0000686", "name"] <- "hsa-miR-34c-5p" 
# mirnames["MI0016760_MIMAT0019232", "name"] <- "hsa-miR-4423-5p"

mircol <- matrix("blue", ncol = 1, nrow = nrow(comb))
rownames(mircol) <- rownames(comb)
colnames(mircol) <- "name"
# mircol["MI0003673_MIMAT0003327", "name"] <- "red" 
# mircol["MI0001648_MIMAT0001541", "name"] <- "red" 
# mircol["MI0003823_MIMAT0010251", "name"] <- "red" 
# mircol["MI0000743_MIMAT0004677", "name"] <- "red" 
# mircol["MI0000742_MIMAT0004676", "name"] <- "red" 
# #mircol["MI0000742_MIMAT0000685", "name"] <- "red" 
# mircol["MI0000743_MIMAT0000686", "name"] <- "red" 


sign_mirs <- read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Preprocess/eQTL_262_output_PC_COPD/Freq_miRNA_COPD_new_fisher1_corr_sortedByFreqCOPD.txt", header=TRUE, sep="\t")

sign_mirs2 <- sign_mirs[rownames(mircol) ,]
mircol[which(sign_mirs2[,"FDR"]<0.2), "name"] <- "red" 


pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Preprocess/eQTL_262_output_PC_COPD/No_COPD_col_paper_fin.pdf")
plot(comb[,"Freq_Control"], comb[,"Freq_COPD"], xlab="Number of connected genes in Control", ylab="Number of connected genes in COPD", pch=19, col=mircol[,"name"], cex.axis=1.3, cex.main=1.3, cex.lab=1.3, las=1, font.lab=2, cex.names=0.5, cex=1.5)
text(comb[,"Freq_Control"], comb[,"Freq_COPD"], labels=mirnames[,"name"], cex= 0.7, pos=4)
dev.off()






###PLOT DC mirs ILD
freq_COPD <- read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Preprocess/eQTL_262_output_PC_COPD/No_miRNA_ILD_new.txt", header=FALSE, sep="\t")

colnames(freq_COPD) <- c("Mir_COPD", "Freq_COPD")
rownames(freq_COPD) <- freq_COPD$Mir_COPD

freq_Control <- read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Preprocess/eQTL_262_output_PC_COPD/No_miRNA_Control_new.txt", header=FALSE, sep="\t")
#freq_Control <- read.table("Freq_miRNA_Control_new.txt", header=FALSE, sep="\t")

colnames(freq_Control) <- c("Mir_Control", "Freq_Control")
rownames(freq_Control) <- freq_Control$Mir_Control

length(intersect(rownames(freq_Control), rownames(freq_COPD)))
length(union(rownames(freq_Control), rownames(freq_COPD)))

comb <- matrix(0, nrow=length(union(rownames(freq_Control), rownames(freq_COPD))), ncol=2)
colnames(comb) <- c("Freq_Control", "Freq_COPD")
rownames(comb) <- union(rownames(freq_Control), rownames(freq_COPD))

comb[rownames(freq_Control),"Freq_Control"] <- freq_Control$Freq_Control
comb[rownames(freq_COPD),"Freq_COPD"] <- freq_COPD$Freq_COPD


mirnames <- matrix("", ncol = 1, nrow = nrow(comb))
rownames(mirnames) <- rownames(comb)
colnames(mirnames) <- "name"
# mirnames["MI0003673_MIMAT0003327", "name"] <- "hsa-miR-449b-5p" 
# mirnames["MI0001648_MIMAT0001541", "name"] <- "hsa-miR-449a" 
# mirnames["MI0003823_MIMAT0010251", "name"] <- "hsa-miR-449c-5p" 
# mirnames["MI0000743_MIMAT0004677", "name"] <- "hsa-miR-34c-3p" 
# mirnames["MI0000742_MIMAT0004676", "name"] <- "hsa-miR-34b-3p" 
# #mirnames["MI0000742_MIMAT0000685", "name"] <- "hsa-miR-34b-5p" 
# mirnames["MI0000743_MIMAT0000686", "name"] <- "hsa-miR-34c-5p" 
# mirnames["MI0016760_MIMAT0019232", "name"] <- "hsa-miR-4423-5p"

mircol <- matrix("blue", ncol = 1, nrow = nrow(comb))
rownames(mircol) <- rownames(comb)
colnames(mircol) <- "name"
# mircol["MI0003673_MIMAT0003327", "name"] <- "red" 
# mircol["MI0001648_MIMAT0001541", "name"] <- "red" 
# mircol["MI0003823_MIMAT0010251", "name"] <- "red" 
# mircol["MI0000743_MIMAT0004677", "name"] <- "red" 
# mircol["MI0000742_MIMAT0004676", "name"] <- "red" 
# #mircol["MI0000742_MIMAT0000685", "name"] <- "red" 
# mircol["MI0000743_MIMAT0000686", "name"] <- "red" 


sign_mirs <- read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Preprocess/eQTL_262_output_PC_ILD/Freq_miRNA_ILD_new_fisher1_corr_sortedByFreqILD.txt", header=TRUE, sep="\t")

sign_mirs2 <- sign_mirs[rownames(mircol) ,]
mircol[which(sign_mirs2[,"FDR"]<0.2), "name"] <- "red" 


pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Preprocess/eQTL_262_output_PC_COPD/No_ILD_col_paper_fin.pdf")
plot(comb[,"Freq_Control"], comb[,"Freq_COPD"], xlab="Number of connected genes in Control", ylab="Number of connected genes in ILD", pch=19, col=mircol[,"name"], cex.axis=1.3, cex.main=1.3, cex.lab=1.3, las=1, font.lab=2, cex.names=0.5, cex=1.5)
text(comb[,"Freq_Control"], comb[,"Freq_COPD"], labels=mirnames[,"name"], cex= 0.7, pos=4)
dev.off()
