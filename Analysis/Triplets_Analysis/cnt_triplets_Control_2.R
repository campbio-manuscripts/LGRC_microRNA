triplets <- read.table('Triplets_eQTL_with_CIT_Control.txt', header=TRUE, sep="\t")
nrow(triplets)

primary.key=paste(as.character(triplets[, "SNP"], as.character(triplets[, "Mir"]), as.character(triplets[, "Gene"])))
filtered.triplets=unique(triplets)
nrow(filtered.triplets)

filtered.triplets.cit <- filtered.triplets[which(filtered.triplets$Causal_Call==1),]
nrow(filtered.triplets.cit)
filtered.triplets.cit.ord <- order(filtered.triplets.cit$P_Causal)
#write.table(filtered.triplets.cit[filtered.triplets.cit.ord,], "triplets_cit_ord_unique_Control_new.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
#write.table(filtered.triplets, "triplets_all_ord_unique_Control_new.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

causal_triplets <- read.table('triplets_cit_ord_unique_Control_new.txt', header=TRUE, sep="\t")
#write.table(causal_triplets[order(causal_triplets$Mir_FDR),'Mir_FDR'], "miR_FDR.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
#write.table(causal_triplets[order(causal_triplets$Gene_FDR),'Gene_FDR'], "Gene_FDR.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
#write.table(causal_triplets[which(causal_triplets$Gene_FDR<0.1 & causal_triplets$Mir_FDR<0.1),], "Gene_miR_FDR0.1.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
#write.table(causal_triplets[which(causal_triplets$Gene_FDR<0.25 & causal_triplets$Mir_FDR<0.25),], "Gene_miR_FDR0.25.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

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
write.table(mir_t_ord_desc, "No_miRNA_Control_new.txt", sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)



triplets <- read.table('triplets_cit_ord_unique_Control_new.txt', header=TRUE, sep="\t")

#mir449b_5p
mir449b_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0003673_MIMAT0003327"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(mir449b_5p) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(mir449b_5p) <- mir449b_5p$'SNP'
write.table(mir449b_5p, "Control_SNP_mir449b_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#mir449b_5p
mir449b_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0003673_MIMAT0003327"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(mir449b_5p) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(mir449b_5p) <- mir449b_5p$'Gene'
write.table(mir449b_5p, "Control_Gene_mir449b_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


#miR449a
miR449a <- unique(triplets[which(triplets[,"Mir"]=="MI0001648_MIMAT0001541"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(miR449a) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(miR449a) <- miR449a$'SNP'
write.table(miR449a, "Control_SNP_miR449a.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR449a
miR449a <- unique(triplets[which(triplets[,"Mir"]=="MI0001648_MIMAT0001541"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(miR449a) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(miR449a) <- miR449a$'Gene'
write.table(miR449a, "Control_Gene_miR449a.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


#miR449c_5p
miR449c_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0003823_MIMAT0010251"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(miR449c_5p) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(miR449c_5p) <- miR449c_5p$'SNP'
write.table(miR449c_5p, "Control_SNP_miR449c_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR449c_5p
miR449c_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0003823_MIMAT0010251"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(miR449c_5p) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(miR449c_5p) <- miR449c_5p$'Gene'
write.table(miR449c_5p, "Control_Gene_miR449c_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR34c_3p
miR34c_3p <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0004677"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(miR34c_3p) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(miR34c_3p) <- miR34c_3p$'SNP'
write.table(miR34c_3p, "Control_SNP_miR34c_3p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR34c_3p
miR34c_3p <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0004677"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(miR34c_3p) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(miR34c_3p) <- miR34c_3p$'Gene'
write.table(miR34c_3p, "Control_Gene_miR34c_3p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR34b_5p
miR34b_3p <- unique(triplets[which(triplets[,"Mir"]=="MI0000742_MIMAT0000685"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(miR34b_3p) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(miR34b_3p) <- miR34b_3p$'SNP'
write.table(miR34b_3p, "Control_SNP_miR34b_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR34b_5p
miR34b_3p <- unique(triplets[which(triplets[,"Mir"]=="MI0000742_MIMAT0000685"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(miR34b_3p) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(miR34b_3p) <- miR34b_3p$'Gene'
write.table(miR34b_3p, "Control_Gene_miR34b_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR34c_5p
miR34c_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0000686"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(miR34c_5p) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(miR34c_5p) <- miR34c_5p$'SNP'
write.table(miR34c_5p, "Control_SNP_miR34c_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR34c_5p
miR34c_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0000686"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(miR34c_5p) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(miR34c_5p) <- miR34c_5p$'Gene'
write.table(miR34c_5p, "Control_Gene_miR34c_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
