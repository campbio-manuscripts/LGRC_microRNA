triplets <- read.table('Triplets_eQTL_with_CIT_COPD.txt', header=TRUE, sep="\t")
nrow(triplets)

primary.key=paste(as.character(triplets[, "SNP"], as.character(triplets[, "Mir"]), as.character(triplets[, "Gene"])))
filtered.triplets=unique(triplets)
nrow(filtered.triplets)

filtered.triplets.cit <- filtered.triplets[which(filtered.triplets$Causal_Call==1),]
nrow(filtered.triplets.cit)
filtered.triplets.cit.ord <- order(filtered.triplets.cit$P_Causal)
write.table(filtered.triplets.cit[filtered.triplets.cit.ord,], "triplets_cit_ord_unique_COPD_new.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
#write.table(filtered.triplets, "triplets_all_ord_unique_COPD_new.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

causal_triplets <- read.table('triplets_cit_ord_unique_COPD_new.txt', header=TRUE, sep="\t")
#write.table(causal_triplets[order(causal_triplets$Mir_FDR),'Mir_FDR'], "miR_FDR.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
#write.table(causal_triplets[order(causal_triplets$Gene_FDR),'Gene_FDR'], "Gene_FDR.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(causal_triplets[which(causal_triplets$Gene_FDR<0.1 & causal_triplets$Mir_FDR<0.1),], "Gene_miR_FDR0.1.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(causal_triplets[which(causal_triplets$Gene_FDR<0.25 & causal_triplets$Mir_FDR<0.25),], "Gene_miR_FDR0.25.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

primary.key=paste(as.character(filtered.triplets.cit[, "Mir"]), as.character(filtered.triplets.cit[, "Gene"]))
mir_gene = subset(filtered.triplets.cit, !duplicated(primary.key))
nrow(mir_gene)

mirs_unique <- unique(as.character(mir_gene$Mir))
mir_t <- matrix(nrow=length(mirs_unique), ncol=1)
colnames(mir_t) <- c("Freq")
rownames(mir_t) <- mirs_unique

#length(unique(as.character(causal_triplets$SNP)))
#2928 COPD; 2945 ILD; 11 common

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
write.table(mir_t_ord_desc, "No_miRNA_COPD_new.txt", sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

freq_COPD <- read.table("No_miRNA_COPD_new.txt", header=FALSE, sep="\t")
#freq_COPD <- read.table("Freq_miRNA_COPD_new.txt", header=FALSE, sep="\t")

colnames(freq_COPD) <- c("Mir_COPD", "Freq_COPD")
rownames(freq_COPD) <- freq_COPD$Mir_COPD
freq_Control <- read.table("No_miRNA_Control_new.txt", header=FALSE, sep="\t")
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
mirnames["MI0003673_MIMAT0003327", "name"] <- "hsa-miR-449b-5p" 
mirnames["MI0001648_MIMAT0001541", "name"] <- "hsa-miR-449a" 
mirnames["MI0003823_MIMAT0010251", "name"] <- "hsa-miR-449c-5p" 
mirnames["MI0000743_MIMAT0004677", "name"] <- "hsa-miR-34c-3p" 
mirnames["MI0000742_MIMAT0004676", "name"] <- "hsa-miR-34b-3p" 
#mirnames["MI0000742_MIMAT0000685", "name"] <- "hsa-miR-34b-5p" 
mirnames["MI0000743_MIMAT0000686", "name"] <- "hsa-miR-34c-5p" 

mircol <- matrix("blue", ncol = 1, nrow = nrow(comb))
rownames(mircol) <- rownames(comb)
colnames(mircol) <- "name"
mircol["MI0003673_MIMAT0003327", "name"] <- "red" 
mircol["MI0001648_MIMAT0001541", "name"] <- "red" 
mircol["MI0003823_MIMAT0010251", "name"] <- "red" 
mircol["MI0000743_MIMAT0004677", "name"] <- "red" 
mircol["MI0000742_MIMAT0004676", "name"] <- "red" 
#mircol["MI0000742_MIMAT0000685", "name"] <- "red" 
mircol["MI0000743_MIMAT0000686", "name"] <- "red" 

jpeg("No_COPD_col.jpeg")
plot(comb[,"Freq_Control"], comb[,"Freq_COPD"], xlab="Control", ylab="COPD", pch=19, col=mircol[,"name"], cex.axis=1.3, cex.main=1.3, cex.lab=1.3, main="Number of predicted mRNA targets")
text(comb[,"Freq_Control"], comb[,"Freq_COPD"], labels=mirnames[,"name"], cex= 0.7, pos=4)
dev.off()



pdf("Barplot_Freq_COPD.pdf", height=15, width=50)
barplot(comb[order(-comb[,"Freq_COPD"]), "Freq_COPD"], col="blue", xaxt='n', cex.lab=2, cex.axis = 3)
dev.off()

comb2 <- comb[order(-comb[,"Freq_COPD"]),]
freq_names <- c("27a-5p", "190b", "449b-5p", "449a", "449c-5p", "4423-5p", "92b-3p", "34c-3p", "205-5p", "23a-5p", "509-3p-2", "509-3p-3","509-3p-1", "30a-3p", "34b-3p", "1185-1-3p", "125b-1-3p", "654-5p", "485-5p", "34c-5p")
pdf("Barplot_Top_Freq_COPD.pdf", height=15, width=50)
barplot(comb2[1:20, "Freq_COPD"],  col="blue", cex.lab=0.1, cex.axis = 3, names.arg=freq_names, cex.names=2.2)
dev.off()

pdf("Barplot_Freq_Control.pdf", height=15, width=40)
barplot(comb[order(-comb[,"Freq_Control"]), "Freq_Control"], col="blue", xaxt='n', cex.lab=2, cex.axis = 3)
dev.off()

comb3 <- comb[order(-comb[,"Freq_Control"]),]
freq_names_ctr <- c("21-5p", "4802-3p", "146a-5p", "378c", "142-3p", "146b-5p", "421", "30a-3p", "378a-5p", "378a-3p", "330-5p", "425-5p", "378i", "26a-5p-1",
"26a-5p-2", "223-5p", "191-5p", "30a-5p", "509-3p-2", "509-3p-3", "509-3p-1", "5571-3p", "301b", "766-3p", "199b-5p", "34a-5p")
pdf("Barplot_Top_Freq_Control.pdf", height=15, width=50)
barplot(comb3[1:26, "Freq_Control"],  col="blue", cex.lab=0.1, cex.axis = 3, names.arg=freq_names_ctr, cex.names=1.9)
dev.off()



# Analyse the SNPS and genes for the interesting miRs
mirs <- c("MI0003673_MIMAT0003327", "MI0001648_MIMAT0001541", "MI0003823_MIMAT0010251", "MI0000743_MIMAT0004677", "MI0000742_MIMAT0004676", "MI0000743_MIMAT0000686")

triplets <- read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/triplets_cit_ord_unique_COPD_new.txt', header=TRUE, sep="\t")

#mir449b_5p
mir449b_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0003673_MIMAT0003327"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(mir449b_5p) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(mir449b_5p) <- mir449b_5p$'SNP'
#write.table(mir449b_5p, "COPD_SNP_mir449b_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#mir449b_5p_2 <- triplets[which(triplets[,"Mir"]=="MI0003673_MIMAT0003327" & triplets[,"SNP"]%in%rownames(mir449b_5p)), c("SNP", "SNP_Chr", "SNP_Start","Gene", "Gene_Chr", "Gene_Start", "Gene_Stop" ) ]
mir449b_5p_all <- triplets[which(triplets[,"Mir"]=="MI0003673_MIMAT0003327"), c("SNP", "SNP_Chr", "SNP_Start","Gene", "Gene_Chr", "Gene_Start", "Gene_Stop" )]

mir_snp <- matrix(nrow=nrow(mir449b_5p), ncol=5)
colnames(mir_snp) <- c("SNP", "SNP_Chr", "SNP_Start", "Number_genes", "List_genes")
rownames(mir_snp) <- rownames(mir449b_5p)

for (i in 1:nrow(mir449b_5p)){
    nr <- 0
    l <- c()
    for (j in 1:nrow(mir449b_5p_all)){
        if (as.character(mir449b_5p_all[j, "SNP"])== as.character(mir449b_5p[i, "SNP"])){
            l <- c(l, as.character(mir449b_5p_all[j, "Gene"]))
            nr <- nr+1
        }
    }
    
    mir_snp[i,"SNP"] <- as.character(mir449b_5p[i,"SNP"])
    mir_snp[i,"SNP_Chr"] <- as.character(mir449b_5p[i,"SNP_Chr"])
    mir_snp[i,"SNP_Start"] <- as.numeric(mir449b_5p[i,"SNP_Start"])
    mir_snp[i,"Number_genes"] <- as.numeric(nr)
    mir_snp[i,"List_genes"] <- toString(l)
}

mir_snp_idx <- order(-as.numeric(mir_snp[,"Number_genes"]))
write.table(mir_snp[mir_snp_idx,], "COPD_SNP_mir449b_5p_all_corr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

pdf("COPD_SNP_gene_distr_mir449b_5p_all.pdf")
plot(density(as.numeric(mir_snp[,"Number_genes"])), main="mir449b_5p")
dev.off()

#mir449b_5p only genes
mir449b_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0003673_MIMAT0003327"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(mir449b_5p) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(mir449b_5p) <- mir449b_5p$'Gene'
write.table(mir449b_5p, "COPD_Gene_mir449b_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


#miR449a
mir449a <- unique(triplets[which(triplets[,"Mir"]=="MI0001648_MIMAT0001541"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(mir449a) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(mir449a) <- mir449a$'SNP'
#write.table(mir449a, "COPD_SNP_miR449a.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mir449a_all <- triplets[which(triplets[,"Mir"]=="MI0001648_MIMAT0001541"), c("SNP", "SNP_Chr", "SNP_Start","Gene", "Gene_Chr", "Gene_Start", "Gene_Stop" ) ]

mir_snp <- matrix(nrow=nrow(mir449a), ncol=5)
colnames(mir_snp) <- c("SNP", "SNP_Chr", "SNP_Start", "Number_genes", "List_genes")
rownames(mir_snp) <- rownames(mir449a)

for (i in 1:nrow(mir449a)){
    nr <- 0
    l <- c()
    for (j in 1:nrow(mir449a_all)){
        if (as.character(mir449a_all[j, "SNP"])== as.character(mir449a[i, "SNP"])){
            l <- c(l, as.character(mir449a_all[j, "Gene"]))
            nr <- nr+1
        }
    }
    
    mir_snp[i,"SNP"] <- as.character(mir449a[i,"SNP"])
    mir_snp[i,"SNP_Chr"] <- as.character(mir449a[i,"SNP_Chr"])
    mir_snp[i,"SNP_Start"] <- as.numeric(mir449a[i,"SNP_Start"])
    mir_snp[i,"Number_genes"] <- as.numeric(nr)
    mir_snp[i,"List_genes"] <- toString(l)
}

mir_snp_idx <- order(-as.numeric(mir_snp[,"Number_genes"]))
write.table(mir_snp[mir_snp_idx,], "COPD_SNP_mir449a_all_corr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

pdf("COPD_SNP_gene_distr_mir449a_all.pdf")
plot(density(as.numeric(mir_snp[,"Number_genes"])), main="mir449a")
dev.off()

#miR449a only genes
mir449a <- unique(triplets[which(triplets[,"Mir"]=="MI0001648_MIMAT0001541"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(mir449a) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(mir449a) <- mir449a$'Gene'
write.table(mir449a, "COPD_Gene_miR449a.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


#miR449c_5p
mir449c_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0003823_MIMAT0010251"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(mir449c_5p) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(mir449c_5p) <- mir449c_5p$'SNP'
#write.table(miR449c_5p, "COPD_SNP_miR449c_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mir449c_5p_all <- triplets[which(triplets[,"Mir"]=="MI0003823_MIMAT0010251"), c("SNP", "SNP_Chr", "SNP_Start","Gene", "Gene_Chr", "Gene_Start", "Gene_Stop" ) ]

mir_snp <- matrix(nrow=nrow(mir449c_5p), ncol=5)
colnames(mir_snp) <- c("SNP", "SNP_Chr", "SNP_Start", "Number_genes", "List_genes")
rownames(mir_snp) <- rownames(mir449c_5p)

for (i in 1:nrow(mir449c_5p)){
    nr <- 0
    l <- c()
    for (j in 1:nrow(mir449c_5p_all)){
        if (as.character(mir449c_5p_all[j, "SNP"])== as.character(mir449c_5p[i, "SNP"])){
            l <- c(l, as.character(mir449c_5p_all[j, "Gene"]))
            nr <- nr+1
        }
    }
    
    mir_snp[i,"SNP"] <- as.character(mir449c_5p[i,"SNP"])
    mir_snp[i,"SNP_Chr"] <- as.character(mir449c_5p[i,"SNP_Chr"])
    mir_snp[i,"SNP_Start"] <- as.numeric(mir449c_5p[i,"SNP_Start"])
    mir_snp[i,"Number_genes"] <- as.numeric(nr)
    mir_snp[i,"List_genes"] <- toString(l)
}

mir_snp_idx <- order(-as.numeric(mir_snp[,"Number_genes"]))
write.table(mir_snp[mir_snp_idx,], "COPD_SNP_mir449c_5p_all_corr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

pdf("COPD_SNP_gene_distr_mir449c_5p_all.pdf")
plot(density(as.numeric(mir_snp[,"Number_genes"])), main="mir449c_5p")
dev.off()


#miR449c_5p only genes
mir449c_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0003823_MIMAT0010251"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(mir449c_5p) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(mir449c_5p) <- mir449c_5p$'Gene'
write.table(mir449c_5p, "COPD_Gene_miR449c_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR34c_3p
mir34c_3p <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0004677"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(mir34c_3p) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(mir34c_3p) <- mir34c_3p$'SNP'
#write.table(mir34c_3p, "COPD_SNP_miR34c_3p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mir34c_3p_all <- triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0004677"), c("SNP", "SNP_Chr", "SNP_Start","Gene", "Gene_Chr", "Gene_Start", "Gene_Stop" ) ]

mir_snp <- matrix(nrow=nrow(mir34c_3p), ncol=5)
colnames(mir_snp) <- c("SNP", "SNP_Chr", "SNP_Start", "Number_genes", "List_genes")
rownames(mir_snp) <- rownames(mir34c_3p)

for (i in 1:nrow(mir34c_3p)){
    nr <- 0
    l <- c()
    for (j in 1:nrow(mir34c_3p_all)){
        if (as.character(mir34c_3p_all[j, "SNP"])== as.character(mir34c_3p[i, "SNP"])){
            l <- c(l, as.character(mir34c_3p_all[j, "Gene"]))
            nr <- nr+1
        }
    }
    
    mir_snp[i,"SNP"] <- as.character(mir34c_3p[i,"SNP"])
    mir_snp[i,"SNP_Chr"] <- as.character(mir34c_3p[i,"SNP_Chr"])
    mir_snp[i,"SNP_Start"] <- as.numeric(mir34c_3p[i,"SNP_Start"])
    mir_snp[i,"Number_genes"] <- as.numeric(nr)
    mir_snp[i,"List_genes"] <- toString(l)
}

mir_snp_idx <- order(-as.numeric(mir_snp[,"Number_genes"]))
write.table(mir_snp[mir_snp_idx,], "COPD_SNP_mir34c_3p_all_corr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR34c_3p only genes
mir34c_3p <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0004677"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(mir34c_3p) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(mir34c_3p) <- mir34c_3p$'Gene'
write.table(mir34c_3p, "COPD_Gene_miR34c_3p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR34b_3p
mir34b_3p <- unique(triplets[which(triplets[,"Mir"]=="MI0000742_MIMAT0004676"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(mir34b_3p) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(mir34b_3p) <- mir34b_3p$'SNP'
#write.table(mir34b_3p, "COPD_SNP_miR34b_3p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mir34b_3p_all <- triplets[which(triplets[,"Mir"]=="MI0000742_MIMAT0004676"), c("SNP", "SNP_Chr", "SNP_Start","Gene", "Gene_Chr", "Gene_Start", "Gene_Stop" ) ]

mir_snp <- matrix(nrow=nrow(mir34b_3p), ncol=5)
colnames(mir_snp) <- c("SNP", "SNP_Chr", "SNP_Start", "Number_genes", "List_genes")
rownames(mir_snp) <- rownames(mir34b_3p)

for (i in 1:nrow(mir34b_3p)){
    nr <- 0
    l <- c()
    for (j in 1:nrow(mir34b_3p_all)){
        if (as.character(mir34b_3p_all[j, "SNP"])== as.character(mir34b_3p[i, "SNP"])){
            l <- c(l, as.character(mir34b_3p_all[j, "Gene"]))
            nr <- nr+1
        }
    }
    
    mir_snp[i,"SNP"] <- as.character(mir34b_3p[i,"SNP"])
    mir_snp[i,"SNP_Chr"] <- as.character(mir34b_3p[i,"SNP_Chr"])
    mir_snp[i,"SNP_Start"] <- as.numeric(mir34b_3p[i,"SNP_Start"])
    mir_snp[i,"Number_genes"] <- as.numeric(nr)
    mir_snp[i,"List_genes"] <- toString(l)
}

mir_snp_idx <- order(-as.numeric(mir_snp[,"Number_genes"]))
write.table(mir_snp[mir_snp_idx,], "COPD_SNP_mir34b_3p_all_corr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

pdf("COPD_SNP_gene_distr_mir34b_3p_all.pdf")
plot(density(as.numeric(mir_snp[,"Number_genes"])), main="mir34b_3p")
dev.off()

#miR34b_3p only genes
mir34b_3p <- unique(triplets[which(triplets[,"Mir"]=="MI0000742_MIMAT0004676"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(mir34b_3p) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(mir34b_3p) <- mir34b_3p$'Gene'
write.table(mir34b_3p, "COPD_Gene_miR34b_3p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)



#miR34c_5p
mir34c_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0000686"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(mir34c_5p) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(mir34c_5p) <- mir34c_5p$'SNP'
#write.table(mir34c_5p, "COPD_SNP_miR34c_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mir34c_5p_all <- triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0000686"), c("SNP", "SNP_Chr", "SNP_Start","Gene", "Gene_Chr", "Gene_Start", "Gene_Stop" ) ]

mir_snp <- matrix(nrow=nrow(mir34c_5p), ncol=5)
colnames(mir_snp) <- c("SNP", "SNP_Chr", "SNP_Start", "Number_genes", "List_genes")
rownames(mir_snp) <- rownames(mir34c_5p)

for (i in 1:nrow(mir34c_5p)){
    nr <- 0
    l <- c()
    for (j in 1:nrow(mir34c_5p_all)){
        if (as.character(mir34c_5p_all[j, "SNP"])== as.character(mir34c_5p[i, "SNP"])){
            l <- c(l, as.character(mir34c_5p_all[j, "Gene"]))
            nr <- nr+1
        }
    }
    
    mir_snp[i,"SNP"] <- as.character(mir34c_5p[i,"SNP"])
    mir_snp[i,"SNP_Chr"] <- as.character(mir34c_5p[i,"SNP_Chr"])
    mir_snp[i,"SNP_Start"] <- as.numeric(mir34c_5p[i,"SNP_Start"])
    mir_snp[i,"Number_genes"] <- as.numeric(nr)
    mir_snp[i,"List_genes"] <- toString(l)
}

mir_snp_idx <- order(-as.numeric(mir_snp[,"Number_genes"]))
write.table(mir_snp[mir_snp_idx,], "COPD_SNP_mir34c_5p_all_corr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR34c_5p only genes
mir34c_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0000686"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(mir34c_5p) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(mir34c_5p) <- mir34c_5p$'Gene'
write.table(mir34c_5p, "COPD_Gene_miR34c_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
