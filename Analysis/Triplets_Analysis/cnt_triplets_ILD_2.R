triplets <- read.table('Triplets_eQTL_with_CIT_ILD.txt', header=TRUE, sep="\t")
nrow(triplets)

primary.key=paste(as.character(triplets[, "SNP"], as.character(triplets[, "Mir"]), as.character(triplets[, "Gene"])))
filtered.triplets=unique(triplets)
nrow(filtered.triplets)

filtered.triplets.cit <- filtered.triplets[which(filtered.triplets$Causal_Call==1),]
nrow(filtered.triplets.cit)
filtered.triplets.cit.ord <- order(filtered.triplets.cit$P_Causal)
#write.table(filtered.triplets.cit[filtered.triplets.cit.ord,], "triplets_cit_ord_unique_ILD_new.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
#write.table(filtered.triplets, "triplets_all_ord_unique_ILD_new.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

causal_triplets <- read.table('triplets_cit_ord_unique_ILD_new.txt', header=TRUE, sep="\t")
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
#write.table(mir_t_ord_desc, "Freq_miRNA_ILD_new.txt", sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

write.table(mir_t_ord_desc, "No_miRNA_ILD_new.txt", sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)



freq_ILD <- read.table("Freq_miRNA_ILD_new.txt", header=FALSE, sep="\t")
colnames(freq_ILD) <- c("Mir_ILD", "Freq_ILD")
rownames(freq_ILD) <- freq_ILD$Mir_ILD
freq_Control <- read.table("Freq_miRNA_Control_new.txt", header=FALSE, sep="\t")
colnames(freq_Control) <- c("Mir_Control", "Freq_Control")
rownames(freq_Control) <- freq_Control$Mir_Control

length(intersect(rownames(freq_Control), rownames(freq_ILD)))
length(union(rownames(freq_Control), rownames(freq_ILD)))

comb <- matrix(0, nrow=length(union(rownames(freq_Control), rownames(freq_ILD))), ncol=2)
colnames(comb) <- c("Freq_Control", "Freq_ILD")
rownames(comb) <- union(rownames(freq_Control), rownames(freq_ILD))

comb[rownames(freq_Control),"Freq_Control"] <- freq_Control$Freq_Control
comb[rownames(freq_ILD),"Freq_ILD"] <- freq_ILD$Freq_ILD

pdf("Barplot_Freq_ILD.pdf", height=15, width=50)
barplot(comb[order(-comb[,"Freq_ILD"]), "Freq_ILD"], col="blue", xaxt='n', cex.lab=2, cex.axis = 3)
dev.off()

comb2 <- comb[order(-comb[,"Freq_ILD"]),]
freq_names <- c("92b-3p", "449a", "200a-5p", "31-5p", "92b-5p", "449c-5p", "200b-3p", "31-3p", "190b", "449b-5p", "511-1", "511-2", "34c-5p", "34c-3p", "146b-5p", "2110", "34b-5p")

pdf("Barplot_Top_Freq_ILD.pdf", height=15, width=50)
barplot(comb2[1:17, "Freq_ILD"],  col="blue", cex.lab=0.1, cex.axis = 3, names.arg=freq_names, cex.names=2.2)
dev.off()

pdf("Freq_ILD.pdf")
plot(comb[,"Freq_Control"], comb[,"Freq_ILD"], xlab="Control", ylab="ILD", pch=19, col="blue")
dev.off()


triplets <- read.table('triplets_cit_ord_unique_ILD_new.txt', header=TRUE, sep="\t")

#mir449b_5p
mir449b_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0003673_MIMAT0003327"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(mir449b_5p) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(mir449b_5p) <- mir449b_5p$'SNP'
#write.table(mir449b_5p, "COPD_SNP_mir449b_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mir449b_5p_all <- triplets[which(triplets[,"Mir"]=="MI0003673_MIMAT0003327"), c("SNP", "SNP_Chr", "SNP_Start","Gene", "Gene_Chr", "Gene_Start", "Gene_Stop" ) ]

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
write.table(mir_snp[mir_snp_idx,], "ILD_SNP_mir449b_5p_all_corr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

pdf("ILD_SNP_gene_distr_mir449b_5p_all.pdf")
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
#write.table(miR449a, "COPD_SNP_miR449a.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

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
write.table(mir_snp[mir_snp_idx,], "ILD_SNP_mir449a_all_corr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

pdf("ILD_SNP_gene_distr_mir449a_all.pdf")
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
#write.table(mir449c_5p, "COPD_SNP_miR449c_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

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
write.table(mir_snp[mir_snp_idx,], "ILD_SNP_mir449c_5p_all_corr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

pdf("ILD_SNP_gene_distr_mir449c_5p_all.pdf")
plot(density(as.numeric(mir_snp[,"Number_genes"])), main="mir449c-5p")
dev.off()


#miR449c_5p
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
write.table(mir_snp[mir_snp_idx,], "ILD_SNP_mir34c_3p_all_corr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)



#miR34c_3p
mir34c_3p <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0004677"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(mir34c_3p) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(mir34c_3p) <- mir34c_3p$'Gene'
write.table(mir34c_3p, "COPD_Gene_miR34c_3p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR34b_5p
mir34b_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0000742_MIMAT0000685"), c("SNP", "SNP_Chr", "SNP_Start")])
colnames(mir34b_5p) <- c("SNP", "SNP_Chr", "SNP_Start")
rownames(mir34b_5p) <- mir34b_5p$'SNP'
#write.table(mir34b_3p, "COPD_SNP_miR34b_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mir34b_5p_all <- triplets[which(triplets[,"Mir"]=="MI0000742_MIMAT0000685"), c("SNP", "SNP_Chr", "SNP_Start","Gene", "Gene_Chr", "Gene_Start", "Gene_Stop" ) ]

mir_snp <- matrix(nrow=nrow(mir34b_5p), ncol=5)
colnames(mir_snp) <- c("SNP", "SNP_Chr", "SNP_Start", "Number_genes", "List_genes")
rownames(mir_snp) <- rownames(mir34b_5p)

for (i in 1:nrow(mir34b_5p)){
    nr <- 0
    l <- c()
    for (j in 1:nrow(mir34b_5p_all)){
        if (as.character(mir34b_5p_all[j, "SNP"])== as.character(mir34b_5p[i, "SNP"])){
            l <- c(l, as.character(mir34b_5p_all[j, "Gene"]))
            nr <- nr+1
        }
    }
    
    mir_snp[i,"SNP"] <- as.character(mir34b_5p[i,"SNP"])
    mir_snp[i,"SNP_Chr"] <- as.character(mir34b_5p[i,"SNP_Chr"])
    mir_snp[i,"SNP_Start"] <- as.numeric(mir34b_5p[i,"SNP_Start"])
    mir_snp[i,"Number_genes"] <- as.numeric(nr)
    mir_snp[i,"List_genes"] <- toString(l)
}

mir_snp_idx <- order(-as.numeric(mir_snp[,"Number_genes"]))
write.table(mir_snp[mir_snp_idx,], "ILD_SNP_mir34b_5p_all_corr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR34b_5p
mir34b_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0000742_MIMAT0000685"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(mir34b_5p) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(mir34b_5p) <- mir34b_5p$'Gene'
write.table(mir34b_3p, "COPD_Gene_miR34b_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#mir34c_5p
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
write.table(mir_snp[mir_snp_idx,], "ILD_SNP_mir34c_5p_all_corr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#miR34c_5p
mir34c_5p <- unique(triplets[which(triplets[,"Mir"]=="MI0000743_MIMAT0000686"), c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")])
colnames(mir34c_5p) <- c("Gene", "Gene_Chr", "Gene_Start", "Gene_Stop")
rownames(mir34c_5p) <- mir34c_5p$'Gene'
write.table(mir34c_5p, "COPD_Gene_miR34c_5p.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
