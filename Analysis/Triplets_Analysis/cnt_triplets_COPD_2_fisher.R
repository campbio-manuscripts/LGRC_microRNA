triplets <- read.table('Triplets_eQTL_with_CIT_COPD.txt', header=TRUE, sep="\t")
nrow(triplets)

triplets_c <- read.table('Triplets_eQTL_with_CIT_Control.txt', header=TRUE, sep="\t")
nrow(triplets_c)

primary.key=paste(as.character(triplets[, "SNP"], as.character(triplets[, "Mir"]), as.character(triplets[, "Gene"])))
filtered.triplets=unique(triplets)
nrow(filtered.triplets)

primary.key_c=paste(as.character(triplets_c[, "SNP"], as.character(triplets_c[, "Mir"]), as.character(triplets_c[, "Gene"])))
filtered.triplets_c=unique(triplets_c)
nrow(filtered.triplets_c)

filtered.triplets.cit <- filtered.triplets[which(filtered.triplets$Causal_Call==1),]
nrow(filtered.triplets.cit)

filtered.triplets.cit_c <- filtered.triplets_c[which(filtered.triplets_c$Causal_Call==1),]
nrow(filtered.triplets.cit_c)

primary.key=paste(as.character(filtered.triplets.cit[, "Mir"]), as.character(filtered.triplets.cit[, "Gene"]))
mir_gene = subset(filtered.triplets.cit, !duplicated(primary.key))
nrow(mir_gene)

primary.key_c=paste(as.character(filtered.triplets.cit_c[, "Mir"]), as.character(filtered.triplets.cit_c[, "Gene"]))
mir_gene_c = subset(filtered.triplets.cit_c, !duplicated(primary.key_c))
nrow(mir_gene_c)

mirs_unique <- unique(as.character(mir_gene$Mir))
mirs_unique_c <- unique(as.character(mir_gene_c$Mir))

mirs_unique_all <- union(mirs_unique, mirs_unique_c)

mir_t <- matrix(nrow=length(mirs_unique_all), ncol=4)
colnames(mir_t) <- c("Freq_COPD", "Freq_Control", "p-value", "FDR")
rownames(mir_t) <- mirs_unique_all

for (i in mirs_unique_all){
    nr <- 0
    for (j in 1:nrow(mir_gene)){
        if (as.character(mir_gene[j,"Mir"])==i){
            nr <- nr + 1
        }
     mir_t[i, "Freq_COPD"] <- nr/nrow(mir_gene)    
    }
    
    nr_c <- 0
    for (j in 1:nrow(mir_gene_c)){
        if (as.character(mir_gene_c[j,"Mir"])==i){
            nr_c <- nr_c + 1
        }
     mir_t[i, "Freq_Control"] <- nr_c/nrow(mir_gene_c)    
    }
    x <- c(nr, nr_c, nrow(mir_gene)-nr , nrow(mir_gene_c)-nr_c)
    mir_t[i, "p-value"] <- as.numeric(fisher.test(t(matrix(x, ncol=2, nrow=2)))[1])
}

mir_t[, "FDR"] <- p.adjust(mir_t[, "p-value"], "fdr")

#mir_t[mirs_unique, "FDR"] <- p.adjust(mir_t[mirs_unique, "p-value"], "fdr")
#mir_t[mirs_unique_c, "FDR"] <- p.adjust(mir_t[mirs_unique_c, "p-value"], "fdr")


mir_t_ord <- order(-mir_t[,1])
mir_t_ord_desc <- mir_t[mir_t_ord,]
write.table(mir_t_ord_desc, "Freq_miRNA_COPD_new_fisher1_corr.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

mir_t_ord <- order(-mir_t[,2])
mir_t_ord_desc <- mir_t[mir_t_ord,]
write.table(mir_t_ord_desc, "Freq_miRNA_COPD_new_fisher2_corr.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
