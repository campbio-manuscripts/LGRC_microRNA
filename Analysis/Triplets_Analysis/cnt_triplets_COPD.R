shared.library.path <- file.path("/unprotected/projects/cbmhive/R_packages", sprintf("R-%s", getRversion()));
.libPaths(c(shared.library.path, .libPaths()));


triplets <- read.table('Triplets_eQTL_with_CIT_COPD.txt', header=TRUE, sep="\t")
nrow(triplets)

triplets_ord <- order(triplets$Mir)
triplets <- triplets[triplets_ord,]

triplets2 <- matrix(ncol=ncol(triplets), nrow=0)
colnames(triplets2) <- colnames(triplets)
triplets2 <- rbind(triplets2, triplets[1,])
triplets3 <- triplets2

for (i in 2 : nrow(triplets)){
nr <- 0
nr_p <- 0
for (j in 1 : i-1){
    if (identical(as.character(triplets[j, "SNP"]), as.character(triplets[i, "SNP"])) &
    identical(as.character(triplets[j, "Mir"]), as.character(triplets[i, "Mir"])) &
    identical(as.character(triplets[j, "Gene"]), as.character(triplets[i, "Gene"])))
    {
      nr <- nr+1  
       
    }
    if (identical(as.character(triplets[j, "Mir"]), as.character(triplets[i, "Mir"])) &
    identical(as.character(triplets[j, "Gene"]), as.character(triplets[i, "Gene"])))
    {
      nr_p <- nr_p+1  
       
    }

}
if (nr == 0 ){
    triplets2 <- rbind(triplets2, triplets[i,])
}

if (nr_p == 0 ){
    triplets3 <- rbind(triplets3, triplets[i,])
}

}
write.table(triplets2, "triplets_all_unique_COPD.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

triplets_cit <- triplets2[which(triplets2$Causal_Call==1),]
nrow(triplets_cit)
triplets_cit_ord <- order(triplets_cit$P_Causal)
write.table(triplets_cit[triplets_cit_ord,], "triplets_cit_ord_unique_COPD.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#unique mir/gene
# triplets2 <- matrix(ncol=ncol(triplets), nrow=0)
# colnames(triplets2) <- colnames(triplets)
# triplets2 <- rbind(triplets2, triplets[1,])

# for (i in 2 : nrow(triplets)){
# nr <- 0
# for (j in 1 : i-1){
    # if (identical(as.character(triplets[j, "Mir"]), as.character(triplets[i, "Mir"])) &
    # identical(as.character(triplets[j, "Gene"]), as.character(triplets[i, "Gene"])))
    # {
      # nr <- nr+1  
       
    # }
# }
# if (nr == 0 ){
    # triplets2 <- rbind(triplets2, triplets[i,])
# }
# }
write.table(triplets3, "triplets_miRNA_gene_unique_COPD.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

triplets_cit <- triplets3[which(triplets3$Causal_Call==1),]
nrow(triplets_cit)

mirs_unique <- unique(as.character(triplets_cit$Mir))
mir_t <- matrix(nrow=length(mirs_unique), ncol=1)
colnames(mir_t) <- c("Freq")
rownames(mir_t) <- mirs_unique

for (i in mirs_unique){
    nr <- 0
    for (j in 1:nrow(triplets_cit)){
        if (as.character(triplets_cit[j,"Mir"])==i){
            nr <- nr + 1
        }
     mir_t[i, "Freq"] <- nr/length(mirs_unique)    
    }
}
mir_t_ord <- order(-mir_t[,1])
mir_t_ord_desc <- mir_t[mir_t_ord,1]
write.table(mir_t_ord_desc, "Freq_miRNA_COPD.txt", sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

