mirs <- c("MI0003673_MIMAT0003327", "MI0001648_MIMAT0001541", "MI0003823_MIMAT0010251", "MI0000743_MIMAT0004677", "MI0000742_MIMAT0004676", "MI0000743_MIMAT0000686", "MI0000742_MIMAT0000685")

triplets <- read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/triplets_cit_ord_unique_COPD_new.txt', header=TRUE, sep="\t")
triplets2 <- read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_ILD/triplets_cit_ord_unique_ILD_new.txt', header=TRUE, sep="\t")
nrow(unique(triplets[,cbind("Mir","Gene")]))

#run from here
mirs <- c("MI0003673_MIMAT0003327", "MI0001648_MIMAT0001541", "MI0003823_MIMAT0010251", "MI0000743_MIMAT0004677", "MI0000742_MIMAT0004676", "MI0000743_MIMAT0000686", "MI0000742_MIMAT0000685")

tr = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/Triplets_eQTL_with_CIT_COPD.txt', header=TRUE, sep="\t")

tr_cit <- tr[which(tr$Causal_Call==1),]

nrow(unique(tr_cit[,cbind("Mir","Gene")]))


copd <- unique(tr_cit[,cbind("Mir","Gene")])

mir.genes = list()
mir.index = 1
different.mir = unique(copd[, "Mir"])
for (mir in different.mir) {
    mir.genes[mir.index] = list(copd[copd[, "Mir"] == mir, "Gene"])
    mir.index = mir.index+1
}
names(mir.genes) = different.mir
mir.matrix = matrix(nrow = length(different.mir), ncol = length(different.mir))
rownames(mir.matrix) = colnames(mir.matrix) = different.mir
for (mir1 in 1:(length(mir.genes) - 1)) {
    mir.matrix[mir1, mir1]=1
    for (mir2 in (mir1 + 1):length(mir.genes)) {
        mir.genes.1 = mir.genes[[mir1]]
        mir.genes.2 = mir.genes[[mir2]]
        mir.matrix[mir2, mir1] = mir.matrix[mir1, mir2] = length(intersect(mir.genes.1, mir.genes.2)) / length(union(mir.genes.1, mir.genes.2))
    }
}
mir.matrix[length(mir.genes), length(mir.genes)] = 1

# names(mir.genes) = different.mir
# valid.mir.indexes = sapply(1:length(mir.genes), function(x) { length(mir.genes[[x]]) > 1 })
# valid.different.mir=different.mir[valid.mir.indexes]
# valid.different.mir.number = length(valid.different.mir)
# mir.matrix = matrix(nrow = valid.different.mir.number, ncol = valid.different.mir.number)
# rownames(mir.matrix) = colnames(mir.matrix) = valid.different.mir
# for (mir1 in 1:(valid.different.mir.number - 1)) {
    # mir.matrix[mir1, mir1]=1
    # for (mir2 in (mir1 + 1):valid.different.mir.number) {
        # mir.genes.1 = mir.genes[[mir1]]
        # mir.genes.2 = mir.genes[[mir2]]
        # mir.matrix[mir2, mir1] = mir.matrix[mir1, mir2] = length(intersect(mir.genes.1, mir.genes.2)) / length(union(mir.genes.1, mir.genes.2))
    # }
# }
# mir.matrix[valid.different.mir.number, valid.different.mir.number] = 1


#hclust.result=hclust(dist(mir.matrix))
#matrix.reordered=mir.matrix[hclust.result$order, hclust.result$order]
library(gplots) 
jpeg(file="/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/heatmap_no_genes_all.COPD.diag.jpeg", width=2000, height=2000)
par(cex.main=3)
orangescale <- colorRampPalette(c("white", "orange", "red"))
heatmap.plot=heatmap.2(mir.matrix,
         main="",
         col=orangescale,
         #col=bluered,
         Rowv=TRUE,
         # Rowv=FALSE,
         Colv=TRUE,
         # Colv=FALSE,
         trace = "none",
         density.info = "none",
         #xlab = "cell lines",
         #ylab = "data types",
         hclustfun = hclust,
         scale=c("none"),
        ColSideColors=FALSE,
        #labRow=c(" "," "),
        cexRow=2, 
        cexCol=2,
        lmat=rbind(c(5, 4, 0), c(0, 1, 0), c(3, 2, 0)),
        lwid=c(2, 10, 3),
        lhei=c(2, 0.15, 7)
)
dev.off()
mir.matrix.carpet=heatmap.plot$carpet
mir.matrix.carpet=mir.matrix[heatmap.plot$rowInd, heatmap.plot$colInd]
mir.matrix.carpet=mir.matrix.carpet[nrow(mir.matrix.carpet):1,]


ids = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/miRNAs_expression_annotation_v19.txt', header=TRUE, sep="\t")
rownames(ids) <- ids[,"V5"]
common.matrix.rownames=colnames(mir.matrix.carpet)[colnames(mir.matrix.carpet)%in% rownames(ids)]
mir.matrix.carpet2 <- mir.matrix.carpet[common.matrix.rownames, common.matrix.rownames]
ids2 <- ids[common.matrix.rownames,]

colnames(mir.matrix.carpet2) = rownames(mir.matrix.carpet2) = ids2[, "V4"]
#write.csv(mir.matrix.carpet2, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/carpet2.copd.csv")
colnames(mir.matrix.carpet2) = rownames(mir.matrix.carpet2) = sapply(ids2[, "V4"], function(x) strsplit(as.character(x), "_")[[1]][2])

library(gplots) 
orangescale <- colorRampPalette(c("white", "orange", "red"))
mir.matrix.heatmap <- mir.matrix.carpet2[1:75,1:75]
diag(mir.matrix.heatmap) <- 0
 jpeg(file="/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/heatmap_no_genes_all.COPD.clustered.mirs.diag.jpeg", width=4000, height=4000)
 heatmap.plot=heatmap.2(mir.matrix.heatmap[75:1,1:75],
          main="",
          col=orangescale,
          #col=bluered,
          Rowv=TRUE,
          #Rowv=FALSE,
          Colv=TRUE,
          #Colv=FALSE,
          trace = "none",
          density.info = "none",
          #xlab = "cell lines",
          #ylab = "data types",
         # #hclustfun = hclust,
          scale=c("none"),
         ColSideColors=FALSE,
         #labRow=c(" "," "),
         cexRow=5, 
         cexCol=5,
         #lmat=rbind(c(5, 1, 6), c(3, 2, 4), c(7,8,9)), lhei=c(2.5, 10,2), lwid=c(2, 10, 2)
         #5 - legend, 4 - column dendrogram, 3 - rowdendrogram, 2 - heatmap, 6,7,1 - padding
         lmat=rbind(c(5, 0,0, 0), c(0,0,4,0),c(0,3, 2, 0), c(0,7,1,6)), lhei=c(2, 2.5,10,2.5), lwid=c(2,2.5, 10, 2.5)
         
         #good key
         #lmat=rbind(c(5, 0,0, 0), c(0,0,4,0),c(0,3, 2, 0), c(0,7,1,6)), lhei=c(2, 2.5,20,2.5), lwid=c(2,2.5, 2, 2.5)
          )
 dev.off()

####OLDER analysis

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

y <- c()
yn <- c()
w <- c()
wn <- c()

y[[1]] <- as.character(mir449b_5p_genes_pos)
yn[[1]] <- as.character(mir449b_5p_genes_neg)
w[[1]] <- as.character(mir449b_5p_genes_pos2)
wn[[1]] <- as.character(mir449b_5p_genes_neg2)

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

y[[2]] <- as.character(mir449a_genes_pos)
yn[[2]] <- as.character(mir449a_genes_neg)
w[[2]] <- as.character(mir449a_genes_pos2)
wn[[2]] <- as.character(mir449a_genes_neg2)

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

y[[3]] <- as.character(mir449c_5p_genes_pos)
yn[[3]] <- as.character(mir449c_5p_genes_neg)
w[[3]] <- as.character(mir449c_5p_genes_pos2)
wn[[3]] <- as.character(mir449c_5p_genes_neg2)

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

y[[4]] <- as.character(mir34c_3p_genes_pos)
yn[[4]] <- as.character(mir34c_3p_genes_neg)
w[[4]] <- as.character(mir34c_3p_genes_pos2)
wn[[4]] <- as.character(mir34c_3p_genes_neg2)

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

y[[5]] <- as.character(mir34b_3p_genes_pos)
yn[[5]] <- as.character(mir34b_3p_genes_neg)
w[[5]] <- as.character(mir34b_3p_genes_pos2)
wn[[5]] <- as.character(mir34b_3p_genes_neg2)


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

y[[6]] <- as.character(mir34c_5p_genes_pos)
yn[[6]] <- as.character(mir34c_5p_genes_neg)
w[[6]] <- as.character(mir34c_5p_genes_pos2)
wn[[6]] <- as.character(mir34c_5p_genes_neg2)

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

y[[7]] <- as.character(mir34b_5p_genes_pos)
yn[[7]] <- as.character(mir34b_5p_genes_neg)
w[[7]] <- as.character(mir34b_5p_genes_pos2)
wn[[7]] <- as.character(mir34b_5p_genes_neg2)


copd <- matrix(ncol=7, nrow=7)
for (i in 1:7){
    for (j in 1:7){
        copd[i, j] <- length(intersect(y[[i]], y[[j]]))
    }
}

copd2 <- copd[1:6, 1:6]

#colnames(copd) <- c("miR-449b-5p", "miR-449a", "miR-449c-5p", "miR-34c-3p", "miR-34b-3p", "miR-34c-5p", "miR-34b-5p")
#rownames(copd) <- c("miR-449b-5p", "miR-449a", "miR-449c-5p", "miR-34c-3p", "miR-34b-3p", "miR-34c-5p", "miR-34b-5p")

colnames(copd2) <- c("miR-449b-5p", "miR-449a", "miR-449c-5p", "miR-34c-3p", "miR-34b-3p", "miR-34c-5p")
rownames(copd2) <- c("miR-449b-5p", "miR-449a", "miR-449c-5p", "miR-34c-3p", "miR-34b-3p", "miR-34c-5p")

library(gplots)    #create a variable with colors that will be used in heat map

#bluered <- colorRampPalette(c("darkblue", "blue", "white","red", "darkred"))


jpeg(file="/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/heatmap_no_genes.COPD.jpeg", width=1000, height=1000)
par(cex.main=3)
orangescale <- colorRampPalette(c("white", "orange", "red", "darkred"))
heatmap.2(copd2,
         main="",
         col=orangescale,
         #col=bluered,
         Rowv=TRUE,
         #Rowv=FALSE,
         Colv=TRUE,
         #Colv=FALSE,
         trace = "none",
         density.info = "none",
         #xlab = "cell lines",
         #ylab = "data types",
         hclustfun = hclust,
         scale=c("none"),
        ColSideColors=FALSE,
        #labRow=c(" "," "),
        cexRow=2, 
        cexCol=2,
        lmat=rbind(c(5, 4, 0), c(0, 1, 0), c(3, 2, 0)),
        lwid=c(2, 10, 3),
        lhei=c(2, 0.15, 7)
)
dev.off()


#ILD

ild <- matrix(ncol=7, nrow=7)
for (i in 1:7){
    for (j in 1:7){
        ild[i, j] <- length(intersect(w[[i]], w[[j]]))
    }
}

ild2 <- ild[c(1,2,3,4,6,7), c(1,2,3,4,6,7)]

colnames(ild2) <- c("miR-449b-5p", "miR-449a", "miR-449c-5p", "miR-34c-3p", "miR-34c-5p", "miR-34b-5p")
rownames(ild2) <- c("miR-449b-5p", "miR-449a", "miR-449c-5p", "miR-34c-3p", "miR-34c-5p", "miR-34b-5p")


library(gplots)    #create a variable with colors that will be used in heat map

#bluered <- colorRampPalette(c("darkblue", "blue", "white","red", "darkred"))


jpeg(file="/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/heatmap_no_genes.ILD.jpeg", width=1000, height=1000)
par(cex.main=3)
orangescale <- colorRampPalette(c("white", "orange", "red", "darkred"))
heatmap.2(ild2,
         main="",
         col=orangescale,
         #col=bluered,
         Rowv=TRUE,
         #Rowv=FALSE,
         Colv=TRUE,
         #Colv=FALSE,
         trace = "none",
         density.info = "none",
         #xlab = "cell lines",
         #ylab = "data types",
         hclustfun = hclust,
         scale=c("none"),
        ColSideColors=FALSE,
        #labRow=c(" "," "),
        cexRow=2, 
        cexCol=2,
        lmat=rbind(c(5, 4, 0), c(0, 1, 0), c(3, 2, 0)),
        lwid=c(2, 10, 3),
        lhei=c(2, 0.15, 7)
)
dev.off()