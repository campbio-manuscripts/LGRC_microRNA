setwd("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Target_Predict/")
integer.matrix=read.table("130605_IsomiR_Targets_All_Refseq_Matched.txt")
setwd("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Preprocess/")
real.matrix=read.table("cortable_ild_fdr_isoforms.txt")
cor.matrix=read.table("cortable_ild_cor_isoforms.txt")
filtered.mirnas=c("MI0000743_MIMAT0000686_111384175", "MI0000743_MIMAT0000686_111384176")
integer.matrix.filtered = integer.matrix[,filtered.mirnas]
real.matrix.filtered = real.matrix[,filtered.mirnas]
cor.matrix.filtered = cor.matrix[,filtered.mirnas]
filter.mirna.indexes=function(mirna.name){
    integer.matrix.filtered[,mirna.name]>0&real.matrix.filtered[,mirna.name]<0.25&cor.matrix.filtered[,mirna.name]<0
}

filter.mirna.indexes.targetscan=function(mirna.name){
    integer.matrix.filtered[,mirna.name]>0
}

first.mirna.indexes.targetscan=filter.mirna.indexes.targetscan("MI0000743_MIMAT0000686_111384175")
second.mirna.indexes.targetscan=filter.mirna.indexes.targetscan("MI0000743_MIMAT0000686_111384176")

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Target_Predict/Targets_cor.pdf")
plot(density(cor.matrix.filtered[first.mirna.indexes.targetscan,"MI0000743_MIMAT0000686_111384175"]), lwd=3, col="blue", main=paste("Distribution of gene correlations with miR-34c","\n","KS-test p<1e-7; t-test p<1e-11"), ylim=c(0, 2.5), cex.lab=1.5, cex.axis=1.5, las=1, font.lab=2, xlab="Spearman correlation", font.axis=2, cex.main=1.5)
lines(density(cor.matrix.filtered[second.mirna.indexes.targetscan,"MI0000743_MIMAT0000686_111384176"]), lwd=3, col="red")
lines(density(c(cor.matrix.filtered[second.mirna.indexes.targetscan & first.mirna.indexes.targetscan,"MI0000743_MIMAT0000686_111384175"], cor.matrix.filtered[second.mirna.indexes.targetscan & first.mirna.indexes.targetscan,"MI0000743_MIMAT0000686_111384176"])), lwd=3, col="yellow")
lines(density(c(cor.matrix.filtered[!(second.mirna.indexes.targetscan & first.mirna.indexes.targetscan),"MI0000743_MIMAT0000686_111384175"], cor.matrix.filtered[!(second.mirna.indexes.targetscan & first.mirna.indexes.targetscan),"MI0000743_MIMAT0000686_111384176"])), lwd=3, col="black")
par(font=2)
legend("topright", legend=c("5' isomiR", "Canonical", "Both", "Neither"), col=c("blue", "red", "yellow", "black"), cex=1.5, lty = 1, lwd=3)
dev.off()

t.test(cor.matrix.filtered[first.mirna.indexes.targetscan,"MI0000743_MIMAT0000686_111384175"], c(cor.matrix.filtered[!(second.mirna.indexes.targetscan & first.mirna.indexes.targetscan),"MI0000743_MIMAT0000686_111384175"], cor.matrix.filtered[!(second.mirna.indexes.targetscan & first.mirna.indexes.targetscan),"MI0000743_MIMAT0000686_111384176"]))
t.test(cor.matrix.filtered[second.mirna.indexes.targetscan,"MI0000743_MIMAT0000686_111384176"], c(cor.matrix.filtered[!(second.mirna.indexes.targetscan & first.mirna.indexes.targetscan),"MI0000743_MIMAT0000686_111384175"], cor.matrix.filtered[!(second.mirna.indexes.targetscan & first.mirna.indexes.targetscan),"MI0000743_MIMAT0000686_111384176"]))
t.test(c(cor.matrix.filtered[second.mirna.indexes.targetscan & first.mirna.indexes.targetscan,"MI0000743_MIMAT0000686_111384175"], cor.matrix.filtered[second.mirna.indexes.targetscan & first.mirna.indexes.targetscan,"MI0000743_MIMAT0000686_111384176"]), c(cor.matrix.filtered[!(second.mirna.indexes.targetscan & first.mirna.indexes.targetscan),"MI0000743_MIMAT0000686_111384175"], cor.matrix.filtered[!(second.mirna.indexes.targetscan & first.mirna.indexes.targetscan),"MI0000743_MIMAT0000686_111384176"]))

first.mirna.indexes=filter.mirna.indexes(filtered.mirnas[1])
second.mirna.indexes=filter.mirna.indexes(filtered.mirnas[2])
first.mirna.valid=rownames(integer.matrix.filtered)[first.mirna.indexes]
second.mirna.valid=rownames(integer.matrix.filtered)[second.mirna.indexes]
intersection.valid=rownames(integer.matrix.filtered)[first.mirna.indexes&second.mirna.indexes]
setwd("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Target_Predict/")
write.table(matrix(nrow=length(first.mirna.valid), ncol=1, data=first.mirna.valid), paste("Valid_",filtered.mirnas[1],".ild.csv", sep=""), row.names=F,col.names=F, quote=F)
write.table(matrix(nrow=length(second.mirna.valid), ncol=1, data=second.mirna.valid), paste("Valid_",filtered.mirnas[2],".ild.csv", sep=""), row.names=F,col.names=F, quote=F)
write.table(matrix(nrow=length(intersection.valid), ncol=1, data=intersection.valid), "Valid_intersection.ild.csv", row.names=F,col.names=F, quote=F)

all(rownames(integer.matrix)==rownames(real.matrix)&rownames(integer.matrix)==rownames(cor.matrix))
selected.intersection=read.csv("Valid_MI0000743_MIMAT0000686_111384176.ild.csv", header=F,stringsAsFactor=F)[,1]

ids <- read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/Preprocess/OtherDatasets/121212_LGRC_Gene_Array_Annotation2.txt", stringsAsFactor=F)
colnames(ids)=ids[1,]
ids=ids[2:nrow(ids),]
rownames(ids)=ids[,"NAME"]
intersection.gene.names = ids[selected.intersection,"GENE_SYMBOL"]
write.table(matrix(nrow=length(intersection.gene.names), ncol=1, data=intersection.gene.names), "Targets.ild.canonical.csv", row.names=F,col.names=F, quote=F)
