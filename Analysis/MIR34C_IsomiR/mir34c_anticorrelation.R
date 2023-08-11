library(affy)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

mir = readRDS("../../Data/MicroRNA/LGRC_MicroRNA_Overlap_Gene.rds")
gene = readRDS("../../Data/Gene/LGRC_Gene_Overlap_MicroRNA.rds")
iso = readRDS("../../Data/MicroRNA/LGRC_MicroRNA_Isoform_RPM_Log2_Batch_Corrected.rds")


targets = readRDS("../Target_Prediction/130826_IsomiR_Targets_Matched.rds")
targets.locus = substring(colnames(targets), 1, 22)
mir.34c.isomir = "MI0000743_MIMAT0000686_111384175"
mir.34c.canonical = "MI0000743_MIMAT0000686_111384176"

matrix.cor = function(m1, m2, method="spearman") {
  all.corr <- t(cor(t(m1), t(m2), method=method))
  all.corr.t <- all.corr / (sqrt( (1 - (all.corr^2)) / (ncol(m1) - 2)))
  all.corr.t.p <- 2*(1 - pt(abs(all.corr.t), ncol(m1) - 2))
  f = p.adjust(all.corr.t.p, 'fdr')
  f.m = matrix(f, nrow=nrow(all.corr.t.p), ncol=ncol(all.corr.t.p))
  colnames(f.m) = rownames(m1)
  rownames(f.m) = rownames(m2)
  return(list(Cor=all.corr, Cor.t=all.corr.t, Cor.p=all.corr.t.p, Cor.fdr=f.m))
}

ild.cor = matrix.cor(exprs(gene[,gene$ILD == 1]), exprs(mir[,mir$ILD == 1]))
copd.cor = matrix.cor(exprs(gene[,gene$COPD == 1]), exprs(mir[,mir$COPD == 1]))
control.cor = matrix.cor(exprs(gene[,gene$Control == 1]), exprs(mir[,mir$Control == 1]))



mir34c.isomir.anti.targets.ind = ild.cor$Cor.fdr["MI0000743_MIMAT0000686",] < 0.25 & ild.cor$Cor["MI0000743_MIMAT0000686",] < 0 & targets[,mir.34c.isomir] > 0
mir34c.isomir.anti.targets.names = fData(gene)[mir34c.isomir.anti.targets.ind,"GENE_SYMBOL"]

mir34c.canonical.anti.targets.ind = ild.cor$Cor.fdr["MI0000743_MIMAT0000686",] < 0.25 & ild.cor$Cor["MI0000743_MIMAT0000686",] < 0 & targets[,mir.34c.canonical] > 0
mir34c.canonical.anti.targets.names = fData(gene)[mir34c.canonical.anti.targets.ind,"GENE_SYMBOL"]

write.table(setdiff(mir34c.isomir.anti.targets.names, mir34c.canonical.anti.targets.names), "MIR34C_ISOMIR_ONLY.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(setdiff(mir34c.canonical.anti.targets.names, mir34c.isomir.anti.targets.names), "MIR34C_CANONICAL_ONLY.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(c(as.character(mir34c.isomir.anti.targets.names), as.character(mir34c.canonical.anti.targets.names)), "MIR34C_UNION.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(intersect(mir34c.isomir.anti.targets.names, mir34c.canonical.anti.targets.names), "MIR34C_INTERSECT.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


i1 = unique(c(as.character(mir34c.isomir.anti.targets.names), as.character(mir34c.canonical.anti.targets.names)))
i2 = setdiff(as.character(mir34c.isomir.anti.targets.names), as.character(mir34c.canonical.anti.targets.names))
i3 = setdiff(as.character(mir34c.canonical.anti.targets.names), as.character(mir34c.isomir.anti.targets.names))
m = matrix(0, ncol=2, nrow=length(i1))
colnames(m) = c("Canonical", "IsomiR")
rownames(m) = i1
m[as.character(mir34c.canonical.anti.targets.names),1] = 1
m[as.character(mir34c.isomir.anti.targets.names),2] = 1

library(limma)
pdf("MIR34C_IsomiR_Target_Overlap.pdf", useDingbats=FALSE)
par(mfrow=c(2,2))

mir34c.isomir.targets.ind = targets[,mir.34c.isomir] > 0
mir34c.canonical.targets.ind = targets[,mir.34c.canonical] > 0
mir34c.all.targets.ind = mir34c.canonical.targets.ind | mir34c.isomir.targets.ind

plot(density(ild.cor$Cor["MI0000743_MIMAT0000686",!mir34c.all.targets.ind]), lwd=2, xlab="Spearman correlation coefficient", main="Distribution of gene correlations to miR-34c")
lines(density(ild.cor$Cor["MI0000743_MIMAT0000686",!mir34c.isomir.targets.ind & mir34c.canonical.targets.ind]), lwd=2, col="red")
lines(density(ild.cor$Cor["MI0000743_MIMAT0000686",mir34c.isomir.targets.ind & !mir34c.canonical.targets.ind]), lwd=2, col="blue")
lines(density(ild.cor$Cor["MI0000743_MIMAT0000686",mir34c.isomir.targets.ind & mir34c.canonical.targets.ind]), lwd=2, col="yellow")
legend("topright", c("IsomiR only", "Canonical only", "Both", "Neither"), title="Predicted targets", col=c("blue", "red", "yellow", "black"), lwd=2)

wilcox.test(ild.cor$Cor["MI0000743_MIMAT0000686",!mir34c.isomir.targets.ind & mir34c.canonical.targets.ind], ild.cor$Cor["MI0000743_MIMAT0000686",!mir34c.all.targets.ind])
wilcox.test(ild.cor$Cor["MI0000743_MIMAT0000686",mir34c.isomir.targets.ind & !mir34c.canonical.targets.ind], ild.cor$Cor["MI0000743_MIMAT0000686",!mir34c.all.targets.ind])
wilcox.test(ild.cor$Cor["MI0000743_MIMAT0000686",mir34c.isomir.targets.ind & mir34c.canonical.targets.ind], ild.cor$Cor["MI0000743_MIMAT0000686",!mir34c.all.targets.ind])

m = cbind(IsomiR=mir34c.isomir.anti.targets.ind+0, Canonical=mir34c.canonical.anti.targets.ind+0)
vc = vennCounts(m)
vennDiagram(vc)

dev.off()




## Summarize Counts/RPM of isoforms

iso.rpm = readRDS("../../Data/MicroRNA/LGRC_MicroRNA_Isoform_RPM_Log2_Present.rds")
iso.counts = readRDS("../../Data/MicroRNA/LGRC_MicroRNA_Isoform_Counts.rds")
fdata = fData(iso.rpm)


forward = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, fdata$Chrom, fdata$Start+1, fdata$Start+7, as.character=F))
reverse = as.character(reverseComplement(getSeq(BSgenome.Hsapiens.UCSC.hg19, fdata$Chrom, fdata$Start-7, fdata$Start-1, as.character=F)))

seed = ifelse(fdata$Strand == "+", forward, reverse)
seed = gsub("T", "U", seed)

fdata = data.frame(fdata, Seed=seed, Avg_RPM=rowMeans(exprs(iso.rpm)), Total_Count=rowSums(exprs(iso.counts)), stringsAsFactors=FALSE)
fdata.select = fdata[grep("miR-34c-5p|miR-34b-5p|miR-449a|miR-449b|miR-449c", fdata$Name),]

write.table(fdata, "Isoform_Seed_Summary.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(fdata.select, "Isoform_Seed_Summary_MIR34_449.txt", sep="\t", quote=FALSE, row.names=FALSE)

sessionInfo()
date()
