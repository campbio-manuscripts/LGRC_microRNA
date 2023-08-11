## Josh Campbell
## 8/28/2013
##
## Performs differential expression and clustering analysis on LGRC miRNA data

## Libraries and custom functions #########################################################################################################################################
library(affy)
library(MASS)
library(heatmap3)
library(flashClust)
library(vioplot)
library(limma)
library(sva)

## Blue, White, Red color scheme for heatmaps
bwr = colorRampPalette(c("blue", "white", "red"))(100)

## Average hierarchical clustering function
ha = function(c) hclust(c, method="average")

## Read ExpressionSets which have data for Differential expression and clustering analysis
data.counts = readRDS("../../Data/MicroRNA/LGRC_MicroRNA_Counts.rds")
data.rpm = readRDS("../../Data/MicroRNA/LGRC_MicroRNA_RPM_Log2_Batch_Corrected.rds")

## Set up covariates
age = data.counts$Age
gender = as.factor(data.counts$Gender)
fev1 = as.numeric(data.counts$FEV1_Percent)
fev1.fvc = as.numeric(data.counts$FEV1_FVC)
dlco = as.numeric(data.counts$DLCO)
bode = as.numeric(data.counts$BODE)
emphysema = as.numeric(data.counts$Emphysema)
reads = data.counts$Mirna_Reads / 1e6
smoke = as.factor(data.counts$Smoke)
protocol = as.factor(data.counts$Protocol)
status = as.factor(data.counts$Disease)
flowcell = as.factor(data.counts$Flow_Cell)
index = as.factor(data.counts$Index)
reads.aligned = as.numeric(data.counts$Reads_aligned) / 1e6
copd = as.factor(data.counts$COPD)
ild = as.factor(data.counts$ILD)
disease = as.factor(as.numeric(data.counts$ILD | data.counts$COPD))

status.col = rep("grey25", length(status))
status.col[status == 1] = "purple"
status.col[status == 2] = "yellow"

fev1.fvc[fev1.fvc > 1.1] = NA  ## Set 2 patients with abnormally high FEV1/FVC to NA

## Perform differential expression analysis using NB linear models
results.anova <- matrix(NA, nrow=nrow(data.counts), ncol=16)
results.summary <- matrix(NA, nrow=nrow(data.counts), ncol=36)
for(i in 1:nrow(data.counts)) {
	m1 <- glm.nb(exprs(data.counts)[i,] ~ reads + protocol + smoke + gender + age + status, maxit=10000)
    m2 <- glm.nb(exprs(data.counts)[i,] ~ reads + protocol + smoke + gender + age, maxit=10000)

	results.anova[i,] = unlist(anova(m1, m2))
	results.summary[i,] = c(summary(m1)$coefficients)
}
results.anova.fdr <- p.adjust(results.anova[,16], "fdr")
results = cbind(results.summary, results.anova[,16], results.anova.fdr)

cnames = c("GLM_Intercept", "GLM_Aligned_Reads", "GLM_Protocol", "GLM_Smoke_Former", "GLM_Smoke_Never", "GLM_Gender", "GLM_Age", "GLM_ILD", "GLM_COPD")
colnames(results) = c(paste(rep(cnames, 4), rep(c("Coef", "Std_Error", "Zvalue", "Pvalue"), each=9), sep="_"), "GLM_ANOVA_Pvalue", "GLM_ANOVA_FDR")

sig.ind = (abs(results[,"GLM_ILD_Coef"]) > log(1.25) | abs(results[,"GLM_COPD_Coef"]) > log(1.25)) & results[,"GLM_ANOVA_FDR"] < 0.1
print(sum(sig.ind))
######################################################################################################################
#### Run Consensus Clustering using ConsensusClusterPlus
######################################################################################################################

source("cc_with_plots.R")
cc_with_plots(t(scale(t(exprs(data.rpm)[sig.ind,]))), maxK=15, output.prefix="CC_Samples", clusterAlg="pam", distance="euclidean", pFeature=1, reps=1000)
cc_with_plots(scale(t(exprs(data.rpm)[sig.ind,])), maxK=15, output.prefix="CC_Mirna", clusterAlg="pam", distance="euclidean", pFeature=1, reps=1000)
cc.sample = readRDS("CC_Samples_results/CC_results.rds")
cc.mir = readRDS("CC_Mirna_results/CC_results.rds")

cc.sample.cluster.temp = cc.sample[[5]][["consensusClass"]]
sample.cluster = rep(NA, length(cc.sample.cluster.temp))
sample.cluster[cc.sample.cluster.temp == 4] = 1
sample.cluster[cc.sample.cluster.temp == 5] = 2
sample.cluster[cc.sample.cluster.temp == 3] = 3
sample.cluster[cc.sample.cluster.temp == 2] = 4
sample.cluster[cc.sample.cluster.temp == 1] = 5

sample.cluster.f = sample.cluster
#sample.cluster.f[sample.cluster == 2] = 1
sample.cluster.f = as.factor(sample.cluster.f)


cc.mir.cluster.temp = cc.mir[[4]][["consensusClass"]]
cc.mir.cluster = rep(NA, length(cc.mir.cluster.temp))
cc.mir.cluster[cc.mir.cluster.temp == 4] = 1
cc.mir.cluster[cc.mir.cluster.temp == 2] = 2
cc.mir.cluster[cc.mir.cluster.temp == 1] = 4
cc.mir.cluster[cc.mir.cluster.temp == 3] = 3

pdf("Cluster_Heatmap.pdf", useDingbats=FALSE)
heatmap3(exprs(data.rpm)[sig.ind,], col=bwr, ColSideColors=cbind(sample.cluster, status.col), RowSideColors=cc.mir.cluster, col.clustering="semi", row.clustering="semi", col.dendrogram=FALSE, row.dendrogram=FALSE, labRow=FALSE, labCol=FALSE, hclustfun=ha)
heatmap3(exprs(data.rpm)[sig.ind,], col=bwr, ColSideColors=cbind(status.col, sample.cluster), RowSideColors=cc.mir.cluster, col.clustering="semi", row.clustering="semi", col.dendrogram=FALSE, row.dendrogram=FALSE, labRow=FALSE, labCol=FALSE, hclustfun=ha)
heatmap3(exprs(data.rpm)[sig.ind,], col=bwr, ColSideColors=cbind(sample.cluster, status.col), RowSideColors=cc.mir.cluster, col.clustering="un", row.clustering="un", col.dendrogram=FALSE, row.dendrogram=FALSE, labRow=FALSE, labCol=FALSE, hclustfun=ha)
dev.off()

pdf("Sample_Barbplot.pdf", useDingbats=FALSE)
ta = table(status, sample.cluster)
ta.norm = sweep(ta, 2, colSums(ta), "/")
par(mfrow=c(2,2))
barplot(ta, col=c("grey25", "purple", "yellow"), xlab="Cluster", ylab="Number of samples")
barplot(ta.norm, col=c("grey25", "purple", "yellow"), xlab="Cluster", ylab="Number of samples")
dev.off()

## Write out new RDS files with batch information
data.rpm.b = data.rpm
pData(data.rpm.b) = data.frame(pData(data.rpm), Sample_Cluster=sample.cluster, check.names=FALSE, stringsAsFactors=FALSE)

all.mir.cluster = rep(NA, nrow(data.rpm))
all.mir.cluster[sig.ind] = cc.mir.cluster
fData(data.rpm.b) = data.frame(fData(data.rpm), results, Mirna_Cluster=all.mir.cluster, check.names=FALSE, stringsAsFactors=FALSE)
saveRDS(data.rpm.b, "LGRC_MicroRNA_Isoform_RPM_Log2_Batch_Corrected_With_Cluster.rds")

## Test to verify that proportions of Disease/COPD/ILD samples are different between clusters
summary(glm(disease ~ reads + protocol + smoke + gender + age + sample.cluster.f, family="binomial"))
summary(glm(copd ~ reads + protocol + smoke + gender + age + sample.cluster.f, family="binomial"))
summary(glm(ild ~ reads + protocol + smoke + gender + age + sample.cluster.f, family="binomial"))

## Test for association of clinical variables with cluster status, both within and across disease types
copd.ind = copd == 1
ild.ind = ild == 1
disease.ind = disease == 1

# FEV1 Percent Predicted
summary(lm(fev1[disease.ind] ~ reads[disease.ind] + protocol[disease.ind] + smoke[disease.ind] + gender[disease.ind] + age[disease.ind] + sample.cluster.f[disease.ind]))
summary(lm(fev1[copd.ind] ~ reads[copd.ind] + protocol[copd.ind] + smoke[copd.ind] + gender[copd.ind] + age[copd.ind] + sample.cluster.f[copd.ind]))
summary(lm(fev1[ild.ind] ~ reads[ild.ind] + protocol[ild.ind] + smoke[ild.ind] + gender[ild.ind] + age[ild.ind] + sample.cluster.f[ild.ind]))
anova(lm(fev1[disease.ind] ~ reads[disease.ind] + protocol[disease.ind] + smoke[disease.ind] + gender[disease.ind] + age[disease.ind] + sample.cluster.f[disease.ind]))
anova(lm(fev1[copd.ind] ~ reads[copd.ind] + protocol[copd.ind] + smoke[copd.ind] + gender[copd.ind] + age[copd.ind] + sample.cluster.f[copd.ind]))
anova(lm(fev1[ild.ind] ~ reads[ild.ind] + protocol[ild.ind] + smoke[ild.ind] + gender[ild.ind] + age[ild.ind] + sample.cluster.f[ild.ind]))
kruskal.test(fev1[disease.ind] ~ sample.cluster.f[disease.ind])
kruskal.test(fev1[copd.ind] ~ sample.cluster.f[copd.ind])
kruskal.test(fev1[ild.ind] ~ sample.cluster.f[ild.ind])


# DLCO Percent Predicted
summary(lm(dlco[disease.ind] ~ reads[disease.ind] + protocol[disease.ind] + smoke[disease.ind] + gender[disease.ind] + age[disease.ind] + sample.cluster.f[disease.ind]))
summary(lm(dlco[copd.ind] ~ reads[copd.ind] + protocol[copd.ind] + smoke[copd.ind] + gender[copd.ind] + age[copd.ind] + sample.cluster.f[copd.ind]))
summary(lm(dlco[ild.ind] ~ reads[ild.ind] + protocol[ild.ind] + smoke[ild.ind] + gender[ild.ind] + age[ild.ind] + sample.cluster.f[ild.ind]))
anova(lm(dlco[disease.ind] ~ reads[disease.ind] + protocol[disease.ind] + smoke[disease.ind] + gender[disease.ind] + age[disease.ind] + sample.cluster.f[disease.ind]))
anova(lm(dlco[copd.ind] ~ reads[copd.ind] + protocol[copd.ind] + smoke[copd.ind] + gender[copd.ind] + age[copd.ind] + sample.cluster.f[copd.ind]))
anova(lm(dlco[ild.ind] ~ reads[ild.ind] + protocol[ild.ind] + smoke[ild.ind] + gender[ild.ind] + age[ild.ind] + sample.cluster.f[ild.ind]))
kruskal.test(dlco[disease.ind] ~ sample.cluster.f[disease.ind])
kruskal.test(dlco[copd.ind] ~ sample.cluster.f[copd.ind])
kruskal.test(dlco[ild.ind] ~ sample.cluster.f[ild.ind])

# BODE 
summary(lm(bode[disease.ind] ~ reads[disease.ind] + protocol[disease.ind] + smoke[disease.ind] + gender[disease.ind] + age[disease.ind] + sample.cluster.f[disease.ind]))
summary(lm(bode[copd.ind] ~ reads[copd.ind] + protocol[copd.ind] + smoke[copd.ind] + gender[copd.ind] + age[copd.ind] + sample.cluster.f[copd.ind]))
summary(lm(bode[ild.ind] ~ reads[ild.ind] + protocol[ild.ind] + smoke[ild.ind] + gender[ild.ind] + age[ild.ind] + sample.cluster.f[ild.ind]))
anova(lm(bode[disease.ind] ~ reads[disease.ind] + protocol[disease.ind] + smoke[disease.ind] + gender[disease.ind] + age[disease.ind] + sample.cluster.f[disease.ind]))
anova(lm(bode[copd.ind] ~ reads[copd.ind] + protocol[copd.ind] + smoke[copd.ind] + gender[copd.ind] + age[copd.ind] + sample.cluster.f[copd.ind]))
anova(lm(bode[ild.ind] ~ reads[ild.ind] + protocol[ild.ind] + smoke[ild.ind] + gender[ild.ind] + age[ild.ind] + sample.cluster.f[ild.ind]))
kruskal.test(bode[disease.ind] ~ sample.cluster.f[disease.ind])
kruskal.test(bode[copd.ind] ~ sample.cluster.f[copd.ind])
kruskal.test(bode[ild.ind] ~ sample.cluster.f[ild.ind])

# FEV1/FVC 
summary(lm(fev1.fvc[disease.ind] ~ reads[disease.ind] + protocol[disease.ind] + smoke[disease.ind] + gender[disease.ind] + age[disease.ind] + sample.cluster.f[disease.ind]))
summary(lm(fev1.fvc[copd.ind] ~ reads[copd.ind] + protocol[copd.ind] + smoke[copd.ind] + gender[copd.ind] + age[copd.ind] + sample.cluster.f[copd.ind]))
summary(lm(fev1.fvc[ild.ind] ~ reads[ild.ind] + protocol[ild.ind] + smoke[ild.ind] + gender[ild.ind] + age[ild.ind] + sample.cluster.f[ild.ind]))
anova(lm(fev1.fvc[disease.ind] ~ reads[disease.ind] + protocol[disease.ind] + smoke[disease.ind] + gender[disease.ind] + age[disease.ind] + sample.cluster.f[disease.ind]))
anova(lm(fev1.fvc[copd.ind] ~ reads[copd.ind] + protocol[copd.ind] + smoke[copd.ind] + gender[copd.ind] + age[copd.ind] + sample.cluster.f[copd.ind]))
anova(lm(fev1.fvc[ild.ind] ~ reads[ild.ind] + protocol[ild.ind] + smoke[ild.ind] + gender[ild.ind] + age[ild.ind] + sample.cluster.f[ild.ind]))
kruskal.test(fev1.fvc[disease.ind] ~ sample.cluster.f[disease.ind])
kruskal.test(fev1.fvc[copd.ind] ~ sample.cluster.f[copd.ind])
kruskal.test(fev1.fvc[ild.ind] ~ sample.cluster.f[ild.ind])

# Percent Emphysema 
summary(lm(emphysema[disease.ind] ~ reads[disease.ind] + protocol[disease.ind] + smoke[disease.ind] + gender[disease.ind] + age[disease.ind] + sample.cluster.f[disease.ind]))
summary(lm(emphysema[copd.ind] ~ reads[copd.ind] + protocol[copd.ind] + smoke[copd.ind] + gender[copd.ind] + age[copd.ind] + sample.cluster.f[copd.ind]))
summary(lm(emphysema[ild.ind] ~ reads[ild.ind] + protocol[ild.ind] + smoke[ild.ind] + gender[ild.ind] + age[ild.ind] + sample.cluster.f[ild.ind]))
anova(lm(emphysema[disease.ind] ~ reads[disease.ind] + protocol[disease.ind] + smoke[disease.ind] + gender[disease.ind] + age[disease.ind] + sample.cluster.f[disease.ind]))
anova(lm(emphysema[copd.ind] ~ reads[copd.ind] + protocol[copd.ind] + smoke[copd.ind] + gender[copd.ind] + age[copd.ind] + sample.cluster.f[copd.ind]))
anova(lm(emphysema[ild.ind] ~ reads[ild.ind] + protocol[ild.ind] + smoke[ild.ind] + gender[ild.ind] + age[ild.ind] + sample.cluster.f[ild.ind]))
kruskal.test(emphysema[disease.ind] ~ sample.cluster.f[disease.ind])
kruskal.test(emphysema[copd.ind] ~ sample.cluster.f[copd.ind])
kruskal.test(emphysema[ild.ind] ~ sample.cluster.f[ild.ind])



# Other QC metrics
summary(lm(reads ~ sample.cluster.f))
anova(lm(reads ~ sample.cluster.f))
kruskal.test(reads ~ sample.cluster.f)

summary(lm(reads.aligned ~ sample.cluster.f))
anova(lm(reads.aligned ~ sample.cluster.f))
kruskal.test(reads.aligned ~ sample.cluster.f)


table(protocol, sample.cluster.f)
fisher.test(protocol, sample.cluster.f)



stripbox = function(e, c, group.names, ...) {
  xlim = c(0.5, length(unique(c))+0.5)
  boxplot(e ~ c, border="black", lwd=2, boxwex=0.75, outline=FALSE, ..., yaxt="n")
  stripchart(e ~ c, add=TRUE, vertical=TRUE, pch=19, xlim=xlim, method="jitter", jitter=0.3, col="grey50", yaxt="n")
  axis(2, las=1)
}  

pdf("Sample_Cluster_Assocations.pdf", useDingbats=FALSE)
par(mfrow=c(3,3), mar=c(2,4,2,2))
stripbox(dlco[ild.ind], sample.cluster.f[ild.ind], group.names=1:5, ylab="DLCO", col="purple")
stripbox(dlco[copd.ind], sample.cluster.f[copd.ind], group.names=1:5, ylab="DLCO", col="yellow")
stripbox(fev1[copd.ind], sample.cluster.f[copd.ind], group.names=1:5, ylab="FEV1 Percent Predicted", col="yellow")
stripbox(fev1.fvc[copd.ind], sample.cluster.f[copd.ind], group.names=1:5, ylab="FEV1/FVC", col="yellow")
stripbox(emphysema[copd.ind], sample.cluster.f[copd.ind], group.names=1:5, ylab="Percent Emphysema", col="yellow")
stripbox(bode[copd.ind], sample.cluster.f[copd.ind], group.names=1:5, ylab="BODE", col="yellow")
dev.off()






## Perform Gene Differential Expression
gene = readRDS("../../Data/Gene/LGRC_Gene_Overlap_MicroRNA.rds")
egene = exprs(gene)
pgene = pData(gene)
fgene = fData(gene)

smoke.na = is.na(pgene$Smoke)
gene.model1 = model.matrix(~ as.factor(Smoke) + as.factor(Gender) + as.numeric(Age), data=pgene[!smoke.na,])
gene.model2 = model.matrix(~ as.factor(Smoke) + as.factor(Gender) + as.numeric(Age) + as.factor(Disease), data=pgene[!smoke.na,])

res.gene.anova = f.pvalue(egene[,!smoke.na], gene.model2, gene.model1)
res.gene.anova.fdr = p.adjust(res.gene.anova, 'fdr')

res.gene.fit = lmFit(egene[,!smoke.na], gene.model2)
res.gene.e = ebayes(res.gene.fit)
res.gene.e.fdr = apply(res.gene.e$p.value, 2, p.adjust, 'fdr')

gene.results = cbind(res.gene.fit$coefficients, res.gene.fit$stdev.unscaled, res.gene.e$t, res.gene.e$p.value, res.gene.anova, res.gene.anova.fdr)

cnames = c("LM_Intercept", "LM_Smoke_Former", "LM_Smoke_Never", "LM_Gender", "LM_Age", "LM_ILD", "LM_COPD")
colnames(gene.results) = c(paste(rep(cnames, 4), rep(c("Coef", "Std_Dev_Unscaled", "Tvalue", "Pvalue"), each=7), sep="_"), "LM_ANOVA_Pvalue", "LM_ANOVA_FDR")

fData(gene) = cbind(fgene, gene.results)
saveRDS(gene, "../../Data/Gene/LGRC_Gene_Overlap_MicroRNA_with_Diff_Exp.rds")


sessionInfo()
date()
