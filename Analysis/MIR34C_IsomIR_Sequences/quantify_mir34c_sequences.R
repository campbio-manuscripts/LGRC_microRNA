library(MASS)
library(affy)
library(ggplot2)

seq.all <- read.table(gzfile("mir34c_sequence_counts.txt.gz"), row.names=1, header=TRUE)
data.counts <- readRDS("../Differential_Expression/LGRC_MicroRNA_Isoform_RPM_Log2_Batch_Corrected_With_Cluster.rds")
seq.qc <- seq.all[,colnames(data.counts)]

# Calculate RPM and total
reads.aligned = data.counts$Mirna_Reads / 1e6
seq.rpm = sweep(seq.qc, 2, reads.aligned, "/")
seq.rpm.log2 = log2(seq.rpm+1)
seq.total <- rowSums(seq.qc)

# Reorder by average RPM and select the top 25
seq.rpm.mean <- apply(seq.rpm.log2, 1, mean)
o <- order(seq.rpm.mean, decreasing = TRUE)[1:25]
seq.counts <- seq.qc[o,]
seq.rpm <- seq.rpm[o,]
seq.rpm.log2 <- seq.rpm.log2[o,]
seq.rpm.mean <- seq.rpm.mean[o]
seq.total <- seq.total[o]

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

fev1.fvc[fev1.fvc > 1.1] = NA  ## Set 2 patients with abnormally high FEV1/FVC to NA

## Perform differential expression analysis using NB linear models
results.anova <- matrix(NA, nrow=nrow(seq.counts), ncol=12)
results.summary <- matrix(NA, nrow=nrow(seq.counts), ncol=36)
for(i in 1:nrow(seq.counts)) {
  try({
	m1 <- glm.nb(as.numeric(seq.counts[i,]) ~ reads + protocol + smoke + gender + age + status, maxit=10000)
    m2 <- glm.nb(as.numeric(seq.counts[i,]) ~ reads + protocol + smoke + gender + age, maxit=10000)

	results.anova[i,] = unlist(anova(m1, m2)[,-c(1, 5)])  # Get rid of 1st and 5th column to make sure it stays numeric
	results.summary[i,] = c(summary(m1)$coefficients)
  }, silent = TRUE)	
}
results.anova.fdr <- p.adjust(results.anova[,12], "fdr")
results = cbind(results.summary, results.anova[,12], results.anova.fdr)

cnames = c("GLM_Intercept", "GLM_Aligned_Reads", "GLM_Protocol", "GLM_Smoke_Former", "GLM_Smoke_Never", "GLM_Gender", "GLM_Age", "GLM_ILD", "GLM_COPD")
colnames(results) = c(paste(rep(cnames, 4), rep(c("Coef", "Std_Error", "Zvalue", "Pvalue"), each=9), sep="_"), "GLM_ANOVA_Pvalue", "GLM_ANOVA_FDR")

# See how many miRNA sequences were significantly associated with ILD
sig.ind = (abs(results[,"GLM_ILD_Coef"]) > log(1.25) | abs(results[,"GLM_COPD_Coef"]) > log(1.25)) & results[,"GLM_ANOVA_FDR"] < 0.1
print(sum(sig.ind, na.rm = TRUE))

# Create results table for supplement
seed <- substring(rownames(seq.counts), 2, 8)
canonical <- rep(NA, length(seed))
canonical[seed == "GGCAGTG"] <- "5p"
canonical[seed == "ATCACTA"] <- "3p"

results.df <- data.frame(Sequence = rownames(seq.counts), Seed = seed, Canonical_Seed = canonical, Average_RPM = seq.rpm.mean, Total = seq.total, results)
write.table(results.df, "mir34c_top25_DE_results.txt", row.names=FALSE, sep="\t", quote=FALSE)


# Plot results in a bar plot
color <- seed
color[seed == "GGCAGTG"] <- "GGCAGTG (5p Canonical)"
color[seed == "ATCACTA"] <- "ATCACTA (3p Canonical)"
color[seed == "AGGCAGT"] <- "AGGCAGT (5p isomiR)"
color[seed == "GCAGTGT"] <- "GCAGTGT (5p isomiR)"
color[seed == "TCACTAA"] <- "TCACTAA (3p isomiR)"
color <- factor(color, levels = c("GGCAGTG (5p Canonical)", "ATCACTA (3p Canonical)", "AGGCAGT (5p isomiR)", "GCAGTGT (5p isomiR)", "TCACTAA (3p isomiR)"))

# Create sequences with appropriate spacing
sequences <- rownames(seq.counts)
sequences[seed != "AGGCAGT"] <- paste0(" ", sequences[seed != "AGGCAGT"])   # Add space to all other sequences except 5' isomiR
sequences[seed == "GCAGTGT"] <- paste0(" ", sequences[seed == "GCAGTGT"])   # Add another space to other 5' isomiR
sequences[seed == "TCACTAA"] <- paste0(" ", sequences[seed == "TCACTAA"])   # Add another space to 3' isomiR
sequences <- factor(sequences, levels = rev(sequences))

pdf("mir34c_top25_expression.pdf", useDingbats=FALSE)
df <- data.frame(Sequence = sequences, Expression=results.df$Average_RPM, Index=nrow(results.df):1, color=color)
ggplot(df, aes(x = Expression, y = Sequence, fill=color)) + geom_col() + theme_bw() +
	scale_fill_brewer(palette="Dark2") +
	guides(color=guide_legend(title="Seed")) +
	ylab("Top 25 hsa-miR-34c sequences") + xlab("Expression (Average RPM)") +
	scale_x_continuous(expand = c(0, 0)) + 
	theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(hjust = 0, family="mono", size = 12), panel.grid=element_blank(), legend.position = c(0.75, 0.15)) +
	guides(fill=guide_legend("Seed"))
dev.off()




sessionInfo()
