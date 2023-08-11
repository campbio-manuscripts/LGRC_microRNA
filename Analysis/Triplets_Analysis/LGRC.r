
module load R/R-3.0.0_gnu-4.4.6

# install
R --vanilla

library(devtools)
install_github("plink2R", username="gabraham", subdir="plink2R")
install_github("flashpca", username="gabraham", subdir="flashpcaR")


# run 
R --vanilla

library(plink2R)

plinkfile <- "/restricted/projectnb/lgrc/SNP/Data/LGRC_SNP_full_cohort"
system.time({r <- read_plink(plinkfile)})

Error in .Call("read_plink", PACKAGE = "plink2R", bedfile, famfile, impute_int,  :
  negative length vectors are not allowed
Timing stopped at: 45.073 14.602 59.806


dim(r$bed)
dim(r$bim)
dim(r$fam)

X = r$bed# TODO: check
library(flashpcaR)
pcas <- flashpca(X, do_loadings=TRUE, verbose=TRUE, stand="binom", ndim=10, nextra=100) # stand = "center" -> compatible with R prcomp; default: binom


####################################

module load boost/1.57.0

in sandbox directory: /restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/

wget http://bitbucket.org/eigen/eigen/get/3.2.4.tar.gz
mv 3.2.4.tar.gz eigen-3.2.4.tar.gz
tar -xvf eigen-3.2.4.tar.gz

git clone git://github.com/gabraham/flashpca
cd flashpca

edit: Makefile -> set the following lines
EIGEN_INC=/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/eigen-eigen-10219c95fe65/
BOOST_INC=/share/pkg/boost/1.57.0/install/include/
BOOST_LIB=/share/pkg/boost/1.57.0/install/lib/

in shell:
make all

#####################################

module load boost/1.57.0
module load plink/1.90b3b

mkdir LGRC_Data
cd LGRC_Data
cp /restricted/projectnb/lgrc/SNP/Data/LGRC_SNP_full_cohort.* ./

cp ../flashpca/exclusion_regions.txt ./
plink --bfile LGRC_SNP_full_cohort --indep-pairwise 1000 50 0.05 --exclude range exclusion_regions.txt
plink --bfile LGRC_SNP_full_cohort --extract plink.prune.in --make-bed --out LGRC_SNP_full_cohort_pruned

#plink --bed LGRC_SNP_Overlap_miR_Seq_Gene_Array.txt --indep-pairwise 1000 50 0.05 --exclude range exclusion_regions.txt
#plink --bed LGRC_SNP_Overlap_miR_Seq_Gene_Array.txt --extract plink.prune.in --make-bed --out LGRC_SNP_Overlap_miR_Seq_Gene_Array.txt
#../flashpca/flashpca --bed LGRC_SNP_Overlap_miR_Seq_Gene_Array.txt --numthreads 8

../flashpca/flashpca --bfile LGRC_SNP_full_cohort_pruned --numthreads 8
Output:
[Mon Jun 15 23:19:32 2015] arguments: flashpca ../flashpca/flashpca --bfile LGRC_SNP_full_cohort_pruned --numthreads 8
[Mon Jun 15 23:19:32 2015] Start flashpca (version 1.2)
[Mon Jun 15 23:19:32 2015] Using 8 OpenMP threads
[Mon Jun 15 23:19:32 2015] seed: 1
[Mon Jun 15 23:19:32 2015] Detected pheno file LGRC_SNP_full_cohort_pruned.fam, 1274 samples
[Mon Jun 15 23:19:32 2015] Analyzing BED file 'LGRC_SNP_full_cohort_pruned.bed', found 65679 SNPs
[Mon Jun 15 23:19:32 2015] Reading BED file 'LGRC_SNP_full_cohort_pruned.bed'
[Mon Jun 15 23:19:32 2015] Detected BED file: LGRC_SNP_full_cohort_pruned.bed with 20951601 bytes, 1274 samples, 65679 SNPs.
[Mon Jun 15 23:19:34 2015] Loaded genotypes: 1274 samples, 65679 SNPs
[Mon Jun 15 23:19:34 2015] PCA begin
[Mon Jun 15 23:19:58 2015] PCA done
[Mon Jun 15 23:19:58 2015] Writing 10 eigenvalues to file eigenvalues.txt
[Mon Jun 15 23:19:58 2015] Writing 10 eigenvectors to file eigenvectors.txt
[Mon Jun 15 23:19:58 2015] Writing 10 PCs to file pcs.txt
[Mon Jun 15 23:19:58 2015] Writing 10 proportion variance explained to file pve.txt
[Mon Jun 15 23:19:58 2015] Goodbye!


module load R/R-3.0.0_gnu-4.4.6
R CMD BATCH ../flashpca/plot.R
Output file: pcs.png



####New analysis 11/28/2015
cd /restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/
module load boost/1.57.0
module load plink/1.90b3b

Extract samples
plink --bfile LGRC_SNP_full_cohort --keep mylist.txt --make-bed --out LGRC_SNP_295_cohort
plink --bfile LGRC_SNP_295_cohort --list

plink --bfile LGRC_SNP_295_cohort --indep-pairwise 1000 50 0.05 --exclude range exclusion_regions.txt
plink --bfile LGRC_SNP_295_cohort --extract plink.prune.in --make-bed --out LGRC_SNP_295_cohort_pruned
../flashpca/flashpca --bfile LGRC_SNP_295_cohort_pruned --numthreads 8
module load R/R-3.0.0_gnu-4.4.6
R CMD BATCH ../flashpca/plot.R
outliers: samples 10, 236

plink --bfile LGRC_SNP_full_cohort --recodeAD


module load boost/1.57.0
module load plink/1.90b3b

plink --bfile LGRC_SNP_full_cohort --keep mylist.txt --make-bed --out LGRC_SNP_293_cohort

plink --bfile LGRC_SNP_293_cohort --indep-pairwise 1000 50 0.05 --exclude range exclusion_regions.txt
plink --bfile LGRC_SNP_293_cohort --extract plink.prune.in --make-bed --out LGRC_SNP_293_cohort_pruned
../flashpca/flashpca --bfile LGRC_SNP_293_cohort_pruned --numthreads 8
module load R/R-3.0.0_gnu-4.4.6
R CMD BATCH ../flashpca/plot.R

The SNP outliers:
236 LT234774LU
10  LT007392RU

plink --bfile LGRC_SNP_full_cohort --keep mylist_264.txt --make-bed --out LGRC_SNP_264_cohort
plink --bfile LGRC_SNP_264_cohort --list

plink --bfile LGRC_SNP_264_cohort --indep-pairwise 1000 50 0.05 --exclude range exclusion_regions.txt
plink --bfile LGRC_SNP_264_cohort --extract plink.prune.in --make-bed --out LGRC_SNP_264_cohort_pruned
../flashpca/flashpca --bfile LGRC_SNP_264_cohort_pruned --numthreads 8
module load R/R-3.0.0_gnu-4.4.6
R CMD BATCH ../flashpca/plot.R


#TODO
/restricted/projectnb/lgrc/small_RNA/lgrc_365/20141014_EQTL/Control/
/restricted/projectnb/lgrc/small_RNA/lgrc_365/Preprocess/
In Preprocess.R filter the 2 outliers and re-run
In MatrixEQTL add the PCA covariates

###MatrixEQTL.R
# Matrix eQTL 
library(affy)
library(MatrixEQTL)

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelANOVA; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = "../../Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt"; 
snps_location_file_name = "../../../Paper_Analysis/Data/SNP/LGRC_SNP_Locations.txt";

# MicroRNA expression file name
mir_expression_file_name = "../../Preprocess/MicroRNA/LGRC_MicroRNA_Overlap_SNP_Control_Samples.txt";
mir_expression_rds_name = "../../Preprocess/MicroRNA/LGRC_MicroRNA_Overlap_SNP_All_Samples.rds";
mir_location_file_name = "../../../Paper_Analysis/Data/MicroRNA/LGRC_MicroRNA_Locations.txt";

# Gene expression file name
gene_expression_file_name = "../../Preprocess/Gene/LGRC_Gene_Overlap_SNP_Control_Samples.txt";
gene_expression_rds_name = "../../Preprocess/Gene/LGRC_Gene_Overlap_SNP_All_Samples.rds";
gene_location_file_name = "../../../Paper_Analysis/Data/Gene/LGRC_Gene_Locations.txt";

# Covariates file name
# Set covariates_file_name to character() for no covariates
pdata = pData(readRDS(mir_expression_rds_name))
pdata$Gender = as.factor(pdata$Gender)
pdata$Age = as.numeric(pdata$Age)
pdata$Smoke = as.factor(pdata$Smoke)
mir.edata = read.table(mir_expression_file_name, header=TRUE, row.names=1, sep="\t")
gene.edata = read.table(gene_expression_file_name, header=TRUE, row.names=1, sep="\t")

demo = t(model.matrix(~ Smoke + Age + Gender, data=pdata[colnames(mir.edata),])[,-1])
write.table(demo, "covar.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
covariates_file_name = "covar.txt";

## Combine Mirna and Gene expression files
combined.expr = rbind(mir.edata, gene.edata)
write.table(combined.expr, "expression.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
expression_file_name = "expression.txt"

# Output file name
output_file_name_cis = "eQTL_Cis.txt";
output_file_name_tra = "eQTL_Trans.txt";

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-2;
pvOutputThreshold_tra = 1e-4;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Transform gene expresssion data to ranks for Kruskal-Wallis
#for( sl in 1:length(gene) ) {
#  mat = gene[[sl]];
#  mat = t(apply(mat, 1, rank, ties.method = "average"));
#  gene[[sl]] = mat;
#}
#rm(sl, mat);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Get location data
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
mirpos = read.table(mir_location_file_name, header = TRUE, stringsAsFactors = FALSE);
exprpos = rbind(mirpos, genepos)
exprpos = exprpos[!is.na(exprpos[,3]),]

## Run the analysis
me = Matrix_eQTL_main(
    snps = snps, 
    gene = gene, 
    cvrt = cvrt,
    output_file_name      = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE, 
    output_file_name.cis  = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos, 
    genepos = exprpos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');

## Make the histogram of local and distant p-values
pdf("eQTL_histogram.pdf", useDingbats=FALSE)
plot(me)
dev.off()

system("gzip eQTL_Cis.txt")
system("gzip eQTL_Trans.txt")
sessionInfo()


/protected/projects/pulmarray/HoggCOPD/A1ATGrifols/experiments/Expmt_0033_PhaseII_Immune_GSVA



/protected/projects/lgrc/small_RNA/lgrc_365/Integrated_Networks2/ILD/run_cit.R
/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/
cov = read.table('covar.txt', header=TRUE, sep="\t")
ncol(cov)

#Count associations
cis = read.table('eQTL_Cis_corr.txt.gz', header=TRUE, sep="\t")
cis_gene <- cis[which(cis$FDR<0.05 & grepl("A", as.character(cis$gene))),]
nrow(cis_gene)
length(unique(as.character(cis_gene$gene)))

#cis_all <- cis[which(cis$FDR<0.05),]

cis_mir <- cis[which(cis$FDR<0.05 & grepl("M", as.character(cis$gene))),]
cis_mir_all <- cis[which(grepl("M", as.character(cis$gene))),]
nrow(cis_mir)
length(unique(as.character(cis_mir$gene)))

length(unique(as.character(cis_gene$SNP)))+length(unique(as.character(cis_mir$SNP)))

trans = read.table('eQTL_Trans_corr.txt.gz', header=TRUE, sep="\t")

triplets_ILD = read.table('Triplets_eQTL_with_CIT_COPD.txt', header=TRUE, sep="\t")


trans_gene <- trans[which(trans$FDR<0.05 & grepl("A", as.character(trans$gene))),]
nrow(trans_gene)
length(unique(as.character(trans_gene$gene)))

#trans_all <- trans[which(trans$FDR<0.05),]

trans_mir <- trans[which(trans$FDR<0.05 & grepl("M", as.character(trans$gene))),]
nrow(trans_mir)
length(unique(as.character(trans_mir$gene)))

length(unique(as.character(trans_gene$SNP)))+length(unique(as.character(trans_mir$SNP)))


pairs <- read.table('Triplets_eQTL_ILD.txt', header=TRUE, sep="\t")
mirgene <- pairs[which(pairs$Mir_Gene_Cor_FDR<0.1),]
triangle <- pairs[which(pairs$Mir_Gene_Cor_FDR<0.1 & pairs$Mir_P<1e-04 & pairs$Gene_P<1e-04),]
nrow(pairs)
nrow(triangle)

triplets <- read.table('Triplets_eQTL_with_CIT_Control.txt', header=TRUE, sep="\t")
nrow(triplets)

triplets2 <- matrix(ncol=ncol(triplets), nrow=0)
colnames(triplets2) <- colnames(triplets)
triplets2 <- rbind(triplets2, triplets[1,])

for (i in 2: nrow(triplets)){
    if ((as.character(triplets[i-1, "SNP"])!=as.character(triplets[i, "SNP"])) |
    (as.character(triplets[i-1, "Mir"])!=as.character(triplets[i, "Mir"])) |
    (as.character(triplets[i-1, "Gene"])!=as.character(triplets[i, "Gene"])))
    {
        triplets2 <- rbind(triplets2, triplets[i,])
        #print(i)
    }
}
write.table(triplets2, "triplets_all_unique_Control.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

triplets_cit <- triplets2[which(triplets2$Causal_Call==1),]
nrow(triplets_cit)
triplets_cit_ord <- order(triplets_cit$P_Causal)
write.table(triplets_cit[triplets_cit_ord,], "triplets_cit_ord_unique_Control.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)



#unique mir/gene
triplets2 <- matrix(ncol=ncol(triplets), nrow=0)
colnames(triplets2) <- colnames(triplets)
triplets2 <- rbind(triplets2, triplets[1,])

for (i in 2: nrow(triplets)){
    if ((as.character(triplets[i-1, "Mir"])!=as.character(triplets[i, "Mir"])) |
    (as.character(triplets[i-1, "Gene"])!=as.character(triplets[i, "Gene"])))
    {
        triplets2 <- rbind(triplets2, triplets[i,])
    }
}

triplets_cit <- triplets2[which(triplets2$Causal_Call==1),]
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
write.table(mir_t_ord_desc, "Freq_miRNA_Control.txt", sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)


sum(cor.network$FDR<0.1)

> "MI0000743_MIMAT0000686" %in% triplets_cit$Mir
[1] TRUE
> "MI0000743_MIMAT0004677" %in% triplets_cit$Mir
[1] TRUE


# check correlation of the triplets
mir.exprs = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/MicroRNA/LGRC_MicroRNA_Overlap_SNP_Control_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))
gene.exprs = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/Gene/LGRC_Gene_Overlap_SNP_Control_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))

snp.overlap.temp = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp)))
snp.overlap = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc)
colnames(snp.overlap["rs4434979_C",])



pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/test1_mir_gene.pdf")
plot(gene.exprs["A_23_P121956",], mir.exprs["MI0017450_MIMAT0019982",], main=c("Cor coef ", round(cor(gene.exprs["A_23_P121956",], mir.exprs["MI0017450_MIMAT0019982",]), 3)))
abline(lm(as.numeric(mir.exprs["MI0017450_MIMAT0019982",])~ as.numeric(gene.exprs["A_23_P121956",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/test1_snp_mir.pdf")
#plot(snp.overlap["rs4434979_C",], mir.exprs["MI0017450_MIMAT0019982",], main=c("Cor coef ", round(cor(snp.overlap["rs4434979_C",], mir.exprs["MI0017450_MIMAT0019982",]), 3)))
plot(mir.exprs["MI0017450_MIMAT0019982",], snp.overlap["rs4434979_C",])
#abline(lm(as.numeric(mir.exprs["MI0017450_MIMAT0019982",])~ as.numeric(snp.overlap["rs4434979_C",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/test1_snp_gene.pdf")
#plot(snp.overlap["rs4434979_C",], mir.exprs["MI0017450_MIMAT0019982",], main=c("Cor coef ", round(cor(snp.overlap["rs4434979_C",], mir.exprs["MI0017450_MIMAT0019982",]), 3)))
plot(gene.exprs["A_23_P121956",], snp.overlap["rs4434979_C",])
#abline(lm(as.numeric(mir.exprs["MI0017450_MIMAT0019982",])~ as.numeric(snp.overlap["rs4434979_C",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/test2_mir_gene.pdf")
plot(gene.exprs["A_23_P428129",], mir.exprs["MI0000477_MIMAT0000449",], main=c("Cor coef ", round(cor(gene.exprs["A_23_P428129",], mir.exprs["MI0000477_MIMAT0000449",]), 3)))
abline(lm(as.numeric(mir.exprs["MI0000477_MIMAT0000449",])~ as.numeric(gene.exprs["A_23_P428129",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/test3_snp_mir.pdf")
boxplot(mir.exprs["MI0000477_MIMAT0000449", which(snp.overlap["rs2302836_G",]==0)], mir.exprs["MI0000477_MIMAT0000449", which(snp.overlap["rs2302836_G",]!=0)], col=c("green", "red"), names=c("no SNP", "SNP"), ylab=c("miRNA expression"))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/test3_snp_gene.pdf")
#plot(gene.exprs["A_23_P428129",], snp.overlap["rs2302836_G",])
boxplot(gene.exprs["A_23_P428129", which(snp.overlap["rs2302836_G",]==0)], gene.exprs["A_23_P428129", which(snp.overlap["rs2302836_G",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "SNP"), ylab=c("gene expression"))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/test4_snp_mir.pdf")
boxplot(mir.exprs["MI0017450_MIMAT0019982", which(snp.overlap["rs4434979_C",]==0)], mir.exprs["MI0017450_MIMAT0019982", which(snp.overlap["rs4434979_C",]!=0)], col=c("green", "red"), names=c("no SNP", "SNP"), ylab=c("miRNA expression"))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/test4_snp_gene.pdf")
#plot(gene.exprs["A_23_P121956",], snp.overlap["rs2302836_G",])
boxplot(gene.exprs["A_23_P121956", which(snp.overlap["rs4434979_C",]==0)], gene.exprs["A_23_P121956", which(snp.overlap["rs4434979_C",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "SNP"), ylab=c("gene expression"))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/test3_snp_mir.pdf")
boxplot(mir.exprs["MI0000263_MIMAT0004553", which(snp.overlap["kgp2233006_T",]==0)], mir.exprs["MI0000263_MIMAT0004553", which(snp.overlap["kgp2233006_T",]!=0)], col=c("green", "red"), names=c("no SNP", "SNP"), ylab=c("miRNA expression"))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/test3_snp_gene.pdf")
boxplot(gene.exprs["A_23_P343398", which(snp.overlap["kgp2233006_T",]==0)], gene.exprs["A_23_P343398", which(snp.overlap["kgp2233006_T",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "SNP"), ylab=c("gene expression"))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/test3_mir_gene.pdf")
plot(gene.exprs["A_23_P343398",], mir.exprs["MI0000263_MIMAT0004553",], main=c("Cor coef ", round(cor(gene.exprs["A_23_P343398",], mir.exprs["MI0000263_MIMAT0004553",]), 3)), col="blue", pch=19)
abline(lm(as.numeric(mir.exprs["MI0000263_MIMAT0004553",])~ as.numeric(gene.exprs["A_23_P343398",])))
dev.off()


pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/test4_mir_gene_ctr.pdf")
plot(gene.exprs["A_23_P106299",], mir.exprs["MI0000458_MIMAT0000434",], main=c("Cor coef ", round(cor(gene.exprs["A_23_P106299",], mir.exprs["MI0000458_MIMAT0000434",]), 3)), col="blue", pch=19)
abline(lm(as.numeric(mir.exprs["MI0000458_MIMAT0000434",])~ as.numeric(gene.exprs["A_23_P106299",])))
dev.off()


pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/test4_snp_mir_ctr.pdf")
boxplot(mir.exprs["MI0000458_MIMAT0000434", which(snp.overlap["rs13288748_A",]==0)], mir.exprs["MI0000458_MIMAT0000434", which(snp.overlap["rs13288748_A",]!=0)], col=c("green", "red"), names=c("no SNP", "SNP"), ylab=c("miRNA expression"))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/test4_snp_gene_ctr.pdf")
boxplot(gene.exprs["A_23_P106299", which(snp.overlap["rs13288748_A",]==0)], gene.exprs["A_23_P106299", which(snp.overlap["rs13288748_A",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "SNP"), ylab=c("gene expression"))
dev.off()


mir.exprs = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/MicroRNA/LGRC_MicroRNA_Overlap_SNP_COPD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))
gene.exprs = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/Gene/LGRC_Gene_Overlap_SNP_COPD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))

snp.overlap.temp = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_COPD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp)))
snp.overlap = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_COPD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc)
colnames(snp.overlap["rs4434979_C",])



pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/test_mir_gene_map4k4.pdf")
plot(gene.exprs["A_23_P348524",], mir.exprs["MI0001648_MIMAT0001541",], main=c("Cor coef ", round(cor(gene.exprs["A_23_P348524",], mir.exprs["MI0001648_MIMAT0001541",]), 3)), col="blue", pch=19)
abline(lm(as.numeric(mir.exprs["MI0001648_MIMAT0001541",])~ as.numeric(gene.exprs["A_23_P348524",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/test_snp_mir_map4k4.pdf")
boxplot(mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs10875366_A",]==0)], mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs10875366_A",]!=0)], col=c("green", "red"), names=c("no SNP", "SNP"), ylab=c("miRNA expression"))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/test_snp_gene_map4k4.pdf")
boxplot(gene.exprs["A_23_P348524", which(snp.overlap["rs10875366_A",]==0)], gene.exprs["A_23_P348524", which(snp.overlap["rs10875366_A",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "SNP"), ylab=c("gene expression"))
dev.off()


pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/test_mir_gene_map4k4_449c.pdf")
plot(gene.exprs["A_23_P348524",], mir.exprs["MI0003823_MIMAT0010251",], main=c("Cor coef ", round(cor(gene.exprs["A_23_P348524",], mir.exprs["MI0003823_MIMAT0010251",]), 3)), col="blue", pch=19)
abline(lm(as.numeric(mir.exprs["MI0003823_MIMAT0010251",])~ as.numeric(gene.exprs["A_23_P348524",])))
dev.off()


pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/test_snp_mir_map4k4_449c.pdf")
boxplot(mir.exprs["MI0003823_MIMAT0010251", which(snp.overlap["rs10875366_A",]==0)], mir.exprs["MI0003823_MIMAT0010251", which(snp.overlap["rs10875366_A",]!=0)], col=c("green", "red"), names=c("no SNP", "SNP"), ylab=c("miRNA expression"))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/test_snp_gene_map4k4_449c.pdf")
boxplot(gene.exprs["A_23_P348524", which(snp.overlap["rs10875366_A",]==0)], gene.exprs["A_23_P348524", which(snp.overlap["rs10875366_A",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "SNP"), ylab=c("gene expression"))
dev.off()


x <- c(501, 105, 693-501, 693-105)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))


x <- c(525, 105, 693-525, 693-105)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))


x <- c(13651, 377, 254708-13651, 17024-377)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))

x <- c(10319, 377, 701188-10319, 17024-377)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))

x <- c(10, 29, 5, 100)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))

x <- c(42, 463, 5, 100)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t), alternative="greater")

https://support.bioconductor.org/p/71702/


/protected/projects/pulmarray/HoggCOPD/A1ATGrifols/experiments/Expmt_0033_PhaseII_Immune_GSVA
/protected/projects/lgrc/small_RNA/lgrc_365/Integrated_Networks2/ILD/run_cit.R
/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/

/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/

#Run GSVA
gene_list = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/COPD_Gene_miR34b_3p.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
ann_gene = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/121212_LGRC_Gene_Array_Annotation2.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(ann_gene) <- ann_gene$NAME
ann_gene[gene_list$Gene, "GENE"]

exp_matrix <- as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Marcet/Marcet_GSE1247_entrezgcdf_v16_noat.txt", sep="\t", stringsAsFactors=FALSE))
final_list <- intersect(rownames(exp_matrix), ann_gene[gene_list$Gene, "GENE"])
final_list_gsva <- list(set1=final_list)
mygsva <- gsva(exp_matrix, final_list_gsva, min.sz=10, max.sz=500, verbose=TRUE, rnaseq=FALSE)
gsva.val <- mygsva$es.obs
write.table(gsva.val, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Marcet/gsva.val.txt", sep="\t")

ann_samples = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Marcet/Marcet_sample_description.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(ann_samples) <- ann_samples$Sample_ID

#cels <- substr(colnames(gsva.val),4,9)
cels <- colnames(gsva.val)
cels[order(gsva.val)]
names_arg <- paste(ann_samples[cels[order(gsva.val)],"Donor"], "_", ann_samples[cels[order(gsva.val)],"Treatment"])
pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/Marcet/TargetsCOPD_miR34b_3p_in_Marcet.pdf", height=5, width=15)
barplot(gsva.val[order(gsva.val)], names.arg=names_arg, cex.names=0.7)
dev.off()


#Hogg
hogg <- readRDS("/protected/projects/pulmarray/HoggCOPD/A1ATGrifols/experiments/Expmt_0033_PhaseII_Immune_GSVA/phaseIIparenchymaNormalizedData.RDS")
gene_exp <- exprs(hogg)
pheno_data <- pData(hogg)
write.table(gene_exp, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/rehoggdata/parenchima.exp.txt", sep="\t", quote=FALSE)
write.table(pheno_data, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/rehoggdata/parenchima.pheno.txt", sep="\t", quote=FALSE)


hogg <- readRDS("/protected/projects/pulmarray/HoggCOPD/A1ATGrifols/experiments/Expmt_0033_PhaseII_Immune_GSVA/phaseIIairwayNormalizedData.RDS")
gene_exp <- exprs(hogg)
pheno_data <- pData(hogg)
write.table(gene_exp, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/rehoggdata/airway.exp.txt", sep="\t", quote=FALSE)
write.table(pheno_data, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/OtherDatasets/rehoggdata/airway.pheno.txt", sep="\t", quote=FALSE)


mirs_rds <- readRDS("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/MicroRNA/LGRC_MicroRNA_Overlap_SNP_All_Samples.rds")
library(affy)
write.table(fData(mirs_rds), "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/LGRC_MicroRNA_Locations_corr.txt", sep="\t", quote=FALSE)
write.table(fData(mirs_rds)[,c("Chrom", "Start", "Stop")], "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/LGRC_MicroRNA_Locations_corr.txt", sep="\t", quote=FALSE)


mir_f_rds <- readRDS("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/MicroRNA/LGRC_MicroRNA_Counts_351.rds")
library(affy)
rownames(pData(mir_f_rds))
write.table(pData(mir_f_rds), "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/Filtered_miRNA_demo.txt", sep="\t", quote=FALSE)
write.table(good.samples, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/Filtered_miRNA_samples.txt", sep="\t", quote=FALSE)

copd <- pheno[which(pheno$Disease==0 | pheno$Disease==2),]
ild <- pheno[which(pheno$Disease==0 | pheno$Disease==1),]
ta = table(copd$Smoke, copd$Disease)
htest$p.value
nrow(pheno[which(pheno$COPD==1 & pheno$Smoke==1),])
nrow(pheno[which(pheno$COPD==1 & pheno$Smoke==2),])
nrow(pheno[which(pheno$COPD==1 & pheno$Smoke==3),])
htest = t.test(pheno$Pack_Years[pheno$ILD == 1], pheno$Pack_Years[pheno$Control == 1])
htest = t.test(pheno$Age[pheno$ILD == 1], pheno$Age[pheno$Control == 1])


#kgp6539253_C (PRKCA) MI0003673_MIMAT0003327 (miR-449b-5p)  A_23_P14035 (AKAP3)
mir.exprs = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/MicroRNA/LGRC_MicroRNA_Overlap_SNP_COPD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))
gene.exprs = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/Gene/LGRC_Gene_Overlap_SNP_COPD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))

snp.overlap.temp = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_COPD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp)))
snp.overlap = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_COPD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc)
colnames(snp.overlap["rs525770_C",])

snp.overlap.temp.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc.ctr = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp.ctr)))
snp.overlap.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc.ctr)
colnames(snp.overlap.ctr["rs525770_C",])




##new good snp
#rs525770_C 25/107 samples in COPD and 0/33 in Control
#A_23_P92562, A_23_P45811, A_23_P77714, A_23_P214950, A_23_P30283
#MI0001648_MIMAT0001541

#A_23_P32629 = HS6ST3


pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_mir_gene_copd6.pdf")
plot(gene.exprs["A_23_P32629",], mir.exprs["MI0001648_MIMAT0001541",], ylab="miR-449a expression", xlab="HS6ST3 expression", main=c(paste("Cor coef ", round(cor(gene.exprs["A_23_P32629",], mir.exprs["MI0001648_MIMAT0001541",], method="pearson"), 3)), paste("p-value ", round(cor.test(gene.exprs["A_23_P32629",], mir.exprs["MI0001648_MIMAT0001541",], method="pearson")$p.value, 5))), col="blue", pch=19, cex.lab=1.3)
abline(lm(as.numeric(mir.exprs["MI0001648_MIMAT0001541",])~ as.numeric(gene.exprs["A_23_P32629",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_snp_mir_copd6.pdf")
boxplot(mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]==0)], mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]!=0)], col=c("green", "red"), names=c("no SNP", "rs525770_C"), ylab=c("miR-449a expression"), xlab="HS6ST3", main=paste("p-value ", round(t.test(mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]==0)], mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_snp_gene_copd6.pdf")
boxplot(gene.exprs["A_23_P32629", which(snp.overlap["rs525770_C",]==0)], gene.exprs["A_23_P32629", which(snp.overlap["rs525770_C",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "rs525770_C"), ylab=c("HS6ST3 expression"), xlab="HS6ST3", main=paste("p-value ", round(t.test(gene.exprs["A_23_P32629", which(snp.overlap["rs525770_C",]==0)], gene.exprs["A_23_P32629", which(snp.overlap["rs525770_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_mir_gene_copd.pdf")
plot(gene.exprs["A_23_P92562",], mir.exprs["MI0001648_MIMAT0001541",], ylab="miR-449a expression", xlab="ADH7 expression", main=c(paste("Cor coef ", round(cor(gene.exprs["A_23_P92562",], mir.exprs["MI0001648_MIMAT0001541",], method="pearson"), 3)), paste("p-value ", round(cor.test(gene.exprs["A_23_P92562",], mir.exprs["MI0001648_MIMAT0001541",], method="pearson")$p.value, 5))), col="blue", pch=19, cex.lab=1.3)
abline(lm(as.numeric(mir.exprs["MI0001648_MIMAT0001541",])~ as.numeric(gene.exprs["A_23_P92562",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_snp_mir_copd.pdf")
boxplot(mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]==0)], mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]!=0)], col=c("green", "red"), names=c("no SNP", "rs525770_C"), ylab=c("miR-449a expression"), xlab="HS6ST3", main=paste("p-value ", round(t.test(mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]==0)], mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_snp_gene_copd.pdf")
boxplot(gene.exprs["A_23_P92562", which(snp.overlap["rs525770_C",]==0)], gene.exprs["A_23_P92562", which(snp.overlap["rs525770_C",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "rs525770_C"), ylab=c("ADH7 expression"), xlab="HS6ST3", main=paste("p-value ", round(t.test(gene.exprs["A_23_P92562", which(snp.overlap["rs525770_C",]==0)], gene.exprs["A_23_P92562", which(snp.overlap["rs525770_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_mir_gene_copd2.pdf")
plot(gene.exprs["A_23_P45811",], mir.exprs["MI0001648_MIMAT0001541",], ylab="miR-449a expression", xlab="DIO1 expression", main=c(paste("Cor coef ", round(cor(gene.exprs["A_23_P45811",], mir.exprs["MI0001648_MIMAT0001541",], method="pearson"), 3)), paste("p-value ", round(cor.test(gene.exprs["A_23_P45811",], mir.exprs["MI0001648_MIMAT0001541",], method="pearson")$p.value, 5))), col="blue", pch=19, cex.lab=1.3)
abline(lm(as.numeric(mir.exprs["MI0001648_MIMAT0001541",])~ as.numeric(gene.exprs["A_23_P45811",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_snp_mir_copd2.pdf")
boxplot(mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]==0)], mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]!=0)], col=c("green", "red"), names=c("no SNP", "rs525770_C"), ylab=c("miR-449a expression"), xlab="HS6ST3", main=paste("p-value ", round(t.test(mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]==0)], mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_snp_gene_copd2.pdf")
boxplot(gene.exprs["A_23_P45811", which(snp.overlap["rs525770_C",]==0)], gene.exprs["A_23_P45811", which(snp.overlap["rs525770_C",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "rs525770_C"), ylab=c("DIO1 expression"), xlab="HS6ST3", main=paste("p-value ", round(t.test(gene.exprs["A_23_P45811", which(snp.overlap["rs525770_C",]==0)], gene.exprs["A_23_P45811", which(snp.overlap["rs525770_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()


pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_mir_gene_copd3.pdf")
plot(gene.exprs["A_23_P77714",], mir.exprs["MI0001648_MIMAT0001541",], ylab="miR-449a expression", xlab="CLUAP1 expression", main=c(paste("Cor coef ", round(cor(gene.exprs["A_23_P77714",], mir.exprs["MI0001648_MIMAT0001541",], method="pearson"), 3)), paste("p-value ", round(cor.test(gene.exprs["A_23_P77714",], mir.exprs["MI0001648_MIMAT0001541",], method="pearson")$p.value, 5))), col="blue", pch=19, cex.lab=1.3)
abline(lm(as.numeric(mir.exprs["MI0001648_MIMAT0001541",])~ as.numeric(gene.exprs["A_23_P77714",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_snp_mir_copd3.pdf")
boxplot(mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]==0)], mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]!=0)], col=c("green", "red"), names=c("no SNP", "rs525770_C"), ylab=c("miR-449a expression"), xlab="HS6ST3", main=paste("p-value ", round(t.test(mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]==0)], mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_snp_gene_copd3.pdf")
boxplot(gene.exprs["A_23_P77714", which(snp.overlap["rs525770_C",]==0)], gene.exprs["A_23_P77714", which(snp.overlap["rs525770_C",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "rs525770_C"), ylab=c("CLUAP1 expression"), xlab="HS6ST3", main=paste("p-value ", round(t.test(gene.exprs["A_23_P77714", which(snp.overlap["rs525770_C",]==0)], gene.exprs["A_23_P77714", which(snp.overlap["rs525770_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()


pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_mir_gene_copd4.pdf")
plot(gene.exprs["A_23_P214950",], mir.exprs["MI0001648_MIMAT0001541",], ylab="miR-449a expression", xlab="PERP expression", main=c(paste("Cor coef ", round(cor(gene.exprs["A_23_P214950",], mir.exprs["MI0001648_MIMAT0001541",], method="pearson"), 3)), paste("p-value ", round(cor.test(gene.exprs["A_23_P214950",], mir.exprs["MI0001648_MIMAT0001541",], method="pearson")$p.value, 5))), col="blue", pch=19, cex.lab=1.3)
abline(lm(as.numeric(mir.exprs["MI0001648_MIMAT0001541",])~ as.numeric(gene.exprs["A_23_P214950",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_snp_mir_copd4.pdf")
boxplot(mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]==0)], mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]!=0)], col=c("green", "red"), names=c("no SNP", "rs525770_C"), ylab=c("miR-449a expression"), xlab="HS6ST3", main=paste("p-value ", round(t.test(mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]==0)], mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_snp_gene_copd4.pdf")
boxplot(gene.exprs["A_23_P214950", which(snp.overlap["rs525770_C",]==0)], gene.exprs["A_23_P214950", which(snp.overlap["rs525770_C",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "rs525770_C"), ylab=c("PERP expression"), xlab="HS6ST3", main=paste("p-value ", round(t.test(gene.exprs["A_23_P214950", which(snp.overlap["rs525770_C",]==0)], gene.exprs["A_23_P214950", which(snp.overlap["rs525770_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()


pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_mir_gene_copd5.pdf")
plot(gene.exprs["A_23_P30283",], mir.exprs["MI0001648_MIMAT0001541",], ylab="miR-449a expression", xlab="FAM174A expression", main=c(paste("Cor coef ", round(cor(gene.exprs["A_23_P30283",], mir.exprs["MI0001648_MIMAT0001541",], method="pearson"), 3)), paste("p-value ", round(cor.test(gene.exprs["A_23_P30283",], mir.exprs["MI0001648_MIMAT0001541",], method="pearson")$p.value, 5))), col="blue", pch=19, cex.lab=1.3)
abline(lm(as.numeric(mir.exprs["MI0001648_MIMAT0001541",])~ as.numeric(gene.exprs["A_23_P30283",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_snp_mir_copd5.pdf")
boxplot(mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]==0)], mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]!=0)], col=c("green", "red"), names=c("no SNP", "rs525770_C"), ylab=c("miR-449a expression"), xlab="HS6ST3", main=paste("p-value ", round(t.test(mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]==0)], mir.exprs["MI0001648_MIMAT0001541", which(snp.overlap["rs525770_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/HS6ST3_snp_gene_copd5.pdf")
boxplot(gene.exprs["A_23_P30283", which(snp.overlap["rs525770_C",]==0)], gene.exprs["A_23_P30283", which(snp.overlap["rs525770_C",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "rs525770_C"), ylab=c("FAM174A expression"), xlab="HS6ST3", main=paste("p-value ", round(t.test(gene.exprs["A_23_P30283", which(snp.overlap["rs525770_C",]==0)], gene.exprs["A_23_P30283", which(snp.overlap["rs525770_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

###############
pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/PRKCA_mir_gene_copd.pdf")
plot(gene.exprs["A_23_P14035",], mir.exprs["MI0003673_MIMAT0003327",], ylab="miR-449b-5p expression", xlab="AKAP3 expression", main=c(paste("Cor coef ", round(cor(gene.exprs["A_23_P14035",], mir.exprs["MI0003673_MIMAT0003327",], method="pearson"), 3)), paste("p-value ", round(cor.test(gene.exprs["A_23_P14035",], mir.exprs["MI0003673_MIMAT0003327",], method="pearson")$p.value, 5))), col="blue", pch=19, cex.lab=1.3)
abline(lm(as.numeric(mir.exprs["MI0003673_MIMAT0003327",])~ as.numeric(gene.exprs["A_23_P14035",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/PRKCA_snp_mir_copd.pdf")
boxplot(mir.exprs["MI0003673_MIMAT0003327", which(snp.overlap["kgp6539253_C",]==0)], mir.exprs["MI0003673_MIMAT0003327", which(snp.overlap["kgp6539253_C",]!=0)], col=c("green", "red"), names=c("no SNP", "kgp6539253_C"), ylab=c("miR-449b-5p expression"), xlab="PRKCA", main=paste("p-value ", round(t.test(mir.exprs["MI0003673_MIMAT0003327", which(snp.overlap["kgp6539253_C",]==0)], mir.exprs["MI0003673_MIMAT0003327", which(snp.overlap["kgp6539253_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/PRKCA_snp_gene_copd.pdf")
boxplot(gene.exprs["A_23_P14035", which(snp.overlap["kgp6539253_C",]==0)], gene.exprs["A_23_P14035", which(snp.overlap["kgp6539253_C",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "kgp6539253_C"), ylab=c("AKAP3 expression"), xlab="PRKCA", main=paste("p-value ", round(t.test(gene.exprs["A_23_P14035", which(snp.overlap["kgp6539253_C",]==0)], gene.exprs["A_23_P14035", which(snp.overlap["kgp6539253_C",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()
#################

#IL1RAP (rs9830737_G, chr3:190537545) MI0003673_MIMAT0003327 (miR-449b-5p) A_23_P312300 (SCGB2A1)

mir.exprs = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/MicroRNA/LGRC_MicroRNA_Overlap_SNP_ILD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))
gene.exprs = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/Gene/LGRC_Gene_Overlap_SNP_ILD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))

snp.overlap.temp = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_COPD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp)))
snp.overlap = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_COPD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc)
#colnames(snp.overlap["rs4434979_C",])

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/IL1RAP_mir_gene_ild.pdf")
plot(gene.exprs["A_23_P312300",], mir.exprs["MI0003673_MIMAT0003327",], ylab="miR-449b-5p expression", xlab="SCGB2A1 expression", main=c(paste("Cor coef ", round(cor(gene.exprs["A_23_P312300",], mir.exprs["MI0003673_MIMAT0003327",], method="pearson"), 3)), paste("p-value ", round(cor.test(gene.exprs["A_23_P312300",], mir.exprs["MI0003673_MIMAT0003327",], method="pearson")$p.value, 5))), col="blue", pch=19, cex.lab=1.3)
abline(lm(as.numeric(mir.exprs["MI0003673_MIMAT0003327",])~ as.numeric(gene.exprs["A_23_P312300",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/IL1RAP_snp_mir_ild.pdf")
boxplot(mir.exprs["MI0003673_MIMAT0003327", which(snp.overlap["rs9830737_G",]==0)], mir.exprs["MI0003673_MIMAT0003327", which(snp.overlap["rs9830737_G",]!=0)], col=c("green", "red"), names=c("no SNP", "rs9830737_G"), ylab=c("miR-449b-5p expression"), xlab="IL1RAP", main=paste("p-value ", round(t.test(mir.exprs["MI0003673_MIMAT0003327", which(snp.overlap["rs9830737_G",]==0)], mir.exprs["MI0003673_MIMAT0003327", which(snp.overlap["rs9830737_G",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/IL1RAP_snp_gene_ild.pdf")
boxplot(gene.exprs["A_23_P312300", which(snp.overlap["rs9830737_G",]==0)], gene.exprs["A_23_P312300", which(snp.overlap["rs9830737_G",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "rs9830737_G"), ylab=c("SCGB2A1 expression"), xlab="IL1RAP", main=paste("p-value ", round(t.test(gene.exprs["A_23_P312300", which(snp.overlap["rs9830737_G",]==0)], gene.exprs["A_23_P312300", which(snp.overlap["rs9830737_G",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

#CLCA1 (rs618555_T, chr1:86481084)  MI0003823_MIMAT0010251 (miR-449b-5p) A_23_P338325 (ELK3), A_24_P48177 (ST3GAL2)
pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/CLCA1_mir_gene_ild.pdf")
plot(gene.exprs["A_23_P338325",], mir.exprs["MI0003823_MIMAT0010251",], ylab="miR-449c-5p expression", xlab="ELK3 expression", main=c(paste("Cor coef ", round(cor(gene.exprs["A_23_P338325",], mir.exprs["MI0003823_MIMAT0010251",], method="pearson"), 3)), paste("p-value ", round(cor.test(gene.exprs["A_23_P338325",], mir.exprs["MI0003823_MIMAT0010251",], method="pearson")$p.value, 5))), col="blue", pch=19, cex.lab=1.3)
abline(lm(as.numeric(mir.exprs["MI0003823_MIMAT0010251",])~ as.numeric(gene.exprs["A_23_P338325",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/CLCA1_snp_mir_ild.pdf")
boxplot(mir.exprs["MI0003823_MIMAT0010251", which(snp.overlap["rs618555_T",]==0)], mir.exprs["MI0003823_MIMAT0010251", which(snp.overlap["rs618555_T",]!=0)], col=c("green", "red"), names=c("no SNP", "rs9830737_G"), ylab=c("miR-449c-5p expression"), xlab="CLCA1", main=paste("p-value ", round(t.test(mir.exprs["MI0003823_MIMAT0010251", which(snp.overlap["rs618555_T",]==0)], mir.exprs["MI0003823_MIMAT0010251", which(snp.overlap["rs618555_T",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/CLCA1_snp_gene_ild.pdf")
boxplot(gene.exprs["A_23_P338325", which(snp.overlap["rs618555_T",]==0)], gene.exprs["A_23_P338325", which(snp.overlap["rs618555_T",]!=0)], col=c("darkgreen", "darkred"), names=c("no SNP", "rs9830737_G"), ylab=c("ELK3 expression"), xlab="CLCA1", main=paste("p-value ", round(t.test(gene.exprs["A_23_P338325", which(snp.overlap["rs618555_T",]==0)], gene.exprs["A_23_P312300", which(snp.overlap["rs618555_T",]!=0)])$p.value,5)), cex.lab=1.3)
dev.off()

trans <- read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/eQTL_Trans_corr.txt.gz', header=TRUE, sep="\t")
y <- trans[which(trans$p<1e-4),]
"A_23_P69908" %in% y$gene
y[which(y[,"gene"] == "A_23_P69908"), ]
#write.table(y[which(y[,"gene"] == "A_23_P69908"), c("SNP", "p.value", "FDR")], "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/SNPS_assoc_GLRX.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

y[which(y[,"gene"] == "A_23_P117082"), ]
#write.table(y[which(y[,"gene"] == "A_23_P117082"), c("SNP", "p.value", "FDR")], "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/SNPS_assoc_HEBP1.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

trans_ctr <- read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/eQTL_Trans_corr.txt.gz', header=TRUE, sep="\t")
y_ctr <- trans_ctr[which(trans_ctr$p<1e-4),]
"A_23_P69908" %in% y$gene
y_ctr[which(y_ctr[,"gene"] == "A_23_P69908"), ]
#write.table(y_ctr[which(y_ctr[,"gene"] == "A_23_P69908"), c("SNP", "p.value", "FDR")], "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/SNPS_assoc_GLRX_Control.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

y_ctr[which(y_ctr[,"gene"] == "A_23_P117082"), ]
#write.table(y[which(y_ctr[,"gene"] == "A_23_P117082"), c("SNP", "p.value", "FDR")], "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/SNPS_assoc_HEBP1_Control.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/proj1_microarray_sol/
library(affy)
x <- readRDS("all.data.rma.rds")
gene_exp <- exprs(x)
write.table(gene_exp, "cc.gene.exp.tab.csv", sep="\t", quote=FALSE)


/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/Biomarker/
/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/Biomarker_Pipeline/rabbit/



copd1 <- read.delim("C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\COPD_SNP_mir34b_3p_all_corr_BIO.txt", sep="\t", check.names=FALSE, stringsAsFactors=FALSE) 
copd2 <- read.delim("C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\COPD_SNP_mir34c_3p_all_corr_BIO.txt", sep="\t", check.names=FALSE, stringsAsFactors=FALSE) 
copd3 <- read.delim("C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\COPD_SNP_mir34c_5p_all_corr_BIO.txt", sep="\t", check.names=FALSE, stringsAsFactors=FALSE) 
copd4 <- read.delim("C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\COPD_SNP_mir449a_all_corr_BIO.txt", sep="\t", check.names=FALSE, stringsAsFactors=FALSE) 
copd5 <- read.delim("C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\COPD_SNP_mir449b_5p_all_corr_BIO.txt", sep="\t", check.names=FALSE, stringsAsFactors=FALSE) 
copd6 <- read.delim("C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\COPD_SNP_mir449c_5p_all_corr_BIO.txt", sep="\t", check.names=FALSE, stringsAsFactors=FALSE) 

for (i in 1:nrow(copd2)){
    if (!(copd2$SNP[i] %in% copd1$SNP)){
        copd1 <- rbind(copd1, copd2[i,])
        print(i)
    }
}

for (i in 1:nrow(copd3)){
    if (!(copd3$SNP[i] %in% copd1$SNP)){
        copd1 <- rbind(copd1, copd3[i,])
        print(i)
    }
}

for (i in 1:nrow(copd4)){
    if (!(copd4$SNP[i] %in% copd1$SNP)){
        copd1 <- rbind(copd1, copd4[i,])
        print(i)
    }
}

for (i in 1:nrow(copd5)){
    if (!(copd5$SNP[i] %in% copd1$SNP)){
        copd1 <- rbind(copd1, copd5[i,])
        print(i)
    }
}

for (i in 1:nrow(copd6)){
    if (!(copd6$SNP[i] %in% copd1$SNP)){
        copd1 <- rbind(copd1, copd6[i,])
        print(i)
    }
}

write.table(copd1[,c("SNP", "SNP_Chr", "SNP_Start", "SNP_BIO")], "C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\COPD_SNPs.txt", sep="\t", quote=FALSE)



####
ild1 <- read.delim("C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\ILD_SNP_mir34b_5p_all_corr_BIO.txt", sep="\t", check.names=FALSE, stringsAsFactors=FALSE) 
ild2 <- read.delim("C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\ILD_SNP_mir34c_3p_all_corr_BIO.txt", sep="\t", check.names=FALSE, stringsAsFactors=FALSE) 
ild3 <- read.delim("C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\ILD_SNP_mir34c_5p_all_corr_BIO.txt", sep="\t", check.names=FALSE, stringsAsFactors=FALSE) 
ild4 <- read.delim("C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\ILD_SNP_mir449a_all_corr_BIO.txt", sep="\t", check.names=FALSE, stringsAsFactors=FALSE) 
ild5 <- read.delim("C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\ILD_SNP_mir449b_5p_all_corr_BIO.txt", sep="\t", check.names=FALSE, stringsAsFactors=FALSE) 
ild6 <- read.delim("C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\ILD_SNP_mir449c_5p_all_corr_BIO.txt", sep="\t", check.names=FALSE, stringsAsFactors=FALSE) 

for (i in 1:nrow(ild2)){
    if (!(ild2$SNP[i] %in% ild1$SNP)){
        ild1 <- rbind(ild1, ild2[i,])
        print(i)
    }
}

for (i in 1:nrow(ild3)){
    if (!(ild3$SNP[i] %in% ild1$SNP)){
        ild1 <- rbind(ild1, ild3[i,])
        print(i)
    }
}

for (i in 1:nrow(ild4)){
    if (!(ild4$SNP[i] %in% ild1$SNP)){
        ild1 <- rbind(ild1, ild4[i,])
        print(i)
    }
}

for (i in 1:nrow(ild5)){
    if (!(ild5$SNP[i] %in% ild1$SNP)){
        ild1 <- rbind(ild1, ild5[i,])
        print(i)
    }
}

for (i in 1:nrow(ild6)){
    if (!(ild6$SNP[i] %in% ild1$SNP)){
        ild1 <- rbind(ild1, ild6[i,])
        print(i)
    }
}

write.table(ild1[,c("SNP", "SNP_Chr", "SNP_Start", "SNP_BIO")], "C:\\Users\\Ana\\Desktop\\LGRC\\SNPs_corr\\ILD_SNPs.txt", sep="\t", quote=FALSE)



###For Jen SNPs in 

#A_23_P154784 BPIFB1
#BPIFA1
#MUC5B A_24_P102650	727897	MUC5B


mir.exprs = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/MicroRNA/LGRC_MicroRNA_Overlap_SNP_COPD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))

gene.exprs.ild = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/Gene/LGRC_Gene_Overlap_SNP_ILD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))

gene.exprs.ctr = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/Gene/LGRC_Gene_Overlap_SNP_Control_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))

gene.exprs.copd = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/Gene/LGRC_Gene_Overlap_SNP_COPD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))



write.table(gene.exprs.ild[c('A_23_P154784', 'A_24_P102650'),], "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/ILD_2genes.txt", sep="\t", quote=FALSE)

write.table(gene.exprs.ctr[c('A_23_P154784', 'A_24_P102650'),], "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/Control_2genes.txt", sep="\t", quote=FALSE)

t.test(gene.exprs.ild['A_23_P154784',], gene.exprs.ctr['A_23_P154784',]) #p-value = 1.217e-07
t.test(gene.exprs.ild['A_24_P102650',], gene.exprs.ctr['A_24_P102650',]) #p-value = 6.057e-06

snp.overlap.temp = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_ILD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp)))
snp.overlap = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_ILD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc)
colnames(snp.overlap["rs525770_C",])

snp.overlap.temp.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc.ctr = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp.ctr)))
snp.overlap.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc.ctr)
colnames(snp.overlap.ctr["rs525770_C",])


length(which(snp.overlap["rs2378252_C",]==1 | snp.overlap["rs2378252_C",]==2)) #79/113

length(which(snp.overlap.ctr["rs2378252_C",]==1 | snp.overlap.ctr["rs2378252_C",]==2)) #32/38


length(which(snp.overlap["kgp3480870_C",]==1 | snp.overlap["kgp3480870_C",]==2)) #78/113

length(which(snp.overlap.ctr["kgp3480870_C",]==1 | snp.overlap.ctr["kgp3480870_C",]==2)) #31/38


length(which(snp.overlap["rs910869_T",]==1 | snp.overlap["rs910869_T",]==2)) #69/113

length(which(snp.overlap.ctr["rs910869_T",]==1 | snp.overlap.ctr["rs910869_T",]==2)) #30/38


length(which(snp.overlap["kgp9581225_G",]==1 | snp.overlap["kgp9581225_G",]==2)) #82/113

length(which(snp.overlap.ctr["kgp9581225_G",]==1 | snp.overlap.ctr["kgp9581225_G",]==2)) #32/38


length(which(snp.overlap["kgp2473673_C",]==1 | snp.overlap["kgp2473673_C",]==2)) #68/113

length(which(snp.overlap.ctr["kgp2473673_C",]==1 | snp.overlap.ctr["kgp2473673_C",]==2)) #30/38


length(which(snp.overlap["kgp5909883_C",]==1 | snp.overlap["kgp5909883_C",]==2)) #68/113

length(which(snp.overlap.ctr["kgp5909883_C",]==1 | snp.overlap.ctr["kgp5909883_C",]==2)) #30/38

rs2378252_C	chr20	33283403
kgp3480870_C	chr20	33285053
rs910869_T	chr20	33292777
kgp9581225_G	chr20	33294353
kgp2473673_C	chr20	33306904
kgp5909883_C	chr20	33276853

x <- c(79, 113, 32, 38)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))

x <- c(78, 113, 31, 38)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))

x <- c(69, 113, 30, 38)
t <- matrix(x, ncol=2, nrow=2)


x <- c(82, 113, 32, 38)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))

x <- c(68, 113, 30, 38)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))


x <- c(77, 113, 21, 38)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))


snp.overlap.temp = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_COPD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp)))
snp.overlap = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_COPD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc)


snp.overlap.temp.ild = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_ILD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp.ild)))
snp.overlap.ild = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_ILD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc)

snp.overlap.temp.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc.ctr = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp.ctr)))
snp.overlap.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc.ctr)


snpname <- "rs1847293_G" #homoz 0.0007 20/31, 0/15; heteroz 0.01 77/31, 18/15

snpname <- "rs1792753_T" #homoz 0.0001 heteroz 0.004

snpname <- "rs525770_C" #heteroz 0.001

snpname <- "rs475882_G" #heteroz 0.001

snpname <- "kgp4823739_A"

snp1.ctr <- snp.overlap.ctr[snpname,which(snp.overlap.ctr[snpname,]==1)]
snp2.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==2)]
snp.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==1 | snp.overlap.ctr[snpname,]==2)]
nosnp.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==0)]


snp1.copd <- snp.overlap[snpname, which(snp.overlap[snpname,]==1)]
snp2.copd <- snp.overlap[snpname,which(snp.overlap[snpname,]==2)]
snp.copd <- snp.overlap[snpname, which(snp.overlap[snpname,]==1 | snp.overlap[snpname,]==2)]
nosnp.copd <- snp.overlap[snpname, which(snp.overlap[snpname,]==0)]


x <- c(length(snp.copd), length(nosnp.copd), length(snp.ctr), length(nosnp.ctr))
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))


x <- c(length(snp2.copd), length(nosnp.copd), length(snp2.ctr), length(nosnp.ctr))
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))



length(which(snp.overlap["kgp1625425_C",]==1 | snp.overlap["kgp1625425_C",]==2)) #79/113
length(which(snp.overlap.ctr["kgp1625425_C",]==1 | snp.overlap.ctr["kgp1625425_C",]==2)) #32/38

x <- c(77, 111, 20, 38)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))
#0.4431

length(which(snp.overlap["rs1483025_A",]==1 | snp.overlap["rs1483025_A",]==2)) #79/113
length(which(snp.overlap.ctr["rs1483025_A",]==1 | snp.overlap.ctr["rs1483025_A",]==2)) #32/38
x <- c(79, 111, 20, 38)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))
#0.3615

length(which(snp.overlap["rs1847293_G",]==1 | snp.overlap["rs1847293_G",]==2)) #79/113
length(which(snp.overlap.ctr["rs1847293_G",]==1 | snp.overlap.ctr["rs1847293_G",]==2)) #32/38
x <- c(77, 111, 15, 38)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))
#0.11


snp.overlap.temp.all = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_All_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc.all = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp.all)))
snp.overlap.all = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_All_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc.ctr)
colnames(snp.overlap.all["kgp1625425_C",])

length(which(snp.overlap.all["rs1847293_G",]==1 | snp.overlap.all["rs1847293_G",]==2)) #172/256=67% freq
length(which(snp.overlap.all["rs1847293_G",]==0)) #84/256
x <- c(77, 111, 172-77, 256-111)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))
#0.84


length(which(snp.overlap["rs525770_C",]==1 | snp.overlap["rs525770_C",]==2)) #25/107 freq

length(which(snp.overlap.all["rs525770_C",]==1 | snp.overlap.all["rs525770_C",]==2)) #65/261 freq = 25%

x <- c(25, 107, 65-25, 261-107)
t <- matrix(x, ncol=2, nrow=2)
fisher.test(t(t))


snp_copd = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/COPD_SNPs.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)

rownames(snp_copd) <- snp_copd$SNP

snp_copd2 <-c()
for (snp in rownames(snp_copd)){
 
    a <- length(which(snp.overlap.all[snp,]==1 | snp.overlap.all[snp,]==2)) #65/261 freq = 25%
    #print(a)
    b <- length(which(!is.na(snp.overlap.all[snp,]))) #65/261 freq = 25%
    #print(b)
    #print(a/b)

    acopd <- length(which(snp.overlap.ild[snp,]==1 | snp.overlap.ild[snp,]==2)) 
    bcopd <- length(which(!is.na(snp.overlap.ild[snp,])))

    actr <- length(which(snp.overlap.ctr[snp,]==1 | snp.overlap.ctr[snp,]==2)) 
    bctr <- length(which(!is.na(snp.overlap.ctr[snp,]))) 
 
    #x <- c(acopd, bcopd-acopd, a, b-a)
    x <- c(acopd, bcopd-acopd, actr, bctr-actr)
    t <- matrix(x, ncol=2, nrow=2)
    ft=fisher.test(t(t))    
    #print(ft$p.value)
    if (ft$p.value<0.1){
       snp_copd2 <- rbind(snp_copd2, c(snp, ft$p.value))
    }
    
    #if ((a/b) > 0.05){
    #    snp_copd2 <- rbind(snp_copd2, snp.overlap.all[snp,])
    #}
}

snp_copd_all <- snp_copd2
#rs1145786_G p=0.07 


snp_copd_ctr <- snp_copd2
#rs525770_C p<0.005
#rs475882_G p<0.005

snp_ild_ctr <- snp_copd2
#none

snp_ild_all <- snp_copd2
#kgp5978297_C p=0.05


#tr = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/Triplets_eQTL_with_CIT_COPD.txt', header=TRUE, sep="\t")

tr = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/Triplets_eQTL_with_CIT_COPD_cis_genes2.txt', header=TRUE, sep="\t")

snp_copd = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/COPD_SNPs.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
rownames(snp_copd) <- snp_copd$SNP

mirs <- c("MI0003673_MIMAT0003327", "MI0001648_MIMAT0001541", "MI0003823_MIMAT0010251", "MI0000743_MIMAT0004677", "MI0000742_MIMAT0004676", "MI0000743_MIMAT0000686", "MI0000742_MIMAT0000685")

tr <- unique(tr)
tr[which(tr[,'SNP'] %in% rownames(snp_copd) & tr[,'Gene'] %in% mirs), ]

#tr[which(tr[,'SNP'] %in% rownames(snp_copd) & tr[,'Mir'] %in% mirs), ]


#scale free
tr = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/Triplets_eQTL_with_CIT_Control.txt', header=TRUE, sep="\t")
t1 = unique(tr[,1:2])
t2 = unique(tr[,2:3])
snips=sort(table(t1[,1]))
mirna=sort(table(c(t1[,2],t2[,1])))
genes=sort(table(t2[,2]))
#plot(hist(log(c(snips,mirna, genes))), col="blue", main="Frequency of node degree \n SNP/miRNA/mRNA CIT network: COPD", xlab="node degree", ylab="log frequency")
#log.frequencies <- log10(c(snips,mirna,genes)) 
#plot(density(log.frequencies), col="blue", main="Frequency of node degree \n SNP/miRNA/mRNA CIT network: COPD", xlab="log node degree", ylab="frequency", xlim=c(0.1, max(log.frequencies)), lwd=2)

log.frequencies <- c(snips,mirna,genes)
table.frequencies=table(log.frequencies)
plot(log10(as.numeric(names(table.frequencies))), log10(table.frequencies), col="blue", main=paste("Frequency of node degree \n SNP/miRNA/mRNA CIT network: Control \n r = ", round(cor(log10(as.numeric(names(table.frequencies))), log10(table.frequencies), method="spearman"), digits=2)), xlab="log node degree", ylab="log frequency", lwd=1, pch=19)#, type="l")#, yaxt='n')
#max.frequencies = max(table.frequencies)
#y.values=seq(0, max.frequencies, by=max.frequencies/10)
#axis(side=2, at=y.values, labels=round(y.values))


tr = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_ILD/Triplets_eQTL_with_CIT_ILD.txt', header=TRUE, sep="\t")
t1 = unique(tr[,1:2])
t2 = unique(tr[,2:3])
snips=sort(table(t1[,1]))
mirna=sort(table(c(t1[,2],t2[,1])))
genes=sort(table(t2[,2]))
plot(hist(log(c(snips,mirna, genes))), col="red", main="Frequency of node degree \n SNP/miRNA/mRNA causal network", xlab="node degree", ylab="log frequency")

tr = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/Triplets_eQTL_with_CIT_Control.txt', header=TRUE, sep="\t")
t1 = unique(tr[,1:2])
t2 = unique(tr[,2:3])
snips=sort(table(t1[,1]))
mirna=sort(table(c(t1[,2],t2[,1])))
genes=sort(table(t2[,2]))
plot(hist(log(c(snips,mirna, genes))), col="red", main="Frequency of node degree \n SNP/miRNA/mRNA causal network", xlab="node degree", ylab="log frequency")


a <- c()
for (i in 1:ncol(cor.network$FDR)){
    for (j in 1:nrow(cor.network$FDR)){
        if (cor.network$FDR[j,i] < 0.1){
            a <- rbind(a, c(rownames(cor.network$FDR)[j], colnames(cor.network$FDR)[i]))
        }
    }
}


mirgene = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/cortable_control.txt', header=TRUE, sep="\t")
cis = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/eQTL_Cis_corr.txt.gz', header=TRUE, sep="\t")
trans = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/eQTL_Trans_corr.txt.gz', header=TRUE, sep="\t")

cis.trans=rbind(cis, trans)
cis.trans.filtered=subset(cis.trans, p.value<1e-4)[,1:2]
colnames(mirgene)=colnames(cis.trans.filtered)
mirgene2 <- mirgene[which((mirgene[,1] %in% cis.trans[,2]) & (mirgene[,2] %in% cis.trans[,2])),]
cis.trans.filtered.edges=rbind(cis.trans.filtered, mirgene2)
cis.trans.frequencies = table(c(cis.trans.filtered.edges[,1], cis.trans.filtered.edges[,2]))
cis.trans.table = table(cis.trans.frequencies)#/nrow(cis.trans.frequencies)

plot(log10(as.numeric(names(cis.trans.table))), log10(cis.trans.table), col="red", main=paste("Frequency of node degree \n SNP/miRNA/mRNA entire network: Control \n r = ", round(cor(log10(as.numeric(names(cis.trans.table))), log10(cis.trans.table), method="pearson"), digits=2)), xlab="log node degree", ylab="log frequency", lwd=1, pch=19)#, type="l")

cis.trans.filtered.edges2=mirgene
cis.trans.frequencies2 = table(c(cis.trans.filtered.edges2[,1], cis.trans.filtered.edges2[,2]))
cis.trans.table2 = table(cis.trans.frequencies2)#/nrow(cis.trans.frequencies)

plot(log10(as.numeric(names(cis.trans.table2))), log10(cis.trans.table2), col="magenta", main=paste("Frequency of node degree \n miRNA/mRNA network: Control \n r = ", round(cor(log10(as.numeric(names(cis.trans.table2))), log10(cis.trans.table2), method="pearson"), digits=2)), xlab="log node degree", ylab="log frequency", lwd=1, pch=19)#, type="l")

plot(as.numeric(names(cis.trans.table)), (cis.trans.table), col="magenta", main="Frequency of node degree \n miRNA/mRNA network: COPD", xlab="node degree", ylab="relative frequency", lwd=1, type="l", yaxt='n')
max.frequencies.entire = max(cis.trans.table)
y.values.entire=seq(0, max.frequencies.entire, by=max.frequencies.entire/10)
axis(side=2, at=y.values.entire, labels=(y.values.entire))



#SNP/miR/mRNA all

mirgene = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/cortable_control.txt', header=TRUE, sep="\t")
cis = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/eQTL_Cis_corr.txt.gz', header=TRUE, sep="\t")
trans = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_Control/eQTL_Trans_corr.txt.gz', header=TRUE, sep="\t")

cis.trans=rbind(cis, trans)
cis.trans.filtered=subset(cis.trans, p.value<1e-4)[,1:2]
colnames(mirgene)=colnames(cis.trans.filtered)
mirgene2 <- mirgene[which((mirgene[,1] %in% cis.trans[,2]) & (mirgene[,2] %in% cis.trans[,2])),]
cis.trans.filtered.edges=rbind(cis.trans.filtered, mirgene2)
cis.trans.frequencies = table(c(cis.trans.filtered.edges[,1], cis.trans.filtered.edges[,2]))
cis.trans.table = table(cis.trans.frequencies)#/nrow(cis.trans.frequencies)

plot(log10(as.numeric(names(cis.trans.table))), log10(cis.trans.table), col="red", main=paste("Frequency of node degree \n SNP/miRNA/mRNA entire network: Control \n r = ", round(cor(log10(as.numeric(names(cis.trans.table))), log10(cis.trans.table), method="pearson"), digits=3)), xlab="log node degree", ylab="log frequency", lwd=1, pch=19)#, type="l")


mirgene = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/cortable_copd.txt', header=TRUE, sep="\t")
cis = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/eQTL_Cis_corr.txt.gz', header=TRUE, sep="\t")
trans = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_COPD/eQTL_Trans_corr.txt.gz', header=TRUE, sep="\t")

cis.trans=rbind(cis, trans)
cis.trans.filtered=subset(cis.trans, p.value<1e-4)[,1:2]
colnames(mirgene)=colnames(cis.trans.filtered)
mirgene2 <- mirgene[which((mirgene[,1] %in% cis.trans[,2]) & (mirgene[,2] %in% cis.trans[,2])),]
cis.trans.filtered.edges=rbind(cis.trans.filtered, mirgene2)
cis.trans.frequencies = table(c(cis.trans.filtered.edges[,1], cis.trans.filtered.edges[,2]))
cis.trans.table = table(cis.trans.frequencies)#/nrow(cis.trans.frequencies)

plot(log10(as.numeric(names(cis.trans.table))), log10(cis.trans.table), col="red", main=paste("Frequency of node degree \n SNP/miRNA/mRNA entire network: COPD \n r = ", round(cor(log10(as.numeric(names(cis.trans.table))), log10(cis.trans.table), method="pearson"), digits=3)), xlab="log node degree", ylab="log frequency", lwd=1, pch=19)#, type="l")


mirgene = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/cortable_ild.txt', header=TRUE, sep="\t")
cis = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_ILD/eQTL_Cis_corr.txt.gz', header=TRUE, sep="\t")
trans = read.table('/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/eQTL_262_output_PC_ILD/eQTL_Trans_corr.txt.gz', header=TRUE, sep="\t")

cis.trans=rbind(cis, trans)
cis.trans.filtered=subset(cis.trans, p.value<1e-4)[,1:2]
colnames(mirgene)=colnames(cis.trans.filtered)
mirgene2 <- mirgene[which((mirgene[,1] %in% cis.trans[,2]) & (mirgene[,2] %in% cis.trans[,2])),]
cis.trans.filtered.edges=rbind(cis.trans.filtered, mirgene2)
cis.trans.frequencies = table(c(cis.trans.filtered.edges[,1], cis.trans.filtered.edges[,2]))
cis.trans.table = table(cis.trans.frequencies)#/nrow(cis.trans.frequencies)

plot(log10(as.numeric(names(cis.trans.table))), log10(cis.trans.table), col="red", main=paste("Frequency of node degree \n SNP/miRNA/mRNA entire network: ILD \n r = ", round(cor(log10(as.numeric(names(cis.trans.table))), log10(cis.trans.table), method="pearson"), digits=3)), xlab="log node degree", ylab="log frequency", lwd=1, pch=19)#, type="l")

