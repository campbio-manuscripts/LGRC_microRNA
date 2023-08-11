
#COPD
snp.overlap.temp = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_COPD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp)))
snp.overlap = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_COPD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc)


snp.overlap.temp.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc.ctr = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp.ctr)))
snp.overlap.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc.ctr)


#cis SNPs copd
snpname <- "kgp2263580_T" #"kgp2263580_T" #"kgp4823739_A"
snp1.ctr <- snp.overlap.ctr[snpname,which(snp.overlap.ctr[snpname,]==1)]
snp2.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==2)]
snp.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==1 | snp.overlap.ctr[snpname,]==2)]
nosnp.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==0)]

snp1.copd <- snp.overlap[snpname, which(snp.overlap[snpname,]==1)]
snp2.copd <- snp.overlap[snpname,which(snp.overlap[snpname,]==2)]
snp.copd <- snp.overlap[snpname, which(snp.overlap[snpname,]==1 | snp.overlap[snpname,]==2)]
nosnp.copd <- snp.overlap[snpname, which(snp.overlap[snpname,]==0)]

x <- c(length(snp.copd), length(nosnp.copd), length(snp.ctr), length(nosnp.ctr))
t1 <- matrix(x, ncol=2, nrow=2)
ft1 <- fisher.test(t(t1), alternative="greater")

y <- c(length(snp2.copd), length(nosnp.copd)+length(snp1.copd), length(snp2.ctr), length(nosnp.ctr)+length(snp1.ctr))
t2 <- matrix(y, ncol=2, nrow=2)
ft2 <- fisher.test(t(t2), alternative="greater")

snpname <- "kgp2263580_T"

mir.copd = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/MicroRNA/LGRC_MicroRNA_Overlap_SNP_COPD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))
gene.copd = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/Gene/LGRC_Gene_Overlap_SNP_COPD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/gene_mir_copd.pdf")
plot(gene.copd["A_23_P64712",], mir.copd["MI0001648_MIMAT0001541",], main=c("Cor coef ", round(cor(gene.copd["A_23_P64712",], mir.copd["MI0001648_MIMAT0001541",]), 3)), xlab="TCTN2", ylab="hsa-miR-449a", pch=19, col="blue")
abline(lm(as.numeric(mir.copd["MI0001648_MIMAT0001541",])~ as.numeric(gene.copd["A_23_P64712",])))
dev.off()


pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/gene_mir_copd4.pdf")
plot(gene.copd["A_23_P64712",], mir.copd["MI0003823_MIMAT0010251",], main=c("Cor coef ", round(cor(gene.copd["A_23_P64712",], mir.copd["MI0003823_MIMAT0010251",]), 3)), xlab="TCTN2", ylab="hsa-miR-449c-5p", pch=19, col="blue")
abline(lm(as.numeric(mir.copd["MI0003823_MIMAT0010251",])~ as.numeric(gene.copd["A_23_P64712",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/snp_mir5.pdf")
boxplot(mir.copd["MI0001648_MIMAT0001541", which(snp.overlap["kgp2263580_T",]!=2)], mir.copd["MI0001648_MIMAT0001541", which(snp.overlap["kgp2263580_T",]==2)], col=c("green", "red"), names=c("common SNP", "homozygous kgp2263580_T"), ylab=c("hsa-miR-449a expression"))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/snp_gene_5.pdf")
boxplot(gene.copd["A_23_P64712", which(snp.overlap["kgp2263580_T",]!=2)], gene.copd["A_23_P64712", which(snp.overlap["kgp2263580_T",]==2)], col=c("darkgreen", "darkred"), names=c("common SNP", "homozygous kgp2263580_T"), ylab=c("TCTN2 expression"))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/snp_mir6.pdf")
boxplot(mir.copd["MI0001648_MIMAT0001541", which(snp.overlap["kgp4823739_A",]!=2)], mir.copd["MI0001648_MIMAT0001541", which(snp.overlap["kgp4823739_A",]==2)], col=c("green", "red"), names=c("common SNP", "homozygous kgp4823739_A"), ylab=c("hsa-miR-449a expression"))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/snp_gene6.pdf")
boxplot(gene.copd["A_23_P64712", which(snp.overlap["kgp4823739_A",]!=2)], gene.copd["A_23_P64712", which(snp.overlap["kgp4823739_A",]==2)], col=c("darkgreen", "darkred"), names=c("common SNP", "homozygous kgp4823739_A"), ylab=c("TCTN2 expression"))
dev.off()



# A_23_P64712
# hsa-miR-449a	MI0001648_MIMAT0001541 #X
# hsa-miR-34b-5p	MI0000742_MIMAT0000685 
# hsa-miR-34b-3p	MI0000742_MIMAT0004676 #X
# hsa-miR-34c-5p	MI0000743_MIMAT0000686 #X
# hsa-miR-449c-5p	MI0003823_MIMAT0010251 #X
# hsa-miR-190b	MI0005545_MIMAT0004929

#All SNPs COPD

snp.overlap.temp = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_COPD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp)))
snp.overlap = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_COPD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc)

snp.overlap.temp.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc.ctr = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp.ctr)))

snp.overlap.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc.ctr)

snpname_copd = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/COPD_SNPs.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)

snp_copd1 <-c()
snp_copd2 <-c()

p1 <- c()
p2 <- c()
p3 <- c()
n <- c()

for (snpname in snpname_copd$SNP){

    snp1.ctr <- snp.overlap.ctr[snpname,which(snp.overlap.ctr[snpname,]==1)]
    snp2.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==2)]
    snp.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==1 | snp.overlap.ctr[snpname,]==2)]
    nosnp.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==0)]

    snp1.copd <- snp.overlap[snpname, which(snp.overlap[snpname,]==1)]
    snp2.copd <- snp.overlap[snpname,which(snp.overlap[snpname,]==2)]
    snp.copd <- snp.overlap[snpname, which(snp.overlap[snpname,]==1 | snp.overlap[snpname,]==2)]
    nosnp.copd <- snp.overlap[snpname, which(snp.overlap[snpname,]==0)]

    x <- c(length(snp.copd), length(nosnp.copd), length(snp.ctr), length(nosnp.ctr))
    t1 <- matrix(x, ncol=2, nrow=2)
    ft1 <- fisher.test(t(t1), alternative="greater")

    p1 <- c(p1, ft1$p.value)
    
    y <- c(length(snp2.copd), length(nosnp.copd)+length(snp1.copd), length(snp2.ctr), length(nosnp.ctr)+length(snp1.ctr))
    t2 <- matrix(y, ncol=2, nrow=2)
    ft2 <- fisher.test(t(t2), alternative="greater")

    p2 <- c(p2, ft2$p.value)
    
    z <- c(length(snp1.copd), length(snp2.copd), length(nosnp.copd), length(snp1.ctr), length(snp2.ctr), length(nosnp.ctr))
    t3 <- matrix(z, ncol=2, nrow=3)
    ft3 <- fisher.test(t(t3))
    n <- rbind(n, z)

    p3 <- c(p3, ft3$p.value)
   
    if (ft1$p.value<0.05){
       snp_copd1 <- rbind(snp_copd1, c(snpname, ft1$p.value, x))
    }

    if (ft2$p.value<0.05){
       snp_copd2 <- rbind(snp_copd2, c(snpname, ft2$p.value, y))
    }
 
}

q1 <- p.adjust(p1, method="fdr")
q2 <- p.adjust(p2, method="fdr")
q3 <- p.adjust(p3, method="fdr")

snps_res <- cbind(snpname_copd$SNP, as.numeric(p1), as.numeric(q1), as.numeric(p2), as.numeric(q2), as.numeric(p3), as.numeric(q3), n)
colnames(snps_res) <- c("SNP", "p_heteroz", "q_heteroz", "p_homoz", "q_homoz", "p_3classes", "q_3classes", "number_heteroz_COPD", "number_homoz_COPD", "number_no_snp_COPD", "number_heteroz_ctr", "number_homoz_ctr", "number_no_snp_ctr") 

write.table(snps_res, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/snps_res_copd.csv", row.names=FALSE, sep=",")


#COPD

#snp_copd1
> snp_copd1
      [,1]           [,2]                   [,3] [,4] [,5] [,6]
 [1,] "kgp1625425_C" "0.0486058200372168"   "77" "34" "20" "18" chr11
 [2,] "rs1483025_A"  "0.0259096763007725"   "79" "31" "20" "18" chr11
 [3,] "rs1847293_G"  "0.00651405225234148"  "77" "31" "15" "18" chr11
 [4,] "rs1785864_A"  "0.0184576270757168"   "81" "30" "20" "18" chr11
 [5,] "rs1699081_G"  "0.0307110881766364"   "79" "32" "20" "18" chr11
 [6,] "rs1792753_T"  "0.00281583252221874"  "81" "25" "16" "17" chr11
 [7,] "rs525770_C"   "0.00056063202884204"  "25" "82" "0"  "33" chr13
 [8,] "rs475882_G"   "0.000646587053859013" "26" "82" "0"  "31" chr13
 [9,] "kgp8017900_C" "0.0473351213131134"   "69" "42" "16" "20" chr3


 #snp_copd2
 > snp_copd2
      [,1]           [,2]                   [,3] [,4]  [,5] [,6]
 [1,] "rs1847293_G"  "0.00304540673621405"  "20" "88"  "0"  "33" chr11
 [2,] "rs4147268_T"  "0.00500168026350384"  "18" "93"  "0"  "35" chr1
 [3,] "rs4436447_A"  "0.0054573657735184"   "18" "92"  "0"  "34" chr1
 [4,] "rs587985_C"   "0.0236663565272242"   "13" "95"  "0"  "34" chr11
 [5,] "rs6796037_T"  "0.0192920648667593"   "14" "97"  "0"  "34" chr3
 [6,] "rs9814167_C"  "0.0192920648667593"   "14" "97"  "0"  "34" chr3
 [7,] "rs9877255_C"  "0.0192920648667593"   "14" "97"  "0"  "34" chr3
 [8,] "kgp7153288_G" "0.0192920648667593"   "14" "97"  "0"  "34" chr3
 [9,] "rs1792753_T"  "0.000733537233998933" "24" "82"  "0"  "33" chr11
[10,] "kgp3700677_G" "0.00570504155055906"  "18" "93"  "0"  "34" chr8
[11,] "kgp5237943_C" "0.0400175319252243"   "11" "100" "0"  "36" chr2
[12,] "kgp4164025_A" "0.0320481911306162"   "12" "96"  "0"  "34" chr8
[13,] "rs727581_G"   "0.0320401807898038"   "12" "99"  "0"  "35" chr8
[14,] "rs4676866_A"  "0.0472166745129795"   "10" "101" "0"  "38" chr3
[15,] "kgp8017900_C" "0.0157809845681123"   "14" "97"  "0"  "36" chr3
[16,] "rs17009077_G" "0.0432542440662351"   "11" "100" "0"  "35" chr3












#ILD

snp.overlap.temp.ild = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_ILD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp.ild)))
snp.overlap.ild = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_ILD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc)

snp.overlap.temp.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc.ctr = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp.ctr)))
snp.overlap.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc.ctr)

# cis SNP ild
snpname <- "rs618555_T" #"rs618555_T" #"kgp7779523_T"
snp1.ctr <- snp.overlap.ctr[snpname,which(snp.overlap.ctr[snpname,]==1)]
snp2.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==2)]
snp.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==1 | snp.overlap.ctr[snpname,]==2)]
nosnp.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==0)]

snp1.ild <- snp.overlap.ild[snpname, which(snp.overlap.ild[snpname,]==1)]
snp2.ild <- snp.overlap.ild[snpname,which(snp.overlap.ild[snpname,]==2)]
snp.ild <- snp.overlap.ild[snpname, which(snp.overlap.ild[snpname,]==1 | snp.overlap.ild[snpname,]==2)]
nosnp.ild <- snp.overlap.ild[snpname, which(snp.overlap.ild[snpname,]==0)]

x <- c(length(snp.ild), length(nosnp.ild), length(snp.ctr), length(nosnp.ctr))
t1 <- matrix(x, ncol=2, nrow=2)
ft1 <- fisher.test(t(t1), alternative="greater")

y <- c(length(snp2.ild), length(nosnp.ild)+length(snp1.ild), length(snp2.ctr), length(nosnp.ctr)+length(snp1.ctr))
t2 <- matrix(y, ncol=2, nrow=2)
ft2 <- fisher.test(t(t2), alternative="greater")

#rs618555_T ft2 p = 0.06

mir.ild = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/MicroRNA/LGRC_MicroRNA_Overlap_SNP_ILD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))
gene.ild = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/Gene/LGRC_Gene_Overlap_SNP_ILD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/gene_mir_ild.pdf")
plot(gene.ild["A_23_P45751",], mir.ild["MI0003823_MIMAT0010251",], main=c("Cor coef ", round(cor(gene.ild["A_23_P45751",], mir.ild["MI0003823_MIMAT0010251",]), 3)), xlab="CLCA4", ylab="hsa-miR-449c-5p", pch=19, col="blue")
abline(lm(as.numeric(mir.ild["MI0003823_MIMAT0010251",])~ as.numeric(gene.ild["A_23_P45751",])))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/gene_mir_ild_mir92b.pdf")
plot(gene.ild["A_23_P45751",], mir.ild["MI0003560_MIMAT0004792",], main=c("Cor coef ", round(cor(gene.ild["A_23_P45751",], mir.ild["MI0003560_MIMAT0004792",]), 3)), xlab="CLCA4", ylab="hsa-miR-92b-5p", pch=19, col="blue")
abline(lm(as.numeric(mir.ild["MI0003560_MIMAT0004792",])~ as.numeric(gene.ild["A_23_P45751",])))
dev.off()

#hsa-miR-449c-5p	MI0003823_MIMAT0010251
#hsa-miR-92b-5p	MI0003560_MIMAT0004792

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/snp_mir_ild2.pdf")
boxplot(mir.ild["MI0003823_MIMAT0010251", which(snp.overlap.ild["rs618555_T",]!=2)], mir.ild["MI0003823_MIMAT0010251", which(snp.overlap.ild["rs618555_T",]==2)], col=c("green", "red"), names=c("common SNP", "homozygous rs618555_T"), ylab=c("hsa-miR-449c-5p expression"))
dev.off()



pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/snp_mir_ild_mir92b.pdf")
boxplot(mir.ild["MI0003560_MIMAT0004792", which(snp.overlap.ild["rs618555_T",]!=2)], mir.ild["MI0003560_MIMAT0004792", which(snp.overlap.ild["rs618555_T",]==2)], col=c("green", "red"), names=c("common SNP", "homozygous rs618555_T"), ylab=c("hsa-miR-92b-5p expression"))
dev.off()


pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/snp_gene_ild2.pdf")
boxplot(gene.ild["A_23_P45751", which(snp.overlap.ild["rs618555_T",]!=2)], gene.ild["A_23_P45751", which(snp.overlap.ild["rs618555_T",]==2)], col=c("darkgreen", "darkred"), names=c("common SNP", "homozygous rs618555_T"), ylab=c("CLCA4 expression"))
dev.off()

pdf("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/snp_gene_ild4.pdf")
boxplot(gene.ild["A_23_P360990", which(snp.overlap.ild["kgp7779523_T",]==0)], gene.ild["A_23_P360990", which(snp.overlap.ild["kgp7779523_T",]!=0)], col=c("darkgreen", "darkred"), names=c("common SNP", "kgp7779523_T"), ylab=c("DYDC1 expression"))
dev.off()

#kgp7779523_T	A_23_P360990


# all SNPs ild

snp.overlap.temp.ild = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_ILD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp.ild)))
snp.overlap.ild = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_ILD_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc)

snp.overlap.temp.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=100)
cc.ctr = rep(c("character", "integer"), c(1, ncol(snp.overlap.temp.ctr)))
snp.overlap.ctr = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP/LGRC_SNP_Control_Samples_GroupMin5.txt", row.names=1, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses=cc.ctr)

snpname_ild = read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/ILD_SNPs.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)

snp_ild1 <-c()
snp_ild2 <-c()
p1 <- c()
p2 <- c()
p3 <- c()
n <- c()

for (snpname in snpname_ild$SNP){

    snp1.ctr <- snp.overlap.ctr[snpname,which(snp.overlap.ctr[snpname,]==1)]
    snp2.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==2)]
    snp.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==1 | snp.overlap.ctr[snpname,]==2)]
    nosnp.ctr <- snp.overlap.ctr[snpname, which(snp.overlap.ctr[snpname,]==0)]

    snp1.ild <- snp.overlap.ild[snpname, which(snp.overlap.ild [snpname,]==1)]
    snp2.ild <- snp.overlap.ild[snpname,which(snp.overlap.ild [snpname,]==2)]
    snp.ild <- snp.overlap.ild[snpname, which(snp.overlap.ild [snpname,]==1 | snp.overlap.ild [snpname,]==2)]
    nosnp.ild <- snp.overlap.ild[snpname, which(snp.overlap.ild [snpname,]==0)]

    
    
    
    x <- c(length(snp.ild), length(nosnp.ild), length(snp.ctr), length(nosnp.ctr))
    t1 <- matrix(x, ncol=2, nrow=2)
    ft1 <- fisher.test(t(t1), alternative="greater")

    p1 <- c(p1, ft1$p.value)
    
    y <- c(length(snp2.ild), length(nosnp.ild)+length(snp1.ild), length(snp2.ctr), length(nosnp.ctr)+length(snp1.ctr))
    t2 <- matrix(y, ncol=2, nrow=2)
    ft2 <- fisher.test(t(t2), alternative="greater")

    p2 <- c(p2, ft2$p.value)
        
    z <- c(length(snp1.ild), length(snp2.ild), length(nosnp.ild), length(snp1.ctr), length(snp2.ctr), length(nosnp.ctr))
    t3 <- matrix(z, ncol=2, nrow=3)
    ft3 <- fisher.test(t(t3))
    n <- rbind(n, z)

    p3 <- c(p3, ft3$p.value)
    
    
    
    if (ft1$p.value<0.05){
       snp_ild1 <- rbind(snp_ild1, c(snpname, ft1$p.value, x))
    }
    #print(ft1$p.value)

    if (ft2$p.value<0.05){
       snp_ild2 <- rbind(snp_ild2, c(snpname, ft2$p.value, y))
    }
    #print(ft2$p.value)
}


q1 <- p.adjust(p1, method="fdr")
q2 <- p.adjust(p2, method="fdr")
q3 <- p.adjust(p3, method="fdr")

snps_res_ild <- cbind(snpname_ild$SNP, as.numeric(p1), as.numeric(q1), as.numeric(p2), as.numeric(q2), as.numeric(p3), as.numeric(q3), n)
colnames(snps_res_ild) <- c("SNP", "p_heteroz", "q_heteroz", "p_homoz", "q_homoz", "p_3classes", "q_3classes", "number_heteroz_ILD", "number_homoz_ILD", "number_no_snp_ILD", "number_heteroz_ctr", "number_homoz_ctr", "number_no_snp_ctr") 

write.table(snps_res_ild, "/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/snps_res_ild.csv", row.names=FALSE, sep=",")



#> snp_ild1
     [,1]            [,2]                 [,3] [,4] [,5] [,6]
[1,] "rs619146_C"    "0.0353678398140796" "61" "52" "13" "24" chr5
[2,] "kgp6113667_C"  "0.0432761073491943" "59" "53" "12" "23" chr17
[3,] "kgp11194278_A" "0.0263773802658788" "81" "32" "20" "18" chr16

#> snp_ild2
     [,1]            [,2]                   [,3] [,4]  [,5] [,6]
     "rs618555_T"   "0.06097"   "10"    "103"   "0" "35" (this SNP is CIS to CLCA4 gene ~ 300kb distance chr1; CLCA4 expression is positively associated to mirR-449c-5p expression)
[1,] "kgp10062335_G" "0.0358850024845802"   "12" "100" "0"  "34" chr3
[2,] "kgp6639144_T"  "0.0338815817531201"   "12" "101" "0"  "35" chr20
[3,] "kgp9719798_T"  "0.0174421408384399"   "14" "97"  "0"  "35" chr18
[4,] "rs7861757_A"   "0.0286606346588138"   "12" "101" "0"  "37" chr9
[5,] "kgp6113667_C"  "0.00976558534677444"  "16" "96"  "0"  "35" chr17
[6,] "rs17530621_G"  "0.00842374919225583"  "17" "96"  "0"  "34" chr4
[7,] "rs7641416_A"   "0.000668063147265113" "24" "89"  "0"  "36" chr3
[8,] "kgp1100673_G"  "0.00622061478812738"  "18" "95"  "0"  "34" chr11



z <- c(24/113, 0/38, 50/113, 20/38, 39/113, 16/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("ILD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_rs7641416_A_ILD.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()

barplot(c(24/113, 0/38, 50/113, 20/38, 39/113, 16/38), col=c("red", "green"), cex.axis=2.5)

50	24	39	20	0	16

number_heteroz_ILD	number_homoz_ILD	number_no_snp_ILD	number_heteroz_ctr	number_homoz_ctr	number_no_snp_ctr

43	16	53	12	0	23

z <- c(16/113, 0/38, 43/113, 12/38, 53/113, 23/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("ILD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_kgp6113667_C.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()


45	17	51	19	0	15


z <- c(17/113, 0/38, 45/113, 19/38, 51/113, 15/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("ILD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_rs17530621_G.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()


49	5	59	10	7	20

z <- c(5/113, 7/38, 49/113, 10/38, 59/113, 20/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("ILD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_rs8061907_G.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()


58	18	37	18	0	16

z <- c(18/113, 0/38, 58/113, 18/38, 37/113, 16/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("ILD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_rkgp1100673_G_ILD.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()


57	24	25	16	0	17

z <- c(24/111, 0/38, 57/111, 16/38, 25/111, 17/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("COPD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_rs1792753_T_COPD.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()

25	0	82	0	0	33

z <- c(0/111, 0/38, 25/111, 0/38, 82/111, 33/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("COPD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_rs525770_C_COPD.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()


26	0	82	0	0	31

z <- c(0/111, 0/38, 26/111, 0/38, 82/111, 31/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("COPD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_rs475882_G_COPD.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()

57	20	31	15	0	18
z <- c(20/111, 0/38, 57/111, 15/38, 31/111, 18/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("COPD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_rs1847293_G_COPD.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()



51	18	42	20	0	15
z <- c(18/111, 0/38, 51/111, 20/38, 42/111, 15/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("COPD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_rs4147268_T_COPD.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()



42	7	62	19	7	12

z <- c(7/111, 7/38, 42/111, 19/38, 62/111, 12/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("COPD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_rs1145786_G_COPD.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()


49	18	43	16	0	18

z <- c(18/111, 0/38, 49/111, 16/38, 43/111, 18/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("COPD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_rs4436447_A_COPD.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()


42	18	51	18	0	16

z <- c(18/111, 0/38, 42/111, 18/38, 51/111, 16/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("COPD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_kgp3700677_G_COPD.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()



55	14	42	16	0	20

z <- c(14/111, 0/38, 55/111, 16/38, 42/111, 20/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("COPD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_kgp8017900_C_COPD.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()



42	11	58	21	0	14

55	14	42	16	0	20

z <- c(11/111, 0/38, 42/111, 21/38, 58/111, 14/38)
tab <- matrix(z, ncol=3, nrow=2)
rownames(tab) <- c("COPD", "CONTROL")
jpeg("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/SNP_rs17009077_G_COPD.jpeg")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(tab,margin=2)
barplot(prop, col=heat.colors(nrow(prop)), width=2,  cex.axis=2.5)
legend("topright",inset=c(-0.28,0), fill=heat.colors(nrow(prop)), legend=rownames(tab))
dev.off()
