copd = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/MicroRNA/LGRC_MicroRNA_Overlap_SNP_COPD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))

ild = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/MicroRNA/LGRC_MicroRNA_Overlap_SNP_ILD_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))

ctr = as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_sandbox/LGRC_Data/Preprocess/MicroRNA/LGRC_MicroRNA_Overlap_SNP_Control_Samples.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))

pheno <- as.matrix(read.table("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/111228_LGRC_demographics.txt", sep="\t", quote="", header=TRUE, row.names=1, stringsAsFactors=FALSE))

pheno_copd <- pheno[colnames(copd),]
length(pheno_copd[which(pheno_copd[,"Cigarette.smoking"]=="2-Ever (&gt;100)"),"Cigarette.smoking"])
length(pheno_copd[which(pheno_copd[,"Cigarette.smoking"]=="3-Never"),"Cigarette.smoking"])
length(pheno_copd[which(pheno_copd[,"Cigarette.smoking"]=="1-Current"),"Cigarette.smoking"])


pheno_ild <- pheno[colnames(ild),]
length(pheno_ild[which(pheno_ild[,"Cigarette.smoking"]=="2-Ever (&gt;100)"),"Cigarette.smoking"])
length(pheno_ild[which(pheno_ild[,"Cigarette.smoking"]=="3-Never"),"Cigarette.smoking"])
length(pheno_ild[which(pheno_ild[,"Cigarette.smoking"]=="1-Current"),"Cigarette.smoking"])

pheno_ctr <- pheno[colnames(ctr),]
length(pheno_ctr[which(pheno_ctr[,"Cigarette.smoking"]=="2-Ever (&gt;100)"),"Cigarette.smoking"])
length(pheno_ctr[which(pheno_ctr[,"Cigarette.smoking"]=="3-Never"),"Cigarette.smoking"])
length(pheno_ctr[which(pheno_ctr[,"Cigarette.smoking"]=="1-Current"),"Cigarette.smoking"])


ta = rbind(c(2, 24, 12), c(4, 71, 38))
htest = fisher.test(ta) 
htest$p.value 

ta = rbind(c(2, 24, 12), c(7, 99, 5))
htest = fisher.test(ta) 
htest$p.value 

ta = rbind(c(4, 71, 38), c(7, 99, 5))
htest = fisher.test(ta) 
htest$p.value 

mean(as.numeric(pheno_copd[, "age"]))
sd(as.numeric(pheno_copd[, "age"]))

mean(as.numeric(pheno_ild[, "age"]))
sd(as.numeric(pheno_ild[, "age"]))

mean(as.numeric(pheno_ctr[, "age"]))
sd(as.numeric(pheno_ctr[, "age"]))

t.test(as.numeric(pheno_ctr[, "age"]), as.numeric(pheno_ild[, "age"]))
t.test(as.numeric(pheno_ctr[, "age"]), as.numeric(pheno_copd[, "age"]))
t.test(as.numeric(pheno_copd[, "age"]), as.numeric(pheno_ild[, "age"]))


mean(as.numeric(pheno_copd[, "Pack.years"]), na.rm=TRUE)
sd(as.numeric(pheno_copd[, "Pack.years"]), na.rm=TRUE)

mean(as.numeric(pheno_ild[, "Pack.years"]), na.rm=TRUE)
sd(as.numeric(pheno_ild[, "Pack.years"]), na.rm=TRUE)

mean(as.numeric(pheno_ctr[, "Pack.years"]), na.rm=TRUE)
sd(as.numeric(pheno_ctr[, "Pack.years"]), na.rm=TRUE)

t.test(as.numeric(pheno_ctr[, "Pack.years"]), as.numeric(pheno_ild[, "Pack.years"]))
t.test(as.numeric(pheno_ctr[, "Pack.years"]), as.numeric(pheno_copd[, "Pack.years"]))
t.test(as.numeric(pheno_copd[, "Pack.years"]), as.numeric(pheno_ild[, "Pack.years"]))


length(pheno_copd[which(pheno_copd[,"Gender"]=="1-Male"),"Gender"])
length(pheno_copd[which(pheno_copd[,"Gender"]=="2-Female"),"Gender"])


length(pheno_ild[which(pheno_ild[,"Gender"]=="1-Male"),"Gender"])
length(pheno_ild[which(pheno_ild[,"Gender"]=="2-Female"),"Gender"])

length(pheno_ctr[which(pheno_ctr[,"Gender"]=="1-Male"),"Gender"])
length(pheno_ctr[which(pheno_ctr[,"Gender"]=="2-Female"),"Gender"])

ta = rbind(c(22, 16), c(65, 46))
htest = fisher.test(ta) 
htest$p.value 

ta = rbind(c(22, 16), c(61, 52))
htest = fisher.test(ta) 
htest$p.value 

ta = rbind(c(65, 46), c(61, 52))
htest = fisher.test(ta)
htest$p.value

mean(as.numeric(pheno_copd[, "Pre.FEV1.FVC.ratio"]), na.rm=TRUE)
sd(as.numeric(pheno_copd[, "Pre.FEV1.FVC.ratio"]), na.rm=TRUE)

mean(as.numeric(pheno_ild[, "Pre.FEV1.FVC.ratio"]), na.rm=TRUE)
sd(as.numeric(pheno_ild[, "Pre.FEV1.FVC.ratio"]), na.rm=TRUE)

mean(as.numeric(pheno_ctr[, "Pre.FEV1.FVC.ratio"]), na.rm=TRUE)
sd(as.numeric(pheno_ctr[, "Pre.FEV1.FVC.ratio"]), na.rm=TRUE)

t.test(as.numeric(pheno_ctr[, "Pre.FEV1.FVC.ratio"]), as.numeric(pheno_ild[, "Pre.FEV1.FVC.ratio"]))
t.test(as.numeric(pheno_ctr[, "Pre.FEV1.FVC.ratio"]), as.numeric(pheno_copd[, "Pre.FEV1.FVC.ratio"]))
t.test(as.numeric(pheno_copd[, "Pre.FEV1.FVC.ratio"]), as.numeric(pheno_ild[, "Pre.FEV1.FVC.ratio"]))


mean(as.numeric(pheno_copd[, "RCL..Emphysema.Percent"]), na.rm=TRUE)
sd(as.numeric(pheno_copd[, "RCL..Emphysema.Percent"]), na.rm=TRUE)

mean(as.numeric(pheno_ild[, "RCL..Emphysema.Percent"]), na.rm=TRUE)
sd(as.numeric(pheno_ild[, "RCL..Emphysema.Percent"]), na.rm=TRUE)

mean(as.numeric(pheno_ctr[, "RCL..Emphysema.Percent"]), na.rm=TRUE)
sd(as.numeric(pheno_ctr[, "RCL..Emphysema.Percent"]), na.rm=TRUE)

t.test(as.numeric(pheno_ctr[, "RCL..Emphysema.Percent"]), as.numeric(pheno_ild[, "RCL..Emphysema.Percent"]))
t.test(as.numeric(pheno_ctr[, "RCL..Emphysema.Percent"]), as.numeric(pheno_copd[, "RCL..Emphysema.Percent"]))
t.test(as.numeric(pheno_copd[, "RCL..Emphysema.Percent"]), as.numeric(pheno_ild[, "RCL..Emphysema.Percent"]))
