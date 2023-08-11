
cc_with_plots = function(data,innerLinkage="average",finalLinkage="average",maxK=15,reps=1000,pItem=0.8,pFeature=0.8,output.prefix="cc_results",clusterAlg="hc",distance="pearson",seed=12345,plot="pngBMP",writeTable=TRUE,verbose=TRUE) {
  require(ConsensusClusterPlus)
  require(flashClust)
  require(heatmap3)
  require(cluster)

  output.dir = paste(output.prefix, "_results", sep="")
  bwr = colorRampPalette(c("blue", "white", "red"))(100)
  
  cc.results = ConsensusClusterPlus(data,innerLinkage=innerLinkage,finalLinkage=finalLinkage,maxK=maxK,reps=reps,pItem=pItem,pFeature=pFeature,title=output.dir,clusterAlg=clusterAlg,distance=distance,seed=seed,plot=plot,writeTable=writeTable,verbose=verbose)
  icl = calcICL(cc.results,title=output.dir,plot="pdf")
  saveRDS(cc.results, paste(output.dir, "/CC_results.rds", sep=""))
  saveRDS(icl, paste(output.dir, "/ICL_results.rds", sep=""))
  
  
  ## Plot heatmaps
  for (i in 2:maxK) {
    cat("Generating heatmap for k = ", i, "\n")
    filename = paste(output.dir, "/heatmap_k", i, ".jpg", sep="")
    
    label = cc.results[[i]][["consensusClass"]] 
    
    jpeg(filename)
    heatmap3(data, col=bwr, ColSideColors=label, col.clustering="semisupervised", labRow=FALSE, labCol=FALSE, hclustfun=function(c) hclust(c, method=finalLinkage))
    dev.off()
  }
  
  ## Plot silhouette width
  sample.sil = c()
  cluster.sil.mean = c()
  filename = paste(output.dir, "/silhouette_width.pdf", sep="")
  
  pdf(filename, useDingbats=FALSE)
  for (i in 2:maxK) {
	label = cc.results[[i]][["consensusClass"]] 
    sil = silhouette(label, as.dist(1-cc.results[[i]][[1]]))
    plot(sil, main=i)
    cluster.sil.mean = c(cluster.sil.mean, mean(sil[,3]))
    sample.sil = cbind(sample.sil, sil[,3])
  }
  plot(2:maxK, cluster.sil.mean, pch=19, lwd=2, xlab="k", ylab="Mean Silhouette Width")
  write.table(sample.sil, paste(output.dir, "/silhouette_width.txt", sep=""), sep="\t", quote=FALSE)
  dev.off()
  
  
  ## Plot Cophenetic cluster
  cop = c()
  for(i in 2:maxK) {
    d = dist(cc.results[[i]][[1]])
    cop = c(cop, cor(d, cophenetic(hclust(d, method=finalLinkage))))
  }

  filename = paste(output.dir, "/cophenetic.pdf", sep="")
  pdf(filename, useDingbats=FALSE)
  plot(2:maxK, cop, type="b", ylab="Cophenetic Correlation", xlab="Cluster", main="Cophenetic Comparison", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  dev.off()

  ## Plot Average Item Consensus
  avg.item = aggregate(icl[[1]][,3], by=list(icl[[1]][,1]), mean, na.rm=T)

  filename = paste(output.dir, "/average_item_consensus.pdf", sep="")
  pdf(filename, useDingbats=FALSE)
  plot(avg.item[,1], avg.item[,2], type="b", ylab="Average Item Consensus", xlab="Cluster", main="Item Consensus Comparison", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0.4, 1))
  dev.off()
}





