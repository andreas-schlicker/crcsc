# This file contains functions for performing a de-novo subtype discovery using
# iNMF. 
# 
# Author: Andreas Schlicker (a.schlicker@nki.nl)
###############################################################################

##' Runs NMF using the given parameters.
##' @param geneLists list with feature vectors to run NMF on
##' @param data data matrix for running NMF
##' @param ranks vector with ranks to test; default: 2:10
##' @param method which NMF method to use; default: nfNMF
##' @param nrun number of NMF runs with random initializations; default: 30
##' @param threads number of concurrent threads; default: 40
##' @return list with the NMF objects
##' @author Andreas Schlicker
runNMF = function(geneLists, data, ranks=2:10, method="nsNMF", nrun=30, threads=40) {
	lapply(1:length(geneLists), 
		   function(x) { nmf(data[geneLists[[x]], ], rank=ranks, method=method, nrun=30, .pbackend=threads) })
}

##' Uses the cophenetic clustering coefficient for suggesting 
##' a rank for NMF. The function locates the first two maxima
##' of the coefficient for each NMF object in the input.
##' @param estranks list with all NMF objects containing the quality measures
##' @return list with the statistics for the first and second maximum
##' @author Andreas Schlicker
getRankStats = function(estranks) {
	n = nrow(nmfEstRanks[[1]]$measures)
	
	# Look for the first and the second maximum
	firstMax = rep(0, times=n)
	names(firstMax) = as.character(nmfEstRanks[[1]]$measures[, "rank"])
	secondMax = rep(0, times=n)
	names(secondMax) = as.character(nmfEstRanks[[1]]$measures[, "rank"])
	
	for (i in 1:length(nmfEstRanks)) {
		# Compute the sign of the difference at each point
		diffSigns = sign(c(0, diff(nmfEstRanks[[1]]$measures[, "cophenetic"])))
		# Check where the values change direction
		diffSignSums = sapply(1:length(diffSigns), function(x) { diffSigns[x] - diffSigns[x+1] })
		
		# k at the first and second maximum
		fMax = which(diffSignSums > 0)[1]
		sMax = which(diffSignSums > 0)[2]
		
		firstMax[fMax] = firstMax[fMax] + 1
		secondMax[sMax] = secondMax[sMax] + 1
	}
	
	list(first=firstMax, second=secondMax)
}

##' Get the most common rank.
##' @param rankStats rank statistics as returned by getRankStats
##' @param which which maximum to select, either "first" or "second"
##' @return the rank as integer
##' @author Andreas Schlicker
selectRank = function(rankStats, which=c("first", "second")) {
	which = match.arg(which)
	
	as.integer(names(rankStats[[which]])[which.max(rankStats[[which]])])
}

##' Compute coclustering across all NMF runs. Will save a heatmap of the 
##' coclustering matrix into the given file if the heatmap parameter is defined.
##' @param nmfs list with all NMF solutions
##' @param heatmap absolute filename for saving the heatmap; default: NULL (no heatmap)
##' @return matrix where each element indicates how often the samples in the row and column
##' were coclustered
##' @author Andreas Schlicker
getCoClustering = function(nmfs, heatmap=NULL) {
	coclust = matrix(0, nrow=ncol(coef(nmfs[[1]])), ncol=ncol(coef(nmfs[[1]])))
	colnames(coclust) = colnames(coef(nmfs[[1]]))
	rownames(coclust) = colnames(coef(nmfs[[1]]))
	
	for (i in 1:length(nmfs)) {
		# Assign each sample to the cluster with the highest coefficient
		tempClust = apply(coef(nmfs[[i]]), 2, which.max)
		# The different clusters
		clusters = split(names(tempClust), tempClust)
		for (j in 1:length(tempClust)) {
			coclust[names(tempClust)[j], clusters[[tempClust[j]]]] = coclust[names(tempClust)[j], clusters[[tempClust[j]]]] + 1
		}
	}
	
	if (!is.null(heatmap)) {
		if (!require(gplots)) {
			warning("Can't plot heatmap without package \"gplots\"!")
		} else {
			png(heatmap, width=4000, height=3000, res=300)
			aheatmap(coclust, scale="none", 
					 color=colorpanel(49, low="white", high="blue4"),
				 	 breaks=seq(min(coclust), max(coclust), length.out=50),
				 	 annCol=dataset[colnames(coclust)], annRow=dataset[rownames(coclust)],
				 	 distfun="correlation")
			dev.off()
		}
	}
	
	coclust
}

##' Assign each sample to one of the clusters
##' @param coclust coclustering matrix
##' @param rank number of clusters
##' @return vector with cluster assignment for each sample
##' @author Andreas Schlicker
getCoClusters = function(coclust, rank) {
	distMat = as.dist(1-cor(coclust))
	subtypes = cutree(hclust(distMat, method="complete"), k=rank)
	sil = silhouette(subtypes, distMat)
	rownames(sil) = colnames(coclust)
	
	list(subtypes=subtypes, silhouette=sil)
}

##' Determines core clusters using the silhouette width. Only samples 
##' exceeding the silhouette width threshold make it to a core cluster.
##' @param sil silhouette object
##' @param cut-off the silhouette width cut-off; default: 0.7
##' The default cut-off come from this site:
##' http://www.unesco.org/webworld/idams/advguide/Chapt7_1_1.htm
##' @return list with the core clusters
##' @author Andreas Schlicker
getCoreClusters = function(sil, cutoff=0.7) {
	samples = sil[, "sil_width"][which(sil[, "sil_width"] > cutoff)]
	clusters = sil[names(samples), "cluster"]
	
	split(names(clusters), clusters)
}

##' Determines subtype signatures based on differential expression between 
##' the clusters.
##' @param clusters list with sample cluster
##' @param data data matrix
##' @param noTestGenes number of genes with highest IQR to be tested for 
##' differential expression
##' @param FDRcutoff FDR threshold for calling a gene differentially expressed
##' @return list with the subtype signatures
##' @author Andreas Schlicker
getSubtypeSignatures = function(clusters, data, noTestGenes=5000, FDRcutoff=0.05) {
	vars = names(sort(apply(data[, unlist(clusters)], 1, IQR), decreasing=TRUE))[1:noTestGenes]
	
	signatures = data.frame()
	for (n in names(clusters)) {
		otherSamples = setdiff(colnames(data), clusters[[n]])
		
		signatures = rbind(signatures,
			data.frame(subtype=n,
					   genes=vars,
					   pvalue=sapply(vars,
							function(y) { wilcox.test(data[y, clusters[[n]]], 
													  data[y, otherSamples],
													  alternative="greater")$p.value })))
	}
	
	signatures = cbind(signatures, FDR=p.adjust(signatures[, "pvalue"], method="BH"))
	signatures = subset(signatures, FDR < FDRcutoff)
	
	list(signatures=split(signatures$genes, signatures$subtype), statistics=signatures)
}

##' Computes the p-value for stratification of datasets into subtypes.
##' @param subtyping vector with subtype assignment for each sample
##' @param dataset vector with dataset assignemnt for each sample
##' @return Fisher exact test p-value
##' @author Andreas Schlicker
testStratification = function(subtyping, dataset) {
	temp = data.frame(subtypes=subtyping, dataset[names(subtyping)])
	
	fisher.test(table(temp), simulate.p.value=TRUE)$p.value
}
