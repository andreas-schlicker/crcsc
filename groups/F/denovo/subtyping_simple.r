# Performs de novo subtyping using a simpler NMF protocol. 
# 
# Author: Andreas Schlicker (a.schlicker@nki.nl)
###############################################################################

library(NMF)
library(fastcluster)
library(gplots)
library(stringr)
library(synapseClient)
library(rGithubClient)
library(MASS)

synapseLogin()

# GitHib repository
crcRepo = getRepo("andreas-schlicker/crcsc")
thisScript = getPermlink(crcRepo, "groups/F/denovo/subtyping_simple.r")

# Merged expression dataset
load(getFileLocation(synGet("syn2431652")))
# Sample name conversion
load(getFileLocation(synGet("syn2502279")))

# This data is not yet in Synapse
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/CRC_SUBTYPES/cms_classification.rdata")

# Transpose
m.ALL = t(m.ALL)
# Make the data positive by subtracting the minimum from each gene 
allExprs = m.ALL - apply(m.ALL, 1, min)

set.seed(123)

# Split samples into discovery and validation set
discoverySet = c()
validationSet = c()
for (n in c("agendia_gse42284", "agendia_ico208", "agendia_vhb70", "amc_ajccii", "french", "petacc3")) {
	samples = grep(n, colnames(allExprs), value=TRUE)
	discoverySet = c(discoverySet, sample(samples, floor(length(samples)/2)))
	validationSet = c(validationSet, setdiff(samples, discoverySet))
}
exprsDiscov = allExprs[, discoverySet]
exprsValid = allExprs[, validationSet]

# Which dataset does a sample belong to?
dataset = sapply(discoverySet, function(x) { str_split(x, "-")[[1]][1] })

### iNMF iterations
# List with all current subtypes
inmfSteps = list(list(samples=discoverySet, history="all"))
# List with all discovered subtypes
inmfSubtypes = list()
counter = 1

while (length(inmfSteps) > 0) {
	# Get the next subtype to split
	current = inmfSteps[[1]]
	inmfSteps[[1]] = NULL
	print(paste("Subtyping ", current$history, " ...", sep="")) 
	
	# Estimate the rank
	print("Estimating rank ...")
	nmfEstRanks = runNMF(list(rownames(exprsDiscov)), exprsDiscov[, current$samples], ranks=2:10, nrun=30, threads=40)
	
	# Get the most common rank(s)
	rankEstimate = getRankStats(nmfEstRanks)
	selectedRank = selectRank(rankEstimate, "first")
	print(paste("Chose rank: ", selectedRank, sep=""))
	
	# Run NMF with this rank
	nmfRun = runNMF(list(rownames(exprsDiscov)), exprsDiscov[, current$samples], ranks=selectedRank, nrun=200, threads=40)
	clustering = apply(coef(nmfRun[[1]]), 2, which.max)
	signatures = extractFeatures(nmfRun[[1]])
	print("Finished clustering ...")
	
	pvalue = testStratification(clustering, dataset)
	subtypes = split(names(clustering), clustering)
	
	inmfSubtypes[[current$history]] = list(nmfEstRanks=nmfEstRanks, nmfRun=nmfRun,
			signatures=signatures, subtypes=subtypes, 
			pvalue=pvalue, history=current$history)
	
	#if (pvalue > 0.05) {
	for (n in names(subtypes)) {
		if (length(subtypes[[n]]) > 50) {
			inmfSteps = c(inmfSteps, list(list(samples=subtypes[[n]], history=paste(current$history, n, sep="-"))))
		}
	}
	#}
	
	save.image(paste("inmf_", counter, ".rdata", sep=""))
	counter = counter + 1
}


### Apply subtyping to all samples
subtypes = list(all=colnames(m.ALL))
subtypeDf = data.frame(sample=colnames(m.ALL), subtype="all", stringsAsFactors=FALSE)
rownames(subtypeDf) = colnames(m.ALL)

for (n in names(inmfSubtypes)) {
	if (length(subtypes[[n]]) > 400) {
		H = ginv(basis(inmfSubtypes[[n]]$nmfRun[[1]])) %*% allExprs[, subtypes[[n]]]
		subs = apply(H, 2, which.max)
		subtypeDf[names(subs), "subtype"] = paste(n, subs, sep="-")
		if (length(unique(subs)) > 1) {
			for (k in unique(subs)) {
				subtypes[[paste(n, k, sep="-")]] = names(subs)[subs==k]
			}
		}
	}
}

rownames(subtypeDf) = all.ids.rev[rownames(subtypeDf)]

subtypeDf = cbind(subtypeDf, 
				  cms4=cms4[rownames(subtypeDf), "CMS"], 
				  cms5=cms5[rownames(subtypeDf), "CMS"],
				  dataset=sapply(rownames(subtypeDf), function(x) { str_split(x, "\\.|-")[[1]][1] } ))

#table(subtypeDf[, c("subtype", "cms4")])
#table(subtypeDf[, c("subtype", "cms5")])
#table(subtypeDf[, c("subtype", "dataset")])

### Save the subtyping result in Synapse
# Write temporary file
filePath = file.path(tempdir(), "iNMF_denovo_subtyping_simple.tsv")
write.table(subtypeDf, file=filePath, sep="\t", quote=FALSE)
# List with used resources
resources = list(list(entity="syn2431652", wasExecuted=F),
				 list(entity="syn2502279", wasExecuted=F),
				 list(url=thisScript, name=basename(thisScript), wasExecuted=T))

# And upload 
synFile = File(path=filePath, parentId="syn2502277")
failed = TRUE
tries = 0
while (failed && (tries < 5)) {
	res = tryCatch(synStore(synFile, used=resources),
			error=function(e) NA)
	if (!is.na(res)) {
		failed=FALSE
	}
	tries = tries + 1
}
unlink(filePath)


### Save the object
filePath = file.path(tempdir(), "iNMF_denovo_subtyping_simple.rdata")
save(inmfSubtypes, subtypes, subtypeDf, file=filePath)
# And upload 
synFile = File(path=filePath, parentId="syn2502277")
failed = TRUE
tries = 0
while (failed && (tries < 5)) {
	res = tryCatch(synStore(synFile, used=resources),
			error=function(e) NA)
	if (!is.na(res)) {
		failed=FALSE
	}
	tries = tries + 1
}
unlink(filePath)
