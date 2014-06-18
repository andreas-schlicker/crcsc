# Uses the de novo subtyping based on the simpler NMF protocol. Before applying
# a subtyping step, this script filters the genes to retain only the most 
# significant for each subtype.
# 
# Author: schlandi
###############################################################################

library(NMF)
library(synapseClient)
library(rGithubClient)
library(MASS)

synapseLogin()

# GitHib repository
crcRepo = getRepo("andreas-schlicker/crcsc")
firstScript = getPermlink(crcRepo, "groups/F/denovo/subtyping_simple.r")
thisScript = getPermlink(crcRepo, "groups/F/denovo/subtyping_simple_filtered.r")

# Merged expression dataset
load(getFileLocation(synGet("syn2431652")))
# Sample name conversion
load(getFileLocation(synGet("syn2502279")))
# De novo subtyping results obtained using all genes
load(getFileLocation(synGet("syn2502633")))

# Transpose
m.ALL = t(m.ALL)
# Make the data positive by subtracting the minimum from each gene 
allExprs = m.ALL - apply(m.ALL, 1, min)


### Apply subtyping to all samples
subtypesFiltered = list(all=colnames(m.ALL))
subtypeSigs = list(all=rownames(m.ALL))
subtypeDf = data.frame(subtypeDf, subtypeFiltered="", stringsAsFactors=FALSE)

rownames(subtypeDf) = all.ids[rownames(subtypeDf)]

for (n in names(inmfSubtypes)) {
	if (length(subtypes[[n]]) > 400) {
		# Get the filtered list of features and translate indexes to feature names
		features = lapply(extractFeatures(inmfSubtypes[[n]]$nmfRun[[1]]), function(x) { featureNames(inmfSubtypes[[n]]$nmfRun[[1]])[x] })
		# Remember the features
		subtypeSigs[[n]] = features
		# Subtype using the significant features only
		H = ginv(basis(inmfSubtypes[[n]]$nmfRun[[1]])[unique(unlist(features)), ]) %*% allExprs[unique(unlist(features)), subtypes[[n]]]
		subs = apply(H, 2, which.max)
		subtypeDf[names(subs), "subtypeFiltered"] = paste(n, subs, sep="-")
		if (length(unique(subs)) > 1) {
			for (k in unique(subs)) {
				subtypesFiltered[[paste(n, k, sep="-")]] = names(subs)[subs==k]
			}
		}
	}
}

rownames(subtypeDf) = all.ids.rev[rownames(subtypeDf)]

#> table(subtypeDf[, "subtype"] == subtypeDf[, "subtypeFiltered"])                                                                                                                                                                                                     [156/648]
#
#FALSE  TRUE
#  260  4235

#> 260 / 4495
#[1] 0.05784205

#> unlist(lapply(subtypeSigs, function(x) { sum(unlist(lapply(x, function(y) { length(y) })))}))
#        all       all-1       all-2     all-1-1     all-2-1     all-2-2   all-1-1-1   all-2-1-1   all-2-1-2   all-2-2-1   all-2-2-2 all-1-1-1-1 all-1-1-1-2 all-2-1-1-1 all-2-1-2-2 all-2-2-2-2
#        354         615         893         862         413         812         891         892         949         720         646         635         764         721         824         797

### Save the subtyping result in Synapse
# Write temporary file
filePath = file.path(tempdir(), "iNMF_denovo_subtyping_simple_filtered.tsv")
write.table(subtypeDf, file=filePath, sep="\t", quote=FALSE)
# List with used resources
resources = list(list(entity="syn2431652", wasExecuted=F),
				 list(entity="syn2502279", wasExecuted=F),
				 list(entity="syn2502633", wasExecuted=F),
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
filePath = file.path(tempdir(), "iNMF_denovo_subtyping_simple_filtered.rdata")
save(inmfSubtypes, subtypes, subtypesFiltered, subtypeSigs, subtypeDf, file=filePath)
# And upload 
synFile = File(path=filePath, parentId="syn2502277")
failed = TRUE
tries = 0
while (failed && (tries < 5)) {
	res = tryCatch(synStore(synFile, used=resources),
			error=function(e) NA)
	if (class(res) != "File") {
		failed=FALSE
	}
	tries = tries + 1
}
unlink(filePath)

