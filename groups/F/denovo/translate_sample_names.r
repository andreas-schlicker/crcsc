# The combined dataset and the Phase 1 consensus clustering use different
# sample names. These scheme have to harmonized. 
# 
# Author: Andreas Schlicker (a.schlicker@nki.nl)
###############################################################################

library(synapseClient)
library(stringr)
library(GEOquery)

synapseLogin()

load(getFileLocation(synGet("syn2431652")))
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/CRC_SUBTYPES/cms_classification.rdata")

# We need the GEO entry for the AMC set
gse33113 = getGEO("GSE33113")

# That's the sample names in the combined dataset
all.ids = rownames(m.ALL)
names(all.ids) = rownames(m.ALL)

names(all.ids) = gsub("agendia_gse42284-", "agendia_gse42284\\.", names(all.ids))
names(all.ids) = gsub("agendia_ico208-", "agendia_ico208\\.", names(all.ids))
names(all.ids) = gsub("agendia_vhb70-", "agendia_vhb70\\.", names(all.ids))
names(all.ids)[grep("amc_ajccii", names(all.ids))] = sapply(grep("amc_ajccii", names(all.ids), value=TRUE), function(x) { paste(str_split(x, "_")[[1]][1:2], collapse="_") })
names(all.ids)[grep("amc_ajccii", names(all.ids))] = gsub("-", "\\.", sapply(grep("amc_ajccii", names(all.ids), value=TRUE), function(x) { y=str_split(x, "-")[[1]]; y=paste(y[1], as.character(pData(gse33113[[1]])[gsub(".CEL", "", y[2]), 1]), sep="."); y }))
names(all.ids)[grep("french", names(all.ids))] = gsub("-", "\\.", sapply(grep("french", names(all.ids), value=TRUE), function(x) { str_split(x, "_")[[1]][1] }))
names(all.ids) = gsub("gse13067-", "gse13067\\.", names(all.ids))
names(all.ids) = gsub("gse13294-", "gse13294\\.", names(all.ids))
names(all.ids) = gsub("gse14333-", "gse14333\\.", names(all.ids))
names(all.ids) = gsub("gse17536-", "gse17536\\.", names(all.ids))
names(all.ids)[grep("gse20916", names(all.ids))] = gsub("-", "\\.", sapply(grep("gse20916", names(all.ids), value=TRUE), function(x) { str_split(x, "_")[[1]][1] }))
names(all.ids) = gsub("gse2109-", "gse2109\\.", names(all.ids))
names(all.ids) = gsub("gse23878-", "gse23878\\.", names(all.ids))
names(all.ids)[grep("gse37892", names(all.ids))] = gsub("-", "\\.", sapply(grep("gse37892", names(all.ids), value=TRUE), function(x) { str_split(x, "_")[[1]][1] }))
names(all.ids) = gsub("kfsyscc-", "kfsyscc\\.", names(all.ids))
names(all.ids) = gsub("mdanderson-", "mdanderson\\.", names(all.ids))
names(all.ids)[grep("nki_az", names(all.ids))] = gsub("-", "\\.", sapply(grep("nki_az", names(all.ids), value=TRUE), function(x) { paste(str_split(x, "_")[[1]][1:2], collapse="_") }))
names(all.ids) = gsub("petacc3-PETACC3", "petacc\\.PETACC3", names(all.ids))
names(all.ids) = gsub("tcgacrc_merged", "tcga_rnaseqAll", names(all.ids))
names(all.ids) = gsub(".CEL", "", names(all.ids))

all.ids.rev = names(all.ids)
names(all.ids.rev) = all.ids

save(all.ids, all.ids.rev, file="/srv/nfs4/medoid-bulk/NKI/a.schlicker/CRC_SUBTYPES/cms_classification_id_translation.rdata")


