require(DESeq2)
require(tidyverse)
require(IHW)

###################################################################################################
#1. load data
###################################################################################################
# read in meta data
metafile <- ("Kim_4970_core_metadata.tsv")
metaData.core <- read.table(metafile, header = T, sep = "\t")

metaData.core <- metaData.core %>%
    mutate(Tissue = factor(Tissue, levels = c("blood", "tumor"))) %>% 
    column_to_rownames("ID")

# read in htseq count table
countfile <- ("Kim_4970_core_featurecount_symbol.tsv")
counttable.core <- read.table(countfile, header = T, sep = "\t")

counttable.core <- counttable.core %>%
    column_to_rownames("symbol")

###################################################################################################
#2. DE analysis
###################################################################################################

ddsHTseq.core <- DESeqDataSetFromMatrix( countData = counttable.core,
                                         colData = metaData.core,
                                         design = ~ 1)

# Reset the levels of colData factor variables
colData(ddsHTseq.core) <- DataFrame(droplevels(data.frame(colData(ddsHTseq.core))))

# filter out low expression genes
ddsHTseq.core <- ddsHTseq.core[rowSums(DESeq2::counts(ddsHTseq.core) > 2) > 2, ]

# Set the DEG design
design(ddsHTseq.core) <- ~ PatientID + Tissue

# Calculate size factor
ddsHTseq.core <- estimateSizeFactors(ddsHTseq.core)

# Run the DEG analysis
dds.core <- DESeq(ddsHTseq.core, parallel = F, betaPrior=FALSE)
rld.core <- rlogTransformation(dds.core, blind=F)
vsd.core <- varianceStabilizingTransformation(dds.core, blind=F)

resultsNames(dds.core)
#tumor vs blood
Res.ihw.tumor_vs_blood.all <- results(dds.core,name="Tissue_tumor_vs_blood",filterFun = ihw)
Res.ihw.tumor_vs_blood.all <- Res.ihw.tumor_vs_blood.all[order(Res.ihw.tumor_vs_blood.all$padj),]
write.table(Res.ihw.tumor_vs_blood.all, file ="Kim_4970_tumor_vs_blood.all.tsv", sep = "\t", quote = F, col.names=NA)

Res.ihw.tumor_vs_blood <- Res.ihw.tumor_vs_blood.all[which(Res.ihw.tumor_vs_blood.all$padj < 0.05),]
Res.ihw.tumor_vs_blood <- Res.ihw.tumor_vs_blood[order(Res.ihw.tumor_vs_blood$padj),]
write.table(Res.ihw.tumor_vs_blood, file ="Kim_4970_tumor_vs_blood.tsv", sep = "\t", quote = F, col.names=NA)

Res.ihw.tumor_vs_blood.rnk <- data.frame(results(dds.core,name="Tissue_tumor_vs_blood",filterFun = ihw))
Res.ihw.tumor_vs_blood.rnk <- Res.ihw.tumor_vs_blood.rnk[,"stat", drop = F]
Res.ihw.tumor_vs_blood.rnk <- Res.ihw.tumor_vs_blood.rnk[order(-Res.ihw.tumor_vs_blood.rnk$stat),,drop = F]
write.table(Res.ihw.tumor_vs_blood.rnk,file="Kim_4970_tumor_vs_blood.rnk",quote=F,sep="\t",col.names="#gene\twald_stat")

