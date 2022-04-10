require(DESeq2)
require(tidyverse)
require(IHW)

##################################################################################
#1. load data
##################################################################################

# read meta data
metafile <- ("Kim_5482_new_metadata.tsv")
metaData.new <- read.table(metafile, header = T, sep = "\t")

metaData.new <- metaData.new %>%
    mutate(Treatment = factor(Treatment, levels = c("untreated", "AC"))) %>% 
    column_to_rownames("ID")

# read featurecount count table
countfile <- ("Kim_5482_new_featurecount_symbol.tsv")
counttable.new <- read.table(countfile, header = T, sep = "\t")

counttable.new <- counttable.new %>%
    column_to_rownames("symbol")

##################################################################################
#2. load data
##################################################################################
ddsHTseq.new <- DESeqDataSetFromMatrix( countData = counttable.new,
                                        colData = metaData.new,
                                        design = ~ 1)

# Reset the levels of colData factor variables
colData(ddsHTseq.new) <- DataFrame(droplevels(data.frame(colData(ddsHTseq.new))))

# dplyr::filter out genes
ddsHTseq.new <- ddsHTseq.new[rowSums(DESeq2::counts(ddsHTseq.new) > 2) > 2, ]

# Set the DEG design
design(ddsHTseq.new) <- ~ MouseID + Treatment

# Calculate size factor
ddsHTseq.new <- estimateSizeFactors(ddsHTseq.new)

# Run the DEG analysis
dds.new <- DESeq(ddsHTseq.new, parallel = F, betaPrior=FALSE)
rld.new <- rlogTransformation(dds.new, blind=F)
vsd.new <- varianceStabilizingTransformation(dds.new, blind=F)


resultsNames(dds.new)
#AC_vs_WT
Res.ihw.AC_vs_WT.new.all <- results(dds.new,name="Treatment_AC_vs_untreated",filterFun = ihw)
Res.ihw.AC_vs_WT.new.all <- Res.ihw.AC_vs_WT.new.all[order(Res.ihw.AC_vs_WT.new.all$padj),]
write.table(Res.ihw.AC_vs_WT.new.all, file ="Kim_5482_AC_vs_WT.new_all.tsv", sep = "\t", quote = F, col.names=NA)

Res.ihw.AC_vs_WT.new <- Res.ihw.AC_vs_WT.new.all[which(Res.ihw.AC_vs_WT.new.all$padj < 0.05),]
Res.ihw.AC_vs_WT.new <- Res.ihw.AC_vs_WT.new[order(Res.ihw.AC_vs_WT.new$padj),]
write.table(Res.ihw.AC_vs_WT.new, file ="Kim_5482_AC_vs_WT.new.tsv", sep = "\t", quote = F, col.names=NA)

Res.ihw.AC_vs_WT.new.rnk <- data.frame(Res.ihw.AC_vs_WT.new.all)
Res.ihw.AC_vs_WT.new.rnk <- Res.ihw.AC_vs_WT.new.rnk[,"stat", drop = F]
Res.ihw.AC_vs_WT.new.rnk <- Res.ihw.AC_vs_WT.new.rnk[order(-Res.ihw.AC_vs_WT.new.rnk$stat),,drop = F]
write.table(Res.ihw.AC_vs_WT.new.rnk,file="Kim_5482_AC_vs_WT.new.rnk",quote=F,sep="\t",col.names="#gene\twald_stat")
