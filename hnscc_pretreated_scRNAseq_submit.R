require(tidyverse)
require(Seurat)
require(SeuratWrappers)
require(cowplot)
require(RColorBrewer)
require(viridis)
require(CellID)
require(DropletUtils)
require(patchwork)
require(celldex)
require(future)
require(Cairo)
options(future.globals.maxSize = 80000 * 1024^2)

wrkdir <- "./"
dtadir <- "/scratch/TBI/Projects/PI_Kim/Kim_3219/data/Count"

setwd(wrkdir)

#################################################################################################################
## 00. read in filtered count data, meta data, prepare for annotation data, combine multiple samples
#################################################################################################################
meta.tbl <- read.delim(paste(wrkdir,"Kim_3219_2853_meta.tsv", sep = "/"), sep="\t", header=TRUE)
meta.tbl

sample_paths <- c()
sample_paths["pat01.t"] <- paste0(dtadir, "/3219-SR-1/filtered_feature_bc_matrix")
sample_paths["pat02.t"] <- paste0(dtadir, "/3219-SR-3/filtered_feature_bc_matrix")
sample_paths["pat03.t"] <- paste0(dtadir, "/3219-SR-5/filtered_feature_bc_matrix")
sample_paths["pat04.t"] <- paste0(dtadir, "/3521-SR-1/filtered_feature_bc_matrix")
sample_paths["pat05.t"] <- paste0(dtadir, "/Kim_2853/Count_1/filtered_feature_bc_matrix")
sample_paths["pat05.b"] <- paste0(dtadir, "/Kim_2853/Count_2/filtered_feature_bc_matrix")
sample_paths["pat06.t"] <- paste0(dtadir, "/Kim_2853/Count_3/filtered_feature_bc_matrix")
sample_paths["pat06.b"] <- paste0(dtadir, "/Kim_2853/Count_4/filtered_feature_bc_matrix")

## prepare the seurat object
all.data <- Read10X(data.dir = sample_paths, gene.column = 2, unique.features = TRUE)

all.seurat <- CreateSeuratObject(counts = all.data, project = names(sample_paths), min.cells = 3, min.features = 200)

## Standard pre-processing workflow
all.seurat$"percent.mt" <- PercentageFeatureSet(all.seurat, pattern = "^MT-")

# set the filtering cutoffs based on plots
all.list <- SplitObject(all.seurat, split.by = "orig.ident")
for (i in 1:length(all.list)) {
  meta <- unlist(strsplit(names(sample_paths)[i], "[.]"))
  all.list[[i]]$patient <- meta[1]
  all.list[[i]]$tissue <- meta[2]
  all.list[[i]]$treatment <- "pre-treated"
  all.list[[i]] <- subset(all.list[[i]], subset = percent.mt < 10 & nCount_RNA > 1500 & nCount_RNA < 15000)
  all.list[[i]] <- SCTransform(all.list[[i]], 
                               vars.to.regress = c("nCount_RNA","percent.mt"),
                               method = 'glmGamPoi',
                               return.only.var.genes = F,
                               verbose = T,
                               variable.features.n = 5000, 
                               conserve.memory = T)
}
all.features <- SelectIntegrationFeatures(object.list = all.list, nfeatures = 5000)
all.list <- PrepSCTIntegration(object.list = all.list, anchor.features = all.features, verbose = T)

all.anchors <- FindIntegrationAnchors(object.list = all.list, normalization.method = "SCT", 
                                      anchor.features = all.features, verbose = T, dims = 1:50)
hnscc.integrated <- IntegrateData(anchorset = all.anchors, normalization.method = "SCT", 
                                  verbose = T, dims = 1:50)

rm(all.anchors,all.data,all.list,all.seurat,all.features)


#################################################################################################################
## 01. extra QC: remove doublets and ambient RNA contamination
#################################################################################################################

pat <- "_doublet.txt$"

doublet.all <- list.files(path = wrkdir, pattern = pat, all.files = TRUE, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE)

doublet.DT <- list()

for (i in 1:length(doublet.all) ) {
  infile = paste(wrkdir, doublet.all[i], sep = "/")
  df <- read.csv(infile, header = T, stringsAsFactors = FALSE)
  samples <- gsub("(.*)\\_doublet.txt", "\\1", doublet.all[i])
  doublet.DT[[samples]] <- df %>%
    dplyr::mutate(barcode = paste(samples, barcode, sep = "_"))
}

doublet.anno <- bind_rows(doublet.DT)
doublet.anno <- doublet.anno %>%
  dplyr::filter(barcode %in% Cells(hnscc.integrated)) %>%
  dplyr::arrange(Cells(hnscc.integrated))

hnscc.integrated$doublet_scores <- doublet.anno$doublet_scores
hnscc.integrated$predicted_doublets <- doublet.anno$predicted_doublets

p.d <- DimPlot(hnscc.integrated, reduction = "tsne",label = F,group.by="predicted_doublets", cols = c("grey", "black"), label.size = 6, label.box = F)
pdf(file = "DimPlot_extra_doublet_filter.pdf", width = 4, height = 4, pointsize = 2, onefile = T)
p.d
dev.off()
hnscc.integrated <- subset(hnscc.integrated, subset = predicted_doublets == "False")

#################################################################################################################
## 02. fine tuning parameters
#################################################################################################################
#cluster integrated data
hnscc.integrated <- RunPCA(hnscc.integrated, features = VariableFeatures(object = hnscc.integrated),npcs = 50)
hnscc.integrated <- RunTSNE(hnscc.integrated, dims = 1:50, verbose = FALSE,tsne.method = "Rtsne")
hnscc.integrated <- RunUMAP(hnscc.integrated, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
hnscc.integrated <- FindNeighbors(hnscc.integrated, dims = 1:50, verbose = FALSE)
hnscc.integrated <- FindClusters(hnscc.integrated, resolution=seq(0.01,1.5,by=0.01), verbose = FALSE,algorithm=1)

best=0.14
hnscc.integrated$seurat_clusters <- hnscc.integrated@meta.data[,which(colnames(hnscc.integrated@meta.data)==paste("integrated_snn_res.",best,sep=""))]
Idents(hnscc.integrated) <- "seurat_clusters"
saveRDS(hnscc.integrated,"hnscc.integrated.rds")

# normalize and scale data in slot "RNA" for downstream visualization purpose
#plan("multiprocess", workers = 1)
hnscc.integrated <- NormalizeData(hnscc.integrated, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(hnscc.integrated)
hnscc.integrated <- ScaleData(hnscc.integrated, features = all.genes)

# markers for each cluster
hnscc.integrated.markers <- FindAllMarkers(hnscc.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST")
write.table(hnscc.integrated.markers, file = "hnscc.cluster.markers.tsv", sep = "\t", quote = F, row.names = F)

#################################################################################################################
## 03. cell type annotation for major cell types
#################################################################################################################
# prepare reference datasets

# GSE127465 (raw data for seurat reference based mapping)
# prepare & FIX meta data and make barcodes.tsv file (error in matrix file, author switched cell and gene columns!!!)
#GSE127465.path <- paste(dtadir, "..", "..", "GSE127465", sep = "/")
#GSE127465.meta <- read.table(paste(GSE127465.path, "GSE127465_human_cell_metadata_54773x25.tsv", sep = "/"), sep = "\t", header = T)
#write(x = paste(GSE127465.meta$Library, GSE127465.meta$Barcode, sep = "_"), file = paste(GSE127465.path, "human", "barcodes.tsv", sep = "/"))
#GSE127465.gene <- read.table(paste(GSE127465.path, "GSE127465_gene_names_human_41861.tsv", sep = "/"), sep = "\t", header = F)
#GSE127465.gene$V2 <- GSE127465.gene$V1
#write.table(GSE127465.gene, "genes.tsv", sep = "\t", col.names = F, row.names = F, quote = F)
#GSE127465.mtx <- read.table(paste(GSE127465.path, "GSE127465_human_counts_normalized_54773x41861.mtx", sep = "/"), sep = " ", header = T, comment.char = "%", check.names = F)
#GSE127465.mtx <- GSE127465.mtx %>%
#  dplyr::arrange(`54773`) %>%
#  dplyr::select(`41861`,`54773`,`44663765`)
#write.table(GSE127465.mtx, paste(GSE127465.path, "human", "matrix.mtx", sep = "/"), sep = " ", quote = F, row.names = F)

GSE127465.data <- Read10X(data.dir = paste(GSE127465.path, "human", sep = "/"), gene.column = 1, unique.features = T)
rownames(GSE127465.meta) <- colnames(GSE127465.data)
GSE127465.seurat <- CreateSeuratObject(counts = GSE127465.data, project = "GSE127465", min.cells = 1, min.features = 1)
GSE127465.seurat <- AddMetaData(GSE127465.seurat, metadata = GSE127465.meta)
GSE127465.list <- SplitObject(GSE127465.seurat, split.by = "Barcoding.emulsion")
for (i in 1:length(GSE127465.list)) {
  GSE127465.list[[i]] <- NormalizeData(GSE127465.list[[i]], verbose = FALSE)
  GSE127465.list[[i]] <- FindVariableFeatures(GSE127465.list[[i]], selection.method = "vst", nfeatures = 5000, 
                                              verbose = FALSE)
}

GSE127465.anchors <- FindIntegrationAnchors(object.list = GSE127465.list, anchor.features = 5000, dims = 1:50)
GSE127465.integrated <- IntegrateData(anchorset = GSE127465.anchors, dims = 1:50)
DefaultAssay(GSE127465.integrated) <- "integrated"
GSE127465.integrated <- ScaleData(GSE127465.integrated, verbose = FALSE)
GSE127465.integrated <- RunPCA(GSE127465.integrated, features = VariableFeatures(object = GSE127465.integrated), npcs = 50, verbose = FALSE)
GSE127465.integrated <- RunTSNE(GSE127465.integrated, dims = 1:50, verbose = FALSE,tsne.method = "Rtsne")
GSE127465.integrated <- RunUMAP(GSE127465.integrated, reduction = "pca", dims = 1:50, verbose = FALSE)

saveRDS(GSE127465.integrated,"GSE127465.integrated.rds")
rm(GSE127465.anchors,GSE127465.data,GSE127465.list)
#GSE127465.integrated <- readRDS("GSE127465.integrated.rds") # load this instead of repeatedly running above steps

DefaultAssay(hnscc.integrated) <- "RNA"
#transfer label with Seurat
hnscc.anchors <- FindTransferAnchors(reference = GSE127465.integrated, query = hnscc.integrated, 
                                     dims = 1:50)
hnscc.predictions <- TransferData(anchorset = hnscc.anchors, refdata = GSE127465.integrated$Minor.subset, # or try Major.cell.type 
                                  dims = 1:50)
hnscc.predictions.major <- TransferData(anchorset = hnscc.anchors, refdata = GSE127465.integrated$Major.cell.type, # or try Major.cell.type 
                                        dims = 1:50)
colnames(hnscc.predictions.major) <- paste0("major.", colnames(hnscc.predictions.major))
hnscc.integrated <- AddMetaData(hnscc.integrated, metadata = hnscc.predictions)
hnscc.integrated <- AddMetaData(hnscc.integrated, metadata = hnscc.predictions.major)
saveRDS(hnscc.integrated,"hnscc.integrated.rds")

# manually refine the major group names
immune.types <- c("bB cells", "bBasophils", "bMonocytes", "bMyeloid precursor-like", "bNeutrophils", "bNK cells", "bpDC", "bT cells", 
                  "tB cells", "tMast cells", "tMoMacDC", "tNeutrophils", "tNK cells", "tpDC", "tPlasma cells", "tT cells")
meta.anno <- hnscc.integrated[[]]
meta.anno$ImmuneOrNot <- ifelse(meta.anno$major.predicted.id %in% immune.types, "Immune Cells", "Non-Immune Cells")
hnscc.integrated <- AddMetaData(hnscc.integrated, metadata = meta.anno[,"ImmuneOrNot", drop = F])
hnscc.integrated$Major.Type <- ifelse(hnscc.integrated$major.predicted.id %in% c("bB cells","tB cells"), "B cells", 
                                      ifelse(hnscc.integrated$major.predicted.id %in% c("bMonocytes","tMoMacDC"), "MoMacDC", 
                                             ifelse(hnscc.integrated$major.predicted.id %in% c("bNK cells","tNK cells"), "NK cells", 
                                                    ifelse(hnscc.integrated$major.predicted.id %in% c("bT cells","tT cells"), "T cells", 
                                                           ifelse(hnscc.integrated$major.predicted.id == "bPlatelets", "Platelets", 
                                                                  ifelse(hnscc.integrated$major.predicted.id == "tMast cells", "Mast cells",
                                                                         ifelse(hnscc.integrated$major.predicted.id == "tNeutrophils", "Neutrophils",
                                                                                ifelse(hnscc.integrated$major.predicted.id == "tpDC", "pDCs",
                                                                                       ifelse(hnscc.integrated$major.predicted.id == "tPlasma cells", "Plasma cells", 
                                                                                              ifelse(hnscc.integrated$major.predicted.id %in% c("bRBC", "tRBC"), "RBCs", 
                                                                                                     ifelse(hnscc.integrated$major.predicted.id %in% c("ND", "Patient1-specific", "Patient3-specific"), "Cancer cells",
                                                                                                            hnscc.integrated$major.predicted.id)))))))))))

hnscc.integrated <- subset(hnscc.integrated, Major.Type %in% setdiff(levels(as.factor(hnscc.integrated$Major.Type)),c("RBCs", "Type II cells", "Platelets")))

hnscc.integrated.v2 <- hnscc.integrated
hnscc.integrated.v2$Major.Type.anno <- hnscc.integrated.v2$Major.Type
major.anno.cols <- c("#008000", "#68ffa9", "#8000f1", "#ff1fc3", "#c2ff82", "#80c600", "#e86c00", "#ff87ab", 
                     "#1b98e0", "#cae9ff", "#ffd53a", "#08fffe", "#e32b00")
names(major.anno.cols) <- levels(as.factor(hnscc.integrated.v2$Major.Type.anno))
Idents(hnscc.integrated.v2) <- hnscc.integrated.v2$Major.Type.anno

hnscc.t.v2 <- subset(hnscc.integrated.v2, subset = tissue == "t")
hnscc.b.v2 <- subset(hnscc.integrated.v2, subset = tissue == "b")
t.cols.v2 <- major.anno.cols[which(names(major.anno.cols) %in% levels(as.factor(hnscc.t.v2$Major.Type.anno)))]
b.cols.v2 <- major.anno.cols[which(names(major.anno.cols) %in% levels(as.factor(hnscc.b.v2$Major.Type.anno)))]
p.t.manual.v2 <- DimPlot(hnscc.t.v2, reduction = "tsne", group.by = "Major.Type.anno", label.size = 3, label = F, label.box = F,  cols = t.cols.v2, repel = TRUE) + 
  NoLegend() + ggtitle("tumor")
p.b.manual.v2 <- DimPlot(hnscc.b.v2, reduction = "tsne", group.by = "Major.Type.anno", label.size = 3, label = F, label.box = F,  cols = b.cols.v2, repel = TRUE) + 
  NoLegend() + ggtitle("blood")
ptb.v2 <- DimPlot(hnscc.integrated.v2, reduction = "tsne",label = F,group.by="tissue", cols = c("#ee6c4d", "#3d5a80"), label.size = 6, label.box = F) + NoLegend()

ptb.v3 <- DimPlot(hnscc.integrated.v2, reduction = "tsne",label = F,group.by="tissue", cols = c("#ee6c4d", "#3d5a80"), label.size = 6, label.box = F) + NoLegend()
p.t.manual.v3 <- DimPlot(hnscc.t.v2, group.by = "Major.Type.anno", pt.size = 0.5, reduction = "tsne") + ggtitle("tumor") + scale_color_manual(values = t.cols.v2, drop = T) + 
  theme(legend.position="right", legend.text = element_text(size =12), aspect.ratio = 1)
p.b.manual.v3 <- DimPlot(hnscc.b.v2, group.by = "Major.Type.anno", pt.size = 0.5, reduction = "tsne") + ggtitle("blood") + scale_color_manual(values = b.cols.v2, drop = T) + 
  theme(legend.position="right", legend.text = element_text(size =12), aspect.ratio = 1)


hnscc.t.markers.v2 <- FindAllMarkers(hnscc.t.v2, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST")
hnscc.t.markers.v2 <- hnscc.t.markers.v2 %>%
  dplyr::filter(p_val_adj < 0.05)
write.table(hnscc.t.markers.v2, file = "hnscc.tumor.major.cell.type.markers.v2.tsv", sep = "\t", quote = F, row.names = F)
t.top20.v2 <- hnscc.t.markers.v2 %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 20, wt = avg_log2FC)

#################################################################################################################
## 05. get subsets of each major type; focus on minor level
#################################################################################################################
#myeloid
hnscc.myeloid <- subset(hnscc.integrated.v2, subset = Major.Type.anno %in% c("Mast cells", "MoMacDC", "pDCs"))
hnscc.myeloid@meta.data <- hnscc.myeloid@meta.data[, which(colnames(hnscc.myeloid@meta.data) %in% c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "patient", "tissue", "treatment", "nCount_SCT", "nFeature_SCT", "seurat_clusters", "predicted.id", "prediction.score.max", "major.predicted.id", "major.prediction.score.max", "first.labels", "labels", "pruned.labels", "reference", "ImmuneOrNot", "Major.Type.anno"))]
DefaultAssay(hnscc.myeloid) <- "integrated"
hnscc.myeloid <- RunPCA(hnscc.myeloid, npcs = 50, assay = "integrated")
hnscc.myeloid <- RunTSNE(hnscc.myeloid, dims = 1:50, verbose = FALSE,tsne.method = "Rtsne")
hnscc.myeloid <- RunUMAP(hnscc.myeloid, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
hnscc.myeloid <- FindNeighbors(hnscc.myeloid, dims = 1:50, verbose = FALSE)
hnscc.myeloid <- FindClusters(hnscc.myeloid, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) # set resolution <1 for fewer clusters (default is 0.5)

best = 0.11
hnscc.myeloid$seurat_clusters=hnscc.myeloid@meta.data[,which(colnames(hnscc.myeloid@meta.data)==paste("integrated_snn_res.",best,sep=""))]
Idents(hnscc.myeloid) <- "seurat_clusters"
p.myeloid.cluster <- DimPlot(hnscc.myeloid, reduction = "tsne",label = T, label.box = F, group.by="seurat_clusters", label.size = 3) + NoLegend() + ggtitle("myeloid cells")

#################################################################################################################
## 06. prepare pancancer myeloid references
#################################################################################################################
#cheng_2021.path <- "/scratch/TBI/Projects/PI_Kim/Kim_3219/reference_data/pancancer_ref/GSE154763"
#novel.cancer.types <- c("ESCA", "KIDNEY", "LYM", "MYE", "OV-FTC", "PAAD", "THCA", "UCEC")
#panc.seurat.list <- list()
#meta.list <- list()
#mtx.list <- list()
#
#for(i in 1:length(novel.cancer.types)){
#  print(novel.cancer.types[i])
#  meta <- read.csv(paste0(cheng_2021.path, "/", "GSE154763_", novel.cancer.types[i], "_metadata.csv.gz"), header = T)
#  meta <- meta %>% 
#    mutate(index = paste0(novel.cancer.types[i], "_", index),
#           Cluster = gsub("^\\w+\\_(\\w+\\_\\w+)$", "\\1", MajorCluster)) %>% 
#    column_to_rownames("index")
#  meta.list[[novel.cancer.types[i]]] <- meta
#  
#  t.norm_mtx <- read.csv(paste0(cheng_2021.path, "/", "GSE154763_", novel.cancer.types[i], "_normalized_expression.csv.gz"), header = T, check.names = F)
#  t.norm_mtx <- t.norm_mtx %>% 
#    mutate(index = paste0(novel.cancer.types[i], "_", index)) %>% 
#    column_to_rownames("index")
#  t.norm_mtx <- t(t.norm_mtx)
#  mtx.list[[novel.cancer.types[i]]] <- as.data.frame(t.norm_mtx) %>% 
#    rownames_to_column("symbol")
#}
#
### organize meta data for LYM and MYE
#meta.list[[3]] <- meta.list[[3]] %>%
#  rownames_to_column("ID") %>%
#  mutate(barcode = gsub("\\w+\\_(\\w+)\\-\\d+", "\\1", ID),
#         library_id = paste(cancer, patient, tissue, sep = "-")) %>%
#  dplyr::select(ID, percent_mito, n_counts, percent_hsp, barcode, batch, library_id, cancer, patient, tissue, n_genes, MajorCluster, source, tech, UMAP1, UMAP2, Cluster) %>%
#  column_to_rownames("ID")
#meta.list[[4]] <- meta.list[[4]] %>%
#  rownames_to_column("ID") %>%
#  mutate(barcode = gsub("\\w+\\_(\\w+)\\-\\d+", "\\1", ID),
#         library_id = paste(cancer, patient, tissue, sep = "-")) %>%
#  dplyr::select(ID, percent_mito, n_counts, percent_hsp, barcode, batch, library_id, cancer, patient, tissue, n_genes, MajorCluster, source, tech, UMAP1, UMAP2, Cluster) %>%
#  column_to_rownames("ID")
#
#panc.meta.all <- do.call(rbind, unname(meta.list))
#
#panc.mtx.all <- Reduce(function(x,y) base::merge(x = x, y = y, by = "symbol", all = T), mtx.list)
#panc.mtx.all <- panc.mtx.all %>% mutate_if(is.numeric , replace_na, replace = 0) %>% column_to_rownames("symbol")
#panc.mtx.all <- as.matrix(panc.mtx.all)
#
#panc.seurat.all <- CreateSeuratObject(counts = panc.mtx.all, 
#                                      project = "pan-cancer", 
#                                      meta.data = panc.meta.all)
#panc.seurat.all$cell.type <- paste0(panc.seurat.all$Cluster, "_", panc.seurat.all$cancer)
##panc.seurat.all <- NormalizeData(panc.seurat.all) #already normalized
#panc.seurat.all <- ScaleData(panc.seurat.all, features = rownames(panc.seurat.all))
#panc.seurat.all <- RunMCA(panc.seurat.all)
#panc.seurat.all <- RunPCA(panc.seurat.all, features = rownames(panc.seurat.all))
#panc.seurat.all <- RunTSNE(panc.seurat.all, dims = 1:50)
#panc.seurat.all <- RunUMAP(panc.seurat.all, dims = 1:50)
##replace UMAP embeddings with published data
#new.embeddings <- as.matrix(panc.seurat.all[[]][,c("UMAP1", "UMAP2")])
#colnames(new.embeddings) <- c("UMAP_1", "UMAP_2")
#panc.seurat.all@reductions$umap@cell.embeddings <- new.embeddings
#UMAP.all.1 <- DimPlot(panc.seurat.all, reduction = "umap", group.by = "Cluster") + ggtitle("pan-cancer", subtitle = "UMAP") + theme(legend.position="right", legend.text = element_text(size =5), aspect.ratio = 1)
#MCA.all <- DimPlot(panc.seurat.all, reduction = "mca", group.by = "Cluster") + ggtitle("pan-cancer", subtitle = "MCA") + theme(legend.position="right", legend.text = element_text(size =5), aspect.ratio = 1)
#UMAP.all <- DimPlot(panc.seurat.all, reduction = "umap", group.by = "cell.type") + ggtitle("pan-cancer", subtitle = "UMAP") + theme(legend.position="right", legend.text = element_text(size =5), aspect.ratio = 1)
#UMAP.all.2 <- DimPlot(panc.seurat.all, reduction = "umap", group.by = "cancer") + ggtitle("pan-cancer", subtitle = "UMAP") + theme(legend.position="right", legend.text = element_text(size =5), aspect.ratio = 1)
#
#pdf(file = paste0("pan-cancer", "_MCA_UMAP_with_legend.pdf"), width = 20, height = 12, onefile = F)
##print(ggarrange(UMAP.all, MCA.all, UMAP.all.1, UMAP.all.2, nrow = 2, common.legend = TRUE, legend = "top"))
#(UMAP.all.1 + MCA.all + UMAP.all + UMAP.all.2)
#dev.off()
#saveRDS(panc.seurat.all, file = "panc.seurat.all.rds")

#################################################################################################################
## 07. annotate myeloid cells with CellID
#################################################################################################################
#get markers for each cell type of pan-cancer data
Idents(panc.seurat.all) <- panc.seurat.all$Cluster
panc.seurat.all.markers <- FindAllMarkers(panc.seurat.all, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,test.use = "wilcox")
write.table(panc.seurat.all.markers, file = "panc.seurat.all.markers.tsv", sep = "\t", quote = F, row.names = F)
panc.seurat.all.markers.top5 <- panc.seurat.all.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 5, wt = avg_log2FC)

###cell type annotation with pan-cancer reference and CellID

DefaultAssay(hnscc.myeloid) <- "RNA"
hnscc.myeloid <- RunMCA(hnscc.myeloid)
original.tSNE <- DimPlot(hnscc.myeloid, reduction = "tsne", group.by = "seurat_clusters")+ ggtitle("tSNE") + theme(legend.text = element_text(size =10), aspect.ratio = 1)
original.MCA <- DimPlot(hnscc.myeloid, reduction = "mca", group.by = "seurat_clusters")  + ggtitle("MCA") + theme(legend.text = element_text(size =10), aspect.ratio = 1)
pdf(file = "hnscc_original_MCA_UMAP.pdf", width = 11, height = 5, onefile = F)
ggarrange(original.MCA, original.tSNE, ncol = 2, common.legend = TRUE, legend = "right")
dev.off()

panc_group_gs <- GetGroupGeneSet(panc.seurat.all, dims = 1:50, n.features = 200, group.by = "Cluster")
HGT_hnscc_group_gs <- RunCellHGT(hnscc.myeloid, pathways = panc_group_gs, dims = 1:50)
hnscc_group_gs_prediction <- rownames(HGT_hnscc_group_gs)[apply(HGT_hnscc_group_gs, 2, which.max)]
hnscc_group_gs_prediction_signif <- ifelse(apply(HGT_hnscc_group_gs, 2, max)>2, yes = hnscc_group_gs_prediction, "Other")
hnscc.myeloid$hnscc_group_gs_prediction <- hnscc_group_gs_prediction_signif
hnscc.myeloid$hnscc_group_gs_prediction <- factor(hnscc.myeloid$hnscc_group_gs_prediction, levels = c("cDC1_CLEC9A", "cDC2_CD1C", "cDC3_LAMP3", "Macro_C1QC", "Macro_FN1", "Macro_GPNMB",
                                                                                                      "Macro_IL1B", "Macro_INHBA", "Macro_ISG15", "Macro_LYVE1", "Macro_NLRP3", "Macro_SPP1",
                                                                                                      "Mono_CD14", "Mono_CD16", "Mono_CD14CD16", "Mast_KIT", "pDC_LILRA4", "Other"))

table(hnscc.myeloid$hnscc_group_gs_prediction)
###refine the subtypes of myeloids based on manually checking of markers overlap
hnscc.myeloid$CellID <- hnscc.myeloid$hnscc_group_gs_prediction
hnscc.myeloid$CellID[which(hnscc.myeloid$CellID == "Macro_FN1")] <- "Mono_CD16" #manually correct improper annotation
hnscc.myeloid$CellID[which(hnscc.myeloid$CellID %in% c("Macro_GPNMB", "Macro_IL1B", "Macro_LYVE1"))] <- "Other" # ignore the cell types with extremely small sample size

hnscc.myeloid$CellID <- droplevels(hnscc.myeloid$CellID)
levels(as.factor(hnscc.myeloid$CellID))


# separate tumor and blood cells to avoid inconsistence
hnscc.t.myeloid <- subset(hnscc.myeloid, subset = tissue == "t")
hnscc.b.myeloid <- subset(hnscc.myeloid, subset = tissue == "b")
hnscc.b.myeloid$CellID[which(hnscc.b.myeloid$CellID %in% c("cDC1_CLEC9A", "Mast_KIT"))] <- "Other" # ignore the cell types with extremely small sample size
hnscc.t.myeloid$CellID <- droplevels(hnscc.t.myeloid$CellID)
hnscc.b.myeloid$CellID <- droplevels(hnscc.b.myeloid$CellID)

ggcolor <- c(
  "Mast_KIT" = "#1688A7", "pDC_LILRA4" = "#7673AE", "cDC1_CLEC9A" = "#b3de69", "cDC2_CD1C" = "#D195F6", "cDC3_LAMP3" = "#7E285E",
  "Macro_C1QC" = "#FD7915", "Macro_INHBA" = "#A443B2", "Macro_ISG15" = "#EB2C1D", "Macro_NLRP3" = "#FFF56A", "Macro_SPP1" = "#FEC718", 
  "Mono_CD14" = "#8197FF", "Mono_CD16" = "#0911E9", "Mono_CD14CD16" = "#1FDBFE", "Other" = "gray"
)
ggcolor1 <- c(
  "Mast_KIT" = "#1688A7", "pDC_LILRA4" = "#7673AE", "cDC1_CLEC9A" = "#b3de69", "cDC2_CD1C" = "#38b000", "cDC3_LAMP3" = "#1a535c",
  "Macro_C1QC" = "#FD7915", "Macro_INHBA" = "#FFE4B5", "Macro_ISG15" = "#EB2C1D", "Macro_NLRP3" = "#FFF56A", "Macro_SPP1" = "#FEC718",
  "Mono_CD14" = "#8197FF", "Mono_CD16" = "#0911E9", "Mono_CD14CD16" = "#1FDBFE", "Other" = "gray"
)

ggPredictions.tsne <- DimPlot(hnscc.myeloid, group.by = "CellID", pt.size = 0.1, reduction = "tsne") + ggtitle("Myeloid cells") + scale_color_manual(values = ggcolor1, drop =FALSE) + 
  theme(legend.position="right", legend.text = element_text(size =12), aspect.ratio = 1)

ggPredictions.umap <- DimPlot(hnscc.myeloid, group.by = "CellID", pt.size = 0.1, reduction = "umap") + NoLegend() + ggtitle("Myeloid cells") + scale_color_manual(values = ggcolor1, drop =FALSE) + 
  theme(legend.position="right", legend.text = element_text(size =12), aspect.ratio = 1)


# markers for subgroups
DefaultAssay(hnscc.t.myeloid) <- "RNA" #very important to change the default assay for DE analysis: "RNA" (the raw data)
#3 annotation: major.predicted.id; predicted.id; labels
Idents(hnscc.t.myeloid) <- hnscc.t.myeloid$CellID
#hnscc.integrated <- NormalizeData(hnscc.integrated, verbose = FALSE)
hnscc.t.myeloid.markers <- FindAllMarkers(hnscc.t.myeloid, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST")
hnscc.t.myeloid.markers <- hnscc.t.myeloid.markers %>%
  dplyr::filter(p_val_adj < 0.05)
write.table(hnscc.t.myeloid.markers, file = "hnscc.t.myeloid.markers_CellID_final.tsv", sep = "\t", quote = F, row.names = F)
t.myeloid.top20 <- hnscc.t.myeloid.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 20, wt = avg_log2FC)
#pdf(file = "HNSCC.t.myeloid.top5.markers.pdf", width = 8, height = 12)
#DoHeatmap(hnscc.t.myeloid, features = t.myeloid.top5$gene, group.by = "myeloid.predicted.id", disp.min = -2, disp.max = 2, raster = F, angle = 90, size = 4) + scale_fill_viridis(option="inferno")
#dev.off()

# markers for subgroups
DefaultAssay(hnscc.b.myeloid) <- "RNA" #very important to change the default assay for DE analysis: "RNA" (the raw data)
#3 annotation: major.predicted.id; predicted.id; labels
Idents(hnscc.b.myeloid) <- hnscc.b.myeloid$CellID
#hnscc.integrated <- NormalizeData(hnscc.integrated, verbose = FALSE)
hnscc.b.myeloid.markers <- FindAllMarkers(hnscc.b.myeloid, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST")
hnscc.b.myeloid.markers <- hnscc.b.myeloid.markers %>%
  dplyr::filter(p_val_adj < 0.05)
write.table(hnscc.b.myeloid.markers, file = "hnscc.b.myeloid.markers_CellID_final.tsv", sep = "\t", quote = F, row.names = F)
b.myeloid.top20 <- hnscc.b.myeloid.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 20, wt = avg_log2FC)


