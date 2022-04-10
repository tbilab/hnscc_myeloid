require(tidyverse)
require(Seurat)
require(SeuratWrappers)
require(cowplot)
require(loomR)
require(SeuratDisk)
require(velocyto.R)


wrkdir <- "/scratch/TBI/Projects/PI_Kim/Kim_3219/update_pancancer"
dtadir <- "/scratch/TBI/Projects/PI_Kim/sc_analysis/velocyto"

setwd(wrkdir)
source("../R/singlecell_viper_fun.R")

# load gene expression seurat object
hnscc.integrated.v2 <- readRDS("hnscc.integrated.updated_07282021.rds")
hnscc.myeloid <- readRDS("hnscc.myeloid.updated_08262021.rds")
hnscc.integrated.v2@meta.data <- hnscc.integrated.v2@meta.data[, which(colnames(hnscc.integrated.v2@meta.data) %in% c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "patient", "tissue", "treatment", "nCount_SCT", "nFeature_SCT", "integrated_snn_res.0.14", "seurat_clusters", "ImmuneOrNot", "doublet_scores", "predicted_doublets", "Major.Type.anno"))]
hnscc.myeloid@meta.data <- hnscc.myeloid@meta.data[, which(colnames(hnscc.myeloid@meta.data) %in% c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "patient", "tissue", "treatment", "nCount_SCT", "nFeature_SCT", "integrated_snn_res.0.14", "seurat_clusters", "ImmuneOrNot", "Major.Type.anno", "CellID"))]
saveRDS(hnscc.integrated.v2,"hnscc_integrated_submit.rds")

#all.meta <- hnscc.integrated.v2@meta.data
#mye.meta <- hnscc.myeloid@meta.data
#mye.meta$CellID <- as.character(mye.meta$CellID)
#CellID <- ifelse(rownames(all.meta) %in% rownames(mye.meta), mye.meta$CellID, "NA")

hnscc.myeloid <- subset(hnscc.myeloid, tissue == "t")

## set the input file path
loom_paths <- c()
loom_paths["pat01.t"] <- paste0(dtadir, "/3219-SR-1/possorted_genome_bam_1L4TR.loom")
#loom_paths["pat01.post"] <- paste0(dtadir, "/3219-SR-2/filtered_feature_bc_matrix")
loom_paths["pat02.t"] <- paste0(dtadir, "/3219-SR-3/possorted_genome_bam_6N45X.loom")
#loom_paths["pat02.post"] <- paste0(dtadir, "/3219-SR-4/filtered_feature_bc_matrix")
loom_paths["pat03.t"] <- paste0(dtadir, "/3219-SR-5/possorted_genome_bam_ZW8JK.loom")
#loom_paths["pat03.post"] <- paste0(dtadir, "/3219-SR-6/filtered_feature_bc_matrix")
loom_paths["pat04.t"] <- paste0(dtadir, "/3521-SR-1/possorted_genome_bam_39PM6.loom")
#loom_paths["pat04.post"] <- paste0(dtadir, "/3521-SR-2/filtered_feature_bc_matrix")
loom_paths["pat05.t"] <- paste0(dtadir, "/Kim_2853_1/possorted_genome_bam_OO5MN.loom")
loom_paths["pat05.b"] <- paste0(dtadir, "/Kim_2853_2/possorted_genome_bam_UF7SI.loom")
loom_paths["pat06.t"] <- paste0(dtadir, "/Kim_2853_3/possorted_genome_bam_BR8P9.loom")
loom_paths["pat06.b"] <- paste0(dtadir, "/Kim_2853_4/possorted_genome_bam_WJXIB.loom")


# read in all loom files, transform into seurat object, rename cell barcodes

velo.seurat.list <- list()

for( i in names(loom_paths) ) {
  print(loom_paths[[i]])
  tmp <- ReadVelocity(file = loom_paths[[i]])
  velo.seurat.list[[i]] <- as.Seurat(tmp)
  velo.seurat.list[[i]] <- RenameCells(velo.seurat.list[[i]], new.names = paste0(i, "_", gsub("^\\w+\\:(\\w+)x$", "\\1", colnames(velo.seurat.list[[i]])), "-1"))  
}

#merge all velo seurat objects
velo.integrate <- merge(x = velo.seurat.list[[1]], y = c(velo.seurat.list[[2]], velo.seurat.list[[3]], velo.seurat.list[[4]], 
                                                         velo.seurat.list[[5]], velo.seurat.list[[6]], velo.seurat.list[[7]], velo.seurat.list[[8]]), merge.data = TRUE)

#filter out cells which are not held in hnscc.integrated
velo.integrate$filtered <- colnames(velo.integrate) %in% colnames(hnscc.integrated.v2)
velo.integrate <- subset(velo.integrate, subset = filtered == "TRUE")

#filter out features which are not held in hnscc.integrated
velo.integrate@assays$spliced <- CreateAssayObject(velo.integrate@assays$spliced[rownames(hnscc.integrated.v2),])
velo.integrate@assays$unspliced <- CreateAssayObject(velo.integrate@assays$unspliced[rownames(hnscc.integrated.v2),])
velo.integrate@assays$ambiguous <- CreateAssayObject(velo.integrate@assays$ambiguous[rownames(hnscc.integrated.v2),])

spliced <- CreateAssayObject(GetAssayData(velo.integrate, assay = "spliced"))
unspliced <- CreateAssayObject(GetAssayData(velo.integrate, assay = "unspliced"))
ambiguous <- CreateAssayObject(GetAssayData(velo.integrate, assay = "ambiguous"))

#get 3 matrices separately
hnscc.integrated.v2[["spliced"]] <- spliced
hnscc.integrated.v2[["unspliced"]] <- unspliced
hnscc.integrated.v2[["ambiguous"]] <- ambiguous

names(hnscc.integrated.v2@meta.data) <- gsub("\\.", "\\_", names(hnscc.integrated.v2@meta.data))

saveRDS(hnscc.integrated.v2,"hnscc_integrated_add_velo_submit.rds")
SaveH5Seurat(hnscc.integrated.v2, filename = "hnscc_velo_submit.h5Seurat")
Convert("hnscc_velo_submit.h5Seurat", dest = "h5ad")


#work on MoMacDC cluster
velo.integrate$myeloids <- colnames(velo.integrate) %in% colnames(hnscc.myeloid)
velo.myeloids <- subset(velo.integrate, subset = myeloids == "TRUE")
velo.myeloids@assays$spliced <- CreateAssayObject(velo.myeloids@assays$spliced[rownames(hnscc.myeloid),])
velo.myeloids@assays$unspliced <- CreateAssayObject(velo.myeloids@assays$unspliced[rownames(hnscc.myeloid),])
velo.myeloids@assays$ambiguous <- CreateAssayObject(velo.myeloids@assays$ambiguous[rownames(hnscc.myeloid),])

myeloids_spliced <- CreateAssayObject(GetAssayData(velo.myeloids, assay = "spliced"))
myeloids_unspliced <- CreateAssayObject(GetAssayData(velo.myeloids, assay = "unspliced"))
myeloids_ambiguous <- CreateAssayObject(GetAssayData(velo.myeloids, assay = "ambiguous"))

hnscc.myeloid[["spliced"]] <- myeloids_spliced
hnscc.myeloid[["unspliced"]] <- myeloids_unspliced
hnscc.myeloid[["ambiguous"]] <- myeloids_ambiguous

names(hnscc.myeloid@meta.data) <- gsub("\\.", "\\_", names(hnscc.myeloid@meta.data))
#hnscc.myeloid$hnscc_group_gs_prediction <- as.character(hnscc.myeloid$hnscc_group_gs_prediction)
hnscc.myeloid$CellID <- as.character(hnscc.myeloid$CellID)
saveRDS(hnscc.myeloid,"hnscc_tumor_myeloids_submit.rds")
SaveH5Seurat(hnscc.myeloid, filename = "hnscc_tumor_myeloids_velo_submit.h5Seurat")
Convert("hnscc_tumor_myeloids_velo_submit.h5Seurat", dest = "h5ad")


