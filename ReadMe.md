# Scripts and input files used in HNSCC myeloid project

## Kim_5482

Mouse BMDM were incubated with Annexin V+/PI- Cal27 cells stained with PKH26+ dye for 1hr. After this, macrophages were washed with ice cold DPBS to remove apoptotic tumor cells that were not engulfed by the BMDM. After additional 6hr, macrophages were harvested and processed for flow sorting to obtain BMDM that had engulfed AC (PKH26+). RNA was isolated using Trizol and used for sequencing.

### Kim_5482_new_metadata.tsv

Sample information file

### Kim_5482_new_featurecount_symbol.tsv

Reads count matrix (FeatureCounts)

### Kim_5482_DE_submit.R

R script for differential gene analysis

## Kim_4970

CD11b+ myeloid cells were flow sorted from previously untreated HNSCC tumors and matched blood (n=7). Briefly, tumors were digested using human tumor digestion buffer following manufacturer’s instructions to obtain single cells. Following this, cells were stained with live dead dye, CD45, CD11b, CD3 and NK1.1. RNA was isolated from sorted CD11b+ cells using Trizol and used for sequencing at VANTAGE next generation sequencing facility at VU.

### Kim_4970_core_metadata.tsv

Sample information file

### Kim_4970_core_featurecount_symbol.tsv

Reads count matrix (FeatureCounts)

### DE_Kim_4970_submit.R

R script for differential gene analysis

## Kim_3219

Head and neck cancer patient samples were weighed and cut into pieces prior to digestion using human tumor digestion mix for 1hr on the tumor dissociator. Cells were then passed through a 70uM cell strainer to remove debris and centrifuged at 300g for 5mins. RBC was lysed with RBC lysis buffer. Cells were then centrifuged again and resuspended at 600cells/uL PBS and dropped off for Chromium Single Cell 3’ Library construction (10X Genomics) and sequencing at VANTAGE sequencing facility at VU following the manufacturer’s instructions. Libraries were sequenced on an Illumina HiSeq4000 and mapped to the human genome (build GRCh38) by CellRanger (Version 3.0).

### Kim_3219_2853_meta.tsv

Sample information file

### hnscc_pretreated_scRNAseq_submit.R

R script for single cell RNAseq data analysis

### hnscc_velocity.R

R script for data format converting and integration of the output matrices of velocity analysis

### hnscc_scVelo_submit.ipynb

Jupyter notebook of scVelo analysis pipeline

## Reference

* Efferocytosis drives myeloid NLRP3 dependent inflammasome signaling and gasdermin D independent secretion of IL-1β to promote tumor growth. Cara Lang*, Sohini Roy*, Yu Wang*, Diana Graves, Yaomin Xu, Henrique Serazani, Michael Korrer, Young J. Kim. 2022, submitted.



