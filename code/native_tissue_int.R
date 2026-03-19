library(tidyverse)
library(Seurat)
library(sctransform)
library(SingleCellExperiment)
library(DropletUtils)
library(tidySingleCellExperiment)
library(scater)
library(scuttle)
library(RColorBrewer)
library(clustree)
library(patchwork)
library(scCustomize)
library(ggalluvial)
library(SeuratWrappers)
library(Azimuth)
library(UCell)
library(scDotPlot)


mtx_files <- list.files("data/GSE228638/", pattern = "mtx.gz", full.names = T)

all_samps <- basename(mtx_files) |>
  gsub(pattern = "_matrix.mtx.gz", replacement = "")

create_sr <- function(x){
  CreateSeuratObject(Seurat::ReadMtx(mtx = x[3], cells = x[1], features = x[2]))
}

all_samps_sr <- vector("list", length(all_samps))

for(i in 1:length(all_samps_sr)){
  print(i)
  all_samps_sr[[i]] <- create_sr(list.files("data/GSE228638", pattern = all_samps[i], full.names = T))
}

fetal_meta <- readr::read_tsv("data/GSE228638/GSE228638_MetaData_HealthyFetalValve.txt.gz")
older_meta <- readr::read_tsv("data/GSE228638/GSE228638_MetaData_HealthyandPAIVS_PulValve.txt.gz")

older_meta <- older_meta |>
  dplyr::rename("Cells" = "x", "CellTypes" = "CellType")


all_samps_sr[[1]] <- AddMetaData(all_samps_sr[[1]], fetal_meta)
all_samps_sr[[2]] <- NULL
pctmt <- function(x) {
  PercentageFeatureSet(x, pattern = "^MT-")
}


for(i in 1:length(all_samps_sr)){
  print(i)
  all_samps_sr[[i]][["percent.mt"]] <- pctmt(all_samps_sr[[i]])
}

#Reverse engineered the thresholds.
#nFeature_RNA > 1000
#nCount_RNA < 100000
#percent.mt < 10
all_samps_sr[[3]] |>
  VlnPlot(features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

for(i in 1:length(all_samps_sr)){
  all_samps_sr[[i]] <- subset(all_samps_sr[[i]], nFeature_RNA > 1000)
  all_samps_sr[[i]] <- subset(all_samps_sr[[i]], nCount_RNA < 100000)
  all_samps_sr[[i]] <- subset(all_samps_sr[[i]], percent.mt < 10)
}


all_samps_sr[[2]] |>
  VlnPlot(features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

all_samps_sr[[1]]$source <- "ctrl1_healthy"
all_samps_sr[[2]]$source <- "ctrl2_healthy"
all_samps_sr[[3]]$source <- "ctrl2_healthy"
all_samps_sr[[4]]$source <- "ctrl2_healthy"
all_samps_sr[[5]]$source <- "ctrl2_healthy"

all_samps_sr[[1]]$orig.ident <- "ctrl1_healthy"
all_samps_sr[[2]]$orig.ident <- "ctrl2_healthy"
all_samps_sr[[3]]$orig.ident <- "ctrl2_healthy"
all_samps_sr[[4]]$orig.ident <- "ctrl2_healthy"
all_samps_sr[[5]]$orig.ident <- "ctrl2_healthy"


all_samps_sr[[2]]$Tissue <- "Aortic Valve"
all_samps_sr[[3]]$Tissue <- "Pulmonary Valve"
all_samps_sr[[4]]$Tissue <- "Mitral Valve"
all_samps_sr[[5]]$Tissue <- "Tricuspid Valve"

all_samps_sr_sct <- lapply(all_samps_sr, SCTransform, vst.flavor = "v2", )

all_samps_sr_sct <- lapply(all_samps_sr_sct, CellCycleScoring, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE, assay = "SCT")

all_samps_sr_sct <- lapply(X = all_samps_sr_sct, FUN = function(x) {
  x <- SCTransform(x, assay="RNA", vars.to.regress = c("S.Score","G2M.Score","percent.mt"), verbose = FALSE)
})


all_samps_sr_sct_features <- SelectIntegrationFeatures(all_samps_sr_sct)

#Prep for integration
all_samps_sr_sct_int <- PrepSCTIntegration(all_samps_sr_sct, anchor.features = all_samps_sr_sct_features)

all_samps_sr_sct_int <- lapply(all_samps_sr_sct_int, FUN = RunPCA, features = all_samps_sr_sct_features)
#Get anchors using our features
all_samps_sr_sct_int_anchors_cca <- FindIntegrationAnchors(object.list = all_samps_sr_sct_int,
                                                           normalization.method = "SCT",
                                                           anchor.features = all_samps_sr_sct_features,
                                                           reduction = "cca",
                                                           dims = 1:30,
                                                           verbose = T)

all_samps_sr_sct_int_cca <- IntegrateData(anchorset = all_samps_sr_sct_int_anchors_cca,
                                          normalization.method = "SCT",
                                          features.to.integrate = unique(c(all_samps_sr_sct_features,
                                                                           cc.genes.updated.2019$s.genes,
                                                                           cc.genes.updated.2019$g2m.genes))
)

# saveRDS(object = all_samps_sr_sct_int_cca, "data/healthy_fetal.RDS")

all_samps_sr_sct_int <- all_samps_sr_sct_int %>%
  RunPCA(dims = 1:50)
#UMAP
all_samps_sr_sct_int <- all_samps_sr_sct_int %>%
  RunUMAP(dims = 1:50)

#Plot and colour by phase
DimPlot(all_samps_sr_sct_int, reduction = "umap", group.by = "CellTypes", order = T)
DimPlot(all_samps_sr_sct_int, reduction = "umap", group.by = "Tissue", order = F)

FeaturePlot(all_samps_sr_sct_int, reduction = "umap", features = "nCount_RNA", order = F)


DimPlot(all_samps_sr_sct_int, reduction = "umap", group.by = "Phase", order = F) +
  DimPlot(all_samps_sr_sct_int, reduction = "umap", group.by = "source", order = F, pt.size = 0.5) +
  FeaturePlot(all_samps_sr_sct_int, reduction = "umap", features = "nCount_RNA", order = F, max.cutoff = 15000)

# Native Adult

native_adult <- readRDS("data/healthy_hv_int.RDS")
DefaultAssay(native_adult) <- "RNA"
native_adult <- DietSeurat(native_adult, layers = "counts", assays = "RNA")

native_adult@meta.data <- native_adult@meta.data |>
  select(orig.ident, percent.mt, nFeature_RNA, nCount_RNA, S.Score, G2M.Score, Phase, samp, VIC_UCell, VEC_UCell, monocyte_UCell, myeloid_UCell, macrophage_UCell, lymphocytes_UCell, mast_UCell, myoFB_UCell)

native_adult@meta.data <- native_adult@meta.data |>
  mutate(orig.ident = if_else(samp == "h1", "adult1_healthy", "adult2_healthy")) |>
  mutate(source = if_else(samp == "h1", "adult1_healthy", "adult2_healthy"))

native_adult$Tissue <- "Aortic Valve"

native_adult_list <- list(
  subset(native_adult, samp == "h1"),
  subset(native_adult, samp == "h2")
)


all_samps_sr[[1]]$source <- "fetal1_healthy"
all_samps_sr[[2]]$source <- "fetal2_healthy"
all_samps_sr[[3]]$source <- "fetal2_healthy"
all_samps_sr[[4]]$source <- "fetal2_healthy"
all_samps_sr[[5]]$source <- "fetal2_healthy"

all_samps_sr[[1]]$orig.ident <- "ctrl1_healthy"
all_samps_sr[[2]]$orig.ident <- "ctrl2_healthy"
all_samps_sr[[3]]$orig.ident <- "ctrl2_healthy"
all_samps_sr[[4]]$orig.ident <- "ctrl2_healthy"
all_samps_sr[[5]]$orig.ident <- "ctrl2_healthy"

all_samps_sr[[2]]$Tissue <- "Aortic Valve"
all_samps_sr[[3]]$Tissue <- "Pulmonary Valve"
all_samps_sr[[4]]$Tissue <- "Mitral Valve"
all_samps_sr[[5]]$Tissue <- "Tricuspid Valve"

all_samp_together <- c(all_samps_sr, native_adult_list)

for(i in 1:length(all_samp_together)){
  print(i)
  all_samp_together[[i]][["percent.mt"]] <- pctmt(all_samp_together[[i]])
}

for(i in 1:length(all_samp_together)){
  all_samp_together[[i]] <- subset(all_samp_together[[i]], nFeature_RNA > 750)
  all_samp_together[[i]] <- subset(all_samp_together[[i]], percent.mt < 5)
}

all_samp_together_sct <- lapply(all_samp_together, SCTransform, vst.flavor = "v2")

all_samp_together_sct <- lapply(all_samp_together_sct, CellCycleScoring, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE, assay = "SCT")

all_samp_together_sct <- lapply(X = all_samp_together_sct, FUN = function(x) {
  x <- SCTransform(x, assay="RNA", vars.to.regress = c("S.Score","G2M.Score","percent.mt"), verbose = FALSE)
})


#Use these features for integration & PCA
all_samp_together_sct_features <- SelectIntegrationFeatures(all_samp_together_sct)

all_samp_together_sct <- lapply(all_samp_together_sct, FUN = RunPCA, features = all_samp_together_sct_features)

all_samp_together_sct_anchors <- vector(mode = "list", length = 7)

for(i in 2:7){
  all_samp_together_sct_anchors[[i]] <- FindTransferAnchors(query = all_samp_together_sct[[i]], reference = all_samp_together_sct[[1]], dims = 1:30, reference.reduction = "pca")
}

predictions <- vector(mode = "list", length = 7)

for(i in 2:7){
  predictions[[i]] <- TransferData(anchorset = all_samp_together_sct_anchors[[i]], refdata = all_samp_together_sct[[1]]$CellTypes, dims = 1:30)
}

for(i in 2:7){
  all_samp_together_sct[[i]] <- AddMetaData(all_samp_together_sct[[i]], metadata = predictions[[i]])
}

for(i in 2:7){
  all_samp_together_sct[[i]]$CellTypes <- all_samp_together_sct[[i]]$predicted.id
}


#Prep for integration
all_samp_together_sct_int <- PrepSCTIntegration(all_samp_together_sct, anchor.features = all_samp_together_sct_features)

#Get anchors using our features
all_samp_together_sct_int_anchors_cca <- FindIntegrationAnchors(object.list = all_samp_together_sct_int,
                                                       normalization.method = "SCT",
                                                       anchor.features = all_samp_together_sct_features,
                                                       reduction = "cca",
                                                       dims = 1:30,
                                                       verbose = T)

all_samp_together_sct_int_cca <- IntegrateData(anchorset = all_samp_together_sct_int_anchors_cca,
                                      normalization.method = "SCT",
                                      features.to.integrate = unique(c(all_samp_together_sct_features,
                                                                       cc.genes.updated.2019$s.genes,
                                                                       cc.genes.updated.2019$g2m.genes))
                                      )

# Save the object for later so we don't have to reprocess from scratch
saveRDS(object = all_samp_together_sct_int_cca, "data/healthy_fetal_and_healthy_adult.RDS")

