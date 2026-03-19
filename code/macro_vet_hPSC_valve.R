library(tidyverse)
library(SingleCellExperiment)
library(Seurat)
library(sctransform)
library(DropletUtils)
library(scater)
library(scuttle)
library(RColorBrewer)
library(clustree)
library(Azimuth)
library(patchwork)
library(scCustomize)
library(UCell)
library(ggalluvial)
library(monocle3)
library(SeuratWrappers)

# Read h5 files as SCE objects

raw_mat_fp <- c("data/2d/cellranger_out/2D-Valve-Cells-GEX/outs/raw_feature_bc_matrix_cellbender_output.h5",
                "data/3d/cellranger_out/3D-VET-GEX/outs/raw_feature_bc_matrix_cellbender_output.h5",
                "data/3d_orig/run_1/outs/raw_feature_bc_matrix_cellbender_output.h5",
                "data/3d_orig/run_2/outs/raw_feature_bc_matrix_cellbender_output.h5",
                "data/3d_orig/run_3/outs/raw_feature_bc_matrix_cellbender_output.h5")

samples_names <- c("GEX_2D",
                   "GEX_3DVET",
                   "ORIG_3D1",
                   "ORIG_3D2",
                   "ORIG_3D3")

all_reads <- lapply(FUN = Read_CellBender_h5_Mat, raw_mat_fp)

all_reads <- lapply(all_reads, FUN = function(x){
  x[,colSums(x) > 0]
})

colnames(all_reads[[1]]) <- paste(colnames(all_reads[[1]]), "_2D", sep = "")
colnames(all_reads[[2]]) <- paste(colnames(all_reads[[2]]), "_3D", sep = "")
colnames(all_reads[[3]]) <- paste(colnames(all_reads[[3]]), "_3D_ORIG1", sep = "")
colnames(all_reads[[4]]) <- paste(colnames(all_reads[[4]]), "_3D_ORIG2", sep = "")
colnames(all_reads[[5]]) <- paste(colnames(all_reads[[5]]), "_3D_ORIG3", sep = "")

all_sce <- vector(mode = "list", length = 5)

all_sce <- lapply(all_reads, FUN = function(x){
  SingleCellExperiment(assays=list("counts"=x))
})

all_sce <- lapply(all_reads, FUN = function(x){
  x$sample <- samples_names[1]
})

all_sce <- mapply(x = all_sce, FUN = function(x, ns){
  x$sample <- ns
  x
}, ns = samples_names)

all_sce_together <- cbind(all_sce[[1]], all_sce[[2]], all_sce[[3]], all_sce[[4]], all_sce[[5]])

all_reads <- all_sce_together

is_mito <- str_detect(pattern = "^MT-", string = rownames(all_reads))
all_reads <- addPerCellQCMetrics(all_reads, subsets=list(Mito=is_mito))

df.qc <- perCellQCMetrics(all_reads, subsets = list(Mito=is_mito))

mt <- df.qc$subsets_Mito_percent > 20
umis <- isOutlier(df.qc$sum, type = "both", log = T, nmads = 3)# %>% table()
detected <- isOutlier(df.qc$detected, type = "lower", log = T, nmads = 3)# %>% table()

discard <- mt | umis | detected

seurat_count_mat <- all_reads@assays@data$counts[,!discard]
rownames(seurat_count_mat) <- rownames(all_reads)

# Now we have a Seurat Object
sce_seurat <- Seurat::CreateSeuratObject(counts = seurat_count_mat,
                                         meta.data = as.data.frame(colData(all_reads))[!discard,])

sce_seurat[["percent.mt"]] <- PercentageFeatureSet(sce_seurat, pattern = "^MT-")

sce_seurat |>
  VlnPlot(features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample")

sce_seurat <- subset(sce_seurat, nFeature_RNA > 1000 & nCount_RNA < 100000 & percent.mt < 10)

sce_seurat |>
  VlnPlot(features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample")

sce_seurat_split <- sce_seurat %>%
  SplitObject(split.by = "sample")

sce_seurat_split <- lapply(sce_seurat_split, FUN = function(x){
  NormalizeData(x, verbose = FALSE)
})

# Perform two pass SCTransform, transform -> score CC -> transform + regress
sce_seurat_split <- lapply(sce_seurat_split, SCTransform, vst.flavor = "v2")

sce_seurat_split <- lapply(sce_seurat_split, FUN = function(x){
  DefaultAssay(x) <- "SCT"
  x
})

sce_seurat_split <- lapply(sce_seurat_split, CellCycleScoring,
                           s.features = cc.genes.updated.2019$s.genes,
                           g2m.features = cc.genes.updated.2019$g2m.genes,
                           set.ident = TRUE, assay = "RNA")

sce_seurat_split <-
  lapply(X = sce_seurat_split, FUN = function(x) {
    x <- SCTransform(x, assay="RNA", vars.to.regress = c("S.Score","G2M.Score"), verbose = FALSE)
  })

# Set up integration
sce_seurat_split_features <- SelectIntegrationFeatures(sce_seurat_split)

sce_seurat_split_int <- PrepSCTIntegration(sce_seurat_split, anchor.features = sce_seurat_split_features)

sce_seurat_split_int <- lapply(sce_seurat_split_int, FUN = RunPCA, features = sce_seurat_split_features)
sce_seurat_split_int_anchors_cca <- FindIntegrationAnchors(object.list = sce_seurat_split_int,
                                                           normalization.method = "SCT",
                                                           anchor.features = sce_seurat_split_features,
                                                           reduction = "cca",
                                                           dims = 1:30,
                                                           verbose = T)

sce_seurat_split_int_cca <- IntegrateData(anchorset = sce_seurat_split_int_anchors_cca,
                                          normalization.method = "SCT",
                                          features.to.integrate = unique(c(sce_seurat_split_features,
                                                                           cc.genes.updated.2019$s.genes,
                                                                           cc.genes.updated.2019$g2m.genes)))

# Prep for find markers
sce_seurat_split_int_cca <- PrepSCTFindMarkers(sce_seurat_split_int_cca)

sce_seurat_split_int_cca <- sce_seurat_split_int_cca %>%
  RunPCA(dims = 1:30)

ElbowPlot(sce_seurat_split_int_cca)

sce_seurat_split_int_cca <- sce_seurat_split_int_cca %>%
  RunUMAP(dims = 1:11)

DimPlot(sce_seurat_split_int_cca, reduction = "umap", group.by = "Phase", order = F) +
  DimPlot(sce_seurat_split_int_cca, reduction = "umap", group.by = "sample", order = F)

sce_seurat_split_int_cca <- RunUMAP(sce_seurat_split_int_cca, dims = 1:14, assay = "SCT", reduction.name = "int_euc", metric = "euclidean", n.neighbors = 20L)

DefaultAssay(sce_seurat_split_int_cca) <- "SCT"

sce_seurat_split_int_cca <- FindNeighbors(sce_seurat_split_int_cca,
                                          reduction = 'pca',
                                          dims = 1:30,
                                          assay = 'SCT', #just to be saaafe
                                          verbose = TRUE
)

sce_seurat_split_int_cca <- FindClusters(sce_seurat_split_int_cca, resolution = 0.2, graph.name = "integrated_snn", verbose = F, algorithm = 3)

saveRDS(object = sce_seurat_split_int_cca, "data/2d3d_3dorig_combined.RDS")
