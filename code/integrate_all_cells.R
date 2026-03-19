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

# Load cached objects and split by sample into lists of their RNA assays

fetal_native_liu <- readRDS("data/healthy_fetal_and_healthy_adult.RDS")

DefaultAssay(fetal_native_liu) <- "RNA"
fetal_native_liu_rna <- DietSeurat(fetal_native_liu, assays = c("RNA"))

fetal_native_liu_rna$orig.ident <- fetal_native_liu_rna$source

fetal_native_liu_rna <- SplitObject(fetal_native_liu_rna, split.by = "source")

PMCC_with_clus_marks <- readRDS("data/2d3d_3dorig_combined.RDS")

DimPlot(PMCC_with_clus_marks, group.by = "cell_type", reduction = "int_euc", label = T)

DefaultAssay(PMCC_with_clus_marks) <- "RNA"
PMCC_with_clus_marks_rna <- DietSeurat(PMCC_with_clus_marks, assays = c("RNA"))

PMCC_with_clus_marks_rna$orig.ident <- PMCC_with_clus_marks_rna$sample

PMCC_with_clus_marks_rna <- SplitObject(PMCC_with_clus_marks_rna, split.by = "sample")

# Perform two step SCTransoform and integrate

all_sr <- lapply(all_sr, FUN = function(x){ #
  SCTransform(x, assay= "RNA", vst.flavor = "v2", verbose = FALSE)
})

all_sr <- lapply(all_sr, CellCycleScoring,
                 s.features = cc.genes.updated.2019$s.genes,
                 g2m.features = cc.genes.updated.2019$g2m.genes,
                 set.ident = TRUE, assay = "RNA")

all_sr <-
  lapply(X = all_sr, FUN = function(x) {
    x <- SCTransform(x, assay= "RNA", vars.to.regress = c("S.Score","G2M.Score"), vst.flavor = "v2", verbose = FALSE)
  })

all_sr_feats <- SelectIntegrationFeatures(all_sr)

#Prep for integration
all_sr_int <- PrepSCTIntegration(all_sr, anchor.features = all_sr_feats)

all_sr_int <- lapply(all_sr_int, FUN = RunPCA, features = all_sr_feats)
#Get anchors using our features
all_sr_int_anchors_cca <- FindIntegrationAnchors(object.list = all_sr_int,
                                                 normalization.method = "SCT",
                                                 anchor.features = all_sr_feats,
                                                 reduction = "cca",
                                                 dims = 1:30,
                                                 verbose = T)

all_sr_int_cca <- IntegrateData(anchorset = all_sr_int_anchors_cca,
                                normalization.method = "SCT",
                                features.to.integrate = unique(c(all_sr_feats,
                                                                 cc.genes.updated.2019$s.genes,
                                                                 cc.genes.updated.2019$g2m.genes))
)


all_sr_int_cca <- FindClusters(all_sr_int_cca, resolution = 0.15, graph.name = "integrated_snn", verbose = F, algorithm = 3)

all_sr_int_cca@meta.data <- all_sr_int_cca@meta.data %>%
  mutate(fixed_ct = case_when(
    integrated_snn_res.0.15 == 0 ~ "VIC",
    integrated_snn_res.0.15 == 1 ~ "VIC",
    integrated_snn_res.0.15 == 2 ~ "VIC",
    integrated_snn_res.0.15 == 3 ~ "Myeloid",
    integrated_snn_res.0.15 == 4 ~ "VEC",
    integrated_snn_res.0.15 == 5 ~ "VEC",
    integrated_snn_res.0.15 == 6 ~ "VIC",
    integrated_snn_res.0.15 == 7 ~ "Lymphocyte",
    integrated_snn_res.0.15 == 8 ~ "Cardiomyocyte"
  ))

all_sr_int_cca@meta.data <- all_sr_int_cca@meta.data %>%
  as.data.frame %>%
  mutate(
    origin = case_when(
      str_detect(orig.ident, pattern = "fetal") ~ "fetal",
      str_detect(orig.ident, pattern = "adult") ~ "adult",
      str_detect(orig.ident, pattern = "ORIG") ~ "3D",
      str_detect(orig.ident, pattern = "3DVET") ~ "3D",
      str_detect(orig.ident, pattern = "GEX_2D") ~ "2D",
    )
  )

saveRDS("data/combined_int.RDS")
