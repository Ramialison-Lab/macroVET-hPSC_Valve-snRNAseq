library(parallel)
library(future)
library(SeuratObject)
library(Seurat)
library(tidyverse)
library(SeuratDisk)
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
library(scran)
library(batchelor)
library(edgeR)
library(scales)

options(future.globals.maxSize = 8000 * 1024^2)

# Read combined object
pmcc_3d_native_sr_obj <- readRDS("data/combined_int.RDS")

DimPlot(pmcc_3d_native_sr_obj, reduction = "int_euc", group.by = "integrated_snn_res.0.15")

c_pal <- brewer.pal(n = length(unique(pmcc_3d_native_sr_obj$fixed_ct)), name = "Set2")
names(c_pal) <- unique(pmcc_3d_native_sr_obj$fixed_ct)


pdf("data/plots/reann_umap.pdf", height = 6, width = 7)
DimPlot(pmcc_3d_native_sr_obj,
        group.by = "fixed_ct",
        reduction = "int_euc")
dev.off()

# umap
pdf("data/plots/fixed_ct_split_origin.pdf", height = 4, width = 12)
DimPlot(pmcc_3d_native_sr_obj, group.by = "fixed_ct", split.by = "origin", reduction = "int_euc", cols = c_pal) & xlab("UMAP 1") & ylab("UMAP 2")
dev.off()

# cell proportions
prop_tbl <- pmcc_3d_native_sr_obj@meta.data |>
  group_by(fixed_ct, origin) |>
  summarise(total = n()) |>
  group_by(origin) |>
  mutate(freq_pct = total / sum(total) * 100)

# Create vector with levels of object@ident
identities <- unique(pmcc_prop_tbl_clus$cell_type_use)

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

# Plot the umap plot with the default ggplot2 colors
DimPlot(object = PMCC_with_clus_marks, reduction = "int_euc") +
  scale_color_manual(values = my_color_palette)
names(my_color_palette) <- c(
  "Cardiomyocyte",
  "Lymphocyte",
  "Myeloid",
  "VEC",
  "VIC"
)

pdf("data/plots/prop_simplified.pdf", height = 5, width = 3.5)
prop_tbl |>
  ggplot(aes(x = origin, y = total, fill = fixed_ct)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c_pal) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

mark_VIC <- c("COL1A1", "COL3A1", "DCN", "VIM")
mark_VEC <- c("PECAM1", "CDH5", "SELE", "SERPINE1")
mark_Myeloid <- c("CD68", "MRC1", "CD163", "CD14")
mark_Mast_cells <- c("CPA3", "KIT", "TPSAB1", "TPSB2")
mark_MyoFB <- c("ACTA2", "MYH11", "MYL9")
mark_Lymphocytes <- c("CD2", "CD3E", "CD3D", "CD7")
mark_CM <- c("TTN", "TNNT2", "MYH7", "MYL2", "MYL9")

Idents(pmcc_3d_native_sr_obj) <- pmcc_3d_native_sr_obj$fixed_ct
pdf("data/plots/dotplot_abridged_markers.pdf", height = 4, width = 12)
DotPlot(pmcc_3d_native_sr_obj, features = list("VIC" = mark_VIC, "VEC" = mark_VEC, "Myeloid" = mark_Myeloid, "Lymphocytes" = mark_Lymphocytes, "Cardiomyocytes" = mark_CM), group.by = "fixed_ct", cluster.idents = F) +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pmcc_3d_native_sr_obj <- PrepSCTFindMarkers(pmcc_3d_native_sr_obj)

pdf("data/plots/native_adult_only_4Genes.pdf", height = 4, width = 5)
DotPlot(subset(pmcc_3d_native_sr_obj_vic_cca, subset = origin %in% c("Adult_native")),
        features = c("EBF2", "MEOX2", "NTRK2", "RHOU"), cluster.idents = F) + #Check with Holly, is this the right gene. Is the one listed in the email a typo
  scale_color_viridis_c(option = "plasma", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_point(aes(size = pct.exp),
             color = "black",
             stroke = 0.3,
             shape = 1,
             show.legend = FALSE)
dev.off()


# some supplementary
psc_only <- subset(pmcc_3d_native_sr_obj_vic_cca, subset = origin %in% c("GEX_2D", "GEX_3D"))
psc_only <- PrepSCTFindMarkers(psc_only)
pdf("data/plots/PSC_4Genes.pdf", height = 4, width = 5)
DotPlot(psc_only,
        features = c("EBF2", "MEOX2", "NTRK2", "RHOU"), cluster.idents = F) +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_point(aes(size = pct.exp),
             color = "black",
             stroke = 0.3,
             shape = 1,
             show.legend = FALSE)
dev.off()


DefaultAssay(psc_only) <- "RNA"

Idents(pmcc_3d_native_sr_obj_vic_cca) <- pmcc_3d_native_sr_obj_vic_cca$origin
pmcc_3d_native_sr_obj_vic_cca <- PrepSCTFindMarkers(pmcc_3d_native_sr_obj_vic_cca)

pdf("data/plots/VIC_only_4gene_dotplot_scaleF.pdf", height = 4, width = 5)
DotPlot(pmcc_3d_native_sr_obj_vic_cca,
        features = c("EBF2", "MEOX2", "NTRK2", "RHOU"), cluster.idents = F, assay = "SCT", scale = F) +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_point(aes(size = pct.exp),
             color = "black",
             stroke = 0.3,
             shape = 1,
             show.legend = FALSE)
dev.off()
