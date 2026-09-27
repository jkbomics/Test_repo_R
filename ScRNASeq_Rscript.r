library(Seurat)
library(scRNAseq)

sce <- ZeiselBrainData()
raw_counts_matrix <- assay(sce, "counts")
obj <- CreateSeuratObject(counts = raw_counts_matrix)
GetAssayData(obj, assay = "RNA", layer = "counts")[1:10, 1:10]


obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

obj$nFeature_RNA
obj$nCount_RNA

obj$percent.mt

VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA
", "percent.mt"), ncol = 3)


FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


obj <- subset(obj, nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt <
                10)


GetAssayData(obj, assay = "RNA", layer = "counts")[1:5, 1:5]

obj <- NormalizeData(obj)

length(VariableFeatures(obj))
head(VariableFeatures(obj))

obj <- FindVariableFeatures(obj)

length(VariableFeatures(obj))
head(VariableFeatures(obj))

gene <- VariableFeatures(obj)[1]
x = FetchData(obj, vars = gene)[,1]

mean(x)
sd(x)
obj <- ScaleData(obj)

scaled_x <- GetAssayData(obj, assay = "RNA", layer = "scale.data")[gene, ]
mean(scaled_x)
sd(scaled_x)

Reductions(obj)
obj <- RunPCA(obj)

obj <- FindNeighbors(obj, dims = 1:10)

obj <- FindClusters(obj)


obj <- RunUMAP(obj, dims = 1:10)


DimPlot(obj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(obj, reduction = "pca")

markers <- FindMarkers(obj, ident.1 = 0)

all_markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(all_markers)

library(dplyr)

top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  )

top_markers

top_markers %>%
  group_by(cluster) %>%
  summarise(
    top_genes = paste(gene, collapse = ", ")
  )

top_markers %>%
  group_by(cluster) %>%
  summarise(
    top_genes = paste(gene, collapse = ", ")
  )

top_genes <- unique(top_markers$gene)
DoHeatmap(
  obj,
  features = top_genes,
  group.by = "seurat_clusters"
) +
  ggtitle("Top marker genes across scRNA-seq clusters")

# this is the code for single cell rna seq
