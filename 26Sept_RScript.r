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
	
library(ggplot2)

top_genes <- unique(top_markers$gene)

DoHeatmap(
  obj,
  features = top_genes,
  group.by = "seurat_clusters"
) +
  ggtitle("Top marker genes across scRNA-seq clusters")
	
	 
MultiAssayExperiment

library(SummarizedExperiment)
library(MultiAssayExperiment)

sample_ids <- paste0("Sample", 1:6)
gene_panel <- c("COL1A1", "ALB", "ACTA2", "TIMP1", "MMP2")

# RNA-seq piece
set.seed(1)
rna_counts <- matrix(
  sample(50:500, length(gene_panel) * 6, replace = TRUE),
  nrow = length(gene_panel),
  dimnames = list(gene_panel, sample_ids)
)
rna_se <- SummarizedExperiment(assays = list(counts = rna_counts))

# methylation piece
set.seed(2)
methyl_values <- matrix(
  round(runif(length(gene_panel) * 6, 0, 1), 2),
  nrow = length(gene_panel),
  dimnames = list(gene_panel, sample_ids)
)
methyl_se <- SummarizedExperiment(assays = list(beta = methyl_values))

# shared sample metadata
shared_sample_metadata <- data.frame(
  condition = rep(c("healthy", "fibrotic"), each = 3),
  row.names = sample_ids
)

mae <- MultiAssayExperiment(
  experiments = list(rnaseq = rna_se, methylation = methyl_se),
  colData = shared_sample_metadata
)

mae