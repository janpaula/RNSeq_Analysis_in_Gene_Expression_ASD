
# =============================================================================
# TRANSCRIPTOMIC ANALYSIS - ASD vs CONTROL
# =============================================================================

library(DESeq2)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tibble)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(STRINGdb)
library(igraph)
library(dplyr)
library(biomaRt)

# =============================================================================
# 1. DATA LOADING
# =============================================================================

counts <- read.table(
  "/home/janicerpaula/Bioinformática/TreinamentoAnaliseTranscriptomicaRNA/GSEs3057364018/GSE64018_countlevel_12asd_12ctl.txt.gz",
  header = TRUE,
  row.names = 1
)

# =============================================================================
# 2.SAMPLES DEFINITION
# =============================================================================

ASD <- c(
  "AN02987_ba41.42.22_8.6", "AN04682_ba41.42.22_8.2",
  "UMB5278_ba41.42.22_7.8", "AN01570_ba41.42.22_7.1",
  "AN00493_ba41.42.22_7.3", "AN12457_ba41.42.22_6",
  "AN08166_ba41.42.22_6.4", "AN08792_ba41.42.22_5.1",
  "AN01971_ba41.42.22_2.9", "AN03632_ba41.42.22_8.1",
  "AN08043_ba41.42.22_7.8", "AN09714_ba41.42.22_7.2"
)

Control <- c(
  "AN17425_ba41.42.22_7.9", "AN07444_ba41.42.22_8.2",
  "UMB4590_ba41.42.22_8.3", "AN10833_ba41.42.22_1.8",
  "AN14757_ba41.42.22_8",   "AN19760_ba41.42.22_6.1",
  "AN12137_ba41.42.22_6.4", "UMB5079_ba41.42.22_7.9",
  "AN08161_ba41.42.22_8.1", "AN08677_ba41.42.22_8",
  "UMB4842_ba41.42.22_8.1", "AN13295_ba41.42.22_6.3"
)

# =============================================================================
# 3. COLDATA AND DESEQ2
# =============================================================================

# Classifica as colunas da matriz counts (contagens de reads) como "ASD" ou "Control" usando ifelse aninhado. 
# Ele verifica se cada nome de coluna está nos vetores ASD ou Control; caso contrário, atribui NA.
condition <- ifelse(colnames(counts) %in% ASD, "ASD",
                    ifelse(colnames(counts) %in% Control, "Control", NA)) 


# Cria um data.frame colData com os nomes das colunas como rownames e o fator condition. 
# Isso é esencial para descrever as amostras experimentalmente.
colData <- data.frame(
  row.names = colnames(counts),
  condition = factor(condition)
)


# O DESeqDataSetFromMatrix constrói o objeto DESeq a partir das contagens e metadados 
# e usando "~ condition" como design (comparação entre grupos). 
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = colData,
  design    = ~ condition
)

# DESeq executa a normalização e estimação de parâmetros, enquanto 
dds <- DESeq(dds)
res <- results(dds)            # Gera a tabela de resultados  
summary(res)                   # Mostra estatísticas como número de genes diferencialmente expressos.

# =============================================================================
# 4. GLOBAL ENSEMBL -> SYMBOL CONVERSION
#    Stage 1: org.Hs.eg.db (fast, offline)
#    Stage 2: biomaRt (covers non-coding genes and recent IDs)
#    Stage 3: fallback to the ENSEMBL ID itself
# =============================================================================

all_ensg <- rownames(counts)

# --- Stage 1: org.Hs.eg.db (fast, offline) ---
symbols_orgdb <- mapIds(
  org.Hs.eg.db,
  keys      = all_ensg,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

missing_after_orgdb <- all_ensg[is.na(symbols_orgdb) | symbols_orgdb == ""]
cat("Genes without symbolo after org.Hs.eg.db:", length(missing_after_orgdb), "\n")

# --- Stage 2: biomaRt (covers non-coding genes and recent IDs) ---
# `useEnsembl() + mirrors` is more robust than `useMart()` and tolerates server failures.
symbols_biomart <- setNames(rep(NA_character_, length(missing_after_orgdb)),
                            missing_after_orgdb)

if (length(missing_after_orgdb) > 0) {
  
  mirrors <- c("www", "useast", "asia")
  mart    <- NULL
  
  for (mirror in mirrors) {
    mart <- tryCatch({
      message("Tentando mirror: ", mirror, ".ensembl.org ...")
      useEnsembl(
        biomart = "ensembl",
        dataset = "hsapiens_gene_ensembl",
        mirror  = mirror
      )
    }, error = function(e) {
      message("  Falhou (", mirror, "): ", conditionMessage(e))
      NULL
    })
    if (!is.null(mart)) {
      message("  Conectado: ", mirror, ".ensembl.org")
      break
    }
  }
  
  if (!is.null(mart)) {
    bm <- tryCatch(
      getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol"),
        filters    = "ensembl_gene_id",
        values     = missing_after_orgdb,
        mart       = mart
      ),
      error = function(e) {
        message("getBM failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(bm)) {
      bm <- bm[bm$hgnc_symbol != "", ]
      bm <- bm[!duplicated(bm$ensembl_gene_id), ]
      symbols_biomart[bm$ensembl_gene_id] <- bm$hgnc_symbol
    }
    cat("Genes still without symbols after biomaRt:",
        sum(is.na(symbols_biomart) | symbols_biomart == ""),
        "(It will be displayed as ENSEMBL ID.)\n")
  } else {
    message("All mirrors failed - skipping stage 2.")
  }
}

# --- Assemble the final unified vector ---
symbols <- symbols_orgdb
symbols[missing_after_orgdb] <- symbols_biomart[missing_after_orgdb]

# Ensure that the vector has correct names.
names(symbols) <- all_ensg

# --- Auxiliar Function (Stage 3: fallback to the ENSEMBL ID itself) ---
ensg_to_sym <- function(ensg_vec) {
  sym <- symbols[ensg_vec]
  # Fallback: NA or empty -> use your own ENSEMBL ID
  bad <- is.na(sym) | sym == ""
  sym[bad] <- ensg_vec[bad]
  unname(sym)
}

gene_table <- data.frame(
  ENSEMBL = all_ensg,
  SYMBOL  = unname(symbols)
)

cat("Full symbol coverage:",
    sum(!is.na(symbols) & symbols != ""), "/", length(symbols), "\n")



# =============================================================================
# 5. VST E PCA GERAL
# =============================================================================

# The vst(dds) function applies Variance Stabilizing Transformation to the normalized 
# data of the DESeq2 object, making the variances more homogeneous for downstream analyses such as PCA.
vsd <- vst(dds)

# The plotPCA automatically generates a principal component plot colored by the experimental condition.
plotPCA(vsd, intgroup = "condition")


# --- Genes Drivers PC1 ---
# Performs manual PCA with prcomp on the transposed matrix of VST values. 
pca      <- prcomp(t(assay(vsd)))

# Loadings (rotations) indicate the contributions of each gene to the PCs.
loadings <- pca$rotation

# Selects the 50 genes with the highest absolute loadings in PC1 (top_PC1), converts 
# Ensembl IDs to gene symbols via the custom function ensg_to_sym, and prepares the matrix for heatmap:

top_PC1      <- head(order(abs(loadings[, 1]), decreasing = TRUE), 50)
genes_PC1    <- rownames(loadings)[top_PC1]        # ENSG (internal use)
symbols_PC1  <- ensg_to_sym(genes_PC1)             # symbol (plots)

# Matrix for heatmap - lines renamed to symbol
mat             <- assay(vsd)[genes_PC1, ]
rownames(mat)   <- symbols_PC1

# Heatmap
pheatmap(mat,
         scale        = "row",
         fontsize_row = 5,
         fontsize_col = 8,
         angle_col    = 45,
         main         = "Top 50 genes - PC1")

# =============================================================================
# 6. PCA (ASD + CLUSTERING)
# =============================================================================

asd_samples <- colData$condition == "ASD"
vsd_asd     <- vsd[, asd_samples]

pca_asd <- prcomp(t(assay(vsd_asd)))
plot(pca_asd$x[, 1], pca_asd$x[, 2],
     xlab = "PC1", ylab = "PC2",
     main = "PCA - only ASD")

dist_mat <- dist(t(assay(vsd_asd)))
hc       <- hclust(dist_mat)
plot(hc)

# Heatmap ASD — lines with símbolos
mat_asd           <- assay(vsd)[genes_PC1, asd_samples]
rownames(mat_asd) <- symbols_PC1

pheatmap(mat_asd,
         scale        = "row",
         fontsize_row = 5,
         fontsize_col = 8,
         angle_col    = 45,
         main         = "Top PC1 genes - ASD samples")

# =============================================================================
# 7. INFLAMATORY SCORE
# =============================================================================

inflam_genes <- c("C1QA","C1QB","C1QC","CD14","ITGB2",
                  "VSIG4","SERPINA1","S100A8","S100A9")

# Convert to ENSEMBL
inflam_ensg       <- names(symbols)[symbols %in% inflam_genes]
inflam_ensg_valid <- inflam_ensg[inflam_ensg %in% rownames(vsd)]

mat_inflam        <- assay(vsd)[inflam_ensg_valid, ]
mat_inflam_scaled <- t(scale(t(mat_inflam)))
inflam_score      <- colMeans(mat_inflam_scaled, na.rm = TRUE)

sample_profile <- ifelse(inflam_score >  0.5, "Inflammatory",
                         ifelse(inflam_score < -0.5, "Neuronal", "Mixed"))

# --- PCA colored by profile ---
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

colors <- ifelse(sample_profile == "Inflammatory", "red",
                 ifelse(sample_profile == "Neuronal", "blue", "gold"))
colors[is.na(colors)] <- "grey"

plot(pca$x[, 1], pca$x[, 2],
     col  = colors, pch = 19, cex = 1.5,
     xlab = paste0("PC1 (", round(var_exp[1] * 100, 1), "%)"),
     ylab = paste0("PC2 (", round(var_exp[2] * 100, 1), "%)"),
     main = "PCA - ASD Cortex")

legend("topright",
       legend = c("Inflammatory","Neuronal","Mixed"),
       col    = c("red","blue","gold"),
       pch    = 19, bty = "n")

# =============================================================================
# 8. DESeq2 BY PROFILE (Inflammatory / Neuronal / Mixed)
# =============================================================================

# Assign to colData of the original dds
colData(dds)$profile <- factor(sample_profile[colnames(dds)])
table(colData(dds)$profile)

design(dds) <- ~ profile
dds         <- DESeq(dds)

# Comparations
res_inf_vs_neu <- results(dds, contrast = c("profile","Inflammatory","Neuronal"))
res_neu_vs_inf <- results(dds, contrast = c("profile","Neuronal","Inflammatory"))
res_mix_vs_inf <- results(dds, contrast = c("profile","Mixed","Inflammatory"))
res_mix_vs_neu <- results(dds, contrast = c("profile","Mixed","Neuronal"))

# DEGs significatives
sig_inf <- res_inf_vs_neu[which(res_inf_vs_neu$padj < 0.05 & res_inf_vs_neu$log2FoldChange >  1), ]
sig_neu <- res_inf_vs_neu[which(res_inf_vs_neu$padj < 0.05 & res_inf_vs_neu$log2FoldChange < -1), ]

# Top genes with symbol
top_inf <- head(sig_inf[order(sig_inf$padj), ], 20)
top_neu <- head(sig_neu[order(sig_neu$padj), ], 20)

top_inf_symbols <- ensg_to_sym(rownames(top_inf))
top_neu_symbols <- ensg_to_sym(rownames(top_neu))

cat("Top inflammatories:\n"); print(top_inf_symbols)
cat("Top neuronals:\n");     print(top_neu_symbols)

# Heatmap of markers - lines with symbols
genes_markers <- c(rownames(top_inf), rownames(top_neu))
mat_markers   <- assay(vsd)[genes_markers, ]
rownames(mat_markers) <- ensg_to_sym(genes_markers)

annotation_col <- data.frame(
  Profile = sample_profile,
  row.names = colnames(vsd)
)

pheatmap(mat_markers,
         scale          = "row",
         annotation_col = annotation_col,
         show_rownames  = TRUE,
         fontsize_row   = 6,
         main           = "Markers by profile")

# =============================================================================
# 9. VOLCANO PLOT
# =============================================================================

volcano_df <- as.data.frame(res_inf_vs_neu)
volcano_df$gene   <- rownames(volcano_df)
volcano_df$symbol <- ensg_to_sym(volcano_df$gene)
volcano_df        <- volcano_df[!is.na(volcano_df$padj), ]

volcano_df$regulation <- "NS"
volcano_df$regulation[volcano_df$padj < 0.05 & volcano_df$log2FoldChange >  1] <- "Up"
volcano_df$regulation[volcano_df$padj < 0.05 & volcano_df$log2FoldChange < -1] <- "Down"

table(volcano_df$regulation)

# label 
genes_up   <- volcano_df[volcano_df$regulation == "Up",   ]
genes_down <- volcano_df[volcano_df$regulation == "Down",  ]
genes_up   <- genes_up[order(-genes_up$log2FoldChange),  ]
genes_down <- genes_down[order(genes_down$log2FoldChange), ]

cat("Top UP:\n");   print(head(genes_up[,   c("symbol","log2FoldChange","padj")], 20))
cat("Top DOWN:\n"); print(head(genes_down[, c("symbol","log2FoldChange","padj")], 20))

# Results
volcano_final <- volcano_df[, c("gene","symbol","log2FoldChange","padj","regulation")]
write.csv(volcano_final, "volcano_results.csv", row.names = FALSE)

# Plot
point_colors <- ifelse(volcano_df$regulation == "Up", "blue",
                       ifelse(volcano_df$regulation == "Down", "red", "gray"))

plot(volcano_df$log2FoldChange,
     -log10(volcano_df$padj),
     pch  = 20, col = point_colors,
     xlab = "Log2 Change of Expression",
     ylab = "-Log10 (P-value adjusted)",
     main = "Volcano Plot: Inflammatory vs Neuronal")

abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-1, 1),     col = "black", lty = 2)

legend("topright",
       legend = c("Genes UP","Genes Down","No Alteration"),
       col    = c("blue","red","gray"),
       pch    = 20)

# =============================================================================
# 10. GO ENRICHMENT
# =============================================================================

entrez_inf <- mapIds(org.Hs.eg.db, keys = rownames(sig_inf),
                     column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
entrez_neu <- mapIds(org.Hs.eg.db, keys = rownames(sig_neu),
                     column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

entrez_inf <- na.omit(entrez_inf)
entrez_neu <- na.omit(entrez_neu)

ego_inf <- enrichGO(gene = entrez_inf, OrgDb = org.Hs.eg.db,
                    ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,
                    readable = TRUE)   # readable=TRUE → símbolos na tabela

ego_neu <- enrichGO(gene = entrez_neu, OrgDb = org.Hs.eg.db,
                    ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,
                    readable = TRUE)

dotplot(ego_inf, showCategory = 15, title = "GO - Inflammatory") +
  theme(
    legend.spacing.y = unit(0.4, "cm"),   # espaço vertical entre itens da legenda
    legend.key.height = unit(0.6, "cm"),  # altura de cada item
    legend.text = element_text(size = 10) # tamanho do texto
  ) +
  guides(
    size  = guide_legend(byrow = TRUE),   # força espaçamento por linha
    color = guide_colorbar(barheight = 5) # altura da barra de cor
  )

dotplot(ego_neu, showCategory = 15, title = "GO - Neuronal") +
  theme(
    legend.spacing.y = unit(0.4, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.text = element_text(size = 10)
  ) +
  guides(
    size  = guide_legend(byrow = TRUE),
    color = guide_colorbar(barheight = 5)
  )

# =============================================================================
# 11. GENES DRIVERS (PCA + DEG)
# =============================================================================

drivers_inf <- intersect(genes_PC1, rownames(sig_inf))
drivers_neu <- intersect(genes_PC1, rownames(sig_neu))

loading_PC1 <- loadings[, 1]

df_inf <- data.frame(
  gene    = drivers_inf,
  symbol  = ensg_to_sym(drivers_inf),
  loading = loading_PC1[drivers_inf],
  padj    = sig_inf[drivers_inf, "padj"]
)
df_inf <- df_inf[order(-abs(df_inf$loading)), ]

df_neu <- data.frame(
  gene    = drivers_neu,
  symbol  = ensg_to_sym(drivers_neu),
  loading = loading_PC1[drivers_neu],
  padj    = sig_neu[drivers_neu, "padj"]
)
df_neu <- df_neu[order(-abs(df_neu$loading)), ]

cat("Top Drivers Inflamatories:\n"); print(head(df_inf, 10))
cat("Top Drivers Neuronals:\n");     print(head(df_neu, 10))

# =============================================================================
# 12. REDE PPI (STRINGdb)
# =============================================================================

genes_all <- unique(c(drivers_inf, drivers_neu))
genes_all <- as.character(genes_all)

df_string <- data.frame(ensembl_gene_id = genes_all)

string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)

mapped <- string_db$map(df_string, "ensembl_gene_id", removeUnmappedRows = TRUE)
mapped <- mapped %>% distinct(STRING_id, .keep_all = TRUE)

interactions <- string_db$get_interactions(mapped$STRING_id)
interactions <- interactions %>%
  filter(from %in% mapped$STRING_id & to %in% mapped$STRING_id)

g <- graph_from_data_frame(interactions[, c("from","to")], directed = FALSE)
g <- delete_vertices(g, degree(g) == 0)

deg <- degree(g)
btw <- betweenness(g, normalized = TRUE)

top_hubs <- sort(deg, decreasing = TRUE)[1:20]
top_btw  <- sort(btw, decreasing = TRUE)[1:20]

# --- STRING -> ENSG -> SÍMBOLO ---
string_to_gene <- dplyr::select(mapped, STRING_id, ensembl_gene_id)

# Create vector labels directly with symbol (without extra biomaRt)
labels_raw <- string_to_gene$ensembl_gene_id[
  match(V(g)$name, string_to_gene$STRING_id)
]
labels <- ensg_to_sym(labels_raw)
labels[is.na(labels) | labels == ""] <- V(g)$name   # fallback STRING_id

# Finals Drivers in symbol
hub_genes_ensg <- string_to_gene$ensembl_gene_id[
  match(names(top_hubs), string_to_gene$STRING_id)
]
hub_genes_ensg       <- unique(hub_genes_ensg[!is.na(hub_genes_ensg)])
drivers_final_ensg   <- intersect(hub_genes_ensg, genes_PC1)
drivers_final_symbols <- ensg_to_sym(drivers_final_ensg)

cat("Final Drivers (symbol):\n"); print(drivers_final_symbols)

# --- Visualization ---
is_hub <- V(g)$name %in% names(top_hubs)[1:10]

plot(g,
     vertex.size        = log(deg + 1) * 3,
     vertex.label       = labels,
     vertex.label.cex   = 0.9,
     vertex.label.color = "black",
     vertex.color       = ifelse(is_hub, "red", "lightblue"),
     main               = "PPI Network - ASD Cortex")

# =============================================================================
# 13. DETECTION OF CLUSTERS 
# =============================================================================

net_clusters    <- cluster_louvain(g)
n_clusters      <- length(net_clusters)
cat("Number of clusters:", n_clusters, "\n")

cluster_list    <- split(V(g)$name, net_clusters$membership)

# Convertion STRING_id -> ENSG -> símbolo
cluster_genes_ensg <- lapply(cluster_list, function(ids) {
  ensg <- string_to_gene$ensembl_gene_id[match(ids, string_to_gene$STRING_id)]
  na.omit(unique(ensg))
})

cluster_genes_sym <- lapply(cluster_genes_ensg, ensg_to_sym)
print(cluster_genes_sym)

# =============================================================================
# 14. ENRICHMENT CLUSTER (GO, keyType ENSEMBL)
# =============================================================================

go_results <- lapply(cluster_genes_ensg, function(genes) {
  if (length(genes) < 3) return(NULL)
  tryCatch(
    enrichGO(
      gene          = genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENSEMBL",
      ont           = "BP",
      pAdjustMethod = "BH",
      readable      = TRUE       # convert to symbol 
    ),
    error = function(e) NULL
  )
})

# --- Rótulo automático ---
cluster_labels <- sapply(go_results, function(res) {
  if (is.null(res) || nrow(res) == 0) return("Unknown")
  terms <- res@result$Description[1:min(5, nrow(res))]
  if (any(grepl("immune|inflammat|complement", terms, ignore.case = TRUE))) return("Inflammatory")
  if (any(grepl("synapse|neuron|neurotransmitter", terms, ignore.case = TRUE))) return("Neuronal")
  if (any(grepl("stress|protein folding|heat shock", terms, ignore.case = TRUE))) return("Stress")
  return("Other")
})

cat("Label of clusters:\n"); print(cluster_labels)

# --- Colorir rede por cluster ---
cluster_colors <- c(
  "Inflammatory" = "red",
  "Neuronal"     = "blue",
  "Stress"       = "orange",
  "Unknown"      = "grey",
  "Other"        = "grey"
)

V(g)$color <- cluster_colors[cluster_labels[net_clusters$membership]]

plot(g,
     vertex.size      = 8,
     vertex.label     = labels,
     vertex.label.cex = 0.6,
     main             = "Clusters annotated for biological function")

# --- Resume ---
cluster_summary <- data.frame(
  Cluster = names(cluster_genes_sym),
  Label   = cluster_labels,
  Size    = sapply(cluster_genes_sym, length),
  TopGenes = sapply(cluster_genes_sym, function(x) paste(head(x, 5), collapse = ", "))
)
print(cluster_summary)
write.csv(cluster_summary, "cluster_summary.csv", row.names = FALSE)
