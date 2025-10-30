# ================================================================
# deseq2_luad_pipeline.R
# Author: [Your Name]
# Description:
# 对 TCGA-LUAD 原始 counts 矩阵进行差异表达分析（Tumor vs Normal）
# 输出：DESeq2 结果表（带 lfcShrink）、火山图、PCA、热图、GO 富集结果
# ================================================================

# -----------------------
# 1) 安装/载入包
# -----------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "DESeq2", "tximport", "apeglm", "ggplot2", "pheatmap", "RColorBrewer",
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "biomaRt",
  "dplyr", "tidyr", "tibble", "EnhancedVolcano"
)

install_if_missing <- function(pkgs) {
  to_install <- pkgs[!(pkgs %in% rownames(installed.packages()))]
  if (length(to_install)) BiocManager::install(to_install, ask = FALSE, update = FALSE)
}
install_if_missing(required_pkgs)

# 加载包
sapply(required_pkgs, library, character.only = TRUE)

# -----------------------
# 2) 路径（修改为你的实际路径）
# -----------------------
counts_path <- "E:/Bioinfo/Projects/Prac2. Differential cancer gene expression analysis/TCGA-LUAD/LUAD_counts_matrix.csv"
sampleinfo_path <- "E:/Bioinfo/Projects/Prac2. Differential cancer gene expression analysis/TCGA-LUAD/LUAD_sample_info.csv"
outdir <- "E:/Bioinfo/Projects/Prac2. Differential cancer gene expression analysis/TCGA-LUAD/DESeq2_results"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------
# 3) 读取数据
# -----------------------
cat("📖 正在读取数据...\n")
counts_df <- read.csv(counts_path, stringsAsFactors = FALSE, check.names = FALSE)
sample_info <- read.csv(sampleinfo_path, stringsAsFactors = FALSE)

if (!colnames(counts_df)[1] %in% c("gene_name", "Gene", "X", "gene")) {
  stop("❌ 请确认 counts CSV 的第一列为 gene_name。当前列名：", colnames(counts_df)[1])
}

counts_df <- counts_df %>%
  rename(gene_name = 1) %>%
  distinct(gene_name, .keep_all = TRUE) %>%
  column_to_rownames("gene_name")

sample_info <- sample_info %>% filter(sample_type %in% c("Tumor", "Normal"))
sample_ids <- sample_info$sample_id

missing_in_counts <- setdiff(sample_ids, colnames(counts_df))
if (length(missing_in_counts) > 0) {
  cat("⚠️ 以下样本在 counts 矩阵中未找到，将移除：\n")
  print(missing_in_counts)
  sample_info <- sample_info %>% filter(sample_id %in% colnames(counts_df))
  sample_ids <- sample_info$sample_id
}

counts_mat <- counts_df[, sample_ids]
counts_mat <- apply(counts_mat, 2, function(x) as.integer(round(as.numeric(x))))
rownames(counts_mat) <- rownames(counts_df)

# -----------------------
# 4) 样本信息
# -----------------------
coldata <- sample_info %>%
  select(sample_id, sample_type, patient_id) %>%
  column_to_rownames("sample_id")
coldata$sample_type <- factor(coldata$sample_type, levels = c("Normal", "Tumor"))

# -----------------------
# 5) 创建 DESeqDataSet
# -----------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = coldata,
  design = ~ sample_type
)

# -----------------------
# 6) 预过滤低表达基因
# -----------------------
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]
cat("✅ 过滤后基因数：", nrow(dds), "\n")

# -----------------------
# 7) 运行 DESeq
# -----------------------
dds <- DESeq(dds)

# -----------------------
# 8) 提取结果 + shrink
# -----------------------
res_names <- resultsNames(dds)
cat("🧮 可用结果名称：", res_names, "\n")

shrink_method <- if ("apeglm" %in% rownames(installed.packages())) "apeglm" else "normal"
res <- lfcShrink(dds, coef = "sample_type_Tumor_vs_Normal", type = shrink_method)

res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "gene_name") %>%
  arrange(padj)

write.csv(res_df, file = file.path(outdir, "DESeq2_results_shrunken.csv"), row.names = FALSE)
cat("💾 差异结果已保存\n")

# -----------------------
# 9) 筛选显著差异基因
# -----------------------
deg <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1)
cat("🔥 显著差异基因数量 (padj<0.05, |log2FC|>=1):", nrow(deg), "\n")
write.csv(deg, file = file.path(outdir, "DEGs_padj0.05_lfc1.csv"), row.names = FALSE)

# -----------------------
# 10) PCA 可视化
# -----------------------
vsd <- vst(dds, blind = FALSE)
pdf(file.path(outdir, "PCA_vst.pdf"), width = 6, height = 5)
p <- plotPCA(vsd, intgroup = "sample_type") + ggtitle("PCA (VST) - Tumor vs Normal")
print(p)
dev.off()
cat("📊 PCA 图已保存\n")

# -----------------------
# 11) 火山图
# -----------------------
pdf(file.path(outdir, "volcano_enhancedvolcano.pdf"), width = 7, height = 7)
print(
  EnhancedVolcano(
    res_df,
    lab = res_df$gene_name,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 3.0,
    title = "Volcano plot: Tumor vs Normal",
    subtitle = paste0("nDEG (padj<0.05 & |log2FC|>=1) = ", nrow(deg))
  )
)
dev.off()
cat("🌋 火山图已保存\n")

# -----------------------
# 12) 热图（Top 50 基因）
# -----------------------
topN <- 50
top_genes <- head(deg$gene_name, topN)
mat <- assay(vsd)[top_genes, , drop = FALSE]
mat_z <- t(scale(t(mat)))
annotation_col <- as.data.frame(coldata["sample_type"])

pdf(file.path(outdir, "heatmap_topDEGs.pdf"), width = 8, height = 10)
pheatmap(
  mat_z, annotation_col = annotation_col,
  show_rownames = TRUE, show_colnames = FALSE,
  clustering_distance_rows = "correlation",
  clustering_method = "complete",
  main = paste0("Top ", topN, " DEGs (z-score)")
)
dev.off()
cat("🔥 热图已保存\n")

# -----------------------
# 13) 基因 ID 映射
# -----------------------
res_df_nodup <- res_df %>% distinct(gene_name, .keep_all = TRUE)
map_res <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = res_df_nodup$gene_name,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "SYMBOL"
)
res_mapped <- left_join(res_df_nodup, map_res, by = c("gene_name" = "SYMBOL"))

# -----------------------
# 14) GO 富集
# -----------------------
deg_entrez <- res_mapped %>%
  filter(!is.na(ENTREZID), padj < 0.05 & abs(log2FoldChange) >= 1) %>%
  pull(ENTREZID) %>%
  unique()

if (length(deg_entrez) >= 10) {
  ego <- enrichGO(
    gene = deg_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  write.csv(as.data.frame(ego), file = file.path(outdir, "GO_BP_enrichment.csv"), row.names = FALSE)
  cat("📈 GO 富集结果已保存\n")
} else {
  cat("⚠️ 差异基因数太少 (<10)，跳过 GO 富集\n")
}

# -----------------------
# 15) 保存对象与会话信息
# -----------------------
saveRDS(dds, file = file.path(outdir, "dds_object.rds"))
saveRDS(vsd, file = file.path(outdir, "vsd_object.rds"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("🎉 All done! 结果保存在：", outdir, "\n")
