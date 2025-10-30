# deseq2_luad_pipeline.R
# 说明：对你生成的 LUAD 原始 counts 矩阵做差异表达分析（Tumor vs Normal）
# 输出：DESeq2 结果表（带 lfcShrink）、火山图、PCA、top gene heatmap、GO 富集结果与可视化

# -----------------------
# 1) 安装/载入包（仅第一次需要 run 下面安装部分）
# -----------------------
required_pkgs <- c(
  "DESeq2", "tximport", "apeglm", "ggplot2", "pheatmap", "RColorBrewer",
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "biomaRt", "dplyr", "tidyr",
  "tibble", "EnhancedVolcano"
)

install_if_missing <- function(pkgs){
  to_install <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if(length(to_install)) BiocManager::install(to_install, ask = FALSE, update = FALSE)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install_if_missing(required_pkgs)

# 加载包
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(apeglm)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(dplyr)
library(tidyr)
library(tibble)
library(EnhancedVolcano)

# -----------------------
# 2) 路径：修改为你本地路径（如果已在同一目录可使用相对路径）
# -----------------------
counts_path <- "E:/Bioinfo/Projects/Prac2. Differential cancer gene expression analysis/TCGA-LUAD/LUAD_counts_matrix.csv"    # 你生成的 counts 矩阵
sampleinfo_path <- "E:/Bioinfo/Projects/Prac2. Differential cancer gene expression analysis/TCGA-LUAD/LUAD_sample_info.csv"  # 你生成的样本信息表
outdir <- "E:/Bioinfo/Projects/Prac2. Differential cancer gene expression analysis/TCGA-LUAD/DESeq2_results"                 # 输出目录
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------
# 3) 读取数据
# -----------------------
cat("Reading data...\n")
counts_df <- read.csv(counts_path, stringsAsFactors = FALSE, check.names = FALSE)
sample_info <- read.csv(sampleinfo_path, stringsAsFactors = FALSE)

# 确认首列是 gene_name
if(colnames(counts_df)[1] %in% c("gene_name", "Gene", "X", "gene")){
  # good
} else {
  stop("请确认 counts CSV 的第一列为 gene_name（基因名）。当前列名：", colnames(counts_df)[1])
}

# 把 gene_name 设为行名，并移除重复基因名（保留第一个）
counts_df <- counts_df %>% 
  rename(gene_name = 1) %>%
  distinct(gene_name, .keep_all = TRUE) %>%
  column_to_rownames("gene_name")

# 保证列顺序与 sample_info 中 sample_id 一致（并只保留 Tumor/Normal）
sample_info <- sample_info %>% filter(sample_type %in% c("Tumor", "Normal"))
sample_ids <- sample_info$sample_id

# 检查所有 sample_id 是否在 counts 列中
missing_in_counts <- setdiff(sample_ids, colnames(counts_df))
if(length(missing_in_counts)>0){
  cat("警告：以下样本在 counts 矩阵中未找到，脚本将移除它们：\n")
  print(missing_in_counts)
  sample_info <- sample_info %>% filter(sample_id %in% colnames(counts_df))
  sample_ids <- sample_info$sample_id
}

# 取出 counts 矩阵对应列
counts_mat <- counts_df[, sample_ids]

# 确保计数为整数（DESeq2 需要原始整数 counts）
counts_mat <- apply(counts_mat, 2, function(x) as.integer(round(as.numeric(x))))
rownames(counts_mat) <- rownames(counts_df)

# -----------------------
# 4) 构建 colData（样本注释）并确保 factor 的对照水平（reference level）
# -----------------------
coldata <- sample_info %>%
  select(sample_id, sample_type, patient_id) %>%
  column_to_rownames("sample_id")

# 把 sample_type 变成 factor，并设定 Normal 为对照（reference）
coldata$sample_type <- factor(coldata$sample_type, levels = c("Normal", "Tumor"))

# -----------------------
# 5) 创建 DESeqDataSet
#    设计式： ~ sample_type
#    若将来需要控制批次可改为 ~ batch + sample_type
# -----------------------
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData = coldata,
                              design = ~ sample_type)

# -----------------------
# 6) 预过滤（去掉低表达基因，常见阈值：至少在 n 样本中 counts >= 10）
# -----------------------
keep <- rowSums(counts(dds) >= 10) >= 5   # 至少 5 个样本表达量 >=10，可调整
dds <- dds[keep, ]
cat("过滤后基因数：", nrow(dds), "\n")

# -----------------------
# 7) 运行 DESeq（估计 size factors、dispersion、拟合模型并做 Wald 检验）
# -----------------------
dds <- DESeq(dds)

# -----------------------
# 8) 提取结果（默认对比： Tumor vs Normal，因为 factor levels 设定）
#    使用 lfcShrink 进行 log2 fold change 收缩，推荐方法 apeglm（若可用）
# -----------------------
res_names <- resultsNames(dds)
cat("结果名：", res_names, "\n")

# 选择 shrink 方法（优先 apeglm）
shrink_method <- if("apeglm" %in% installed.packages()) "apeglm" else "normal"
res <- lfcShrink(dds=dds, coef="sample_type_Tumor_vs_Normal", type=shrink_method)

# 加入基因名列
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "gene_name") %>%
  arrange(padj)

# 保存差异结果（全部）
write.csv(res_df, file = file.path(outdir, "DESeq2_results_shrunken.csv"), row.names = FALSE)
cat("差异结果已保存：", file.path(outdir, "DESeq2_results_shrunken.csv"), "\n")

# -----------------------
# 9) 标记显著差异基因（例：padj < 0.05 & abs(log2FoldChange) >= 1）
# -----------------------
deg <- res_df %>% filter(!is.na(padj)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1)

cat("显著差异基因数量 (padj<0.05, |log2FC|>=1):", nrow(deg), "\n")
write.csv(deg, file = file.path(outdir, "DEGs_padj0.05_lfc1.csv"), row.names = FALSE)

# -----------------------
# 10) PCA 可视化（基于 vst 变换）
# -----------------------
vsd <- vst(dds, blind = FALSE)  # variance stabilizing transform
pdf(file.path(outdir, "PCA_vst.pdf"), width = 6, height = 5)
plotPCA(vsd, intgroup = "sample_type") + ggtitle("PCA (VST) - Tumor vs Normal")
dev.off()
cat("PCA 图已保存\n")

# -----------------------
# 11) 火山图（使用 EnhancedVolcano 或 ggplot2）
# -----------------------
pdf(file.path(outdir, "volcano_enhancedvolcano.pdf"), width = 7, height = 7)
EnhancedVolcano(res_df,
                lab = res_df$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0,
                title = "Volcano plot: Tumor vs Normal",
                subtitle = paste0("nDEG (padj<0.05 & |log2FC|>=1) = ", nrow(deg)))
dev.off()
cat("火山图已保存\n")

# -----------------------
# 12) 热图：取 top N 差异基因的表达
# -----------------------
topN <- 50
top_genes <- head(deg$gene_name, topN)
mat <- assay(vsd)[top_genes, , drop=FALSE]
# 标准化每基因 (row) 为 z-score
mat_z <- t(scale(t(mat)))
annotation_col <- as.data.frame(coldata["sample_type"])
rownames(annotation_col) <- rownames(coldata)

pdf(file.path(outdir, "heatmap_topDEGs.pdf"), width = 8, height = 10)
pheatmap(mat_z, annotation_col = annotation_col,
         show_rownames = TRUE, show_colnames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_method = "complete",
         main = paste0("Top ", topN, " DEGs (z-score)"))
dev.off()
cat("热图已保存\n")

# -----------------------
# 13) 基因ID转换（从基因名 -> ENTREZ ID）以便做 clusterProfiler 富集
# -----------------------
res_df_nodup <- res_df %>% distinct(gene_name, .keep_all = TRUE)
map_res <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = res_df_nodup$gene_name,
                                 columns = c("ENTREZID","SYMBOL"),
                                 keytype = "SYMBOL")
res_mapped <- left_join(res_df_nodup, map_res, by = c("gene_name" = "SYMBOL"))

# -----------------------
# 14) GO 富集（以 DEG 的 ENTREZID 列表）
# -----------------------
deg_entrez <- res_mapped %>%
  filter(!is.na(ENTREZID), padj < 0.05 & abs(log2FoldChange) >= 1) %>%
  pull(ENTREZID) %>% unique()

if(length(deg_entrez) >= 10){
  ego <- enrichGO(gene = deg_entrez,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  write.csv(as.data.frame(ego), file = file.path(outdir, "GO_BP_enrichment.csv"), row.names = FALSE)
  cat("GO 富集结果已保存\n")
  
  # GO 富集结果可视化
  pdf(file.path(outdir, "GO_enrichment_BP.pdf"), width = 7, height = 7)
  barplot(ego, showCategory = 20)
  dev.off()
  cat("GO 富集图已保存\n")
} else {
  cat("差异基因数太少，跳过 GO 富集（需要 >=10 个 ENTREZID）\n")
}

# -----------------------
# 15) 保存会话信息与结束
# -----------------------
saveRDS(dds, file = file.path(outdir, "dds_object.rds"))
saveRDS(vsd, file = file.path(outdir, "vsd_object.rds"))
cat("All done. Results in:", outdir, "\n")
