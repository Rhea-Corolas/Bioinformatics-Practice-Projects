import pandas as pd
from pathlib import Path
import os

# =====================
# 路径设置
# =====================
map_path = Path(r"E:\Bioinfo\Projects\Prac2. Differential cancer gene expression analysis\TCGA-LUAD\file_sample_map_with_type.csv")
data_dir = Path(r"E:\Bioinfo\Projects\Prac2. Differential cancer gene expression analysis\TCGA-LUAD\Gene expression")

# =====================
# 读取映射表
# =====================
meta_df = pd.read_csv(map_path)
print(f"共有 {len(meta_df)} 个样本映射信息")
print(meta_df["sample_type"].value_counts())

# =====================
# 逐个读取并合并
# =====================
expr_dict = {}
gene_index = None
read_count = 0
missing_files = []

for i, row in meta_df.iterrows():
    file_path = data_dir / row["file"]
    sample_id = row["sample_id"]

    if not file_path.exists():
        missing_files.append(str(file_path))
        continue

    # 跳过注释行 "# ..." 并指定 tab 分隔
    df = pd.read_csv(file_path, sep="\t", comment="#")

    # 保留关键列
    if "gene_name" not in df.columns or "unstranded" not in df.columns:
        print(f"警告！ 文件格式异常: {file_path}")
        continue

    # 第一次读取时保存基因名索引
    if gene_index is None:
        gene_index = df["gene_name"]
    
    expr_dict[sample_id] = df["unstranded"].values
    read_count += 1

print(f"完成 成功读取 {read_count} 个样本，共缺失 {len(missing_files)} 个")

if missing_files:
    print("部分缺失文件示例：")
    for f in missing_files[:5]:
        print("  -", f)

# =====================
# 生成表达矩阵 DataFrame
# =====================
expr_df = pd.DataFrame(expr_dict, index=gene_index)
expr_df.index.name = "gene_name"

# =====================
# 保存表达矩阵
# =====================
out_matrix = Path(r"E:\Bioinfo\Projects\Prac2. Differential cancer gene expression analysis\TCGA-LUAD\LUAD_counts_matrix.csv")
expr_df.to_csv(out_matrix)
print(f"完成 表达矩阵已保存至: {out_matrix}")

# =====================
# 输出样本信息文件
# =====================
sample_info = meta_df[["sample_id", "sample_type", "patient_id"]]
out_info = Path(r"E:\Bioinfo\Projects\Prac2. Differential cancer gene expression analysis\TCGA-LUAD\LUAD_sample_info.csv")
sample_info.to_csv(out_info, index=False)
print(f"完成 样本信息已保存至: {out_info}")

print("恭喜 整合完成！矩阵维度：", expr_df.shape)
