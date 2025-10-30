import pandas as pd

# 读取你的映射表
meta_df = pd.read_csv("E:\\Bioinfo\\Projects\\Prac2. Differential cancer gene expression analysis\\TCGA-LUAD\\file_sample_map.csv")

# 定义一个函数，用于提取样本类型
def get_sample_type(sample_id):
    """
    提取TCGA样本类型编码的前两位
    01 -> Tumor, 11 -> Normal, 其他 -> Other
    """
    try:
        code = sample_id.split('-')[3][:2]
        if code == '01':
            return 'Tumor'
        elif code == '11':
            return 'Normal'
        else:
            return 'Other'
    except:
        return 'Unknown'

# 提取病人编号（前3段），方便后续配对分析
meta_df["patient_id"] = meta_df["sample_id"].apply(lambda x: "-".join(x.split("-")[:3]))

# 添加样本类型列
meta_df["sample_type"] = meta_df["sample_id"].apply(get_sample_type)

# 保存
meta_df.to_csv("E:\\Bioinfo\\Projects\\Prac2. Differential cancer gene expression analysis\\TCGA-LUAD\\file_sample_map_with_type.csv", index=False)

print(meta_df.head())
