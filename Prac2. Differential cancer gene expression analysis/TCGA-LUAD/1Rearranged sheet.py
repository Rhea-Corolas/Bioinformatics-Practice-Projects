import json
import pandas as pd
from pathlib import Path

# 读取 metadata.json
# 映射表制作
with open('E:\\Bioinfo\\Projects\\Prac2. Differential cancer gene expression analysis\\TCGA-LUAD\\metadata.cart.2025-10-30.json') as f:
    metadata = json.load(f)

records = []
for item in metadata:
    file_name = item['file_name']
    sample_id = item['associated_entities'][0]['entity_submitter_id']
    records.append({'file': file_name, 'sample_id': sample_id})

meta_df = pd.DataFrame(records)
meta_df.to_csv('E:\\Bioinfo\\Projects\\Prac2. Differential cancer gene expression analysis\\TCGA-LUAD\\file_sample_map.csv', index=False)
meta_df.head()