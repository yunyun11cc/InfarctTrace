import pandas as pd
import os

# 1. 路径设置
DATASET_DIR = r"E:\project\mi-dna-choice\datasets\GSE252027"
COUNTS_FILE = os.path.join(DATASET_DIR, "FB_raw_counts_matrix.csv")
TPM_FILE = os.path.join(DATASET_DIR, "FB_tpm_matrix.csv")
METADATA_FILE = r"E:\project\mi-dna-choice\datasets\GSE252027\metadata\samples_FB.csv" # 注意这里用了你指定的路径

print("正在读取表达矩阵和元数据...")

# 2. 读取数据
# 因为保存时用了 index=True，所以第一列（基因名）作为 index_col=0 读入
counts_df = pd.read_csv(COUNTS_FILE, index_col=0)
tpm_df = pd.read_csv(TPM_FILE, index_col=0)
metadata_df = pd.read_csv(METADATA_FILE)

# 3. 确保元数据和矩阵的样本顺序完全一致
# 提取表达矩阵里的所有 GSM 编号（即列名）
matrix_samples = counts_df.columns.tolist()
# 提取元数据里的 GSM 编号
meta_samples = metadata_df['gsm'].tolist()

# 检查是否有遗漏
missing_in_meta = [s for s in matrix_samples if s not in meta_samples]
if missing_in_meta:
    print(f"⚠️ 警告：以下样本在矩阵中存在，但在元数据中缺失: {missing_in_meta}")

# 将元数据按照表达矩阵的列名顺序进行重新排序（这是为了让它们完美对齐）
metadata_df = metadata_df.set_index('gsm')
metadata_df = metadata_df.loc[matrix_samples] 

# 4. 数据清洗：剔除你之前发现的那些“坏样本”（可选但强烈建议）
# 从之前的记录看，GSM7992335, GSM7992336 的 reads 极低，GSM7992334/2337 属于 D3。
# 你之前图4画了 22 个点，说明剔除了 2 个。
# 我们假设要剔除异常值，你可以通过设定列表来过滤：
outliers = ['GSM7992335', 'GSM7992336'] 
valid_samples = [s for s in matrix_samples if s not in outliers]

print(f"剔除 {len(outliers)} 个异常样本后，保留 {len(valid_samples)} 个样本。")

# 清理最终版本的数据
counts_final = counts_df[valid_samples]
tpm_final = tpm_df[valid_samples]
metadata_final = metadata_df.loc[valid_samples]

# 5. 为 PyDESeq2 准备数据
# PyDESeq2 要求输入是：行是样本，列是基因！
# 所以我们要把它们转置 (.T)
counts_for_deseq = counts_final.T
counts_for_deseq.index.name = 'SampleID'

# 6. 保存整合后的干净数据
counts_out = os.path.join(DATASET_DIR, "FB_counts_clean_for_deseq.csv")
tpm_out = os.path.join(DATASET_DIR, "FB_tpm_clean_for_ml.csv")
meta_out = os.path.join(DATASET_DIR, "FB_metadata_clean.csv")

counts_for_deseq.to_csv(counts_out)
tpm_final.to_csv(tpm_out)
metadata_final.to_csv(meta_out)

print("=" * 50)
print("✅ 数据整合完成！")
print(f"1. 供 PyDESeq2 输入的 Counts 矩阵（行=样本, 列=基因）已保存: {counts_out}")
print(f"2. 供机器学习使用的高清 TPM 矩阵（行=基因, 列=样本）已保存: {tpm_out}")
print(f"3. 完美的元数据（包含处理时间、分组等）已保存: {meta_out}")
print("\n元数据预览：")
print(metadata_final[['title', 'day', 'cell_type_abbreviation']].head())