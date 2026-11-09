import os
import tarfile
import pandas as pd
import io
import gzip # 新增：专门处理 .gz 解压

# 1. 设置路径
TAR_FILE_PATH = r"E:\project\mi-dna-choice\datasets\GSE252027\original\GSE252027_RAW.tar"
OUTPUT_DIR = r"E:\project\mi-dna-choice\datasets\GSE252027"

print(f"正在打开 Tar 包: {TAR_FILE_PATH}")

counts_list = []
tpm_list = []

with tarfile.open(TAR_FILE_PATH, "r:*") as tar:
    members = tar.getmembers()
    print(f"Tar 包中共有 {len(members)} 个文件。")
    
    for member in members:
        # 关键过滤：只处理 FB (成纤维细胞) 的 .sf.gz 文件
        if member.name.endswith('.sf.gz') and '_FB_' in member.name:
            filename = os.path.basename(member.name)
            gsm_id = filename.split('_')[0]
            
            # 从 tar 包中提取该文件的内存流
            f = tar.extractfile(member)
            if f is not None:
                # 【核心修改】：先手动用 gzip 解压二进制流，然后再用 pandas 读取
                try:
                    with gzip.open(f, 'rt') as gz_file:
                        df = pd.read_csv(gz_file, sep='\t', usecols=['Name', 'NumReads', 'TPM'])
                        
                        df = df.set_index('Name')
                        counts_list.append(df['NumReads'].rename(gsm_id))
                        tpm_list.append(df['TPM'].rename(gsm_id))
                        print(f"✅ 成功读取: {gsm_id}")
                except Exception as e:
                    print(f"❌ 读取 {gsm_id} 失败: {e}")

# 3. 合并矩阵
print(f"\n成功读取 {len(counts_list)} 个 FB 样本，正在合并矩阵...")
counts_matrix = pd.concat(counts_list, axis=1)
tpm_matrix = pd.concat(tpm_list, axis=1)

# 4. 数据清洗
counts_matrix = counts_matrix.fillna(0).astype(int)
tpm_matrix = tpm_matrix.fillna(0)

# 5. 保存
os.makedirs(OUTPUT_DIR, exist_ok=True)
counts_out_path = os.path.join(OUTPUT_DIR, "FB_raw_counts_matrix.csv")
tpm_out_path = os.path.join(OUTPUT_DIR, "FB_tpm_matrix.csv")

print("正在保存大矩阵文件...")
counts_matrix.to_csv(counts_out_path, index=True)
tpm_matrix.to_csv(tpm_out_path, index=True)

print("=" * 50)
print("✅ 处理完成！")
print(f"Raw Counts 矩阵形状 (基因数, 样本数): {counts_matrix.shape}")
print(f"TPM 矩阵形状 (基因数, 样本数): {tpm_matrix.shape}")
print(f"文件已保存至: {OUTPUT_DIR}")