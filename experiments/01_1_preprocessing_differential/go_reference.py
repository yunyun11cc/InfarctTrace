"""从 Bioconductor 发布包读取小鼠基因和 GO 注释数据库。

``org.Mm.eg.db`` 提供 Ensembl ↔ Entrez 映射与基因信息，``GO.db`` 提供 GO 条目
名称。两者都以 SQLite 文件的形式打包在源码包 ``inst/extdata/`` 里。

本模块只负责“取到并留痕”：下载指定版本的 tar.gz、校验归档哈希、解出 SQLite
文件，并把版本、URL 与 SHA-256 写进 ``reports/go_sources.json``。分析阶段
（``enrichment.py``）以只读方式打开这些 SQLite 文件，不会再次联网。
"""

import hashlib
import json
import shutil
import tarfile
from urllib.request import urlopen

from paths import ROOT, REPORTS


# Bioconductor 3.22 的注释源码包地址；版本固定，保证结果可复现
SOURCE = 'https://bioconductor.org/packages/3.22/data/annotation/src/contrib/'

# 包名 → 包内 SQLite 文件名
PACKAGES = {'org.Mm.eg.db': 'org.Mm.eg.sqlite', 'GO.db': 'GO.sqlite'}
VERSION = '3.22.0'

# 下载与解压缓存目录，位于项目根的 .tools/
CACHE = ROOT / '.tools'


# ---------------------------------------------------------------------------
# 1. 校验工具
# ---------------------------------------------------------------------------

def sha256(path):
    """返回文件的 SHA-256 十六进制摘要，用于记录注释快照。"""
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


# ---------------------------------------------------------------------------
# 2. 入口
# ---------------------------------------------------------------------------

def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    records = {}
    for package, database in PACKAGES.items():
        archive = CACHE / f'{package}_{VERSION}.tar.gz'
        url = SOURCE + archive.name

        # 已有归档就复用；下载时先写 .part，成功后再原子替换，避免半截文件
        if not archive.exists():
            temporary = archive.with_suffix('.gz.part')
            with urlopen(url, timeout=120) as response, temporary.open('wb') as target:
                shutil.copyfileobj(response, target)
            temporary.replace(archive)

        output = CACHE / database
        temporary = output.with_suffix('.sqlite.part')

        # 只提取已知的 SQLite 文件，分析时直接读取注释表。
        # 注意：解压路径写死在包内结构里，避免整包解压到工作目录。
        with tarfile.open(archive, mode='r:gz') as bundle:
            name = f'{package}/inst/extdata/{database}'
            with bundle.extractfile(name) as source, temporary.open('wb') as target:
                shutil.copyfileobj(source, target)
        temporary.replace(output)

        # 归档与解出的 SQLite 都记哈希：归档证明下载内容，SQLite 证明实际读到的表
        records[package] = {'version': VERSION, 'url': url, 'archive_sha256': sha256(archive),
                            'sqlite_file': database, 'sqlite_sha256': sha256(output)}

    (REPORTS / 'go_sources.json').write_text(json.dumps(records, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print('小鼠基因与 GO 注释数据库已准备。')


if __name__ == '__main__':
    main()
