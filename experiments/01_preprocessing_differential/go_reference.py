"""从 Bioconductor 发布包读取小鼠基因和 GO 注释数据库。"""
import hashlib
import json
import shutil
import tarfile
from urllib.request import urlopen

from paths import ROOT, REPORTS

SOURCE = 'https://bioconductor.org/packages/3.22/data/annotation/src/contrib/'
PACKAGES = {'org.Mm.eg.db': 'org.Mm.eg.sqlite', 'GO.db': 'GO.sqlite'}
VERSION = '3.22.0'
CACHE = ROOT / '.tools'


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    CACHE.mkdir(parents=True, exist_ok=True) 
    REPORTS.mkdir(parents=True, exist_ok=True)
    records = {}
    for package, database in PACKAGES.items():
        archive = CACHE / f'{package}_{VERSION}.tar.gz'
        url = SOURCE + archive.name
        if not archive.exists():
            temporary = archive.with_suffix('.gz.part')
            with urlopen(url, timeout=120) as response, temporary.open('wb') as target:
                shutil.copyfileobj(response, target)
            temporary.replace(archive)
        output = CACHE / database
        temporary = output.with_suffix('.sqlite.part')
        # 只提取已知的 SQLite 文件，分析时直接读取注释表。
        with tarfile.open(archive, mode='r:gz') as bundle:
            name = f'{package}/inst/extdata/{database}'
            with bundle.extractfile(name) as source, temporary.open('wb') as target:
                shutil.copyfileobj(source, target)
        temporary.replace(output)
        records[package] = {'version': VERSION, 'url': url, 'archive_sha256': sha256(archive),
                            'sqlite_file': database, 'sqlite_sha256': sha256(output)}
    (REPORTS / 'go_sources.json').write_text(json.dumps(records, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print('小鼠基因与 GO 注释数据库已准备。')


if __name__ == '__main__':
    main()
