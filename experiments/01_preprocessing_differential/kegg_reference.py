"""获取 KEGG 小鼠通路名称与基因映射，记录分析使用的快照。

KEGG 是持续更新的数据库，同一 URL 在不同时间可能返回不同内容。为了让
富集结果可复现，本模块在首次运行时把三个接口的原始文件保存到
``datasets/annotation/kegg_mmu/``，并把 URL、SHA-256 与 release 信息写进
``manifest.json``；再次运行只校验哈希，不再重新下载。
"""

import hashlib
import json
from datetime import datetime, timezone
from urllib.request import urlopen

from paths import ROOT, REPORTS


# 三个 REST 端点：通路名称、通路成员、数据库版本
SOURCES = {
    'pathway_names.tsv': 'https://rest.kegg.jp/list/pathway/mmu',
    'pathway_genes.tsv': 'https://rest.kegg.jp/link/pathway/mmu',
    'release.txt': 'https://rest.kegg.jp/info/pathway',
}


def main():
    directory = ROOT / 'datasets/annotation/kegg_mmu'
    directory.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    manifest = directory / 'manifest.json'

    # 已存在快照时只做完整性校验：任一文件内容变化即报错，而不是静默换版
    if manifest.exists():
        records = json.loads(manifest.read_text(encoding='utf-8'))
        for name, url in SOURCES.items():
            actual = hashlib.sha256((directory / name).read_bytes()).hexdigest()
            if records['files'][name] != {'url': url, 'sha256': actual}:
                raise ValueError(f'KEGG 快照校验失败：{name}')

    # 首次运行：逐个下载、做基本内容检查，再记录哈希与 release
    else:
        records = {'downloaded_at_utc': datetime.now(timezone.utc).isoformat(), 'files': {}}
        for name, url in SOURCES.items():
            with urlopen(url, timeout=45) as response:
                raw = response.read()

            # pathway 接口应返回以 pathway 开头的文本，mmu 接口应返回 mmu 编号
            expected = b'pathway' if name == 'release.txt' else b'mmu'
            if len(raw) < 100 or not raw.startswith(expected):
                raise ValueError(f'KEGG 返回内容与预期不符：{name}')

            (directory / name).write_bytes(raw)
            records['files'][name] = {'url': url, 'sha256': hashlib.sha256(raw).hexdigest()}

        # release.txt 记录了本次快照对应的 KEGG 版本号
        records['release'] = (directory / 'release.txt').read_text(encoding='utf-8').strip()
        manifest.write_text(json.dumps(records, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

    # 把同一份记录复制到阶段报告，输出目录即可独立追溯注释来源
    (REPORTS / 'kegg_sources.json').write_text(json.dumps(records, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print('KEGG 小鼠通路快照已校验。')


if __name__ == '__main__':
    main()
