"""获取 KEGG 小鼠通路名称与基因映射，记录分析使用的快照。"""
import hashlib
import json
from datetime import datetime, timezone
from urllib.request import urlopen

from paths import ROOT, REPORTS

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
    if manifest.exists():
        records = json.loads(manifest.read_text(encoding='utf-8'))
        for name, url in SOURCES.items():
            actual = hashlib.sha256((directory / name).read_bytes()).hexdigest()
            if records['files'][name] != {'url': url, 'sha256': actual}:
                raise ValueError(f'KEGG 快照校验失败：{name}')
    else:
        records = {'downloaded_at_utc': datetime.now(timezone.utc).isoformat(), 'files': {}}
        for name, url in SOURCES.items():
            with urlopen(url, timeout=45) as response:
                raw = response.read()
            expected = b'pathway' if name == 'release.txt' else b'mmu'
            if len(raw) < 100 or not raw.startswith(expected):
                raise ValueError(f'KEGG 返回内容与预期不符：{name}')
            (directory / name).write_bytes(raw)
            records['files'][name] = {'url': url, 'sha256': hashlib.sha256(raw).hexdigest()}
        records['release'] = (directory / 'release.txt').read_text(encoding='utf-8').strip()
        manifest.write_text(json.dumps(records, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    (REPORTS / 'kegg_sources.json').write_text(json.dumps(records, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print('KEGG 小鼠通路快照已校验。')


if __name__ == '__main__':
    main()
