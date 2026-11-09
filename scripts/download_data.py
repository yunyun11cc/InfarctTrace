"""恢复并校验 GSE252027 源文件及分析所需的基因注释。"""
import concurrent.futures
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def download(row):
    # 文件完整通过历史哈希校验后才将临时下载移入正式位置。
    target = (ROOT / row['path']).resolve()
    if not target.is_relative_to(ROOT):
        raise ValueError(f'Download path escapes project: {row["path"]!r} -> {target}')
    target.parent.mkdir(parents=True, exist_ok=True)
    expected = row['sha256']
    if target.exists():
        actual = sha256(target)
        if actual != expected:
            raise ValueError(f'Existing file differs from manifest; preserved: {target}')
        state = 'existing_verified'
    else:
        partial = target.with_name(target.name + '.part')
        result = subprocess.run([
            'curl.exe', '--fail', '--location', '--silent', '--show-error',
            '--retry', '3', '--connect-timeout', '30', '--max-time', '1200',
            row['source_url'], '--output', str(partial),
        ], capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(f"{row['path']}: {result.stderr.strip()}")
        actual = sha256(partial)
        if partial.stat().st_size != int(row['bytes']) or actual != expected:
            raise ValueError(f'New download differs from historical checksum: {partial}')
        partial.replace(target)
        state = 'downloaded_verified'
    print(f"{state}: {row['path']} ({target.stat().st_size:,} bytes)", flush=True)
    return dict(path=row['path'], source_url=row['source_url'],
                bytes=target.stat().st_size, sha256=actual, status=state)


def main():
    with (ROOT / 'provenance/download_manifest.csv').open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.DictReader(stream))
    # 论文 PDF 等历史记录保留用于追溯，这里只恢复实验输入及相关说明。
    wanted = [r for r in rows if r['geo'] == 'GSE252027' and (
        '/original/' in r['path'] or r['path'].endswith('_family.soft.gz') or
        r['path'].endswith('AMI_revised_submission.zip') or
        (r['path'].endswith('.xml') and 'fulltext' in r['path'].lower()))]
    annotation_record = ROOT / 'provenance/annotation_download_2026-11-08.json'
    if annotation_record.exists():
        wanted.append(json.loads(annotation_record.read_text())['analysis_annotation'])
    results, failures = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        pending = {pool.submit(download, row): row for row in wanted}
        for future in concurrent.futures.as_completed(pending):
            try:
                results.append(future.result())
            except Exception as exc:
                failures.append(dict(path=pending[future]['path'], error=str(exc)))
                print(f'FAILED: {exc}', flush=True)
    report = dict(checked_at_utc=datetime.now(timezone.utc).isoformat(),
                  results=sorted(results, key=lambda r: r['path']), failures=failures)
    out = ROOT / 'outputs/01_preprocessing_differential/reports/download_verification.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    if failures:
        raise SystemExit(f'{len(failures)} downloads failed; see {out}')


if __name__ == '__main__':
    main()
