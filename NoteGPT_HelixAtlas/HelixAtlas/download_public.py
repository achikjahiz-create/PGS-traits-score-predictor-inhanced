import argparse, hashlib, shutil
from pathlib import Path
import requests
BASE = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
def fetch(name, directory):
    target = directory / name
    temporary = directory / (name + '.part')
    if target.exists():
        print('Already present:', target)
        return
    total = 0
    h = hashlib.sha256()
    try:
        with requests.get(BASE + name, stream=True, timeout=(15, 120), allow_redirects=False) as r:
            r.raise_for_status()
            if r.is_redirect: raise RuntimeError('Unexpected redirect')
            expected = int(r.headers.get('Content-Length', 0))
            if expected > 3 * 1024**3: raise RuntimeError('File too large')
            with temporary.open('wb') as f:
                for block in r.iter_content(1024**2):
                    if not block: continue
                    total += len(block)
                    if total > 3 * 1024**3: raise RuntimeError('Size limit exceeded')
                    if shutil.disk_usage(directory).free < 1024**3: raise RuntimeError('Low disk space')
                    f.write(block); h.update(block)
            if expected and total != expected: raise RuntimeError('Truncated download')
            if total == 0: raise RuntimeError('Empty download')
        temporary.replace(target)
        target.with_name(name + '.sha256').write_text(h.hexdigest())
        print('Downloaded:', target, 'bytes:', total)
    finally:
        temporary.unlink(missing_ok=True)
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--chromosomes', type=int, nargs='+', required=True)
    p.add_argument('--directory', default='datasets')
    args = p.parse_args()
    if any(c not in range(1,23) for c in args.chromosomes): p.error('Only chromosomes 1-22')
    directory = Path(args.directory); directory.mkdir(parents=True, exist_ok=True)
    for chrom in sorted(set(args.chromosomes)):
        name = f'ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'
        fetch(name, directory); fetch(name + '.tbi', directory)
