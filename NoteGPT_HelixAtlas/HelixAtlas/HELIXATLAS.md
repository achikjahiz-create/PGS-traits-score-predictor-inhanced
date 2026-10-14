# HelixAtlas

Run `sh start.sh` in an online IDE and forward port 8501. Docker is optional.

Use `python download_public.py --chromosomes 22` to download a real public VCF and its index. This command downloads files only; VCF normalization and automatic genotype-to-PGS analysis are not implemented. Do not mistake the synthetic demo for real analysis.

## Academic dataset broker

`datasets.py` adds a verified catalog of public ancient/modern resources with a
bounded, MD5-checked metadata downloader and an honest `"Group (N)"` sample
browser (counts parsed from disk, never invented).

```bash
python datasets.py catalog
python datasets.py download-metadata AADR_V66_1240K --out datasets
python datasets.py list AADR_V66_1240K --search yamnaya
python datasets.py list AADR_V66_1240K --era modern --min-n 50
python datasets.py backend AADR_V66_1240K
```

AADR comes from **Harvard Dataverse** (`doi:10.7910/DVN/FFIDCW`, v66.p1; the old
`reichdata.hms.harvard.edu` links now 404). Controlled datasets (EGA, NIAGADS)
are catalogued but never auto-downloaded. Multi-GB `.geno` matrices are left to
the pgsc_calc backend, not streamed into the web app. The sidebar expander
"Academic data catalog (metadata only)" exposes all of this in the UI.

The original README documents the scoring prototype. No authentication is included. Use private access for personal data. PDF exports contain tables, not charts. Ancient PGS is disabled. Requirements are compatibility ranges, not a tested lock.
