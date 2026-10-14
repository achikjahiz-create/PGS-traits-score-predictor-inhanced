# Genome Research Studio

Browser-based research prototype for normalized SNP data. Not clinically validated; generated code and tests have not been executed in this environment.

## Downloaded files

| Download | Save locally as | Purpose |
|---|---|---|
| engine.txt | engine.py | Genotype validation, exact-allele scoring, quality flags |
| Latest app text file | app.py | Streamlit interface, downloads, charts, exports |
| Test text file | test_engine.py | Regression tests |
| README.md | README.md | Instructions |

Keep these files together. No guardrails.py file is required: the app implements its warnings directly. Do not use any hard-coded trait IDs or labels from an ancillary file; published labels come from the downloaded scoring-file metadata.

## Install and run

Requires Python 3.11 or later. Use a private deployment for real human genotypes.

```bash
python -m venv .venv
source .venv/bin/activate
# Windows activation: .venv\Scripts\activate
pip install 'streamlit>=1.40,<2' 'pandas>=2.2,<3' 'numpy>=1.26,<3' 'plotly>=5.24,<7' 'requests>=2.32,<3' 'reportlab>=4.2,<5' 'pytest>=8,<10'
pytest -q
streamlit run app.py
```

These are compatibility ranges, not a tested dependency lock. Run the tests and smoke-test the UI; then freeze the tested environment before deployment.

For a private container, add this Dockerfile in the same folder:

```dockerfile
FROM python:3.11-slim
WORKDIR /app
COPY . .
RUN pip install --no-cache-dir 'streamlit>=1.40,<2' 'pandas>=2.2,<3' 'numpy>=1.26,<3' 'plotly>=5.24,<7' 'requests>=2.32,<3' 'reportlab>=4.2,<5' 'pytest>=8,<10'
EXPOSE 8501
CMD ["streamlit", "run", "app.py", "--server.address=0.0.0.0", "--server.port=8501"]
```

```bash
docker build -t genome-studio .
docker run --rm -p 8501:8501 genome-studio
```

Open http://localhost:8501. For internet access, deploy the container behind HTTPS and authenticated access. This delivery is source code, not a running hosted service.

## Input formats

Metadata CSV:

```csv
sample_id,dataset,group,era,build,ploidy,call_rate,contamination
EXAMPLE_1,study_1,cohort_A,modern,GRCh38,2,0.99,0.001
```

Genotypes CSV (one row per sample/site):

```csv
sample_id,chrom,pos,ref,alt,dosage,ploidy
EXAMPLE_1,1,1001,A,G,1,2
```

The illustrative quality values and variant are not real measurements or trait markers. Use measured quality metrics for real data. ALT dosage is in [0,ploidy]; blank dosage means missing. REF/ALT must be aligned to the declared assembly. Only biallelic autosomal A/C/G/T SNPs are accepted. Sample IDs must be globally unique. Mixed-ploidy individuals are unsupported.

Metadata include dataset, cohort label, ancient/modern era, genome build, measured call rate, and measured contamination. For modern data where contamination is unknown, obtain a defensible QC estimate upstream; do not invent zero. Rename/cohort-edit source metadata before uploading. In-app edits affect inclusion and thresholds only.

## Scoring policy

Only additive beta-weighted models are supported. Review model documentation: a beta header alone does not establish that a score lacks interactions, centering, scaling, or other special rules.

For effect-allele dosage d and published weight w:

$$PGS_i^{partial}=\sum_{j\in S}w_jd_{ij}$$

S contains non-palindromic exact-allele sites observed in every selected sample. REF-effect dosage is 2 minus ALT dosage. No strand complements, liftover, imputation, LD adjustment, or missingness rescaling are performed. Duplicate sites, build mismatches, unsupported weights, and inadequate coverage stop scoring.

The coverage denominator is eligible weights after unsupported variants are excluded; rejected counts are reported separately. This can still be a substantially partial model. Changing sample selection can change S; do not compare runs without verifying identical variant sets. Raw scores are not calibrated phenotype estimates, probabilities, or reference percentiles.

## Quality and averages

Absolute call-rate and contamination thresholds are editable. Directional MAD quality flags use technical/cohort strata of at least ten samples; zero-MAD strata are not flagged. No trait-score outlier exclusion is performed. Relatedness, batch effects, damage, authenticity and contamination must be assessed upstream. The app summarizes scored individuals; it does not create synthetic average genomes. Mean SE assumes independent samples and is not prediction uncertainty.

## Scope and limitations

| Capability | Current status |
|---|---|
| Modern normalized SNP scoring | Implemented for eligible additive models |
| Ancient sample browsing/filtering | Implemented; ancient scoring disabled |
| Online score weights | PGS Catalog harmonized files, bounded downloads |
| Online individual genotype datasets | Not automatically imported; provide authorized normalized files |
| Height | Raw published-score computation, not centimetres |
| Educational attainment | Research raw score only; not IQ or individual school prediction |
| Skin color | Raw compatible score only; no pigmentation categories |
| Hair color | Raw compatible score only; no color categories |
| Eye color | Raw compatible score only; no color categories |
| Additional traits | Compatible reviewed models, one score per run |
| Sample selection and cohort averages | Implemented |
| Interactive charts and tables | Implemented |
| PDF | Tables/provenance; charts not embedded; 2,000-sample limit |
| CSV and manifest export | Implemented |
| Native VCF/PLINK/EIGENSTRAT | Requires upstream conversion/backend |
| Unlimited dataset scale | Not supported; 128 MiB per input, memory-bound processing |

## Public dataset broker (`datasets.py`)

`datasets.py` is the academic bridge between the app and public ancient/modern
genotype resources. It holds a **verified catalog**, downloads only the small
index/annotation files (with MD5 + size checks, atomic writes and a byte cap),
and lists samples as `"Group (N)"` from data actually parsed off disk. It never
fabricates counts, never fetches controlled-access data, and never streams
multi-GB genotype matrices into the web app.

| Key | Resource | Access | Build |
|---|---|---|---|
| `AADR_V66_1240K` | Allen Ancient DNA Resource v66.p1, 1240K capture | public | GRCh37 |
| `AADR_V66_HO` | AADR v66.p1 on the Human Origins (~597K) array | public | GRCh37 |
| `AADR_ZARR_HF` | AADR **v62.0** Zarr mirror (third-party, Hugging Face) | public | GRCh37 |
| `KG_1000G_HG38` | 1000 Genomes, hg38 normalized PLINK 2 (Zenodo 15596480) | public | GRCh38 |
| `HUMAN_ORIGINS_ANGOLA_EGA` | Human Origins array, 208 Angolan individuals (EGA) | controlled | GRCh37 |
| `LASI_DAD_NIAGADS` | LASI-DAD WGS, 2,679 QC-passed (NIAGADS DSS NG00067) | controlled | GRCh38 |

Verified facts and corrected claims (checked 2026-10-14):

- The AADR is distributed via **Harvard Dataverse** (`doi:10.7910/DVN/FFIDCW`,
  release v66.0/v66.p1). The older `reichdata.hms.harvard.edu/.../V54.1.p1/`
  directory links **now return HTTP 404** and are not used here.
- AADR v66.p1 1240K contains **23,089 individuals / 3,897 group labels**; the
  released `.geno` is a **patch** file to be merged per the bundled README.
- AADR targets **hg19 (GRCh37)**; liftover is required to score GRCh38 weights.
- The Hugging Face Zarr mirror is a **third-party conversion of v62.0** (older,
  ~17,629 individuals), not an official Reich Lab release.
- The Zenodo 1000G file is a **third-party preprocessing** (v1 = record
  15596480, v2 = 16924674); its `.psam` carries no population column, so labels
  are joined from the official 1000G 3202-sample population file.
- **LASI-DAD is not on Harvard Dataverse and is not public**: controlled access
  is via **NIAGADS DSS `NG00067`** (sample set `snd10033`), related to dbGaP
  `phs000572` — not EGA.
- `EGAD00010002458` is the **Angolan Namib** Human Origins dataset (208
  individuals), controlled, non-commercial use.

Command line:

```bash
python datasets.py catalog
python datasets.py download-metadata AADR_V66_1240K --out datasets
python datasets.py list AADR_V66_1240K --search yamnaya
python datasets.py list AADR_V66_1240K --era modern --min-n 50
python datasets.py backend AADR_V66_1240K     # how to process the .geno
```

In-app, the sidebar expander **“Academic data catalog (metadata only)”** shows
the catalog, lists samples for whatever metadata is on disk, and offers the
bounded metadata download. Ancient/haploid samples remain **browsable but
unscorable** by design; `datasets.to_metadata_frame` maps AADR samples onto the
engine schema for QC review only.


## Scientific and deployment safeguards

Ancient samples and haploid genotypes are refused for scoring. Ancient diploid data also require temporal, ancestry, coverage and damage validation absent here. Educational-attainment scores must not be interpreted as intelligence, innate ability, predicted schooling, or evidence of population superiority. Cross-ancestry portability and confounding require study-specific evaluation.

Validated hair/eye/skin categorical prediction needs separately licensed and implemented models such as HIrisPlex-S, including required markers, coefficients and calibration; these are not fabricated here.

Public genotype releases have usage conditions; controlled biobanks need authorization. Use public ancient datasets responsibly and respect relevant community considerations.

This prototype has no built-in authentication, encryption-at-rest, persistent job queue, automatic retention controls, or multi-user authorization. Uploads are processed in server memory; hosting logs and infrastructure policies still matter. Do not upload real DNA to an untrusted public instance. For large-scale production, use private object storage, chunked validated ingestion, worker jobs, a validated scoring backend, access controls, and independently reviewed reference-cohort calibration.

## Sources checked 2026-10-14

| Source | Link |
|---|---|
| PGS Catalog formats | https://www.pgscatalog.org/downloads/ |
| pgsc_calc workflow | https://pgsc-calc.readthedocs.io/en/latest/ |
| AADR resource description | https://www.nature.com/articles/s41597-024-03031-7 |
| HIrisPlex-S | https://hirisplex.erasmusmc.nl/ |
| Educational-attainment prediction review | https://eric.ed.gov/?ff1=subHeredity&id=EJ1438186&q=descriptor%3A%22Heredity%22 |
