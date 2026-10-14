"""HelixAtlas academic dataset broker.

A verified catalog of public ancient and modern human genotype resources, plus
an honest sample browser and a *bounded, checksum-verified* downloader that
fetches the metadata (SNP index / individual index / annotation) needed to
describe samples and to intersect them with PGS Catalog scoring files.

Design principles
-----------------
1. No fabricated numbers. Every ``"Group (N)"`` count is produced by parsing a
   file that was actually downloaded. There are no hard-coded demo counts.
2. No silent network access. Downloads are explicit, bounded by a byte cap,
   written atomically, and verified against a published MD5 where one exists.
3. No unlicensed harvesting. Controlled-access archives (EGA, NIAGADS DSS,
   dbGaP) are catalogued for discovery but never downloaded automatically; the
   module raises :class:`AccessRestricted` and prints the access route.
4. Genotype matrices stay out of the web app. AADR ``.geno`` files are
   multi-GB; this module hands them to an external backend (pgsc_calc /
   ADMIXTOOLS 2 / PLINK 2) instead of streaming them into Streamlit.
5. Traits are NOT a property of genotype resources. A genotype resource tells
   you *who* was sequenced; which phenotypes are predictable comes from the
   PGS Catalog scoring file, not from the AADR ``.anno``. The former
   ``traits_supported`` / ``trait_filter`` promise is therefore removed.

Catalog provenence
------------------
Facts below were checked against the live sources on 2026-10-14:

* AADR - Harvard Dataverse ``doi:10.7910/DVN/FFIDCW`` (version 14.0, release
  v66.0 / v66.p1). The legacy ``reichdata.hms.harvard.edu/.../V54.1.p1/``
  paths now return HTTP 404 and are NOT used.
* AADR Zarr mirror - Hugging Face ``gebregl/allen_ancient_dna`` (third-party
  conversion, AADR **v62.0**, hg19, CC0).
* 1000 Genomes hg38 normalized - Zenodo record ``15596480`` (v1) /
  ``16924674`` (v2), author David Laub, CC BY 4.0.
* Human Origins, Angola - EGA ``EGAD00010002458`` (208 individuals, 209
  samples; study ``EGAS00001007011``; controlled, non-commercial).
* LASI-DAD - NIAGADS DSS ``NG00067`` (ADSP Umbrella, sample set ``snd10033``),
  2,686 sequenced / 2,679 QC-passed; controlled. Related dbGaP ``phs000572``.

References: AADR paper https://www.nature.com/articles/s41597-024-03031-7 ;
pgsc_calc https://pgsc-calc.readthedocs.io/ ; PGS Catalog
https://www.pgscatalog.org/ .
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import sys
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #

#: Default ceiling for a single metadata file fetched inside the web app.
MAX_METADATA_BYTES = 256 * 1024**2

#: Hard ceiling even when a caller raises the cap, to stop accidents.
ABSOLUTE_BYTES = 4 * 1024**3

#: Hosts this module is willing to talk to (guards against URL tampering).
ALLOWED_HOSTS = {
    'dataverse.harvard.edu',
    'huggingface.co',
    'cdn-lfs.huggingface.co',
    'cdn-lfs-us-1.hf.co',
    'zenodo.org',
    'ftp.ebi.ac.uk',
    'ftp.1000genomes.ebi.ac.uk',
}

#: AADR 1240K Eigenstrat ``.snp`` column layout (no header row).
AADR_SNP_COLUMNS = ['snp_id', 'chrom', 'genetic_pos', 'pos', 'ref', 'alt']


class DatasetError(RuntimeError):
    """Base class for dataset-broker failures."""


class AccessRestricted(DatasetError):
    """Raised when a catalog entry is controlled-access and cannot be auto-fetched."""


class DatasetNotFetched(DatasetError):
    """Raised when sample listing needs a metadata file that is not on disk yet."""


# --------------------------------------------------------------------------- #
# Catalog
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class RemoteFile:
    """A single downloadable artifact with an optional integrity check."""

    filename: str
    url: str
    kind: str                       # anno | ind | snp | psam | geno | index | readme
    bytes: Optional[int] = None
    md5: Optional[str] = None
    metadata: bool = False          # safe/small enough for in-app download


@dataclass(frozen=True)
class Dataset:
    key: str
    name: str
    provider: str
    description: str
    landing_url: str
    access: str                     # public | controlled
    build: str
    license: str
    citation: str
    notes: tuple = ()
    files: tuple = field(default_factory=tuple)

    def file(self, kind: str) -> RemoteFile:
        for f in self.files:
            if f.kind == kind:
                return f
        raise DatasetError(f'{self.key}: no {kind!r} file is catalogued')


def _dataverse(datafile_id: int, filename: str, kind: str, size, md5, metadata):
    return RemoteFile(
        filename=filename,
        url=f'https://dataverse.harvard.edu/api/access/datafile/{datafile_id}',
        kind=kind,
        bytes=size,
        md5=md5,
        metadata=metadata,
    )


_AADR_NOTES = (
    'Distributed by Harvard Dataverse; the old reichdata.hms.harvard.edu '
    'V54/V62 directory links now 404.',
    'The released .geno is a "patch" file; follow the bundled README '
    '(aadr_v66.p1__README.docx) for the exact merge order before use.',
    'AADR targets hg19 (GRCh37). Liftover to GRCh38 is required to score '
    'GRCh38-harmonised PGS weights.',
    'AADR is largely pseudo-haploid; the app intentionally refuses to score '
    'ancient/haploid samples. Use a dedicated ancient-DNA model.',
)

DATASET_CATALOG: dict = {
    # ---------------------------------------------------------------- AADR --
    'AADR_V66_1240K': Dataset(
        key='AADR_V66_1240K',
        name='Allen Ancient DNA Resource (AADR) v66.p1 - 1240K',
        provider='David Reich Lab / Harvard University (Harvard Dataverse)',
        description=('Curated ancient + present-day human genotypes at the ~1.15M '
                     'site "1240K" capture set, with rich annotation (dates, '
                     'locations, coverage, QC). Latest release 2026-04-13.'),
        landing_url='https://doi.org/10.7910/DVN/FFIDCW',
        access='public',
        build='GRCh37',
        license='See Dataverse terms of use; cite AADR (doi:10.7910/DVN/FFIDCW).',
        citation=('Mallick et al. 2024, Scientific Data 11:182 '
                  '(https://doi.org/10.1038/s41597-024-03031-7).'),
        notes=_AADR_NOTES,
        files=(
            _dataverse(13994515, 'v66.p1_1240K.aadr.PUB.anno', 'anno', 13443726,
                       'a2db1ac16f0f3558ed66fb251e1d5c7d', True),
            _dataverse(13994513, 'v66.p1_1240K.aadr.patch.PUB.ind', 'ind', 1017460,
                       '19a434ac954bcd10dbb8dba1d1188a09', True),
            _dataverse(13994514, 'v66.p1_1240K.aadr.patch.PUB.snp', 'snp', 77679819,
                       '50f66178fc81b8aa087cc4b135317e59', True),
            _dataverse(13995422, 'v66.p1__files.md5sum', 'readme', 1769,
                       'ef9e9d78750808e28d33b4074098c5e8', True),
            _dataverse(13994829, 'v66.p1_1240K.aadr.patch.PUB.geno', 'geno', 7117276654,
                       '5ea1d2675a271c81e55b8f8b08b3ff3b', False),
        ),
    ),
    'AADR_V66_HO': Dataset(
        key='AADR_V66_HO',
        name='Allen Ancient DNA Resource (AADR) v66.p1 - Human Origins array',
        provider='David Reich Lab / Harvard University (Harvard Dataverse)',
        description=('Same sample collection typed on the Affymetrix Human '
                     'Origins (~597K site) array; useful where 1240K capture is '
                     'too sparse.'),
        landing_url='https://doi.org/10.7910/DVN/FFIDCW',
        access='public',
        build='GRCh37',
        license='See Dataverse terms of use; cite AADR.',
        citation='Mallick et al. 2024, Scientific Data 11:182.',
        notes=_AADR_NOTES,
        files=(
            _dataverse(13994528, 'v66.p1_HO.aadr.PUB.anno', 'anno', 15465906,
                       '548b011b4adbcd46f0ff66b957c59c3f', True),
            _dataverse(13994526, 'v66.p1_HO.aadr.patch.PUB.ind', 'ind', 1176253,
                       'a14985ffae6dd33a00f4fd5bfefe58c5', True),
            _dataverse(13994527, 'v66.p1_HO.aadr.patch.PUB.snp', 'snp', 36800253,
                       'fbd2cfe09ec6bad794f7607768b57a3a', True),
            _dataverse(13994808, 'v66.p1_HO.aadr.patch.PUB.geno', 'geno', 4029634650,
                       '9f9f0629f803835985bc38856e989995', False),
        ),
    ),
    'AADR_ZARR_HF': Dataset(
        key='AADR_ZARR_HF',
        name='AADR v62.0 (Zarr mirror)',
        provider='Third-party mirror on Hugging Face (not the Reich Lab)',
        description=('AADR converted from Eigenstrat to chunked Zarr, plus the '
                     'v62.0 .anno/.ind/.snp. Convenient for lazy/array loading; '
                     'this is an OLDER release (v62.0) than Dataverse v66.'),
        landing_url='https://huggingface.co/datasets/gebregl/allen_ancient_dna',
        access='public',
        build='GRCh37',
        license='CC0-1.0 mirror; underlying AADR reuse policy still applies.',
        citation=('Converted by the dataset maintainer; underlying data: '
                  'Mallick et al. 2024, Scientific Data 11:182.'),
        notes=(
            'Third-party conversion, hg19, ~17,629 individuals (v62.0).',
            'Zarr values: 0=hom-ref, 1=het, 2=hom-alt, 3=no read.',
            'Prefer the Dataverse v66 release for reproducibility.',
        ),
        files=(
            RemoteFile(
                filename='v62.0_1240k_public.anno',
                url='https://huggingface.co/datasets/gebregl/allen_ancient_dna/resolve/main/v62.0_1240k_public.anno',
                kind='anno', metadata=True),
            RemoteFile(
                filename='v62.0_1240k_public.ind',
                url='https://huggingface.co/datasets/gebregl/allen_ancient_dna/resolve/main/v62.0_1240k_public.ind',
                kind='ind', metadata=True),
            RemoteFile(
                filename='v62.0_1240k_public.snp',
                url='https://huggingface.co/datasets/gebregl/allen_ancient_dna/resolve/main/v62.0_1240k_public.snp',
                kind='snp', metadata=True),
        ),
    ),
    # ------------------------------------------------------- 1000 Genomes --
    'KG_1000G_HG38': Dataset(
        key='KG_1000G_HG38',
        name='1000 Genomes Project - hg38, phased, biallelic, normalized',
        provider='David Laub (Zenodo)',
        description=('1000G 2504 high-coverage samples, left-aligned/normalised '
                     'to hg38, split multiallelics, emitted as PLINK 2 '
                     '(pgen/pvar/psam). Good modern reference panel.'),
        landing_url='https://zenodo.org/records/15596480',
        access='public',
        build='GRCh38',
        license='CC BY 4.0',
        citation=('Laub, D. (2025). 1000 Genomes Project hg38, phased, biallelic, '
                  'normalized [Dataset]. Zenodo. doi:10.5281/zenodo.15596480'),
        notes=(
            'A newer v2 (record 16924674) exists; pin the exact record you cite.',
            'This is a third-party preprocessing of 1000G, not the EBI primary FTP.',
            'The Zenodo .psam carries only #IID/PAT/MAT/SEX (3202 samples, no '
            'labels); population labels come from the official 1000G population '
            'file catalogued here.',
        ),
        files=(
            RemoteFile(
                filename='1kGP.snp_indel.split_multiallelics.psam',
                url='https://zenodo.org/records/15596480/files/1kGP.snp_indel.split_multiallelics.psam?download=1',
                kind='psam', bytes=52105,
                md5='cb8ea444a41cceb27483df26bb76ad1b', metadata=True),
            RemoteFile(
                filename='20130606_g1k_3202_samples_ped_population.txt',
                url=('https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/'
                     '1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt'),
                kind='population', bytes=97483,
                md5='bac176a443f941890ac7b51fb1cabe0d', metadata=True),
            RemoteFile(
                filename='1kGP.snp_indel.split_multiallelics.pvar.zst',
                url='https://zenodo.org/records/15596480/files/1kGP.snp_indel.split_multiallelics.pvar.zst?download=1',
                kind='snp', bytes=2791728742,
                md5='cb00b271504713dbbdbbf17ee0b1825e', metadata=False),
            RemoteFile(
                filename='1kGP.snp_indel.split_multiallelics.pgen.xz',
                url='https://zenodo.org/records/15596480/files/1kGP.snp_indel.split_multiallelics.pgen.xz?download=1',
                kind='geno', bytes=3328599654,
                md5='3e07b63b48e5a205a2b2ae5022e01871', metadata=False),
        ),
    ),
    # --------------------------------------------------- Controlled access --
    'HUMAN_ORIGINS_ANGOLA_EGA': Dataset(
        key='HUMAN_ORIGINS_ANGOLA_EGA',
        name='Human Origins array - Angolan Namib (EGA)',
        provider='European Genome-phenome Archive (EGA)',
        description=('Affymetrix Axiom Human Origins array data for 208 Angolan '
                     'individuals (209 samples), study EGAS00001007011.'),
        landing_url='https://ega-archive.org/datasets/EGAD00010002458',
        access='controlled',
        build='GRCh37',
        license='Non-commercial use only; DAC-approved access required.',
        citation='EGA dataset EGAD00010002458 / study EGAS00001007011.',
        notes=('Requires a Data Access Committee request; cannot be auto-downloaded.',
               'PLINK binaries are delivered after approval, via pyega3.'),
        files=(),
    ),
    'LASI_DAD_NIAGADS': Dataset(
        key='LASI_DAD_NIAGADS',
        name='LASI-DAD whole-genome sequencing (NIAGADS DSS)',
        provider='NIAGADS Data Sharing Service (ADSP Umbrella NG00067)',
        description=('~30x WGS for 2,686 participants from 18 Indian states '
                     '(2,679 QC-passed; ~73.2M autosomal biallelic variants). '
                     'Sample set snd10033.'),
        landing_url='https://dss.niagads.org/datasets/ng00067/',
        access='controlled',
        build='GRCh38',
        license='Controlled access; IRB/data-use approval required.',
        citation='LASI-DAD WGS via NIAGADS DSS NG00067 (snd10033); dbGaP phs000572.',
        notes=('NOT distributed by the Reich Lab and NOT on Harvard Dataverse.',
               'Access is via NIAGADS DSS, not EGA; request individually.'),
        files=(),
    ),
}


def catalog_rows() -> list:
    """Return the catalog as a list of flat dicts (for a DataFrame/table)."""
    rows = []
    for ds in DATASET_CATALOG.values():
        meta = [f.filename for f in ds.files if f.metadata]
        big = [f.filename for f in ds.files if not f.metadata]
        rows.append({
            'key': ds.key,
            'name': ds.name,
            'provider': ds.provider,
            'access': ds.access,
            'build': ds.build,
            'landing_url': ds.landing_url,
            'metadata_files': ', '.join(meta) or '-',
            'large_files': ', '.join(big) or '-',
        })
    return rows


def get(dataset_key: str) -> Dataset:
    if dataset_key not in DATASET_CATALOG:
        raise DatasetError(
            f'Unknown dataset {dataset_key!r}. Known: {sorted(DATASET_CATALOG)}')
    return DATASET_CATALOG[dataset_key]


# --------------------------------------------------------------------------- #
# Download
# --------------------------------------------------------------------------- #


def _check_host(url: str) -> None:
    host = urllib.parse.urlsplit(url).hostname or ''
    if host not in ALLOWED_HOSTS:
        raise DatasetError(f'Refusing to fetch from untrusted host: {host!r}')


def download_file(url: str, dest, *, expected_md5: Optional[str] = None,
                  expected_bytes: Optional[int] = None,
                  max_bytes: int = MAX_METADATA_BYTES,
                  timeout=(15, 300), chunk: int = 1024**2) -> Path:
    """Stream ``url`` to ``dest`` atomically, then verify size and MD5.

    Nothing is left on disk unless the whole file arrived and passed its
    integrity check; partial downloads are removed on failure.
    """
    _check_host(url)
    if max_bytes > ABSOLUTE_BYTES:
        raise DatasetError('Requested byte cap exceeds the absolute ceiling')
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        size = dest.stat().st_size
        if (expected_bytes is None or size == expected_bytes) and (
                expected_md5 is None or _md5(dest) == expected_md5):
            return dest
        raise DatasetError(f'{dest} exists but does not match the manifest')

    part = dest.with_name(dest.name + '.part')
    digest = hashlib.md5()
    total = 0
    try:
        request = urllib.request.Request(url, headers={'User-Agent': 'HelixAtlas/1.0'})
        with urllib.request.urlopen(request, timeout=timeout[1]) as response:
            declared = int(response.headers.get('Content-Length') or 0)
            if declared and declared > max_bytes:
                raise DatasetError(
                    f'{dest.name}: {declared} bytes exceeds cap {max_bytes}')
            if expected_bytes and declared and declared != expected_bytes:
                raise DatasetError(
                    f'{dest.name}: server size {declared} != manifest {expected_bytes}')
            with part.open('wb') as handle:
                while True:
                    block = response.read(chunk)
                    if not block:
                        break
                    total += len(block)
                    if total > max_bytes:
                        raise DatasetError(f'{dest.name}: exceeded byte cap')
                    if shutil.disk_usage(dest.parent).free < 256 * 1024**2:
                        raise DatasetError('Low disk space; aborting download')
                    handle.write(block)
                    digest.update(block)
        if total == 0:
            raise DatasetError('Empty download')
        if expected_bytes and total != expected_bytes:
            raise DatasetError(f'{dest.name}: got {total} bytes, expected {expected_bytes}')
        if expected_md5 and digest.hexdigest() != expected_md5:
            raise DatasetError(f'{dest.name}: MD5 {digest.hexdigest()} != {expected_md5}')
        part.replace(dest)
        return dest
    finally:
        part.unlink(missing_ok=True)


def _md5(path: Path, chunk: int = 1024**2) -> str:
    h = hashlib.md5()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(chunk), b''):
            h.update(block)
    return h.hexdigest()


def download_metadata(dataset_key: str, out_dir='datasets',
                      kinds: Iterable[str] = ('anno', 'ind', 'snp', 'psam',
                                              'population', 'readme'),
                      max_bytes: int = MAX_METADATA_BYTES) -> dict:
    """Download only the small metadata files of a dataset.

    Large genotype matrices are refused here by design; fetch them with the
    external backend (see :func:`describe_backend`).
    """
    ds = get(dataset_key)
    if ds.access != 'public':
        raise AccessRestricted(
            f'{ds.key} is controlled-access. Request it at {ds.landing_url} '
            f'({"; ".join(ds.notes)}).')
    wanted = set(kinds)
    target = Path(out_dir) / ds.key
    fetched = {}
    for remote in ds.files:
        if remote.kind not in wanted or not remote.metadata:
            continue
        path = download_file(remote.url, target / remote.filename,
                             expected_md5=remote.md5,
                             expected_bytes=remote.bytes,
                             max_bytes=max_bytes)
        fetched[remote.kind] = path
    if not fetched:
        raise DatasetError(f'{ds.key}: no matching metadata files to download')
    return fetched


def describe_backend(dataset_key: str) -> str:
    """Return instructions for turning large genotype files into the app's input."""
    ds = get(dataset_key)
    if ds.access != 'public':
        return (f'{ds.key}: controlled access. Obtain approval at {ds.landing_url}, '
                f'then run the pgsc_calc backend on the delivered files.')
    big = [f for f in ds.files if not f.metadata]
    names = ', '.join(f.filename for f in big) or 'n/a'
    return (
        f'{ds.key}: large matrices [{names}] are intentionally NOT pulled into the '
        'web app. Use pgsc_calc (https://pgsc-calc.readthedocs.io/) or ADMIXTOOLS 2 '
        f'to convert/normalise them to hg38 and emit the CSV/TSV the app expects. '
        f'Assembly declared by this catalog entry: {ds.build}.')


# --------------------------------------------------------------------------- #
# Parsers - real counts, no fabrication
# --------------------------------------------------------------------------- #

# Positional indices into the AADR ``.anno`` (its header has duplicate names,
# so positional access is the only safe approach).
AADR_ANNO_INDEX = {
    'genetic_id': 0,
    'persistent_id': 1,
    'individual_id': 2,
    'date_mean_bp': 10,
    'group': 14,
    'locality': 15,
    'political_entity': 16,
    'latitude': 17,
    'longitude': 18,
    'data_type': 21,
    'mean_coverage': 23,
    'snps_2m': 25,
    'snps_1240k': 26,
    'snps_ho': 27,
    'molecular_sex': 30,
    'hapconx_ci': 43,
    'assessment': 47,
}


def load_aadr_anno(path) -> list:
    """Parse an AADR ``.anno`` file into a list of dicts with named fields."""
    with Path(path).open(encoding='utf-8', newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        if len(header) < 49:
            raise DatasetError(f'{path}: unexpected .anno column count {len(header)}')
        rows = []
        for raw in reader:
            if not raw or not raw[0].strip():
                continue
            rows.append({name: raw[i] if i < len(raw) else ''
                         for name, i in AADR_ANNO_INDEX.items()})
    if not rows:
        raise DatasetError(f'{path}: no individuals parsed')
    return rows


def load_aadr_ind(path) -> list:
    """Parse an Eigenstrat/AADR ``.ind`` file -> [{sample_id, sex, group}]."""
    records = []
    with Path(path).open(encoding='utf-8') as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 3:
                continue
            records.append({'sample_id': parts[0], 'sex': parts[1], 'group': parts[2]})
    return records


def load_1000g_psam(path) -> list:
    """Parse a PLINK 2 ``.psam`` -> [{sample_id, population, superpopulation, sex}]."""
    records = []
    with Path(path).open(encoding='utf-8') as handle:
        header = handle.readline().lstrip('#').split()
        for line in handle:
            parts = line.split()
            if len(parts) < len(header):
                continue
            records.append(dict(zip(header, parts)))
    return records


def _gnorm(value: float) -> float:
    """Guard against AADR's float-in-an-int-string quirk (e.g. '1.0')."""
    return float(value)


# --------------------------------------------------------------------------- #
# Sample browser
# --------------------------------------------------------------------------- #


def group_counts(records: Iterable[dict], group_key: str = 'group') -> list:
    """Return ``(group, n)`` pairs sorted by descending n then name."""
    counts: dict = {}
    for record in records:
        counts[record[group_key]] = counts.get(record[group_key], 0) + 1
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))


def list_samples(dataset_key: str, *, data_dir='datasets', group_contains=None,
                 era=None, min_n=1) -> list:
    """List samples grouped by population as ``'Group (N)'``.

    Counts come from a downloaded index/annotation file; if it is missing the
    function raises :class:`DatasetNotFetched` rather than inventing numbers.

    Parameters
    ----------
    data_dir       : where :func:`download_metadata` stored the files.
    group_contains : case-insensitive substring filter on the group label.
    era            : for AADR, ``'ancient'`` or ``'modern'`` (date-lived proxy).
    min_n          : drop groups smaller than this.
    """
    ds = get(dataset_key)
    if ds.access != 'public':
        raise AccessRestricted(
            f'{ds.key} is controlled-access; no sample index is available '
            f'without approval at {ds.landing_url}.')

    records = _load_records(ds, Path(data_dir) / ds.key, era=era)
    pairs = [(g, n) for g, n in group_counts(records) if n >= min_n]
    if group_contains:
        needle = group_contains.lower()
        pairs = [(g, n) for g, n in pairs if needle in g.lower()]
    return [f'{g} ({n})' for g, n in pairs]


def _load_records(ds: Dataset, folder: Path, *, era=None) -> list:
    if ds.key.startswith('AADR'):
        candidates = sorted(folder.glob('*.anno'))
        if not candidates:
            raise DatasetNotFetched(
                f'{ds.key}: no .anno in {folder}. Run '
                f'download_metadata({ds.key!r}) first.')
        records = load_aadr_anno(candidates[-1])
        if era in ('ancient', 'modern'):
            records = [r for r in records if _aadr_era(r) == era]
        return records
    if ds.key == 'KG_1000G_HG38':
        candidates = sorted(folder.glob('*.psam'))
        if not candidates:
            raise DatasetNotFetched(
                f'{ds.key}: no .psam in {folder}. Run '
                f"download_metadata('{ds.key}') first.")
        pop_file = next(iter(sorted(folder.glob('*population*.txt'))), None)
        if pop_file is None:
            raise DatasetNotFetched(
                f'{ds.key}: the .psam has no population column; download the '
                f'1000G population file too (download_metadata, kinds includes '
                f"'population').")
        labels = {}
        for line in Path(pop_file).read_text(encoding='utf-8').splitlines():
            parts = line.split()
            if len(parts) >= 6 and parts[1] != 'SampleID':
                labels[parts[1]] = {'group': parts[5], 'superpopulation': parts[6]}
        records = []
        for row in load_1000g_psam(candidates[-1]):
            sid = row.get('IID') or row.get('sample_id') or row.get('sample')
            if sid in labels:
                records.append({'sample_id': sid, **labels[sid]})
        if not records:
            raise DatasetError(f'{ds.key}: could not join .psam IDs to the panel')
        return records


def _aadr_era(record: dict) -> str:
    """Provisional ancient/modern proxy: modern rows carry date 0 or '..'.

    Documented heuristic, not an authority - verify against the AADR README.
    """
    raw = (record.get('date_mean_bp') or '').strip()
    if raw in ('', '..'):
        return 'ancient'          # undated ancient context
    try:
        return 'modern' if _gnorm(raw) <= 0 else 'ancient'
    except ValueError:
        return 'ancient'


def _aadr_ploidy(record: dict) -> int:
    """Provisional: '.DG'/'diploid' data type -> 2, otherwise pseudo-haploid 1."""
    dtype = (record.get('data_type') or '').lower()
    return 2 if 'diploid' in dtype else 1


def search_samples(dataset_key: str, query: str, **kwargs) -> list:
    """Substring search over the group listing produced by :func:`list_samples`."""
    return [s for s in list_samples(dataset_key, **kwargs) if query.lower() in s.lower()]


def sample_table(dataset_key: str, *, data_dir='datasets', era=None) -> list:
    """Return per-individual rows (id, group, era, sex) for the app's selector."""
    ds = get(dataset_key)
    if ds.access != 'public':
        raise AccessRestricted(f'{ds.key}: controlled-access; no table available.')
    records = _load_records(ds, Path(data_dir) / ds.key, era=era)
    if ds.key.startswith('AADR'):
        return [{'sample_id': r['genetic_id'], 'group': r['group'],
                 'era': _aadr_era(r), 'sex': r.get('molecular_sex', '')}
                for r in records]
    return [{'sample_id': r.get('IID') or r.get('sample_id') or r.get('sample'),
             'group': r['group'],
             'era': 'modern',
             'sex': r.get('SEX', r.get('sex', ''))} for r in records]


# --------------------------------------------------------------------------- #
# Bridge to the scoring engine's metadata schema
# --------------------------------------------------------------------------- #


def to_metadata_frame(dataset_key: str, *, data_dir='datasets', dataset_label=None,
                      snp_site_count: Optional[int] = None):
    """Map AADR samples onto the columns ``engine.validate_metadata`` expects.

    The returned frame is deliberately *provisional*: ``call_rate`` is the
    measured share of 1240K sites hit (only when ``snp_site_count`` is known),
    and ``contamination`` is the upper bound of the male X-based hapConX CI or
    NaN. The app forbids inventing quality values, and it refuses to score
    ancient/haploid samples, so this frame is for browsing and QC review.

    Returns a ``pandas.DataFrame``.
    """
    import pandas as pd

    ds = get(dataset_key)
    folder = Path(data_dir) / ds.key
    if not ds.key.startswith('AADR'):
        raise DatasetError('to_metadata_frame currently supports AADR entries only.')
    candidates = sorted(folder.glob('*.anno'))
    if not candidates:
        raise DatasetNotFetched(
            f'{ds.key}: no .anno in {folder}. Run download_metadata first.')
    if snp_site_count is None:
        snp = sorted(folder.glob('*.snp'))
        snp_site_count = _count_lines(snp[-1]) if snp else None

    rows = []
    for record in load_aadr_anno(candidates[-1]):
        hit = _as_int(record.get('snps_1240k'))
        call_rate = (hit / snp_site_count) if (hit and snp_site_count) else float('nan')
        rows.append({
            'sample_id': record['genetic_id'],
            'dataset': dataset_label or ds.key,
            'group': record['group'],
            'era': _aadr_era(record),
            'build': ds.build,
            'ploidy': _aadr_ploidy(record),
            'call_rate': call_rate,
            'contamination': _aadr_contamination(record),
            'locality': record.get('locality', ''),
            'date_mean_bp': record.get('date_mean_bp', ''),
            'assessment': record.get('assessment', ''),
        })
    return pd.DataFrame(rows)


def _as_int(value) -> Optional[int]:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None


def _aadr_contamination(record: dict) -> float:
    """Upper bound of the male X-based hapConX CI, else NaN (never invented)."""
    raw = (record.get('hapconx_ci') or '').strip()
    if raw.startswith('[') and raw.endswith(']'):
        try:
            high = raw.strip('[]').split(',')[1]
            return float(high)
        except (IndexError, ValueError):
            return float('nan')
    return float('nan')


def _count_lines(path: Path) -> int:
    total = 0
    with Path(path).open('rb') as handle:
        for _ in handle:
            total += 1
    return total


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def _cli(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description='HelixAtlas academic dataset broker (verified catalog + '
                    'bounded metadata download + honest sample browser).')
    sub = parser.add_subparsers(dest='command', required=True)

    sub.add_parser('catalog', help='print the verified dataset catalog')

    p_list = sub.add_parser('list', help='list samples as "Group (N)"')
    p_list.add_argument('dataset')
    p_list.add_argument('--data-dir', default='datasets')
    p_list.add_argument('--search')
    p_list.add_argument('--era', choices=['ancient', 'modern'])
    p_list.add_argument('--min-n', type=int, default=1)

    p_get = sub.add_parser('download-metadata', help='fetch small index/annotation files')
    p_get.add_argument('dataset')
    p_get.add_argument('--out', default='datasets')

    p_be = sub.add_parser('backend', help='show how to process the large genotype files')
    p_be.add_argument('dataset')

    args = parser.parse_args(argv)
    try:
        if args.command == 'catalog':
            for row in catalog_rows():
                print(f"{row['key']:<24} {row['access']:<11} {row['build']:<8} "
                      f"{row['name']}")
            return 0
        if args.command == 'list':
            for line in list_samples(args.dataset, data_dir=args.data_dir,
                                     group_contains=args.search, era=args.era,
                                     min_n=args.min_n):
                print(line)
            return 0
        if args.command == 'download-metadata':
            for kind, path in download_metadata(args.dataset, args.out).items():
                print(f'{kind}: {path}')
            return 0
        if args.command == 'backend':
            print(describe_backend(args.dataset))
            return 0
    except DatasetError as exc:
        print(f'error: {exc}', file=sys.stderr)
        return 2
    return 1


if __name__ == '__main__':
    raise SystemExit(_cli())
