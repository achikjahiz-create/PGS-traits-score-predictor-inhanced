"""Offline tests for the HelixAtlas dataset broker.

These tests never touch the network: they build miniature ``.anno``/``.ind``/
``.psam`` fixtures and exercise the parsers, the sample browser, the era/ploidy
heuristics and the access/host guards.
"""
import csv
import textwrap

import pytest

import datasets as D


# 49-column AADR-style header, matching the real positional layout.
def _anno_header():
    header = ['x'] * 49
    header[0] = 'Genetic ID'
    header[1] = 'Persistent Genetic ID'
    header[2] = 'Individual ID'
    header[10] = 'Date mean in BP'
    header[14] = 'Group ID'
    header[15] = 'Locality'
    header[21] = 'Data type'
    header[23] = 'Mean coverage'
    header[25] = 'SNPs hit 2M'
    header[26] = 'SNPs hit 1240k'
    header[30] = 'Molecular Sex'
    header[43] = 'hapConX 95% CI'
    header[47] = 'ASSESSMENT'
    return header


def _anno_row(gid, pid, group, date, dtype, snps, hapconx, sex='M'):
    row = [''] * 49
    row[0], row[1], row[2] = gid, pid, gid.split('.')[0]
    row[10], row[14], row[15] = date, group, 'Loc'
    row[21], row[26], row[30], row[43], row[47] = dtype, snps, sex, hapconx, 'Pass'
    return row


@pytest.fixture()
def aadr_dir(tmp_path):
    folder = tmp_path / 'AADR_V66_1240K'
    folder.mkdir()
    rows = [
        _anno_row('A.AG', '1', 'Austria_Avar', '1200', '.AG', '700000', '[0.005,0.012]'),
        _anno_row('B.AG', '2', 'Austria_Avar', '1250', '.AG', '650000', '..'),
        _anno_row('C.SG', '3', 'Han', '0', 'Shotgun', '1100000', 'n/a (female)', 'F'),
        _anno_row('D.DG', '4', 'Han', '0', 'Shotgun.diploid', '1100000', '..', 'M'),
        _anno_row('E.AG', '5', 'Russia_Yamnaya', '4500', '.AG', '800000', '..'),
    ]
    with (folder / 'v66.p1_1240K.aadr.PUB.anno').open('w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(_anno_header())
        writer.writerows(rows)
    (folder / 'v66.p1_1240K.aadr.patch.PUB.ind').write_text(textwrap.dedent(
        'A.AG M Austria_Avar\nB.AG M Austria_Avar\nC.SG F Han\n'))
    (folder / 'v66.p1_1240K.aadr.patch.PUB.snp').write_text('rs1 1 0.0 100 A G\n')
    return tmp_path


def test_catalog_entry_lookup():
    assert D.get('AADR_V66_1240K').access == 'public'
    with pytest.raises(D.DatasetError):
        D.get('NOPE')


def test_all_catalog_urls_are_allowlisted():
    for ds in D.DATASET_CATALOG.values():
        for remote in ds.files:
            D._check_host(remote.url)  # raises if host is not trusted


def test_anno_parse_and_group_counts(aadr_dir):
    records = D.load_aadr_anno(aadr_dir / 'AADR_V66_1240K' / 'v66.p1_1240K.aadr.PUB.anno')
    assert len(records) == 5
    assert D.group_counts(records) == [('Austria_Avar', 2), ('Han', 2), ('Russia_Yamnaya', 1)]


def test_list_samples_from_disk(aadr_dir):
    listing = D.list_samples('AADR_V66_1240K', data_dir=aadr_dir)
    assert 'Austria_Avar (2)' in listing
    assert 'Han (2)' in listing


def test_search_and_filters(aadr_dir):
    assert D.search_samples('AADR_V66_1240K', 'yam', data_dir=aadr_dir) == ['Russia_Yamnaya (1)']
    assert D.list_samples('AADR_V66_1240K', data_dir=aadr_dir, min_n=2) == [
        'Austria_Avar (2)', 'Han (2)']
    assert D.list_samples('AADR_V66_1240K', data_dir=aadr_dir, era='modern') == [
        'Han (2)']


def test_era_and_ploidy_heuristics(aadr_dir):
    frame = D.to_metadata_frame('AADR_V66_1240K', data_dir=aadr_dir, snp_site_count=1000000)
    by_id = frame.set_index('sample_id')
    assert by_id.loc['A.AG', 'era'] == 'ancient'
    assert by_id.loc['C.SG', 'era'] == 'modern'
    assert by_id.loc['A.AG', 'ploidy'] == 1
    assert by_id.loc['D.DG', 'ploidy'] == 2
    # call rate = measured sites / snp_site_count
    assert by_id.loc['A.AG', 'call_rate'] == pytest.approx(0.7)
    # contamination from the hapConX upper bound; never invented elsewhere
    assert by_id.loc['A.AG', 'contamination'] == pytest.approx(0.012)
    assert by_id.loc['B.AG', 'contamination'] != by_id.loc['B.AG', 'contamination']  # NaN


def test_missing_metadata_raises_not_invents(tmp_path):
    with pytest.raises(D.DatasetNotFetched):
        D.list_samples('AADR_V66_1240K', data_dir=tmp_path)


def test_controlled_access_is_refused():
    for key in ('LASI_DAD_NIAGADS', 'HUMAN_ORIGINS_ANGOLA_EGA'):
        with pytest.raises(D.AccessRestricted):
            D.list_samples(key)
        with pytest.raises(D.AccessRestricted):
            D.download_metadata(key, out_dir='/tmp/never')


def test_untrusted_host_refused():
    with pytest.raises(D.DatasetError):
        D._check_host('https://evil.example.com/data.geno')


def test_backend_never_returns_geno_for_in_app(monkeypatch):
    text = D.describe_backend('AADR_V66_1240K')
    assert 'pgsc_calc' in text
    assert 'intentionally NOT pulled' in text
