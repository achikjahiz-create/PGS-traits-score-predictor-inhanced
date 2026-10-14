import pandas as pd
import pytest
from engine import score, validate_metadata, validate_genotypes, qc_flags


def fixture():
    m = pd.DataFrame({'sample_id': ['s1', 's2'], 'dataset': ['d', 'd'],
        'group': ['g', 'g'], 'era': ['modern', 'modern'],
        'build': ['GRCh38', 'GRCh38'], 'ploidy': [2, 2],
        'call_rate': [1., 1.], 'contamination': [0., 0.]})
    g = pd.DataFrame({'sample_id': ['s1', 's2'], 'chrom': ['1', '1'],
        'pos': [101, 101], 'ref': ['A', 'A'], 'alt': ['G', 'G'],
        'dosage': [0., 2.], 'ploidy': [2, 2]})
    w = pd.DataFrame({'chrom': ['1'], 'pos': [101], 'effect_allele': ['G'],
        'other_allele': ['A'], 'effect_weight': [.5]})
    return m, g, w


def test_alt_effect():
    m, g, w = fixture()
    result, details = score(g, m, w, 'GRCh38')
    assert result.raw_partial_pgs.tolist() == [0., 1.]
    assert details['coverage'] == 1.


def test_ref_effect():
    m, g, w = fixture()
    w['effect_allele'], w['other_allele'] = 'A', 'G'
    result, _ = score(g, m, w, 'GRCh38')
    assert result.raw_partial_pgs.tolist() == [1., 0.]


def test_missing_common_site_refused():
    m, g, w = fixture()
    g.loc[1, 'dosage'] = float('nan')
    with pytest.raises(ValueError, match='coverage'):
        score(g, m, w, 'GRCh38')


def test_build_refused():
    m, g, w = fixture()
    with pytest.raises(ValueError, match='Build mismatch'):
        score(g, m, w, 'GRCh37')


def test_haploid_refused():
    m, g, w = fixture()
    m['ploidy'] = 1
    with pytest.raises(ValueError, match='haploid'):
        score(g, m, w, 'GRCh38')


def test_palindromic_refused():
    m, g, w = fixture()
    w['effect_allele'], w['other_allele'] = 'A', 'T'
    with pytest.raises(ValueError, match='non-palindromic'):
        score(g, m, w, 'GRCh38')


def test_duplicate_genotypes_refused():
    m, g, _ = fixture()
    with pytest.raises(ValueError, match='Duplicate'):
        validate_genotypes(pd.concat([g, g]), m)


def test_bad_dosage_refused():
    m, g, _ = fixture()
    g.loc[0, 'dosage'] = 3
    with pytest.raises(ValueError, match='dosage'):
        validate_genotypes(g, m)


def test_duplicate_sample_refused():
    m, _, _ = fixture()
    with pytest.raises(ValueError, match='unique'):
        validate_metadata(pd.concat([m, m]))


def test_small_qc_stratum_not_flagged():
    m, _, _ = fixture()
    assert not qc_flags(m).any()
