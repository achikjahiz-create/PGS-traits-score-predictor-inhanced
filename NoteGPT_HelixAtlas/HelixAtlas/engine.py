"""Conservative, autosomal biallelic SNP research scoring engine."""
import numpy as np
import pandas as pd

KEY = ['chrom', 'pos']


def require(df, columns, name):
    missing = set(columns) - set(df.columns)
    if missing:
        raise ValueError(f'{name}: missing columns {sorted(missing)}')


def normalize_variants(df):
    df = df.copy()
    df['chrom'] = df.chrom.astype(str).str.replace(r'^chr', '', regex=True)
    p = pd.to_numeric(df.pos, errors='raise')
    if p.isna().any() or (p < 1).any() or (p % 1 != 0).any():
        raise ValueError('Positions must be positive integers')
    df['pos'] = p.astype('int64')
    return df


def validate_metadata(meta):
    require(meta, ['sample_id', 'dataset', 'group', 'era', 'build', 'ploidy',
                   'call_rate', 'contamination'], 'metadata')
    meta = meta.copy()
    for col in ['sample_id', 'dataset', 'group', 'era', 'build']:
        if meta[col].isna().any():
            raise ValueError(f'Missing {col}')
        meta[col] = meta[col].astype(str)
        if meta[col].str.strip().eq('').any():
            raise ValueError(f'Empty {col}')
    if meta.sample_id.duplicated().any():
        raise ValueError('sample_id must be globally unique')
    if not meta.era.isin(['ancient', 'modern']).all():
        raise ValueError('era must be ancient or modern')
    if not meta.build.isin(['GRCh37', 'GRCh38']).all():
        raise ValueError('build must be GRCh37 or GRCh38')
    meta['ploidy'] = pd.to_numeric(meta.ploidy, errors='raise')
    if not meta.ploidy.isin([1, 2]).all():
        raise ValueError('ploidy must be 1 or 2')
    for col in ['call_rate', 'contamination']:
        meta[col] = pd.to_numeric(meta[col], errors='raise')
        if meta[col].isna().any() or not meta[col].between(0, 1).all():
            raise ValueError(f'{col} must be a measured fraction in [0,1]')
    return meta


def validate_genotypes(g, meta):
    require(g, ['sample_id', 'chrom', 'pos', 'ref', 'alt', 'dosage', 'ploidy'], 'genotypes')
    g = normalize_variants(g)
    g['sample_id'] = g.sample_id.astype(str)
    if not g.sample_id.isin(meta.sample_id).all():
        raise ValueError('Unknown genotype sample_id')
    if g.duplicated(['sample_id'] + KEY).any():
        raise ValueError('Duplicate sample/site: resolve upstream')
    for col in ['ref', 'alt']:
        g[col] = g[col].astype(str).str.upper()
    if not g.ref.isin(list('ACGT')).all() or not g.alt.isin(list('ACGT')).all():
        raise ValueError('Only single-base A/C/G/T SNPs supported')
    if g.ref.eq(g.alt).any():
        raise ValueError('REF equals ALT')
    if not g.chrom.isin([str(i) for i in range(1, 23)]).all():
        raise ValueError('Only autosomes 1–22 supported')
    for col in ['dosage', 'ploidy']:
        g[col] = pd.to_numeric(g[col], errors='raise')
    if not g.ploidy.isin([1, 2]).all():
        raise ValueError('Invalid genotype ploidy')
    if (g.ploidy != g.sample_id.map(meta.set_index('sample_id').ploidy)).any():
        raise ValueError('Genotype/metadata ploidy mismatch')
    observed = g.dosage.notna()
    if ((g.loc[observed, 'dosage'] < 0) | (g.loc[observed, 'dosage'] > g.loc[observed, 'ploidy'])).any():
        raise ValueError('ALT dosage outside [0,ploidy]')
    return g


def qc_flags(meta, cutoff=4.0):
    """Quality-only directional MAD flags within technical/cohort strata."""
    result = pd.Series(False, index=meta.index)
    for _, group in meta.groupby(['dataset', 'group', 'era', 'ploidy', 'build'], dropna=False):
        if len(group) < 10:
            continue
        for col, direction in [('call_rate', -1), ('contamination', 1)]:
            x = group[col]
            median = x.median()
            mad = (x - median).abs().median()
            if mad > 0:
                z = 0.67448975 * (x - median) / mad
                result.loc[group.index] |= direction * z > cutoff
    return result


def prepare_weights(w):
    require(w, ['chrom', 'pos', 'effect_allele', 'other_allele', 'effect_weight'], 'weights')
    w = normalize_variants(w)
    for col in ['effect_allele', 'other_allele']:
        w[col] = w[col].astype(str).str.upper()
    w['effect_weight'] = pd.to_numeric(w.effect_weight, errors='raise')
    if not np.isfinite(w.effect_weight).all():
        raise ValueError('Nonfinite weights')
    if w.duplicated(KEY).any():
        raise ValueError('Duplicate weight positions unsupported')
    valid = w.chrom.isin([str(i) for i in range(1, 23)])
    valid &= w.effect_allele.isin(list('ACGT')) & w.other_allele.isin(list('ACGT'))
    valid &= w.effect_allele.ne(w.other_allele)
    pairs = w.effect_allele + w.other_allele
    valid &= ~pairs.isin(['AT', 'TA', 'CG', 'GC'])
    clean = w.loc[valid].copy()
    if clean.empty:
        raise ValueError('No supported non-palindromic weights')
    return clean, int((~valid).sum())


def score(g, meta, weights, build, min_coverage=0.9):
    if not 0 < min_coverage <= 1:
        raise ValueError('Invalid coverage threshold')
    if meta.empty:
        raise ValueError('No selected samples')
    if not meta.build.eq(build).all():
        raise ValueError('Build mismatch; liftover must be performed upstream')
    if not meta.ploidy.eq(2).all():
        raise ValueError('Pseudo-haploid/haploid scoring is disabled')
    w, rejected = prepare_weights(weights)
    ids = meta.sample_id
    joined = g[g.sample_id.isin(ids)].merge(w, on=KEY, validate='many_to_one')
    ea, oa = joined.effect_allele, joined.other_allele
    forward = ea.eq(joined.alt) & oa.eq(joined.ref)
    reverse = ea.eq(joined.ref) & oa.eq(joined.alt)
    matched = joined[(forward | reverse) & joined.dosage.notna()].copy()
    matched['effect_dosage'] = np.where(ea.loc[matched.index].eq(matched.alt),
                                        matched.dosage, 2 - matched.dosage)
    counts = matched.groupby(KEY).sample_id.nunique()
    common = counts[counts.eq(len(ids))].reset_index()[KEY]
    coverage = len(common) / len(w)
    if coverage < min_coverage:
        raise ValueError(f'Common-site coverage {coverage:.2%} below {min_coverage:.2%}; no score emitted')
    used = matched.merge(common, on=KEY, validate='many_to_one')
    used['contribution'] = used.effect_dosage * used.effect_weight
    totals = used.groupby('sample_id').contribution.sum(min_count=1)
    out = meta.merge(totals.rename('raw_partial_pgs'), left_on='sample_id', right_index=True)
    out['matched_variants'] = len(common)
    out['eligible_variants'] = len(w)
    out['coverage'] = coverage
    return out, {'rejected_variants': rejected, 'eligible_variants': len(w),
                 'common_variants': len(common), 'coverage': coverage}
