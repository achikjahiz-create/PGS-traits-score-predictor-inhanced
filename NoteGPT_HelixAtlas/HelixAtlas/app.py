import gzip
import hashlib
import io
import json
import re
from datetime import datetime, timezone
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd
import plotly.express as px
import requests
import streamlit as st
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from engine import validate_metadata, validate_genotypes, qc_flags, score
import datasets as data_broker

st.set_page_config(page_title='Genome Research Studio', page_icon='🧬', layout='wide')
MAX_BYTES = 128 * 1024 * 1024
WARNING = ('Research only. Raw partial scores are not predicted height, skin tone, hair or eye categories, '
           'intelligence, school performance, or medical risk. No population rankings. '
           'Ancient/pseudo-haploid samples require separate validation; scoring is disabled here. '
           'Changing the cohort can change the common variant set and therefore the score.')


def read_upload(upload):
    if upload.size > MAX_BYTES:
        raise ValueError('Upload exceeds 128 MiB; preprocess with a batch backend')
    blob = upload.getvalue()
    frame = pd.read_csv(io.BytesIO(blob), sep='\t' if upload.name.endswith('.tsv') else ',')
    return frame, hashlib.sha256(blob).hexdigest()


@st.cache_data(ttl=86400, max_entries=8, show_spinner=False)
def fetch_weights(pgs_id, build):
    if not re.fullmatch(r'PGS\d{6}', pgs_id) or build not in ['GRCh37', 'GRCh38']:
        raise ValueError('Invalid PGS identifier/build')
    url = (f'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs_id}/'
           f'ScoringFiles/Harmonized/{pgs_id}_hmPOS_{build}.txt.gz')
    with requests.get(url, stream=True, timeout=(10, 60), allow_redirects=False) as response:
        response.raise_for_status()
        if response.is_redirect:
            raise ValueError('Redirect refused')
        buffer = io.BytesIO()
        for chunk in response.iter_content(1024 * 1024):
            buffer.write(chunk)
            if buffer.tell() > MAX_BYTES:
                raise ValueError('Weight file exceeds web-app limit')
    blob = buffer.getvalue()
    with gzip.GzipFile(fileobj=io.BytesIO(blob)) as handle:
        raw = handle.read(MAX_BYTES + 1)
    if len(raw) > MAX_BYTES:
        raise ValueError('Decompressed weight file exceeds limit')
    text = raw.decode('utf-8-sig')
    headers = {}
    for line in text.splitlines():
        if line.startswith('#') and '=' in line:
            key, value = line.lstrip('#').split('=', 1)
            headers[key.strip()] = value.strip()
    if headers.get('HmPOS_build') != build:
        raise ValueError('Harmonized build header missing/mismatched')
    if headers.get('pgs_id') != pgs_id:
        raise ValueError('PGS file identifier mismatch')
    if headers.get('weight_type', '').lower() != 'beta':
        raise ValueError('Only additive beta weights supported; OR/HR/nonlinear scores refused')
    w = pd.read_csv(io.StringIO(text), sep='\t', comment='#', low_memory=False)
    w = w.rename(columns={'hm_chr': 'chrom', 'hm_pos': 'pos'})
    if w.columns.duplicated().any():
        raise ValueError('Ambiguous weight columns')
    return w, headers, {'url': url, 'sha256': hashlib.sha256(blob).hexdigest(),
                         'retrieved_utc': datetime.now(timezone.utc).isoformat()}


def synthetic_data():
    rng = np.random.default_rng(42)
    ids = [f'DEMO_{i:03}' for i in range(30)]
    meta = pd.DataFrame({'sample_id': ids, 'dataset': 'SYNTHETIC',
                         'group': ['A'] * 15 + ['B'] * 15, 'era': 'modern',
                         'build': 'GRCh38', 'ploidy': 2,
                         'call_rate': rng.uniform(.97, 1, 30), 'contamination': 0.0})
    g = pd.DataFrame([{'sample_id': sid, 'chrom': '1', 'pos': pos, 'ref': 'A',
                       'alt': 'G', 'dosage': rng.binomial(2, .35), 'ploidy': 2}
                      for sid in ids for pos in range(1001, 1101)])
    w = pd.DataFrame({'chrom': '1', 'pos': range(1001, 1101), 'effect_allele': 'G',
                       'other_allele': 'A', 'effect_weight': rng.normal(0, .03, 100)})
    return meta, g, w


def make_pdf(results, summary, manifest):
    buffer = io.BytesIO()
    styles = getSampleStyleSheet()
    story = [Paragraph('Genome Research Studio', styles['Title']),
             Paragraph(escape(WARNING), styles['BodyText']), Spacer(1, 12),
             Paragraph('Provenance and parameters', styles['Heading2'])]
    for key, value in manifest.items():
        story.append(Paragraph(escape(f'{key}: {value}'), styles['BodyText']))
    for title, frame in [('Cohort summaries', summary), ('Individual raw partial scores', results)]:
        story += [Spacer(1, 12), Paragraph(title, styles['Heading2'])]
        columns = ['dataset', 'group', 'n', 'mean', 'sd'] if title == 'Cohort summaries' else [
            'sample_id', 'raw_partial_pgs', 'matched_variants', 'coverage']
        view = frame[columns].copy()
        rows = [[Paragraph(escape(str(v)), styles['BodyText']) for v in columns]]
        for row in view.itertuples(index=False, name=None):
            rows.append([Paragraph(escape(f'{v:.5g}' if isinstance(v, float) else str(v)),
                                   styles['BodyText']) for v in row])
        table = Table(rows, repeatRows=1, hAlign='LEFT')
        table.setStyle(TableStyle([('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#E2E8F0')),
                                   ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                                   ('LINEBELOW', (0, 0), (-1, 0), .5, colors.grey)]))
        story.append(table)
    SimpleDocTemplate(buffer).build(story)
    return buffer.getvalue()


st.title('🧬 Genome Research Studio')
st.caption('Auditable SNP scoring • quality filters • cohort summaries • PDF export')
st.warning(WARNING)
demo = st.sidebar.toggle('Use clearly labeled synthetic demo', value=True)
st.sidebar.caption('Genotypes stay in this server session; use a private deployment for human DNA.')


def render_dataset_catalog():
    """Academic data catalog + honest sample browser (metadata files only)."""
    with st.sidebar.expander('Academic data catalog (metadata only)'):
        st.caption('Verified public ancient/modern genotype resources. Only small '
                   'index/annotation files are fetched here; multi-GB genotype '
                   'matrices stay on the pgsc_calc backend. Counts are parsed from '
                   'the downloaded files, never invented.')
        rows = data_broker.catalog_rows()
        st.dataframe(pd.DataFrame(rows)[['key', 'access', 'build', 'metadata_files']],
                     hide_index=True, use_container_width=True)
        key = st.selectbox('Browse a dataset', list(data_broker.DATASET_CATALOG))
        entry = data_broker.get(key)
        st.caption(f'{entry.name} • access: {entry.access} • build: {entry.build}')
        st.markdown(f'[Landing page / citation]({entry.landing_url})')
        for note in entry.notes:
            st.caption(f'• {note}')
        st.caption(f'Cite: {entry.citation}')
        if entry.access != 'public':
            st.info('Controlled access: request through the landing page. This app '
                    'will not auto-download controlled data.')
            return
        if st.button('Download small metadata files', key=f'dl_{key}'):
            try:
                with st.spinner('Downloading and verifying metadata'):
                    got = data_broker.download_metadata(key, 'datasets')
                st.success('Verified: ' + ', '.join(f'{k}={v.name}' for k, v in got.items()))
            except data_broker.DatasetError as exc:
                st.error(str(exc))
        try:
            listing = data_broker.list_samples(key, data_dir='datasets')
        except data_broker.DatasetNotFetched:
            st.caption('No local metadata yet — download the small index files to '
                       'list samples as "Group (N)".')
            return
        except data_broker.DatasetError as exc:
            st.caption(f'Catalog note: {exc}')
            return
        query = st.text_input('Filter groups', key=f'q_{key}')
        if query:
            listing = [s for s in listing if query.lower() in s.lower()]
        st.write(f'{len(listing)} groups with local genotypes available:')
        st.dataframe(pd.DataFrame({'group (n)': listing}), hide_index=True, height=220)


try:
    render_dataset_catalog()
except Exception as exc:  # catalog is a convenience, never block the scorer
    st.sidebar.caption(f'Catalog unavailable: {exc}')

try:
    if demo:
        meta, genotypes, weights = synthetic_data()
        hashes = {'metadata': 'synthetic-seed-42', 'genotypes': 'synthetic-seed-42'}
        st.info('SYNTHETIC DATA: no real individuals or biological predictions.')
    else:
        st.subheader('Import normalized files')
        st.caption('Metadata and long-form ALT-dosage CSV/TSV; formats are in the README.')
        m_upload = st.file_uploader('Sample metadata', type=['csv', 'tsv'])
        g_upload = st.file_uploader('Genotypes', type=['csv', 'tsv'])
        if m_upload is None or g_upload is None:
            st.stop()
        meta, mh = read_upload(m_upload)
        genotypes, gh = read_upload(g_upload)
        hashes = {'metadata': mh, 'genotypes': gh}
    meta = validate_metadata(meta)
    genotypes = validate_genotypes(genotypes, meta)
    st.sidebar.header('Quality rules')
    minimum_call = st.sidebar.slider('Minimum measured call rate', 0.0, 1.0, .95, .01)
    maximum_contamination = st.sidebar.slider('Maximum measured contamination', 0.0, 1.0, .05, .01)
    mad_cut = st.sidebar.slider('Directional quality MAD threshold', 2.0, 8.0, 4.0, .5)
    exclude_mad = st.sidebar.checkbox('Exclude flagged quality outliers', True)
    meta['quality_outlier'] = qc_flags(meta, mad_cut)
    meta['quality_pass'] = (meta.call_rate >= minimum_call) & (meta.contamination <= maximum_contamination)
    if exclude_mad:
        meta['quality_pass'] &= ~meta.quality_outlier
    st.caption('MAD flags use dataset/group/era/ploidy/build strata with at least 10 samples. '
               'Zero-MAD strata are not flagged; absolute quality thresholds still apply. '
               'Quality inputs must be measured upstream, not invented.')
    datasets = st.multiselect('Datasets', sorted(meta.dataset.unique()), default=sorted(meta.dataset.unique()))
    groups = st.multiselect('Groups', sorted(meta.group.unique()), default=sorted(meta.group.unique()))
    eras = st.multiselect('Ancient or modern', ['ancient', 'modern'], default=['modern'])
    subset = meta[meta.dataset.isin(datasets) & meta.group.isin(groups) & meta.era.isin(eras)]
    st.subheader('Editable sample selection; metadata remain immutable')
    selection = subset.copy()
    selection.insert(0, 'include', selection.quality_pass)
    edited = st.data_editor(selection, hide_index=True, disabled=list(subset.columns),
                            column_config={'include': st.column_config.CheckboxColumn('Include')},
                            key='sample_selection')
    chosen = edited[edited.include & edited.quality_pass].drop(columns='include')
    st.write(f'{len(chosen)} selected quality-passing samples of {len(meta)} imported.')
    if chosen.empty:
        st.stop()
    if (chosen.era.eq('ancient') | chosen.ploidy.ne(2)).any():
        st.error('This version browses ancient metadata but refuses ancient/haploid PGS scoring. '
                 'Select validated modern diploid samples. A dedicated ancient-DNA model is required.')
        st.stop()
    build = st.selectbox('Scoring genome build', ['GRCh38', 'GRCh37'])
    threshold = st.slider('Minimum common-site coverage of eligible weights', .5, 1.0, .9, .01)
    if demo:
        headers = {'pgs_id': 'SYNTHETIC', 'trait_reported': 'Synthetic demonstration', 'weight_type': 'beta'}
        provenance = {'source': 'synthetic seed 42'}
    else:
        st.markdown('Find a published score at [PGS Catalog](https://www.pgscatalog.org/). '
                    'This engine accepts only additive beta-weighted, biallelic autosomal models.')
        pgs_id = st.text_input('Published score ID (PGS followed by six digits)', '').strip().upper()
        if not pgs_id:
            st.stop()
        with st.spinner('Loading harmonized scoring file'):
            weights, headers, provenance = fetch_weights(pgs_id, build)
        st.json(headers, expanded=False)
        st.caption('Check model validation, phenotype, ancestry and license in the Catalog before use.')
        if not st.checkbox('I verified the score uses a simple additive model and accept its license and research limitations'):
            st.stop()
    results, details = score(genotypes, chosen, weights, build, threshold)
    summary = results.groupby(['dataset', 'group', 'build'], sort=False).raw_partial_pgs.agg(
        n='size', mean='mean', sd='std', median='median').reset_index()
    summary['mean_se'] = summary.sd / np.sqrt(summary.n)
    st.subheader('Cohort summaries—not phenotype estimates or population rankings')
    st.dataframe(summary, hide_index=True, use_container_width=True)
    st.caption('SE describes between-sample variation assuming independent samples; '
               'it does not capture model error, relatedness, missing variants, or prediction uncertainty.')
    t1, t2, t3 = st.tabs(['Distribution', 'Samples', 'Quality'])
    with t1:
        st.plotly_chart(px.histogram(results, x='raw_partial_pgs', nbins=25,
                                    title='Raw partial score distribution'), use_container_width=True)
        st.plotly_chart(px.scatter(results, x='call_rate', y='raw_partial_pgs',
                                  hover_name='sample_id', title='Quality versus raw partial score'),
                       use_container_width=True)
    with t2:
        st.dataframe(results, hide_index=True, use_container_width=True)
        selected_id = st.selectbox('Inspect a sample', results.sample_id.tolist())
        st.dataframe(results[results.sample_id.eq(selected_id)].T)
    with t3:
        st.dataframe(meta, hide_index=True, use_container_width=True)
        st.json(details)
    manifest = {'generated_utc': datetime.now(timezone.utc).isoformat(),
                'score_id': headers.get('pgs_id'), 'trait': headers.get('trait_reported'),
                'build': build, 'minimum_call_rate': minimum_call,
                'maximum_contamination': maximum_contamination, 'mad_cutoff': mad_cut,
                'exclude_quality_outliers': exclude_mad, 'minimum_common_coverage': threshold,
                'selected_samples': results.sample_id.tolist(), 'input_sha256': hashes,
                'weight_source': provenance, 'weight_headers': headers, 'matching': details,
                'missingness_policy': 'common observed non-palindromic exact-allele sites only',
                'warning': WARNING}
    st.download_button('Export sample CSV', results.to_csv(index=False).encode(), 'samples.csv', 'text/csv')
    st.download_button('Export cohort CSV', summary.to_csv(index=False).encode(), 'cohorts.csv', 'text/csv')
    st.download_button('Export reproducibility manifest', json.dumps(manifest, indent=2).encode(),
                       'manifest.json', 'application/json')
    st.caption('PDF contains the complete score tables and provenance; interactive charts are web-only in this version.')
    if len(results) <= 2000:
        pdf = make_pdf(results, summary, manifest)
        st.download_button('Export PDF report', pdf, 'research-report.pdf', 'application/pdf')
    else:
        st.info('PDF limited to 2,000 samples to prevent memory exhaustion. Export CSV or narrow selection.')
except (ValueError, KeyError, pd.errors.ParserError, requests.RequestException,
        UnicodeDecodeError, OSError) as exc:
    st.error(f'Analysis stopped safely: {exc}')
    st.stop()
