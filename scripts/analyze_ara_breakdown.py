import argparse
import glob
import json
import os
import pandas as pd
import pyarrow.parquet as pq

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEFAULT_DATA_ROOT = os.path.join(os.path.expanduser("~"), "AraPPLM_Data", "phase2")
DATA_ROOT = os.environ.get("ARAPPLM_DATA_ROOT", DEFAULT_DATA_ROOT)

parser = argparse.ArgumentParser(description="Analyze Arabidopsis source breakdown from Parquet warehouse")
parser.add_argument("--data-dir", default=os.path.join(DATA_ROOT, "evidence_dataset", "species_code=ARA"),
                    help="Path to species_code=ARA partitioned Parquet directory")
parser.add_argument("--output", default=os.path.join(WORKSPACE_ROOT, "data", "interim", "phase2", "audits", "ara_source_breakdown.json"),
                    help="Path to output JSON summary file")
args = parser.parse_args()

data_dir = args.data_dir

results = {}

for res in ['INTACT', 'BIOGRID', 'XLMS', 'STRING']:
    res_path = os.path.join(data_dir, f'source_resource={res}')
    files = glob.glob(os.path.join(res_path, '*.parquet'))
    cols = [
        'evidence_id', 'unordered_pair_key_provisional', 'participant_a_id_value', 'participant_b_id_value',
        'interaction_semantics', 'assay_family', 'detection_method_psi_mi_name',
        'pmid', 'throughput_class', 'intact_miscore', 'string_experimental_score', 'string_combined_score',
        'xlms_interprotein_flag', 'xlms_intraprotein_flag', 'xlms_score', 'xlms_evalue'
    ]
    dfs = []
    for f in files:
        avail_cols = [c for c in cols if c in pq.read_schema(f).names]
        dfs.append(pq.read_table(f, columns=avail_cols).to_pandas())
    df = pd.concat(dfs, ignore_index=True)
    
    pairs = df['unordered_pair_key_provisional']
    unique_pairs = set(pairs)
    
    def is_homo(p):
        parts = str(p).split('__')
        return len(parts) == 2 and parts[0] == parts[1]
    
    homo_count = sum(1 for p in unique_pairs if is_homo(p))
    hetero_count = len(unique_pairs) - homo_count
    
    proteins = set(df['participant_a_id_value'].dropna().astype(str)).union(set(df['participant_b_id_value'].dropna().astype(str)))
    
    res_dict = {
        'total_records': len(df),
        'unique_pairs_total': len(unique_pairs),
        'unique_pairs_hetero': hetero_count,
        'unique_pairs_homo': homo_count,
        'unique_proteins': len(proteins),
        'interaction_semantics': df['interaction_semantics'].value_counts().to_dict(),
        'assay_family': df['assay_family'].value_counts().to_dict(),
        'detection_methods_top10': df['detection_method_psi_mi_name'].value_counts().head(10).to_dict(),
        'crosstab_semantics_assay': {
            idx: {col: int(val) for col, val in row.items() if val > 0}
            for idx, row in pd.crosstab(df['interaction_semantics'], df['assay_family']).iterrows()
        }
    }
    
    if 'pmid' in df.columns:
        valid_pmids = df['pmid'].dropna().loc[lambda x: (x != '') & (x != 'unassigned') & (x != 'None') & (x != 'NA')]
        res_dict['unique_pmids'] = int(valid_pmids.nunique())
        res_dict['top_pmids'] = {str(k): int(v) for k, v in valid_pmids.value_counts().head(10).items()}
    
    if 'throughput_class' in df.columns:
        res_dict['throughput_class'] = {str(k): int(v) for k, v in df['throughput_class'].value_counts(dropna=False).items()}
        
    if res == 'INTACT':
        scores = pd.to_numeric(df['intact_miscore'], errors='coerce').dropna()
        res_dict['miscore_stats'] = {
            'count': int(len(scores)),
            'mean': round(float(scores.mean()), 4),
            'median': round(float(scores.median()), 4),
            'ge_045_records': int((scores >= 0.45).sum()),
            'lt_045_records': int((scores < 0.45).sum())
        }
        df_45 = df[pd.to_numeric(df['intact_miscore'], errors='coerce') >= 0.45]
        res_dict['miscore_stats']['ge_045_unique_pairs'] = int(df_45['unordered_pair_key_provisional'].nunique())
        
    if res == 'STRING':
        exp_scores = pd.to_numeric(df['string_experimental_score'], errors='coerce').dropna()
        comb_scores = pd.to_numeric(df['string_combined_score'], errors='coerce').dropna()
        res_dict['string_stats'] = {
            'exp_gt_0': int((exp_scores > 0).sum()),
            'exp_ge_700': int((exp_scores >= 700).sum()),
            'comb_ge_700': int((comb_scores >= 700).sum()),
        }
        df_exp700 = df[pd.to_numeric(df['string_experimental_score'], errors='coerce') >= 700]
        res_dict['string_stats']['exp_ge_700_unique_pairs'] = int(df_exp700['unordered_pair_key_provisional'].nunique())
        
    if res == 'XLMS':
        res_dict['xlms_stats'] = {
            'interprotein_records': int((df['xlms_interprotein_flag'] == True).sum()),
            'intraprotein_records': int((df['xlms_intraprotein_flag'] == True).sum()),
        }
        df_inter = df[df['xlms_interprotein_flag'] == True]
        res_dict['xlms_stats']['interprotein_unique_pairs'] = int(df_inter['unordered_pair_key_provisional'].nunique())
        
    results[res] = res_dict

os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
with open(args.output, 'w', encoding='utf-8') as f:
    json.dump(results, f, indent=2)

print(f'SUCCESS: Results written to {args.output}')
