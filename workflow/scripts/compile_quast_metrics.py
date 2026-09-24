import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print("Calculating genome sizes, expected contigs, N50, and auN from FASTA indices...")
genome_sizes = {}
expected_contigs = {}
ref_n50s = {}
ref_auns = {}

# Ensure it's treated as a list even if there's only one sample
fai_files = snakemake.input.fai if isinstance(snakemake.input.fai, list) else [snakemake.input.fai]

for fai_file in fai_files:
    fai_path = Path(fai_file)
    # Extracts the sample name (everything before the first '.')
    sample_name = fai_path.name.split('.')[0]
    
    try:
        with open(fai_path, "r") as f:
            lines = f.readlines()
            
            # The second column (index 1) in a .fai file is the sequence length
            lengths = [int(line.split("\t")[1]) for line in lines]
            total_len = sum(lengths)
            
            genome_sizes[sample_name] = total_len
            expected_contigs[sample_name] = len(lengths)
            
            if total_len > 0:
                # Calculate reference auN: sum of squared lengths / total length
                ref_auns[sample_name] = sum(l * l for l in lengths) / total_len
                
                # Calculate reference N50
                lengths.sort(reverse=True)
                cumsum = 0
                for l in lengths:
                    cumsum += l
                    if cumsum >= total_len / 2.0:
                        ref_n50s[sample_name] = l
                        break
            else:
                ref_auns[sample_name] = np.nan
                ref_n50s[sample_name] = np.nan

    except FileNotFoundError:
        print(f"Warning: .fai file not found: {fai_file}")

assembly_metrics = []
for file_path in snakemake.input.reports:
    try:
        p = Path(file_path)
        
        # Path: ../quast/<combo>/<depth>x/<model>/<sample>.<combo>.report.tsv
        model = p.parts[-2]
        depth = p.parts[-3]
        combo = p.parts[-4]
        
        # Extract the sample name from the file name
        sample = p.name.split('.')[0]

        df_quast = pd.read_csv(p, sep='\t', index_col=0)
        col_name = df_quast.columns[0]

        mismatches = float(df_quast.loc['# mismatches per 100 kbp', col_name]) if '# mismatches per 100 kbp' in df_quast.index else np.nan
        indels = float(df_quast.loc['# indels per 100 kbp', col_name]) if '# indels per 100 kbp' in df_quast.index else np.nan
        nga50 = float(df_quast.loc['NGA50', col_name]) if 'NGA50' in df_quast.index else np.nan
        misassemblies = int(df_quast.loc['# misassemblies', col_name]) if '# misassemblies' in df_quast.index else np.nan
        contigs_obs = int(df_quast.loc['# contigs', col_name]) if '# contigs' in df_quast.index else np.nan
        
        aunga = float(df_quast.loc['auNGA', col_name]) if 'auNGA' in df_quast.index else np.nan
        dup_ratio = float(df_quast.loc['Duplication ratio', col_name]) if 'Duplication ratio' in df_quast.index else np.nan

        # Denominators from FASTA index
        ref_n50 = ref_n50s.get(sample, np.nan)
        ref_aun = ref_auns.get(sample, np.nan)
        
        # Calculate the ratio of the selected contiguity metrics
        nga50_norm = nga50 / ref_n50 if pd.notna(ref_n50) and ref_n50 > 0 else np.nan
        aunga_ratio = aunga / ref_aun if pd.notna(ref_aun) and ref_aun > 0 else np.nan

        c_exp = expected_contigs.get(sample, np.nan)
        excess_contigs = contigs_obs - c_exp if pd.notna(c_exp) else np.nan

        assembly_metrics.append({
            'combo': combo, 'depth': depth, 'sample': sample, 'model': model,
            'Mismatches per 100kbp': mismatches,
            'Indels per 100kbp': indels,
            'NGA50': nga50,
            'NGA50_norm': nga50_norm,
            'auNGA': aunga,
            'auNGA_ratio': aunga_ratio, 
            'Duplication_ratio': dup_ratio,
            'misassemblies': misassemblies,
            'excess_contigs': excess_contigs
        })

    except Exception as e:
        print(f"Skipping {file_path} due to error: {e}")

df = pd.DataFrame(assembly_metrics)
df.sort_values(by=["combo", "sample", "model", "depth"], inplace=True)
df.to_csv(snakemake.output.csv, index=False)
print(f"Saved compiled metrics to {snakemake.output.csv}")
