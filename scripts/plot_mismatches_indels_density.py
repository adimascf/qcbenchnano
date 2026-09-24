import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

file_path = 'quast-BPH/contigs_reports/minimap_output/BPH2947__202310-filtlong-barbell-20x-assembly.used_snps'
cols = ['ref_name', 'contig_name', 'ref_pos', 'ref_base', 'asm_base', 'contig_pos']
df = pd.read_csv(file_path, sep='\t', header=None, names=cols)

# A '.' in either the ref or asm column indicates an Indel.
def classify_error(row):
    if '.' in str(row['ref_base']) or '.' in str(row['asm_base']):
        return 'Indel'
    else:
        return 'Mismatch'

df['error_type'] = df.apply(classify_error, axis=1)

# Filter for just the main chromosome to keep the x-axis linear
df_chrom = df[df['ref_name'] == 'chromosome']

plt.figure(figsize=(14, 6))

sns.histplot(
    data=df_chrom,
    x='ref_pos',
    hue='error_type',
    binwidth=100000,
    element='poly',
    fill=False,
    palette={'Mismatch': '#0072b2', 'Indel': '#d55e00'} 
)

plt.title('Error Density along Chromosome (100 kbp bins)', fontsize=14)
plt.xlabel('Genomic Position (bp)', fontsize=12)
plt.ylabel('Errors per 100 kbp', fontsize=12)

ticks = plt.gca().get_xticks()
plt.gca().set_xticklabels([f'{x/1000000:.1f} Mb' for x in ticks])

sns.despine()
plt.tight_layout()

plt.savefig('BPH2947__202310.filtlong-barbell_error_density_10kbp.png', dpi=300)
