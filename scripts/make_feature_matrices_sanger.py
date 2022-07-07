import collections
import pandas as pd
import subprocess
from taigapy import TaigaClient
from tqdm import *

# Window size in bp
WINDOW_SIZE = 2000

# Create output directory
subprocess.run("rm -rf feature_matrices_sanger", shell=True, check=True)
subprocess.run("mkdir -p feature_matrices_sanger", shell=True, check=True)

tc = TaigaClient()

# Read in raw data
df = pd.read_feather('sanger_methylation_data.feather')
df = df.set_index('index')
df = df.T

CCLE_expression = tc.get(name='depmap-a0ab', version=116, file='CCLE_expression')
all_genes_with_expression_data = set([element.split(' ')[0] for element in CCLE_expression.columns])

# Read in 450k probe to genomic position mapping
methylation_probemap_gene_split = tc.get(name='tcga-metadata-d555', version=2, file='methylation_probemap_gene_split')
methylation_probemap_gene_split = methylation_probemap_gene_split.drop_duplicates(['#id'])

# Read in RefSeq gene to TSS data and subset it to only genes found in the CCLE Expression data
refseq_tss_data = pd.read_csv('', sep='\t')
refseq_tss_data = refseq_tss_data[refseq_tss_data['name2'].isin(all_genes_with_expression_data)]

# Map each base pair within window of a gene's TSS to that gene
pos_to_genes = collections.defaultdict(lambda : collections.defaultdict(set))
for index, row in tqdm(refseq_tss_data.iterrows(), total=refseq_tss_data.shape[0]):
    chrom = row['chrom']
    gene = row['name2']
    txStart = row['txStart']
    txEnd = row['txEnd']
    strand = row['strand']
    
    for i in range(-1 * WINDOW_SIZE, WINDOW_SIZE + 1):
        if strand == '+':
            pos_to_genes[chrom][txStart + i].add(gene)
        elif strand == '-':
            pos_to_genes[chrom][txEnd + i].add(gene)
        else:
            assert False
        
# Associate each gene to its proper set of probes
gene_to_probes = collections.defaultdict(set)
for index, row in tqdm(methylation_probemap_gene_split.iterrows(), total=methylation_probemap_gene_split.shape[0]):
    probe = row['#id']
    chrom = row['chrom']
    pos_start = row['chromStart']
    pos_end = row['chromEnd']
    
    for gene in pos_to_genes[chrom][pos_start]:
        gene_to_probes[gene].add(probe)
        
    for gene in pos_to_genes[chrom][pos_end]:
        gene_to_probes[gene].add(probe)
        
# Write each gene's feature matrix to output directory
for gene in tqdm(gene_to_probes):
    gene_probes = gene_to_probes[gene]
    feature_matrix = df[gene_probes]
    feature_matrix.to_csv(f"feature_matrices_sanger/{gene}.csv")
    