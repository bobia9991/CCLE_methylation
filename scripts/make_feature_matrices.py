import collections
import itertools
import numpy as np
import pandas as pd
import subprocess 
from taigapy import TaigaClient
from tqdm import *

# Window size in bp
WINDOW_SIZE = 2000

tc = TaigaClient()

# Delete and create temp directory to store terra output files
subprocess.run("rm -rf methylkit_outputs", shell=True, check=True)
subprocess.run("mkdir -p methylkit_outputs", shell=True, check=True)

# Read in terra data tsv file to get cloud filepaths of output
sample_info_terra = pd.read_csv('', sep='\t')
all_methylkit_output_files = sample_info_terra['methylkit_output'].dropna().tolist()

# Output file containing all google cloud files to download
filelist = open('methylkit_outputs_to_download.txt', 'w+')
for file in all_methylkit_output_files:
    filelist.write(f'{file}\n')
filelist.close()

# Copy all methylkit data down with gsutil. Make sure you have logged in to the right gsutil account or the download will fail
subprocess.run("cat methylkit_outputs_to_download.txt | gsutil -m cp -I methylkit_outputs", shell=True, check=True)
            
# Download CCLE Expression and read in RefSeq gene to TSS data
CCLE_expression = tc.get(name='depmap-a0ab', version=116, file='CCLE_expression')
all_genes_with_expression_data = set([element.split(' ')[0] for element in CCLE_expression.columns])
refseq_tss_data = pd.read_csv('', sep='\t')

# Subset RefSeq gene to TSS data to only genes found in the CCLE Expression data
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

# Read through all input files generated from MethylKit and populate the two dictionaries that store the feature matrices of each gene 
data_dir = 'methylkit_outputs'
gene_to_methylation_df_dict = collections.defaultdict(lambda : collections.defaultdict(dict))
gene_to_coverage_df_dict = collections.defaultdict(lambda : collections.defaultdict(dict))
for file in tqdm(all_methylkit_output_files):
    file = file.split('/')[-1]
    cell_line_name = file.split('_')[0]
    filepath = f'{data_dir}/{file}'
    with open(filepath) as f:
        header_line = f.readline()
        for line in f:
            tokens = line.split('\t')
            cpg_locus_tokens = tokens[0].split('.')
            chrom = cpg_locus_tokens[0]
            pos = int(cpg_locus_tokens[1])
            
            if chrom in pos_to_genes:
                if pos in pos_to_genes[chrom]:
                    cpg_locus_string = f"{chrom}:{pos}"
                    
                    coverage = int(tokens[4])
                    freqC = float(tokens[5])
                    freqT = float(tokens[6])
                    methylation = freqC / (freqC + freqT)
                    for gene in pos_to_genes[chrom][pos]:
                        gene_to_methylation_df_dict[gene][cell_line_name][cpg_locus_string] = methylation
                        gene_to_coverage_df_dict[gene][cell_line_name][cpg_locus_string] = coverage
         
# Create output directory
subprocess.run("rm -rf feature_matrices", shell=True, check=True)
subprocess.run("mkdir -p feature_matrices", shell=True, check=True)
          
# Write each gene's feature matrix to output directory
for gene in tqdm(gene_to_methylation_df_dict):
    df = pd.DataFrame(gene_to_methylation_df_dict[gene]).T
    df.to_csv(f'feature_matrices/{gene}.csv')
    
# Write each gene's coverage of CpG loci to output directory
for gene in tqdm(gene_to_coverage_df_dict):
    df = pd.DataFrame(gene_to_coverage_df_dict[gene]).T
    df.to_csv(f'feature_matrices/{gene}_coverage.csv')
    
# Remove temp directory
subprocess.run("rm -rf methylkit_outputs", shell=True, check=True)
