#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
import gzip
import os
import random

def gene_set_to_bed(args):
    print('making gene set bed file')
    GeneSet = pd.read_csv(args.gene_set_file, header = None, names = ['GENE'])
    all_genes = pd.read_csv(args.gene_coord_file, delim_whitespace = True)
    df = pd.merge(GeneSet, all_genes, on = 'GENE', how = 'inner')
    df['START'] = np.maximum(1, df['START'] - args.windowsize)
    df['END'] = df['END'] + args.windowsize
    iter_df = [['chr'+(str(x1).lstrip('chr')), x2 - 1, x3] for (x1,x2,x3) in np.array(df[['CHR', 'START', 'END']])]
    return BedTool(iter_df).sort().merge()

def replace_with_median(value):
    if isinstance(value, str) and ',' in value:
        values = [float(x) for x in value.split(',')]
        median = np.median(values)
        return median
    else:
        return value
    
def make_annot_files(args, bed_for_annot):
    print('making annot file')
    df_bim = pd.read_csv(args.bimfile, delim_whitespace=True, usecols=[0,1,2,3], names=['CHR','SNP','CM','BP'])
    df_bim['CHR'] = df_bim['CHR'].astype(str)
    df_bim['CHR'] = 'chr' + df_bim['CHR']
    
    random_sequence = ''.join(random.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789') for _ in range(12))
    tmp_bed = args.bed_file + random_sequence + "temp.bed"
    bed_for_annot.saveas(tmp_bed)

    df_annot= pd.read_csv(tmp_bed, sep='\t', header =None)

    num_cols = bed_for_annot.field_count()
    if num_cols==4:
        print("This is a continuous annotation")
        df_annot.columns = ['CHR', 'BP', 'END', 'ANNOT']
        df_annot['BP'] += 1

        #os.remove(tmp_bed)

        #print("This is the first ten lines of df_annot before merging with bim")
        #print(df_annot.head(10))

        df_annot = df_annot[df_annot['CHR'].isin(df_bim['CHR']) & df_annot['BP'].isin(df_bim['BP'])]
        df_annot = pd.merge(df_bim, df_annot.drop('END', axis=1), on=['CHR', 'BP'], how='left')

        df_annot.fillna(0, inplace=True)
        df_annot['ANNOT'] = df_annot['ANNOT'].apply(replace_with_median)
        
        #print("This is after merging with bim")
        print(df_annot.head(10))

        non_zero_scores = (df_annot['ANNOT'] != 0).sum()
        zero_scores = (df_annot['ANNOT'] == 0).sum()
        print('Number of non-zero scores:', non_zero_scores)
        print('Number of zero scores:', zero_scores)

        output_file = args.annot_file
        if args.annot_file.endswith('.gz'):
            with gzip.open(args.annot_file, 'wb') as f:
                df_annot.to_csv(output_file, sep='\t', index=False, float_format='%.6g')
        else:
            df_annot.to_csv(output_file, sep='\t', index=False, float_format='%.6g')

    elif num_cols==3:
        print("This is a binary annotation")

        df_annot.columns = ['CHR', 'BP', 'END']
        df_annot['BP'] += 1
        df_annot['ANNOT'] = 1
   
        #print("This is the first ten lines of df_annot before merging with bim")
        #print(df_annot.head(10))

        df_annot = df_annot[df_annot['CHR'].isin(df_bim['CHR']) & df_annot['BP'].isin(df_bim['BP'])]
        df_annot = pd.merge(df_bim, df_annot.drop('END', axis=1), on=['CHR', 'BP'], how='left')
        
        df_annot.fillna(0, inplace=True)
        df_annot['ANNOT'] = df_annot['ANNOT'].astype(int)
        
        print("This is after merging with bim")
        print(df_annot.head(10))

        non_zero_scores = (df_annot['ANNOT'] != 0).sum()
        zero_scores = (df_annot['ANNOT'] == 0).sum()
        print('Number of non-zero scores:', non_zero_scores)
        print('Number of zero scores:', zero_scores)
        
        output_file = args.annot_file
        if args.annot_file.endswith('.gz'):
            with gzip.open(args.annot_file, 'wb') as f:
                df_annot.to_csv(output_file, sep='\t', index=False, float_format='%.6g')
        else:
            df_annot.to_csv(output_file, sep='\t', index=False, float_format='%.6g')
        
    else:
        print("Unexpected number of columns in bed file")
    # Delete the tmp_bed file
    os.remove(tmp_bed)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene-set-file', type=str, help='a file of gene names, one line per gene.')
    parser.add_argument('--gene-coord-file', type=str, default='ENSG_coord.txt', help='a file with columns GENE, CHR, START, and END, where START and END are base pair coordinates of TSS and TES. This file can contain more genes than are in the gene set. We provide ENSG_coord.txt as a default.')
    parser.add_argument('--windowsize', type=int, help='how many base pairs to add around the transcribed region to make the annotation?')
    parser.add_argument('--bed-file', type=str, help='the UCSC bed file with the regions that make up your annotation')
    parser.add_argument('--nomerge', action='store_true', default=False, help='don\'t merge the bed file; make an annot file wi    th values proportional to the number of intervals in the bedfile overlapping the SNP.')
    parser.add_argument('--bimfile', type=str, help='plink bim file for the dataset you will use to compute LD scores.')
    parser.add_argument('--annot-file', type=str, help='the name of the annot file to output.')

    args = parser.parse_args()

    if args.gene_set_file is not None:
        bed_for_annot = gene_set_to_bed(args)
    else:
        bed_for_annot = BedTool(args.bed_file).sort()
        if not args.nomerge:
            # Check if a 4th column exists in the original BED file
            has_fourth_column = len(bed_for_annot[0].fields) >= 4
            if has_fourth_column:
            # Merge overlapping intervals while retaining the fourth column
                bed_for_annot = bed_for_annot.merge(c=[4], o=['collapse'], delim=',')
                print("bed for annot merge")
                print(bed_for_annot.head(10))
            else:
                bed_for_annot = bed_for_annot.merge()
    make_annot_files(args, bed_for_annot)
