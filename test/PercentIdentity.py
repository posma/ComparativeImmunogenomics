import pandas as pd
import numpy as np
from Bio import Align
from Bio import SeqIO
import argparse

aligner = Align.PairwiseAligner()
aligner.mode = 'global' 

def pi_matrix(seqs, seqs2, flag=True):
        
    pi_matrix = np.zeros((len(seqs), len(seqs2)))
    
    for i in range(len(seqs)):
        for j in range(len(seqs2)):
            
            if i == j and flag:
                pi_matrix[i, j] = 0
            
            else:
                alignment = aligner.align(seqs[i], seqs2[j])[0]
                len_alignment = min(len(alignment[0].strip('-')), len(alignment[1].strip('-')))
                
                pi_matrix[i, j] =  alignment.score / len_alignment
            
    return pi_matrix


def PI_finder(matrix, flag=True):
    
    closest_genes = list()
        
    for row in matrix:
        closest_genes.append(max(row))
        
    if not flag:
        num_columns = len(matrix[0])
        for col in range(num_columns):
            column_values = [matrix[row][col] for row in range(len(matrix))]
            closest_genes.append(max(column_values))

    len_all = len(closest_genes)
    closest_genes_95 = [x for x in closest_genes if x >= 0.95]
    sim_95 = len(closest_genes_95) / len_all

    return np.mean(closest_genes), sim_95



def calculate_gene_statistics(vgene_path, vgene_path2=None):
    
    ## flag is True if self-Vgene

    tmp = pd.read_csv(vgene_path, sep='\t')
    seqs = list(tmp['Sequence'])
    
    if vgene_path2:
        tmp2 = pd.read_csv(vgene_path2, sep='\t')
        seqs2 = list(tmp2['Sequence'])
        flag = False
    else:
        seqs2 = seqs
        flag = True
        
    matrix = pi_matrix(seqs, seqs2)
    pi_mean, sim_95 = PI_finder(matrix, flag)
    
    return pi_mean, sim_95


def calculate(**kwargs):
    pi_mean, sim_95 = calculate_gene_statistics(**kwargs)
    return pi_mean, sim_95


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate average percent identity')
    
    parser.add_argument('-v', '--vgene', type=str, required=True, help='path to tsv file with genes')
    parser.add_argument('-v2', '--vgene2', type=str, required=False, help='path to second tsv file with genes')

    args = parser.parse_args()
    
    kwargs = {'vgene_path': args.vgene}  

    # optional
    if args.vgene2 is not None:
        kwargs['vgene_path2'] = args.vgene2
   
    # call
    pi_mean, sim_95 = calculate(**kwargs)
    print(pi_mean)
    print(sim_95)
    
