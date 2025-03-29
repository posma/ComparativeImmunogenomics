import pandas as pd
import numpy as np
import argparse

def calculate_area_under_line(row, cols):
    areas = []
    for col in cols:
        y_values = row[col].values
        x_values = np.arange(len(y_values))
        area = np.trapz(y_values, x_values)
        areas.append(area)
    return areas

def calculate_area_between_lines(row, q_cols, s_cols):
    q_values = row[q_cols].values
    s_values = row[s_cols].values
    area = np.trapz(np.abs(q_values - s_values), np.arange(len(q_values)))
    return area

def calculate_haplotype_similarity(yass_path, q_contig_length, s_contig_length):
    
    itters = dict()
    q_contig_length = int(q_contig_length)
    s_contig_length = int(s_contig_length)
    tmp = pd.read_csv(yass_path, sep='\t')
    tmp['original_index'] = tmp.index
    final_indeces = list()

    # sorting and finding % covered
    q = [0] * q_contig_length
    s = [0] * s_contig_length
              
    tmp = tmp.sort_values(by=['q. length,'], ascending=False)
    tmp = tmp.reset_index(drop=True)
    tmp = tmp.reset_index()
    final_indeces += list(tmp['original_index'][:101])
        
    for i in range(101):
        q_start = tmp[tmp['index']==i]['# q. start,'].values[0]
        q_end = tmp[tmp['index']==i]['q. end,'].values[0] 
        q[q_start-1:q_end] = [1 for x in q[q_start-1:q_end]]
        itters[f'q_top{i}'] = [sum(q) / q_contig_length]

    tmp = tmp.sort_values(by=['s.length,'], ascending=False)
    tmp = tmp.drop(columns=['index'])
    tmp = tmp.reset_index(drop=True)
    tmp = tmp.reset_index()
    final_indeces += list(tmp['original_index'][:101])

    for i in range(101):
        s_start = tmp[tmp['index']==i]['s. start,'].values[0]
        s_end = tmp[tmp['index']==i]['s. end,'].values[0]
        s[s_start-1:s_end] = [1 for x in s[s_start-1:s_end]]
        itters[f's_top{i}'] = [sum(s) / s_contig_length]
        
    itters = pd.DataFrame(itters)
    q_columns = [f'q_top{i}' for i in range(101)]
    s_columns = [f's_top{i}' for i in range(101)]

    itters['q_area'] = itters.apply(lambda row: np.trapz(row[q_columns], dx=1), axis=1)
    itters['s_area'] = itters.apply(lambda row: np.trapz(row[s_columns], dx=1), axis=1)
    itters['dif_area'] = itters.apply(lambda row: calculate_area_between_lines(row, q_columns, s_columns), axis=1)

    itters['exp_area'] = itters[['q_area', 's_area']].mean(axis=1)
    itters['exp_area'] = itters['exp_area'] - itters['dif_area']
    
    return itters


def calculate(**kwargs):
    result = calculate_haplotype_similarity(**kwargs)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate haplotype similarity')
    
    parser.add_argument('-y', '--yass_csv', type=str, required=True, help='path to yass csv file')
    parser.add_argument('-q', '--q_contig_length', type=str, required=True, help='q contig length')
    parser.add_argument('-s', '--s_contig_length', type=str, required=True, help='s contig length')

    args = parser.parse_args()
    
    kwargs = {'yass_path': args.yass_csv,
             'q_contig_length': args.q_contig_length,
             's_contig_length': args.s_contig_length}  
   
    # call
    result = calculate(**kwargs)
    print(result['exp_area'].values[0])
    

