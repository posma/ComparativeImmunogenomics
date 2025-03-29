import pandas as pd
import argparse

def count_repetitiveness(yass_csv, contig_length, rep=15000, coverage=5, normalization=False, Vstart=0, Vend=0):

    to_int = [int(arg) for arg in [contig_length, rep, coverage, Vstart, Vend]]
    contig_length, rep, coverage, Vstart, Vend = to_int

    # get repeats
    df = pd.read_csv(yass_csv, sep = '\t')
    df = df[df['q. length,']>=int(rep)]
    df = df[df['s.length,']>=int(rep)]
    df = df[df['q. length,']!=contig_length]
    df = df[df['s.length,']!=contig_length]
    q_starts = list(df['# q. start,'])
    q_ends = list(df['q. end,'])
    s_starts = list(df['s. start,'])
    s_ends = list(df['s. end,'])

    starts = q_starts + s_starts
    ends = q_ends + s_ends

    # find array
    array = [0] * contig_length
    for i in range(len(starts)):
        start = starts[i]
        end = ends[i]

        array[start-1:end] = [x + 1 for x in array[start-1:end]]

    # find rep more than coverage

    if normalization:
        count = sum(1 for i in array[int(Vstart-1):int(Vend)] if coverage <= i)
        return count/(Vend-Vstart)
        
    else:
        count = sum(1 for i in array if coverage <= i)
        return count/len(array)



def calculate(**kwargs):
    result = count_repetitiveness(**kwargs)
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate tandomness')
    
    parser.add_argument('-y', '--yass_csv', type=str, required=True, help='path to yass csv file')
    parser.add_argument('-c', '--contig_length', type=str, required=True, help='contig length')
    parser.add_argument('-L', '--L_argument', type=str, required=False, help='smallest length of repeat')
    parser.add_argument('-N', '--N_argument', type=str, required=False, help='at least n repeats')
    parser.add_argument('-n', '--normalization', type=str, required=False, help='if normilization is not on the contig length')
    parser.add_argument('-s', '--start', type=str, required=False, help='start of new normalization array')
    parser.add_argument('-e', '--end', type=str, required=False, help='end of new normalization array')

    args = parser.parse_args()
    
    kwargs = {'yass_csv': args.yass_csv,
             'contig_length': args.contig_length}  

    # optional
    if args.L_argument is not None:
        kwargs['reps'] = args.L_argument
    if args.N_argument is not None:
        kwargs['coverage'] = args.N_argument
    
    if args.normalization is not None:
        kwargs['normalization'] = args.normalization
    if args.start is not None:
        kwargs['Vstart'] = args.start
    if args.end is not None:
        kwargs['Vend'] = args.end
   
    # call
    result = calculate(**kwargs)
    print(result)
    

