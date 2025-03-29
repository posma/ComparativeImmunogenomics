# Comparative analysis of mammalian immunoglobulin (IG) and T-cell receptor (TR) loci

## Scripts for computing IG/TR loci characteristics
All scripts are command line and notebook executable. Jupyter execution can be found in file *test/test.ipynb*. Command line executions are present below:

### Computing locus tandemness
*test/Repeats.py* <br />
usage: python3 Repeats.py -y {path/to/yass_csv/output} -c {int contig length}

### Computing gene statistics within the same locus
*test/PercentIdentity.py* <br />
usage: python3 PercentIdentity.py -v {path/to/vgene/df} 

### Computing haplotype similarity
*test/HapSimilarity.py* <br />
usage: python3 HapSimilarity.py -y {path/to/yass_csv/output} -q {q_contig_length} -s {s_contig_length}

### Computing gene statistics within a haplotype pair
*test/PercentIdentity.py* <br />
usage: python3 PercentIdentity.py -v {path/to/vgene/df} -v2 {path/to/vgene2/df}

## Data with IG/TR loci characteristics
### Curated IG/TR locus sequences
*data/loci/{IGH|IGK|IGL|TRA|TRB}*

### Annotated germline IG/TR genes
*data/Vgenes/{IGH|IGK|IGL|TRA|TRB}* <br />
contains only productive V genes

### IG/TR locus lengths and gene counts
Yana

### Pairwise locus similarities
Yana

### Summary characteristics of IG/TR loci
*data/summary.tsv* <br />
Contains columns: <br />
- R15_5: repetiteveness metric ~ fractions of an input sequence positions covered by repeats of length at least 15kbp at least 5 times
- tandem: tandemness category (non-tandem, semi-tandem, tandem)
- pi_mean: average percent identity 
- sim_95: average percentage of V genes that have at least 95% similarity to another V gene within the same locus 
- Simple_repeat: percentage of locus covered by simple repeats
- LINE/L1: percentage of locus covered by LINE/L1
- LTR/ERVL-MaLR: percentage of locus covered by LTR/ERVL-MaLR
- LTR/ERVL: percentage of locus covered by LTR/ERVL
- LTR/ERV1: percentage of locus covered by LTR/ERV1
-Low_complexity: percentage of locus covered by low complexity repeats

### Effective population sizes
Yana

### Summary characteristics of IG/TR haplotypes
*data/summary_haplotypes.tsv* <br />
Contains columns: <br />
- exp_area: haplotype similarity metrics
- pi_mean: average percent identity between haplotypes genes
- size2: number of pairs as the smallest cut of hierarchical tree
- pair12: number of pairs, where genes come from haplotype 1 AND haplotype 2
- pairPN: number of pairs, where genes are productive AND non-productive
- #genes: number of genes
- 12%pairs: percent of pairs containing haplotype 1 AND haplotype 2 (to # of pairs)
- 12%genes: percent of genes that form pair12 (to # of genes)
- PN%pairs: percent of pairs containing productive AND non-productive (to # of pairs)
- PN%genes: percent of genes that form pairPN (to # of genes)

## Scripts for computing correlations and association P-values 
*data/computations.ipynb*
