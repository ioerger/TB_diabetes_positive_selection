# TB_diabetes_positive_selection

**Calculating pN/pS for a gene** 
input args: a fasta file with aligned sequences, and  a file with reference sequences for each ORF

example:
> python3 pNpS.py Rv0006.seqs.noLC.txt H37Rv3_orfs.fna
Rv3508: num seqs=178, seq len=5703, num codons=1901
1       ATG     M       178     0       9       0       0 {'ATG': 178}
2       TCG     S       178     3       6       0       0 {'TCG': 178}
3       TTC     F       178     1       8       0       0 {'TTC': 178}
...
1898    ACC     T       167     3       6       0       1 {'GCC': 5, 'ACC': 162}
1899    GAC     D       173     1       8       0       0 {'GAC': 173}
1900    GGC     G       173     3       6       0       0 {'GGC': 173}
1901    AGC     S       176     1       8       0       0 {'AGC': 176}
Rv3508 summary: codons=1901, totObsNS=44, totObsS=19, totSitesNS=12073, totSitesSS=5038, NS/S=2.396387, pN=0.003727, pS=0.003969, pN/pS=0.939022
