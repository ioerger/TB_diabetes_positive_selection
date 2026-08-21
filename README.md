# TB_diabetes_positive_selection

**Calculating pN/pS for a gene** 
input args: a fasta file with aligned sequences, and  a file with reference sequences for each ORF

example:
```
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
Rv3508 summary: codons=1901, totObsNS=44, totObsS=19, totSitesNS=12073, totSitesSS=5038, NS/S=2.396387, pN=0.003727, pS=0.003969, pN/pS=0.939022`
```
**Collapsing identical isolates for PAML**
input args: a fasta file with aligned sequences, and a Newick tree of those isolates

identical sequences are collapsed to unique haplotypes, the tree is pruned to those tips, and the files codeml needs are written to --outdir:
align.phy, tree.nwk
(tree.labelled.nwk is also written if the input tree has $1 / $2 labels on the two root clades)

example:

```
> python3 collapse_isolates.py --fasta Rv0006.seqs.noLC.txt --tree diabetes.DB_NDB12.CDS.split.noblen.bintree --outdir .
input    : 178 isolates, 838 codons
collapsed: 36 distinct sequences (haplotypes)  [labels kept separate: $1, $2]
tree     : pruned to 36 tips
group sizes (isolates per haplotype): [115, 16, 9, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
wrote align.phy, tree.nwk, tree.labelled.nwk to .
> cat tree.nwk
(s8,((s30,(s19,s4)),(((((s16,s20),((s5,s25),((s2,s18),(s36,s22)))),((s12,s35),(((s31,(s9,s3)),s17),(s28,s6)))),(((s34,s32),(s15,((s33,(s7,s14)),s1))),(s21,s29))),(((s26,((s13,s24),s11)),(s10,s23)),s27))));
> cat tree.labelled.nwk
(s8 $1,((s30,(s19,s4)),(((((s16,s20),((s5,s25),((s2,s18),(s36,s22)))),((s12,s35),(((s31,(s9,s3)),s17),(s28,s6)))),(((s34,s32),(s15,((s33,(s7,s14)),s1))),(s21,s29))),(((s26,((s13,s24),s11)),(s10,s23)),s27))) $1);
```
