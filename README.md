# TB_diabetes_positive_selection

**Calculating pN/pS for a gene** 

input args: a fasta file with aligned sequences for a given ORF, a file with reference sequences for each ORF, and the orfid

notes: 
 * instead of doing all pairwise comparisons between seqs, this script compares each seq to a reference sequence
 * input seqs should be aligned to reference sequence and all same length
 * the refseqs file should also be a fasta file, with the orfid as the first symbol in the header

`usage: python pNpS.py <seqs_fasta> --refseqs <fasta> --orfid XXXX`

example:
```
> python3 pNpS.py Rv0006.seqs.noLC.txt --refseqs H37Rv3_orfs.fna --orfid Rv0006
Rv0006: num seqs=178, seq len=2514, num codons=838
Rv0006: num seqs=178, seq len=2514, num codons=838
resnum	codon	aa	Nseqs	possS	possNS	obsS	obsNS	alleles
1	ATG	M	178	0	9	0	0	ATG(M):178
2	ACA	T	178	3	6	0	0	ACA(T):178
3	GAC	D	178	1	8	0	0	GAC(D):178
...
834	GAC	D	178	1	8	0	1	GAC(D):177,AAC(N):1
835	CAG	Q	178	1	8	0	0	CAG(Q):178
836	ACG	T	178	3	6	0	0	ACG(T):178
837	GGC	G	178	3	6	0	0	GGC(G):178
838	AAT	N	178	1	8	0	0	AAT(N):178
Rv0006 summary: codons=838, totObsNS=27, totObsS=11, totSitesNS=1877.67, totSitesS=636.33, NS/S=2.950760, pN=0.014904, pS=0.018828, pN/pS=0.791578
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
