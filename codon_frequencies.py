import sys

def read_fasta(filename):
  headers,seqs = [],[]
  seq = ""
  for line in open(filename):
    line = line.rstrip()
    if len(line)==0: continue
    if line[0]==">":
      headers.append(line)
      if seq!="": seqs.append(seq)
      seq = ""
    else: seq += line
  seqs.append(seq)
  return headers,seqs

codon_table = { 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
             'CTA': 'L', 'CTG': 'L', 'CTN': 'L', 'TGG': 'W',
             'TAA': '*', 'TAG': '*', 'TGA': '*', 'ATG': 'M',
             'TTT': 'F', 'TTC': 'F', 'TAT': 'Y', 'TAC': 'Y',
             'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 
             'TCN': 'S', 'AGT': 'S', 'AGC': 'S', 'CCT': 'P', 
             'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CCN': 'P',
             'TGT': 'C', 'TGC': 'C', 'CAT': 'H', 'CAC': 'H',
             'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
             'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
             'CGN': 'R', 'AGA': 'R', 'AGG': 'R', 'ATT': 'I', 
             'ATC': 'I', 'ATA': 'I', 'AAA': 'K', 'AAG': 'K',
             'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
             'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
             'ACN': 'T', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
             'GTG': 'V', 'GTN': 'V', 'GCT': 'A', 'GCC': 'A',
             'GCA': 'A', 'GCG': 'A', 'GCN': 'A', 'GGT': 'G', 
             'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GGN': 'G',
             'TAN': 'X', 'TTN': 'X', 'TGN': 'X', 'CAN': 'X', 
             'ATN': 'X', 'AAN': 'X', 'GAN': 'X', 'AGN': 'X',
             'ANA': 'X', 'ANT': 'X', 'ANG': 'X', 'ANC': 'X', 
             'TNA': 'X', 'TNT': 'X', 'TNG': 'X', 'TNC': 'X', 
             'GNA': 'X', 'GNT': 'X', 'GNG': 'X', 'GNC': 'X', 
             'CNA': 'X', 'CNT': 'X', 'CNG': 'X', 'CNC': 'X', 
             'NAA': 'X', 'NAT': 'X', 'NAG': 'X', 'NAC': 'X', 
             'NTA': 'X', 'NTT': 'X', 'NTG': 'X', 'NTC': 'X', 
             'NGA': 'X', 'NGT': 'X', 'NGG': 'X', 'NGC': 'X', 
             'NCA': 'X', 'NCT': 'X', 'NCG': 'X', 'NCC': 'X', 
             'NNA': 'X', 'NNT': 'X', 'NNG': 'X', 'NNC': 'X', 
             'ANN': 'X', 'TNN': 'X', 'GNN': 'X', 'NNC': 'X', 
             'NAN': 'X', 'NTN': 'X', 'NGN': 'X', 'NCN': 'X', 
             'NNN': 'X'}

Codons,Translation = [],[]
for codon,aa in codon_table.items():
  if 'N' not in codon and aa in "ACDEFGHIKLMNPQRSTVWY": Codons.append(codon); Translation.append(aa)

##########################

if __name__=="__main__":
  seqfile = sys.argv[1]
  H,S = read_fasta(seqfile)

  Nseqs = len(H)
  Nsites = int(len(S[0])/3) # assume it is a multiple, and all seqs the same legnth

  Counts = {}
  for i in range(Nseqs):
    for j in range(Nsites):
      codon = S[i][j*3:(j+1)*3]
      if codon in Codons:
        if codon not in Counts: Counts[codon] = 0
        Counts[codon] += 1

  tot = sum(Counts.values())
  sys.stderr.write("%s codons\n" % tot)
  for codon,count in Counts.items():
    print(codon,count,"%0.6f" % (count/float(tot)))
