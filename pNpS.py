import sys
from codon_frequencies import *

# return list of (item,cnt) sorted by counts; most abundant first

def popularity(lst):
  hash = {}
  for x in lst:
    if x not in hash: hash[x] = 0
    hash[x] += 1
  data = [(hash[x],x) for x in hash.keys()]
  data.sort(reverse=True)
  data = [(y,x) for (x,y) in data]
  return data

# i is index of codon (0-based) - not a multiple of 3

def majority_codon(seqs,i):
  codons = [seq[3*i:3*(i+1)] for seq in seqs]
  pop = popularity(codons)
  return pop[0][0]

# all codons that differ by 1 nuc

def codon_neighbors(codon):
  alleles = []
  for i in range(3):
    for nuc in "AGCT":
      if codon[i]!=nuc:
        mut = [x for x in codon]
        mut[i] = nuc
        mut = ''.join(mut)
        alleles.append(mut)
  return alleles

def trans(codon): return codon_table.get(codon,'?')

####################################
# read refseqs

RefSeqs = {}
refHdrs,refSeqs = read_fasta(sys.argv[2])
for hdr,seq in zip(refHdrs,refSeqs):
  orfid = hdr[1:].split()[0] # strip off leading '>'
  RefSeqs[orfid] = seq
  
####################################
  
orfid = sys.argv[1]
if '/' in orfid: orfid = orfid[orfid.rfind('/')+1:] # extract just the base name of the file
orfid = orfid[:orfid.find('.')]
refseq = RefSeqs[orfid]

Hdrs,Seqs = read_fasta(sys.argv[1])
Nseqs = len(Seqs)
Ncodons = int(len(Seqs[0])/3)
print("%s: num seqs=%s, seq len=%s, num codons=%s" % (orfid,Nseqs,len(Seqs[0]),Ncodons))

goodcodons = {}
for codon,aa in codon_table.items(): # defined in codon_frequencies.py
  if 'N' not in codon: goodcodons[codon] = 1
#print(len(goodcodons)) # 64

totSitesS,totSitesNS = 0,0 # sum of possible alleles
totObsS,totObsNS = 0,0 # observed alleles

for i in range(Ncodons):
  #ref = majority_codon(Seqs,i)
  ref = refseq[3*i:3*(i+1)]
  if ref not in goodcodons: print("warning: skipping aa %d, ref codon %s not recognized" % (i+1,ref)); continue
  refaa = trans(ref)
  obsS,obsNS,locS,locNS = 0,0,0,0

  # count all possible single-nuc changes
  SNVs = codon_neighbors(ref)
  for allele in SNVs:
    mutaa = trans(allele)
    if mutaa==refaa: locS += 1
    else: locNS += 1

  alleles = {}
  good = 0
  for j in range(Nseqs):
    codon = Seqs[j][3*i:3*(i+1)]
    if codon not in goodcodons: continue # skip dashes, N's, etc
    good += 1
    if codon not in alleles: alleles[codon] = 0
    alleles[codon] += 1

  for mut,cnt in alleles.items(): # only good codons, including self
    if mut==ref: continue 
    mutaa = trans(mut)
    if mut in SNVs:
      if mutaa==refaa: obsS += 1
      else: obsNS += 1
    else: obsNS += 1; locNS += 1; print("warning: multiple changes in aa %s: %s -> %s" % (i+1,ref,mut)) # for codons with multiple nuc substitutions

  vals = [i+1,ref,refaa,good,locS,locNS,obsS,obsNS]
  print('\t'.join([str(x) for x in vals]),alleles)

  totObsS += obsS
  totObsNS += obsNS
  totSitesS += locS
  totSitesNS += locNS

pN = (totObsNS+1)/float(totSitesNS+1) # use pseudocounts to avoid div-by-zero
pS = (totObsS+1)/float(totSitesS+1)
pNpS = pN/pS

print("%s summary: codons=%s, totObsNS=%s, totObsS=%s, totSitesNS=%s, totSitesSS=%s, NS/S=%0.6f, pN=%0.6f, pS=%0.6f, pN/pS=%0.6f" % (orfid,Ncodons,totObsNS,totObsS,totSitesNS,totSitesS,totSitesNS/float(totSitesS),pN,pS,pNpS))
