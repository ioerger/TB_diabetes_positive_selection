#!/usr/bin/env python3
import argparse
import re
import sys
from pathlib import Path

sys.setrecursionlimit(100000)   # deep (caterpillar) trees can recurse a few hundred levels


####################################
# read fasta -> {isolate_id: sequence}

def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []   # take the id, ignore any description after a space
        elif line:
            buf.append(line.upper())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


####################################
# minimal Newick tree handling
# topology only: codeml re-estimates branch lengths, so input lengths are skipped
# #1/#2 leaf/branch tags are stripped; $1/$2 on root subtrees are kept as node names

def parse_newick(text):
    s = re.sub(r"#\d+", "", text).strip()          # drop #N branch labels; $N names are kept
    if s.endswith(";"):
        s = s[:-1]
    pos = 0

    def clade():
        nonlocal pos
        node = {"name": None, "children": []}
        if pos < len(s) and s[pos] == "(":         # an internal node: (child, child, ...)
            pos += 1
            while True:
                node["children"].append(clade())
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
        name = []                                  # a leaf name (or internal label) up to :,()
        while pos < len(s) and s[pos] not in ":,()":
            name.append(s[pos])
            pos += 1
        if name:
            node["name"] = "".join(name).strip()
        if pos < len(s) and s[pos] == ":":         # skip ":branchlength"
            pos += 1
            while pos < len(s) and s[pos] not in ",()":
                pos += 1
        return node

    return clade()

# set of leaf names under a node

def leaves(node):
    if not node["children"]:
        return {node["name"]} if node["name"] else set()
    out = set()
    for c in node["children"]:
        out |= leaves(c)
    return out

# {leaf: label} for a tree shaped ( <subtree>$1 , <subtree>$2 )
# each root child's name ($1 / $2) is copied to every leaf under it

def split_labels(root):
    out = {}
    for child in root["children"]:
        lab = child.get("name")
        for leaf in leaves(child):
            out[leaf] = lab
    return out

# copy of the tree containing only the leaves in `keep`
# single-child internal nodes are collapsed

def prune(node, keep):
    if not node["children"]:
        return {"name": node["name"], "children": []} if node["name"] in keep else None
    kids = [p for p in (prune(c, keep) for c in node["children"]) if p is not None]
    if not kids:
        return None
    if len(kids) == 1:          # a degree-2 node after pruning -> splice it out
        return kids[0]
    return {"name": None, "children": kids}

# make the tree strictly binary
# codeml refuses a node with >2 children ('too many daughter nodes, raise MAXNSONS')
# split into a nested chain; codeml re-estimates the new branches to ~0

def resolve_multifurcations(node):
    kids = [resolve_multifurcations(c) for c in node["children"]]
    if len(kids) <= 2:
        return {"name": node["name"], "children": kids}
    right = kids[-1]
    for c in reversed(kids[1:-1]):
        right = {"name": None, "children": [c, right]}
    return {"name": node["name"], "children": [kids[0], right]}

# topology-only Newick (no branch lengths), leaf names replaced via the `alias` map

def to_newick(node, alias):
    def rec(n):
        if not n["children"]:
            return alias[n["name"]]
        return "(" + ",".join(rec(c) for c in n["children"]) + ")"
    return rec(node) + ";"

# like to_newick, but restore $1/$2 on each root subtree from the representatives under it
# use this tree with codeml model=2; model=0 uses plain to_newick (no labels)

def labelled_clades_newick(root, alias, rep_label):
    def rec(n):
        if not n["children"]:
            return alias[n["name"]]
        return "(" + ",".join(rec(c) for c in n["children"]) + ")"
    parts = []
    for child in root["children"]:
        lab = next((rep_label[m] for m in leaves(child) if rep_label.get(m)), "")
        parts.append(rec(child) + (" " + lab if lab else ""))
    return "(" + ",".join(parts) + ");"


####################################
# collapse identical sequences -> haplotypes -> s1..sN aliases
# group by sequence (or (sequence, label) if label_of is set, so $1 and $2 stay separate)
# most abundant first; alias s1 is the most shared sequence
# rep = alphabetically-first member that is also in the tree

def collapse(seqs, tree_leaves, label_of=None):
    by_key = {}
    for name, s in seqs.items():
        if s:                                       # ignore empty records
            key = (s, label_of.get(name)) if label_of else s
            by_key.setdefault(key, []).append(name)

    haps = []
    # order: largest group first, then by representative id so the numbering is reproducible
    for key, members in sorted(by_key.items(), key=lambda kv: (-len(kv[1]), sorted(kv[1])[0])):
        seq = key[0] if label_of else key
        in_tree = [m for m in sorted(members) if m in tree_leaves]
        if not in_tree:                             # no isolate of this haplotype is in the tree
            continue
        haps.append({"rep": in_tree[0], "seq": seq, "members": sorted(members)})
    for i, h in enumerate(haps, 1):
        h["alias"] = f"s{i}"
    return haps


####################################
# write the output files

def write_outputs(haps, tree_root, outdir, label_of=None):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    alias_of = {h["rep"]: h["alias"] for h in haps}
    L = len(haps[0]["seq"])

    # align.phy -- PHYLIP: header "<n> <ncols>", then "<alias padded to 10><sequence>"
    with open(outdir / "align.phy", "w") as fh:
        fh.write(f" {len(haps)} {L}\n")
        for h in haps:
            fh.write(f"{h['alias']:<10}{h['seq']}\n")

    # tree.nwk -- prune the input tree to the representatives, force binary, relabel to aliases
    pruned = prune(tree_root, set(alias_of))
    pruned = resolve_multifurcations(pruned)
    (outdir / "tree.nwk").write_text(to_newick(pruned, alias_of) + "\n")

    # tree.labelled.nwk -- same pruned tree with $1/$2 restored (for codeml model=2)
    if label_of and len(pruned["children"]) == 2:
        rep_label = {h["rep"]: label_of.get(h["rep"], "") for h in haps}
        (outdir / "tree.labelled.nwk").write_text(labelled_clades_newick(pruned, alias_of, rep_label) + "\n")

    return len(haps), L // 3, len(leaves(pruned))


####################################

def main():
    ap = argparse.ArgumentParser(description="Collapse identical isolate sequences for codeml.")
    ap.add_argument("--fasta", default="Rv0290.seqs.noLC.txt", help="per-isolate codon FASTA")
    ap.add_argument("--tree", default="diabetes.CDS.DB.bintree", help="isolate tree (Newick)")
    ap.add_argument("--outdir", default=".", help="where to write the outputs")
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    tree_root = parse_newick(Path(args.tree).read_text())
    tree_leaves = leaves(tree_root)

    lengths = {len(s) for s in seqs.values() if s}
    if len(lengths) != 1:
        sys.exit(f"sequences are not all the same length: {sorted(lengths)}")

    # $-labelled split tree -> keep different labels apart even when the sequence matches
    label_of = split_labels(tree_root)
    labelled = any(label_of.values())
    haps = collapse(seqs, tree_leaves, label_of=label_of if labelled else None)
    n_hap, n_codon, n_tip = write_outputs(haps, tree_root, args.outdir, label_of if labelled else None)

    print(f"input    : {len(seqs)} isolates, {n_codon} codons")
    print(f"collapsed: {n_hap} distinct sequences (haplotypes)"
          + (f"  [labels kept separate: {', '.join(sorted(set(label_of.values())))}]" if labelled else ""))
    print(f"tree     : pruned to {n_tip} tips")
    print("group sizes (isolates per haplotype):", [len(h["members"]) for h in haps])
    print(f"wrote align.phy, tree.nwk{', tree.labelled.nwk' if labelled else ''} to {args.outdir}")


if __name__ == "__main__":
    main()
