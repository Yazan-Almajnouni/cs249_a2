#!/usr/bin/env python3
"""
OLC assembler using a suffix‐automaton overlap finder + string‐graph reduction.

1) Parse FASTQ reads into `reads`.
2) Build a suffix automaton over S = r0$0 r1$1 … r{n-1}${n-1}.
3) Find ALL overlaps u→v of length ≥ min_overlap.
4) Build the full overlap graph (all edges).
5) Transitively reduce the graph: remove any direct edge u→w if there
   is a two‐edge path u→v→w.
6) Extract non‐branching unitigs + cycles from the reduced graph.
7) Layout reads along each path and compute a column‐wise consensus.
8) Write consensus contigs to FASTA.
"""
import argparse
import sys
from collections import defaultdict

class SuffixAutomaton:
    def __init__(self):
        self.next = [dict()]
        self.link = [-1]
        self.len  = [0]
        self.last = 0

    def extend(self, c):
        p, cur = self.last, len(self.next)
        self.next.append({})
        self.len.append(self.len[p] + 1)
        self.link.append(0)
        while p >= 0 and c not in self.next[p]:
            self.next[p][c] = cur
            p = self.link[p]
        if p == -1:
            self.link[cur] = 0
        else:
            q = self.next[p][c]
            if self.len[p] + 1 == self.len[q]:
                self.link[cur] = q
            else:
                clone = len(self.next)
                self.next.append(self.next[q].copy())
                self.len.append(self.len[p] + 1)
                self.link.append(self.link[q])
                while p >= 0 and self.next[p].get(c) == q:
                    self.next[p][c] = clone
                    p = self.link[p]
                self.link[q] = self.link[cur] = clone
        self.last = cur

def find_overlaps_linear(reads, min_ov):
    n = len(reads)
    seps = [chr(0x10FFFF - i) for i in range(n)]
    sep2id = {seps[i]: i for i in range(n)}
    # build concatenated string
    S = "".join(r + seps[i] for i, r in enumerate(reads))
    sa = SuffixAutomaton()
    for c in S:
        sa.extend(c)
    # collect separator‐transitions
    sep_trans = []
    for st in range(len(sa.next)):
        lst = []
        for c, _ in sa.next[st].items():
            if c in sep2id:
                lst.append(sep2id[c])
        sep_trans.append(lst)
    # walk each read
    overlaps = {u: {} for u in range(n)}
    for u, r in enumerate(reads):
        p = l = 0
        for c in r:
            while p > 0 and c not in sa.next[p]:
                p = sa.link[p]
                l = sa.len[p]
            if c in sa.next[p]:
                p = sa.next[p][c]
                l += 1
            else:
                p = l = 0
            if l >= min_ov:
                for v in sep_trans[p]:
                    if v == u:
                        continue
                    prev = overlaps[u].get(v, 0)
                    if l > prev:
                        overlaps[u][v] = l
    return overlaps

def build_full_graph(overlaps):
    adj    = defaultdict(list)
    indeg  = defaultdict(int)
    outdeg = defaultdict(int)
    for u, hits in overlaps.items():
        for v, ov in hits.items():
            adj[u].append((v, ov))
            outdeg[u] += 1
            indeg[v]  += 1
    for u in overlaps:
        adj.setdefault(u, [])
        indeg[u]  += 0
        outdeg[u] += 0
    return adj, indeg, outdeg

def transitively_reduce(adj):
    """
    Remove any edge u->w if there exists v such that u->v and v->w.
    """
    new_adj = {}
    # For each u, collect direct successors as a dict
    i = 0
    for u, outs in adj.items():
        direct = {v: ov for v,ov in outs}
        to_remove = set()
        for v, ov_uv in outs:
            for w, ov_vw in adj.get(v, ()):
                if w != u and w in direct:
                    to_remove.add(w)
        # build reduced list
        new_adj[u] = [(v, direct[v]) for v in direct if v not in to_remove]
        i += 1
        if i % 100 == 0:
            print(f"finished {i}")
    return new_adj

def recompute_degrees(adj):
    indeg  = defaultdict(int)
    outdeg = defaultdict(int)
    for u, outs in adj.items():
        outdeg[u] = len(outs)
        for v, _ in outs:
            indeg[v] += 1
    for u in adj:
        indeg[u]  += 0
        outdeg[u] += 0
    return indeg, outdeg

def extract_contig_paths(adj, indeg, outdeg):
    used = set()
    contigs = []
    def walk(u, v, ov):
        path = [u, v]
        ovs  = [ov]
        used.add((u,v))
        cur = v
        while indeg[cur] == 1 and outdeg[cur] == 1:
            nxt, nov = adj[cur][0]
            if (cur, nxt) in used:
                break
            path.append(nxt)
            ovs.append(nov)
            used.add((cur, nxt))
            cur = nxt
        return path, ovs
    # starts
    for u, outs in adj.items():
        if outdeg[u] and indeg[u] != 1:
            for v, ov in outs:
                if (u,v) not in used:
                    contigs.append(walk(u,v,ov))
    # cycles
    for u, outs in adj.items():
        for v, ov in outs:
            if (u,v) not in used:
                contigs.append(walk(u,v,ov))
    return contigs

def compute_consensus_for_path(path, ovs, reads):
    length = len(reads[path[0]])
    for ov, rid in zip(ovs, path[1:]):
        length += len(reads[rid]) - ov
    cov = [defaultdict(int) for _ in range(length)]
    pos = 0
    for i, rid in enumerate(path):
        seq = reads[rid]
        offset = 0 if i==0 else pos - ovs[i-1]
        for j, b in enumerate(seq):
            cov[offset+j][b] += 1
        pos = offset + len(seq)
    cons = []
    for d in cov:
        if not d:
            cons.append("N")
        else:
            # pick most frequent base
            base = max(sorted(d.keys()), key=lambda x: d[x])
            cons.append(base)
    return "".join(cons)

def write_fasta_consensus(contig_paths, reads, out_fn):
    with open(out_fn, "w") as fh:
        for i, (path, ovs) in enumerate(contig_paths, 1):
            seq = compute_consensus_for_path(path, ovs, reads)
            fh.write(f">contig_{i}\n{seq}\n")

def parse_fastq(fqs):
    for fn in fqs:
        with open(fn) as fh:
            while True:
                hdr = fh.readline()
                if not hdr:
                    break
                seq = fh.readline().strip()
                fh.readline()
                fh.readline()
                yield seq

def main():
    p = argparse.ArgumentParser(
        description="OLC assembler with suffix‐automaton + string‐graph"
    )
    p.add_argument("-i","--input", nargs="+", required=True,
                   help="Input FASTQ file(s)")
    p.add_argument("-n","--min-overlap", type=int, required=True,
                   help="Minimum suffix→prefix overlap")
    p.add_argument("-o","--output", required=True,
                   help="Output FASTA for consensus contigs")
    args = p.parse_args()

    reads = list(parse_fastq(args.input))
    if not reads:
        sys.exit("No reads found in input FASTQ.")

    print(f"number of reads = {len(reads)}")

    overlaps = find_overlaps_linear(reads, args.min_overlap)
    full_adj, indeg, outdeg = build_full_graph(overlaps)

    # 5) string‐graph (transitive) reduction
    sg_adj = transitively_reduce(full_adj)

    # now recompute degrees so we can skip dead‐ends in the 3-step pass
    indeg, outdeg = recompute_degrees(sg_adj)

    # 6) unitig extraction
    contig_paths = extract_contig_paths(sg_adj, indeg, outdeg)
    if not contig_paths:
        sys.exit("No contigs assembled (no overlaps ≥ min-overlap).")

    # 7) consensus & 8) output
    write_fasta_consensus(contig_paths, reads, args.output)
    print(f"Wrote {len(contig_paths)} contig(s) to {args.output}")

if __name__ == "__main__":
    main()

