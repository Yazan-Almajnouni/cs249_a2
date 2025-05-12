#!/usr/bin/env python3
"""
Overlap–Layout–Consensus assembler (Task 1.2)

Takes FASTQ input, finds all‐vs‐all overlaps ≥ n, builds an overlap graph,
extracts non‐branching paths (unitigs), lays out reads, and writes contigs.
"""

import argparse
import sys
from collections import defaultdict
from collections import Counter


def parse_fastq(fastq_files):
    """Yield (read_id, seq) from one or more FASTQ files."""
    rid = 0
    for fn in fastq_files:
        with open(fn) as fh:
            while True:
                hdr = fh.readline()
                if not hdr:
                    break
                seq = fh.readline().strip()
                fh.readline()  # '+'
                fh.readline()  # quals
                yield rid, seq
                rid += 1


def overlap(a, b, min_len):
    """
    Return length of longest suffix of a matching prefix of b
    with length ≥ min_len, else 0.
    """
    max_l = min(len(a), len(b))
    for l in range(max_l, min_len - 1, -1):
        if a.endswith(b[:l]):
            return l
    return 0


def build_overlap_graph(reads, min_ov):
    """
    Build directed overlap graph. Returns
      adj: u -> list of (v, ov_len)
      indeg, outdeg maps
    """
    print("building overlap graph")
    adj = defaultdict(list)
    indeg = defaultdict(int)
    outdeg = defaultdict(int)
    ids = list(reads.keys())
    num_reads = len(ids)
    print(f"{num_reads = }")
    for i, u in enumerate(ids):
        for v in ids:
            if u == v:
                continue
            ov = overlap(reads[u], reads[v], min_ov)
            if ov:
                adj[u].append((v, ov))
                outdeg[u] += 1
                indeg[v] += 1

        if (i + 1) % 100 == 0:
            print(f"finished {i + 1} / {num_reads}")
    # ensure all nodes present
    for u in ids:
        indeg[u] += 0
        outdeg[u] += 0
    return adj, indeg, outdeg


def extract_contigs(adj, indeg, outdeg, reads):
    """
    Find non‐branching paths and simple cycles, then compute a
    column‐wise consensus over each laid‐out read path.
    Returns list of consensus contig sequences.
    """
    print("extracting contigs")
    used = set()   # edges we’ve consumed
    contigs = []   # list of (path, ovs) as before

    def walk(u, v, ov):
        path = [u, v]
        ovs = [ov]
        used.add((u, v))
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

    # 1) non‐branching starts
    for u, outs in adj.items():
        if outdeg[u] and indeg[u] != 1:
            for v, ov in outs:
                if (u, v) not in used:
                    contigs.append(walk(u, v, ov))

    # 2) leftover cycles
    for u, outs in adj.items():
        for v, ov in outs:
            if (u, v) not in used:
                contigs.append(walk(u, v, ov))

    # helper to pick majority base (tie → lexicographically min)
    def majority_vote(bases):
        cnt = Counter(bases)
        maxc = max(cnt.values())
        # all bases with top count
        winners = [b for b, c in cnt.items() if c == maxc]
        return min(winners)

    # now build a consensus string for each (path, ovs)
    seqs = []
    for path, ovs in contigs:
        # initialize columns with the first read
        first = reads[path[0]]
        columns = [[b] for b in first]
        pos = len(columns)

        # for each subsequent read, merge the overlap
        for rid, ov in zip(path[1:], ovs):
            seq = reads[rid]
            start = pos - ov
            # absorb overlap columns
            for i in range(ov):
                columns[start + i].append(seq[i])
            # append the non‐overlapping suffix as new columns
            for b in seq[ov:]:
                columns.append([b])
            pos = len(columns)

        # build consensus by majority vote on each column
        cons = "".join(majority_vote(col) for col in columns)
        seqs.append(cons)

    return seqs


def write_fasta(contigs, out_fn):
    """Write a list of sequences to FASTA."""
    with open(out_fn, "w") as fh:
        for i, seq in enumerate(contigs, 1):
            fh.write(f">contig_{i}\n{seq}\n")


def main():
    p = argparse.ArgumentParser(
        description="OLC assembler (Task 1.2)"
    )
    p.add_argument(
        "-i", "--input", nargs="+", required=True,
        help="Input FASTQ file(s)"
    )
    p.add_argument(
        "-n", "--min-overlap", type=int, required=True,
        help="Minimum overlap length"
    )
    p.add_argument(
        "-o", "--output", required=True,
        help="Output FASTA for contigs"
    )
    args = p.parse_args()

    # load reads
    reads = {}
    for rid, seq in parse_fastq(args.input):
        reads[rid] = seq
    if not reads:
        sys.exit("No reads found in input FASTQ.")

    # build graph
    adj, indeg, outdeg = build_overlap_graph(
        reads, args.min_overlap
    )
    # extract contigs
    contigs = extract_contigs(adj, indeg, outdeg, reads)
    if not contigs:
        sys.exit("No contigs assembled (no overlaps ≥ min-overlap).")

    write_fasta(contigs, args.output)
    print(f"Wrote {len(contigs)} contig(s) to {args.output}")


if __name__ == "__main__":
    main()

