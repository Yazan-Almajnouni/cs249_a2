"""
De Bruijn Graph assembler (Task 1.1)
Takes FASTQ input, builds k-mer graph, finds Eulerian trails,
and writes contigs to FASTA.
"""

import argparse
import sys
from collections import defaultdict
from copy import deepcopy


def parse_fastq(fastq_files):
    """
    Yield sequences from one or more FASTQ files.
    """
    for fname in fastq_files:
        with open(fname) as fh:
            while True:
                name = fh.readline()
                if not name:
                    break
                seq = fh.readline().strip()
                fh.readline()  # plus line
                fh.readline()  # quality line
                yield seq


def build_debruijn_graph(seqs, k, min_count=1):
    """
    Build a De Bruijn graph from input sequences.
    Returns:
      adj: dict node -> list of successor nodes
      indeg, outdeg: dicts of in- and out-degrees
    """
    kmer_counts = defaultdict(int)
    for seq in seqs:
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            kmer_counts[kmer] += 1

    adj = defaultdict(list)
    indeg = defaultdict(int)
    outdeg = defaultdict(int)

    for kmer, cnt in kmer_counts.items():
        if cnt < min_count:
            continue
        prefix = kmer[:-1]
        suffix = kmer[1:]
        adj[prefix].append(suffix)
        outdeg[prefix] += 1
        indeg[suffix] += 1
        # ensure nodes exist
        indeg[prefix] += 0
        outdeg[suffix] += 0

    return adj, indeg, outdeg


def find_eulerian_trails(adj, indeg, outdeg):
    """
    Find all Eulerian trails in the graph (one per component).
    Returns list of contig strings.
    """
    # Collect all nodes: those with outgoing edges (keys of adj)
    # and those that appear as targets (values in adj).
    all_nodes = set(adj.keys())
    for vs in adj.values():
        all_nodes.update(vs)

    # Build a mutable copy of adjacency for every node,
    # defaulting to an empty list if no outgoing edges.
    graph = {u: list(adj.get(u, [])) for u in all_nodes}

    contigs = []

    def pick_start():
        # Prefer nodes where outdeg>indeg (starts of trails)
        for u in all_nodes:
            if graph[u] and outdeg[u] > indeg[u]:
                return u
        # Fallback: any node with remaining edges
        for u in all_nodes:
            if graph[u]:
                return u
        return None

    while True:
        start = pick_start()
        if start is None:
            break

        stack = [start]
        path = []
        while stack:
            v = stack[-1]
            if graph[v]:
                w = graph[v].pop()
                stack.append(w)
            else:
                path.append(stack.pop())
        path.reverse()

        if len(path) < 2:
            continue
        # Reconstruct contig: first node in full, then append last char of each
        seq = path[0]
        for node in path[1:]:
            seq += node[-1]
        contigs.append(seq)

    return contigs



def write_fasta(contigs, out_file):
    """
    Write contigs to a FASTA file.
    """
    with open(out_file, "w") as fh:
        for i, seq in enumerate(contigs, start=1):
            fh.write(f">contig_{i}\n{seq}\n")

def write_gfa(adj, k, gfa_file):
    """
    Write the de Bruijn graph in GFA1 format.
    Nodes are (k-1)-mers; edges link prefix->suffix with (k-2)-base overlap.
    """
    # collect all nodes
    nodes = set(adj.keys())
    for nbrs in adj.values():
        nodes.update(nbrs)

    with open(gfa_file, "w") as fh:
        # header
        fh.write("H\tVN:Z:1.0\n")
        # segments: S <node_id> <sequence>
        for node in sorted(nodes):
            fh.write(f"S\t{node}\t{node}\n")
        # links: L <from> + <to> + <overlap>CIGAR
        ov = k - 2
        if ov < 1:
            ov = 1
        for u in sorted(adj):
            seen = set()
            for v in adj[u]:
                if (u, v) in seen:
                    continue
                seen.add((u, v))
                fh.write(f"L\t{u}\t+\t{v}\t+\t{ov}M\n")


def main():
    p = argparse.ArgumentParser(
        description="De Bruijn Graph Assembler (Task 1.1)"
    )
    p.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Input FASTQ file(s)",
    )
    p.add_argument(
        "-k",
        "--kmer",
        type=int,
        required=True,
        help="k-mer size",
    )
    p.add_argument(
        "--min-count",
        type=int,
        default=1,
        help="Minimum k-mer count to include",
    )
    p.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output FASTA file for contigs",
    )
    p.add_argument(
        "--gfa",
        help="If given, write the DBG as GFA1 to this file",
    )
    args = p.parse_args()

    seqs = list(parse_fastq(args.input))
    if not seqs:
        sys.exit("No reads found in input FASTQ.")

    adj, indeg, outdeg = build_debruijn_graph(
        seqs, args.kmer, args.min_count
    )
    contigs = find_eulerian_trails(adj, indeg, outdeg)
    if not contigs:
        sys.exit("No contigs assembled (graph may be empty).")

    if args.gfa:
        write_gfa(adj, args.kmer, args.gfa)
        print(f"Wrote DBG graph in GFA format to {args.gfa}")

    write_fasta(contigs, args.output)
    print(f"Wrote {len(contigs)} contig(s) to {args.output}")


if __name__ == "__main__":
    main()