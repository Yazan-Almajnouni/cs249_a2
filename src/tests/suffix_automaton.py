class SuffixAutomaton:
    def __init__(self):
        # state 0 = the initial state
        self.next = [dict()]   # list of dict: state -> (char -> next_state)
        self.link = [-1]       # suffix‐link
        self.len  = [0]        # length of longest string in class
        self.last = 0          # the state for the whole string so far

    def extend(self, c):
        p = self.last
        cur = len(self.next)
        self.next.append({})
        self.len.append(self.len[p] + 1)
        self.link.append(0)
        while p >= 0 and c not in self.next[p]:
            self.next[p][c] = cur
            p = self.link[p]
        if p == -1:
            # no p => link to root
            self.link[cur] = 0
        else:
            q = self.next[p][c]
            if self.len[p] + 1 == self.len[q]:
                self.link[cur] = q
            else:
                # split q
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
    """
    reads: list of strings r_0,...,r_{n-1}
    min_ov: minimum overlap length
    returns overlaps[u][v] = max overlap length of suffix(read u)→prefix(read v)
    in O(N + a) time.
    """
    n = len(reads)
    # 1) pick n distinct end‐markers from the Unicode range
    #    (assume n < 0x10FFFF)
    seps = [chr(0x10FFFF - i) for i in range(n)]
    sep2id = {seps[i]: i for i in range(n)}

    # 2) build the big string S = r_0$0 r_1$1 ... r_{n-1}$(n-1)
    S = []
    for i, r in enumerate(reads):
        S.append(r)
        S.append(seps[i])
    S = "".join(S)

    # 3) build suffix-automaton of S in O(|S|)
    sa = SuffixAutomaton()
    for c in S:
        sa.extend(c)

    # 4) extract just the “separator” transitions for each state
    sep_trans = []
    for st in range(len(sa.next)):
        lst = []
        for c, t in sa.next[st].items():
            if c in sep2id:
                lst.append((c, sep2id[c]))
        sep_trans.append(lst)

    # 5) now walk each read once, recording overlaps
    overlaps = {u: {} for u in range(n)}
    for v, r in enumerate(reads):
        p = 0     # current SA state
        l = 0     # current match‐length (should be position+1)
        for c in r:
            # standard “feed c into SA” with suffix‐link fallback
            while p > 0 and c not in sa.next[p]:
                p = sa.link[p]
                l = sa.len[p]
            if c in sa.next[p]:
                p = sa.next[p][c]
                l += 1
            else:
                p = 0
                l = 0
            # if we have at least min_ov, look for outgoing $-edges
            if l >= min_ov:
                for sep_c, u in sep_trans[p]:
                    if u == v:
                        continue
                    prev = overlaps[u].get(v, 0)
                    if l > prev:
                        overlaps[u][v] = l

    return overlaps

strings = ["GCAADAGCA", "CACCGDDAC", "DACGCGA"]

print(find_overlaps_linear(strings, 4))