import igraph
import networkx as nx

from typing import Callable, Iterable

Walk = Callable[[int], Iterable[int]]


def longest_path(b: int, k: int, DS: Iterable[int]):
    """
        Returns the longest remaining path given a decycling set for the De Bruijn graph
            of base b and order k.

        Raises a networkx.NetworkXUnfeasible error if DS is not a decycling set
    """

    G = nx.DiGraph(
        igraph.Graph.De_Bruijn(b, k).to_networkx()
    )

    G.remove_nodes_from(DS)

    return nx.dag_longest_path_length(G) + 1


def record_decycling_set(b: int, k: int, DS: Iterable[int]):
    """
        Write out a decycling set as b-k-l.txt, where
            l is the length of the longest remaining path.
    """
    
    l = longest_path(b, k, DS)

    kmer = kmer_formatter(b, k)
    with open(f'{b}-{k}-{l}.txt', 'w') as out:
        out.write(
            '\n'.join(map(kmer, DS))
        )


def PCR(b: int, k: int):
    """
        Produces an iterator over the PCRs of the (b, k) De Bruijn graph.
        The PCR follows xS --> Sx, continuing until returning to XS.
    """
    
    kmer = kmer_formatter(b, k)

    seen = set()
    for x in range(b ** k):
        if x in seen:
            continue

        x = kmer(x)
        y, cycle = x, [int(x, b)]
        while (y := y[1:] + y[0]) != x:
            cycle.append(int(y, b))

        yield cycle
        seen.update(cycle)


def PCCR(b: int, k: int):
    """
        Produces an iterator over the PCCR cycles of the DBG with base 2 and order k.
        The PCCR cycles follows xS -> S(¬x), continuing until returning to xS.
    """

    kmer = kmer_formatter(b, k)

    swap = {  str(x): str(b - 1 - x) for x in range(b) }

    seen = set()
    for x in range(b ** k):
        if x in seen:
            continue

        x = kmer(x)
        y, cycle = x, [int(x, b)]
        while (y := y[1:] + swap[y[0]]) != x:
            cycle.append(int(y, b))

        yield cycle
        seen.update(cycle)

class Complementer:
    def __init__(self, b: int, k: int):
        self.b = b
        self.k = k
        self.bases = ''.join([str(x) for x in range(self.b)])
        self.comps = self.bases[::-1]
        self.swap = str.maketrans(self.bases, self.comps)
        self.d = self.b ** self.k

    def kmer(self, x: int):
        digits = []
        d = self.d
        for i in range(self.k):
            d //= self.b
            xi, x = divmod(x, d)
            digits.append(str(xi))

        return ''.join(digits)

    def reversal(self, x: int):
        return int(self.kmer(x)[::-1], self.b)

    def complement(self, x: int):
        return int(self.kmer(x).translate(self.swap), self.b)

    def rc(self, x: int):
        return int(self.kmer(x)[::-1].translate(self.swap), self.b)

def complement(x: int, b: int, k: int):
    kmer = kmer_formatter(b, k)
    bases = ''.join(range(b))
    swap = str.maketrans('01', '10')

    return int(kmer(x).translate(swap), 2)


def reversal(x: int, b: int, k: int):
    kmer = kmer_formatter(b, k)
    return int(kmer(x)[::-1], b)


def kmer_formatter(b: int, k: int):
    """
        Returns a string formatting function for base b kmers.
    """
    
    def kmer(x: int):
        d = b ** k

        digits = []
        for i in range(k):
            d //= b
            xi, x = divmod(x, d)
            digits.append(str(xi))

        return ''.join(digits)

    return kmer


def DBG_walk_maps(b: int, k: int):
    """
        Returns the In and Out maps for the base b, order k De Bruijn graph.
    """

    p = b ** (k - 1)
    s = b **  k

    def In(x):
        base = (x // b) % p
        return [
            base + i * p for i in range(b)
        ]

    def Out(x):
        base = (x * b) % s
        return [
            base + i for i in range(b)
        ]

    return In, Out
