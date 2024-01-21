def transposition_table(*pairs):
    X = ''.join(pairs)
    Y = ''.join(pair[::-1] for pair in pairs)

    return str.maketrans(X, Y)


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
