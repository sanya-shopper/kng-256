"""Exact inverses of SHA-256's diffusion maps Sigma_0 and Sigma_1.

Companion script to "SHA-256 as a Mathematical Object" (appendix).
Words are 32-bit integers; bit i is the coefficient of x^i, so that
ROTR^r is an F_2-linear map of the word space F_2^32.

The script (1) builds each map's 32x32 matrix over F_2, (2) inverts it
by Gaussian elimination, (3) reads the inverse off as an xor of
rotations -- the inverse of a circulant is circulant -- and (4) verifies
everything functionally on random words, including the order facts
Sigma_0^8 = ROTR^8 and Sigma_1^16 = id proved in the text.
"""
import random

W = 32
MASK = (1 << W) - 1


def rotr(x, r):
    r %= W
    return ((x >> r) | (x << (W - r))) & MASK


def Sigma0(x):
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)


def Sigma1(x):
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)


def matrix_of(f):
    """Column i (as a bitmask of rows) is f applied to the basis word 2^i."""
    return [f(1 << i) for i in range(W)]


def invert_matrix(cols):
    """Invert an F_2 matrix given by column bitmasks; None if singular."""
    rows = []
    for r in range(W):  # rows of [M | I] as bitmasks over columns
        m = sum(((cols[c] >> r) & 1) << c for c in range(W))
        rows.append([m, 1 << r])
    piv = 0
    for c in range(W):
        sel = next((r for r in range(piv, W) if (rows[r][0] >> c) & 1), None)
        if sel is None:
            return None
        rows[piv], rows[sel] = rows[sel], rows[piv]
        for r in range(W):
            if r != piv and (rows[r][0] >> c) & 1:
                rows[r][0] ^= rows[piv][0]
                rows[r][1] ^= rows[piv][1]
        piv += 1
    return [sum(((rows[r][1] >> c) & 1) << r for r in range(W))
            for c in range(W)]


def inverse_rotation_set(f):
    """Rotation offsets S with f^{-1} = XOR of ROTR^r over r in S."""
    inv = invert_matrix(matrix_of(f))
    if inv is None:
        raise ValueError("map is not invertible")
    y = inv[0]  # image of the basis word 1; ROTR^r(1) sets bit (32-r) % 32
    return sorted((32 - b) % 32 for b in range(W) if (y >> b) & 1)


def xor_of_rotations(S):
    def f(x):
        acc = 0
        for r in S:
            acc ^= rotr(x, r)
        return acc
    return f


def iterate(f, x, n):
    for _ in range(n):
        x = f(x)
    return x


if __name__ == "__main__":
    random.seed(20260716)
    words = [0, 1, MASK] + [random.getrandbits(W) for _ in range(5000)]
    for name, f, order in [("Sigma_0", Sigma0, 32), ("Sigma_1", Sigma1, 16)]:
        S = inverse_rotation_set(f)
        finv = xor_of_rotations(S)
        assert all(finv(f(x)) == x and f(finv(x)) == x for x in words)
        assert all(iterate(f, x, order) == x for x in words)          # f^order = id
        assert not all(iterate(f, x, order // 2) == x for x in words)  # exact order
        assert all(iterate(f, f(x), order - 1) == x for x in words)   # f^-1 = f^(order-1)
        print(f"{name}^-1 = XOR of ROTR^r, r in {S}  ({len(S)} terms, "
              f"order of {name} is {order})")
    assert all(iterate(Sigma0, x, 8) == rotr(x, 8) for x in words)     # Sigma_0^8 = ROTR^8
    print("all identities verified on", len(words), "words")
