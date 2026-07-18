"""Centralizer census for word-narrowed SHA-256 components.

Companion script to "SHA-256 as a Mathematical Object" (the
nonarchimedean-dynamics passage of the geometry section).  The text
argues that SHA-256's designers are in the business of manufacturing
1-Lipschitz dynamics that commute with nothing; this script makes that
claim quantitative at word width w = 8, where a single word has 256
values and every question below is exactly decidable.

For each narrowed component f (rotation amounts of the standard's
Sigma/sigma maps reduced mod 8, round constants narrowed to their top
bytes) the census measures:

  1. the exact number of transformations (arbitrary self-maps of
     Z/2^8) commuting with f -- by a cycle-type closed form when f is
     a permutation, by an endomorphism-counting dynamic program over
     its functional graph when it is not;
  2. the exact number of permutations commuting with f (the
     centralizer of f in Sym(256), again from the cycle type);
  3. commuting maps inside enumerable structured classes: additive
     translations x + c, xor translations x ^ c, and the
     affine-circulant family x -> (xor of rotations of x) ^ c;
  4. for F_2-linear f, the exact dimension of the full linear
     commutant {M : MF = FM}, by Gaussian elimination over F_2.

Every statistic is compared against the same statistic for uniformly
random permutations / random maps of 256 points, so that "commutes
with nothing" becomes "its centralizer statistics sit where the
random baseline puts them" -- a measurement, not a slogan.

A final section runs the Anashin-flavored ablation from the adjacent
project box: the narrowed add-constant map and the Klimov--Shamir
T-function x -> x + (x^2 | 5) are verified to be single 2^k-cycles
for every k <= 16 (transitivity mod 2^k, the finite shadow of 2-adic
ergodicity), then a single rotation is composed in and the cycle
statistics collapse onto the Flajolet--Odlyzko random-permutation
profile.

All counts are exact integers; the counting routines self-test against
brute-force enumeration on small domains before anything is reported.
Pure stdlib; runs in about ten seconds.
"""
import math
import random
from itertools import product

W = 8
N = 1 << W
MASK = N - 1

# Narrowed round constants: top bytes of K_2 = b5c0fbcf and
# K_1 = 71374491 (both odd, so x -> x + K is a full 256-cycle).
KA, KB = 0xB5, 0x71


def rotr(x, r, w=W):
    r %= w
    m = (1 << w) - 1
    return ((x >> r) | (x << (w - r))) & m


# --- narrowed components ------------------------------------------------
# Rotation amounts are the standard's, reduced mod 8; shift amounts are
# reduced mod 8 as well (a shift by >= 8 would annihilate the word).
# The list is chosen to probe the paper's three predicted regimes:
# the pure F_2-linear diffusion layers (expected: over-symmetric), the
# pure 2-adic odometer (expected: under-symmetric -- it commutes with
# all of its own powers and nothing else), and their alternation
# (expected: symmetry statistics of a random permutation).  That is
# the ARX design thesis in miniature, made measurable.
def Sigma0(x):  # ROTR 2,13,22 -> 2,5,6.  F_2-linear, wild 2-adically.
    return rotr(x, 2) ^ rotr(x, 5) ^ rotr(x, 6)


def Sigma1(x):  # ROTR 6,11,25 -> 6,3,1.  Ditto.
    return rotr(x, 6) ^ rotr(x, 3) ^ rotr(x, 1)


def sigma0(x):  # ROTR 7,18, SHR 3 -> ROTR 7,2, SHR 3.  The shift makes
    return rotr(x, 7) ^ rotr(x, 2) ^ (x >> 3)  # it non-circulant.


def sigma1(x):  # ROTR 17,19, SHR 10 -> ROTR 1,3, SHR 2.  Rank-deficient
    return rotr(x, 1) ^ rotr(x, 3) ^ (x >> 2)  # at w=8: not a bijection.


def addK(x):  # the bare odometer: 1-Lipschitz, single 256-cycle (K odd)
    return (x + KA) & MASK


def arx1(x):  # ONE ARX alternation: add (2-adic), then rotate-xor (F_2)
    return Sigma0((x + KA) & MASK)


def arx2(x):  # two alternations: does the generic regime persist?
    return Sigma1((arx1(x) + KB) & MASK)


def chdiag(x):  # the round's chooser Ch(e,f,g) = (e&f)^(~e&g), made a
    # self-map by feeding it rotated copies of one word: the census's
    # nonlinear NON-bijection (image size < 256), probing the regime
    # where commutant counts become functional-graph statistics.
    ch = (x & rotr(x, 2)) ^ (~x & MASK & rotr(x, 5))
    return (ch + KA) & MASK


COMPONENTS = [
    ("Sigma0", Sigma0), ("Sigma1", Sigma1),
    ("sigma0", sigma0), ("sigma1", sigma1),
    ("addK", addK), ("arx1", arx1), ("arx2", arx2),
    ("chdiag", chdiag),
]


# --- cycle-type closed forms (f a permutation) --------------------------
# For a permutation, every centralizer question below is a function of
# the cycle type alone -- the "symmetry mass" of f is readable off its
# orbit structure, which is why the paper can compare components to
# random permutations by comparing cycle types.
def cycle_type(tab):
    n = len(tab)
    seen = [False] * n
    ct = {}
    for x in range(n):
        if not seen[x]:
            length, y = 0, x
            while not seen[y]:
                seen[y] = True
                y = tab[y]
                length += 1
            ct[length] = ct.get(length, 0) + 1
    return ct


def perm_centralizer_size(ct):
    """|C_Sym(f)| = prod k^{c_k} c_k! over the cycle type."""
    out = 1
    for k, c in ct.items():
        out *= k ** c * math.factorial(c)
    return out


def perm_commuting_maps(ct):
    """Number of ALL self-maps g with g f = f g, f a permutation.

    g is free on one point per cycle; on a k-cycle the image point y
    must satisfy f^k(y) = y, and #Fix(f^k) = sum_{d | k} d c_d.
    """
    out = 1
    for k, c in ct.items():
        fix = sum(d * cd for d, cd in ct.items() if k % d == 0)
        out *= fix ** c
    return out


# --- exact commutant count for arbitrary maps ---------------------------
def count_commuting_maps(tab):
    """Exact #{g : g o f = f o g}, f = tab an arbitrary self-map.

    g commutes with f exactly when g is an endomorphism of f's
    functional graph (the digraph x -> f(x)): g(x) -> g(f(x)) must be
    an edge, i.e. g(f(x)) = f(g(x)).  Endomorphisms are counted by
    dynamic programming, without enumeration (the count is often
    astronomical -- a constant map f has n^(n-1) of them):

      * each weak component contributes independently (multiply);
      * a component is a cycle of length k with trees hanging off it;
        the image of the cycle is forced by the image w of one cycle
        node, and w may be any point with f^k(w) = w;
      * for a tree node u whose parent's image is v, the image y of u
        ranges over the f-preimages of v, and u's subtree then
        contributes independently given y:
            cnt[u][v] = sum_{y in f^-1(v)} prod_{c child of u} cnt[c][y],
        computed leaves-first, so the total per component is
            sum_{w : f^k(w) = w} prod_j prod_{c child of z_j} cnt[c][f^j(w)]
        over the cycle nodes z_j.

    Interpretive caveat used in the text: for a NON-bijection this
    count is dominated by collapsing endomorphisms (maps folding trees
    onto each other), so it measures the shape of f's functional graph
    -- a Flajolet--Odlyzko-style statistic -- more than any hidden
    algebra.  The sharp structurelessness claims live in the
    permutation rows and the structured-class censuses.
    """
    n = len(tab)
    preim = [[] for _ in range(n)]
    for x in range(n):
        preim[tab[x]].append(x)

    # cycle detection (color: 0 unvisited, 1 on current path, 2 done)
    on_cycle = [False] * n
    color = [0] * n
    for s in range(n):
        if color[s]:
            continue
        path, x = [], s
        while color[x] == 0:
            color[x] = 1
            path.append(x)
            x = tab[x]
        if color[x] == 1:  # new cycle discovered; mark it
            y = x
            while True:
                on_cycle[y] = True
                y = tab[y]
                if y == x:
                    break
        for y in path:
            color[y] = 2

    # weak components
    comp = [-1] * n
    comps = []
    for s in range(n):
        if comp[s] == -1:
            ci = len(comps)
            comp[s] = ci
            stack, members = [s], [s]
            while stack:
                x = stack.pop()
                for y in preim[x] + [tab[x]]:
                    if comp[y] == -1:
                        comp[y] = ci
                        stack.append(y)
                        members.append(y)
            comps.append(members)

    cnt = [None] * n
    total = 1
    for members in comps:
        z0 = next(x for x in members if on_cycle[x])
        cycle_order = [z0]
        while tab[cycle_order[-1]] != z0:
            cycle_order.append(tab[cycle_order[-1]])
        k = len(cycle_order)

        # tree nodes in BFS order from the cycle outward
        bfs = [p for z in cycle_order for p in preim[z] if not on_cycle[p]]
        i = 0
        while i < len(bfs):
            bfs.extend(preim[bfs[i]])  # preimages of a tree node are tree nodes
            i += 1

        for u in reversed(bfs):  # leaves first
            kids = preim[u]
            row = [0] * n
            for v in range(n):
                tot = 0
                for y in preim[v]:
                    pr = 1
                    for c in kids:
                        pr *= cnt[c][y]
                        if pr == 0:
                            break
                    tot += pr
                row[v] = tot
            cnt[u] = row

        comp_total = 0
        for w in range(n):
            x = w
            for _ in range(k):
                x = tab[x]
            if x != w:  # w must be periodic with period dividing k
                continue
            ways, phi = 1, w
            for z in cycle_order:
                for c in preim[z]:
                    if not on_cycle[c]:
                        ways *= cnt[c][phi]
                if ways == 0:
                    break
                phi = tab[phi]
            comp_total += ways
        total *= comp_total
    return total


def brute_commuting_maps(tab):
    n = len(tab)
    return sum(
        all(gg[tab[x]] == tab[gg[x]] for x in range(n))
        for gg in product(range(n), repeat=n)
    )


# --- structured-class censuses ------------------------------------------
# The raw centralizer counts above answer "how much symmetry"; these
# censuses answer "any symmetry an analyst could USE?".  The classes are
# the maps simple in the state space's two rival geometries (the paper's
# sec on the two group laws): translations x+c (2-adic side), xor
# translations x^c and xors-of-rotations (F_2 side), and their affine
# combination.  A nontrivial commuting member g would be a genuine
# self-reduction of f -- Lubin's theorem is the warning, one regularity
# class up, that such a symmetry drags hidden algebraic structure with
# it.  Finding none, against a calibrated baseline, is the measurement
# of structurelessness the text asks for.
def commuting_in_class(tab, gens):
    """Count g in a finite class of candidate maps with g f = f g."""
    n = len(tab)
    return sum(
        all(tab[gt[x]] == gt[tab[x]] for x in range(n)) for gt in gens
    )


def add_translations():
    return [[(x + c) & MASK for x in range(N)] for c in range(N)]


def xor_translations():
    return [[x ^ c for x in range(N)] for c in range(N)]


def affine_circulant_commutant(tab):
    """Count pairs (S, c) with x -> (xor_{r in S} ROTR^r x) ^ c
    commuting with tab.  2^8 subsets x 2^8 constants = 65536 maps."""
    rots = [[rotr(x, r) for x in range(N)] for r in range(W)]
    circ = {0: [0] * N}
    count = 0
    for S in range(1 << W):
        if S:
            low = (S & -S).bit_length() - 1
            prev = circ[S ^ (1 << low)]
            circ[S] = [prev[x] ^ rots[low][x] for x in range(N)]
        LS = circ[S]
        for c in range(N):
            if all(tab[LS[x] ^ c] == LS[tab[x]] ^ c for x in range(N)):
                count += 1
    return count


# --- F_2 linear algebra -------------------------------------------------
# For the F_2-linear components the FULL linear commutant {M : MF = FM}
# is exactly computable: MF = FM is a linear (Sylvester-type) system in
# M's w^2 entries over F_2.  This is how the census certifies, e.g.,
# that narrowed Sigma_0's linear symmetries are exactly the 2^8
# xors-of-rotations (the circulant algebra F_2[ROTR]) and nothing more.
def linear_cols(tab):
    """Column bitmasks of tab's matrix if tab is F_2-linear, else None."""
    if tab[0] != 0:
        return None
    cols = [tab[1 << i] for i in range(W)]
    for x in range(N):
        y = 0
        for i in range(W):
            if (x >> i) & 1:
                y ^= cols[i]
        if y != tab[x]:
            return None
    return cols


def f2_rank(rows, ncols):
    rows = list(rows)
    rank = 0
    for c in range(ncols):
        piv = next(
            (i for i in range(rank, len(rows)) if (rows[i] >> c) & 1), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(len(rows)):
            if i != rank and (rows[i] >> c) & 1:
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def linear_commutant_dim(cols):
    """dim over F_2 of {M : MF = FM}; the commutant has 2^dim elements."""
    F = [[(cols[j] >> i) & 1 for j in range(W)] for i in range(W)]
    rows = []
    for i in range(W):
        for j in range(W):
            row = 0
            for k in range(W):
                if F[k][j]:
                    row ^= 1 << (i * W + k)   # M[i,k] F[k,j]
                if F[i][k]:
                    row ^= 1 << (k * W + j)   # F[i,k] M[k,j]
            rows.append(row)
    return W * W - f2_rank(rows, W * W)


def matrix_rank(cols):
    return f2_rank(cols, W)  # columns as bitmasks = rows of transpose


# --- Anashin / ablation -------------------------------------------------
# Anashin's theorem: a 1-Lipschitz map of Z_2 is ergodic iff it is
# transitive on Z/2^k for every k -- i.e. iff it is a single 2^k-cycle
# at every truncation.  That infinite hierarchy has a finite shadow we
# can check outright; then composing in ONE rotation (the operation
# that is F_2-linear but 2-adically discontinuous) destroys the
# hierarchy, and the cycle statistics collapse onto the random-
# permutation profile of Flajolet--Odlyzko (expected #cycles ~ ln N +
# gamma, expected longest-cycle fraction the Golomb--Dickman constant
# 0.6243...).  This is the project box's "break the hypotheses and
# measure what dies", run to completion.
def single_cycle_lengths(step, kmax):
    """Largest k <= kmax such that the map is a single 2^k-cycle for
    every k' <= k (transitivity mod 2^k -- ergodicity's finite shadow)."""
    good = 0
    for k in range(1, kmax + 1):
        n, m = 1 << k, (1 << k) - 1
        x, length = step(0) & m, 1
        while x != 0:
            x = step(x) & m
            length += 1
            if length > n:
                break
        if length == n:
            good = k
        else:
            break
    return good


def cycle_stats(tab):
    ct = cycle_type(tab)
    ncyc = sum(ct.values())
    longest = max(ct)
    return ncyc, longest / len(tab)


def truncation_spread(step, k, low):
    """Average #distinct low-`low`-bit outputs over each residue class
    mod 2^low (1 for a T-function; ~2^low (1 - 1/e) for random).

    Measures what remains of 1-Lipschitz causality (output bit i
    depends only on input bits 0..i) after ablation.  The measured
    value for the rotated map is exactly 128 = 2^7: seven of the low
    byte's bits scatter fully, but one bit -- bit 7 of the T-function
    output, carried into the low byte by ROTR^7 -- stays frozen when
    the input's low byte is fixed.  The substrate leaves a scar."""
    n, lm = 1 << k, (1 << low) - 1
    spread = 0
    for a in range(1 << low):
        vals = {step(a + (t << low)) & lm for t in range(1 << (k - low))}
        spread += len(vals)
    return spread / (1 << low)


# --- reporting ----------------------------------------------------------
def fmt(x):
    if x < 10 ** 6:
        return str(x)
    return f"{x:.3e}".replace("e+0", "e").replace("e+", "e")


def percentile_of(value, samples):
    return 100.0 * sum(s < value for s in samples) / len(samples)


def selftest():
    random.seed(20260717)
    for n in (4, 5):
        for _ in range(6):
            tab = [random.randrange(n) for _ in range(n)]
            assert count_commuting_maps(tab) == brute_commuting_maps(tab)
        perm = list(range(n))
        random.shuffle(perm)
        ct = cycle_type(perm)
        assert count_commuting_maps(perm) == perm_commuting_maps(ct)
        assert count_commuting_maps(perm) == brute_commuting_maps(perm)
    print("self-tests passed (search == brute force == closed form)\n")


def main():
    selftest()
    random.seed(20260717)

    tabs = {name: [f(x) for x in range(N)] for name, f in COMPONENTS}

    # random baselines --------------------------------------------------
    # The counts only mean something against the null hypothesis "f is
    # a uniformly random (permutation | map) of 256 points".  Random
    # permutations are cheap to sample in bulk (their statistics are
    # cycle-type functions); the affine-circulant census over random
    # permutations calibrates how many "accidental" structured
    # symmetries a generic permutation carries -- in practice the
    # constant map at each fixed point, Poisson(1) many on average.
    print("=== random baselines on 256 points ===")
    perm_cent, perm_comm = [], []
    base = list(range(N))
    for _ in range(10_000):
        p = base[:]
        random.shuffle(p)
        ct = cycle_type(p)
        perm_cent.append(perm_centralizer_size(ct))
        perm_comm.append(perm_commuting_maps(ct))
    for label, xs in (("|C_Sym|", perm_cent), ("#commuting maps", perm_comm)):
        logs = sorted(math.log10(x) for x in xs)
        print(f"  random permutation {label}: log10 median "
              f"{logs[len(logs)//2]:.2f}, 5-95% "
              f"[{logs[500]:.2f}, {logs[9500]:.2f}]")
    map_comm = []
    for _ in range(100):
        m = [random.randrange(N) for _ in range(N)]
        map_comm.append(count_commuting_maps(m))
    map_comm.sort()
    mlogs = [math.log10(x) for x in map_comm]
    print(f"  random map #commuting maps: log10 median {mlogs[50]:.2f}, "
          f"5-95% [{mlogs[5]:.2f}, {mlogs[95]:.2f}] over 100 samples")
    ntriv = 0
    for _ in range(200):
        p = base[:]
        random.shuffle(p)
        ntriv += affine_circulant_commutant(p) - 1
    print(f"  random permutation nontrivial affine-circulant "
          f"commutants: {ntriv} across 200 samples\n")

    # component census --------------------------------------------------
    adds, xors = add_translations(), xor_translations()
    for name, _ in COMPONENTS:
        tab = tabs[name]
        bij = len(set(tab)) == N
        cols = linear_cols(tab)
        print(f"=== {name} ({'bijective' if bij else 'non-bijective'}"
              f"{', F2-linear' if cols else ''}) ===")
        if cols:
            print(f"  rank {matrix_rank(cols)}/{W}, linear commutant "
                  f"dim {linear_commutant_dim(cols)} "
                  f"(2^{linear_commutant_dim(cols)} commuting linear maps)")
        if bij:
            ct = cycle_type(tab)
            cent, comm = perm_centralizer_size(ct), perm_commuting_maps(ct)
            check = count_commuting_maps(tab)
            assert check == comm, (name, check, comm)
            print(f"  cycle type {{len: count}} = "
                  f"{dict(sorted(ct.items()))}")
            print(f"  |C_Sym| = {fmt(cent)} "
                  f"(random-perm percentile {percentile_of(cent, perm_cent):.0f})")
            print(f"  #commuting maps = {fmt(comm)} "
                  f"(random-perm percentile {percentile_of(comm, perm_comm):.0f})")
        else:
            comm = count_commuting_maps(tab)
            print(f"  #commuting maps = {fmt(comm)} "
                  f"(log10 = {math.log10(comm):.2f}; random-map percentile "
                  f"{percentile_of(comm, map_comm):.0f})")
        print(f"  commuting additive translations x+c: "
              f"{commuting_in_class(tab, adds)}/256")
        print(f"  commuting xor translations x^c: "
              f"{commuting_in_class(tab, xors)}/256")
        print(f"  commuting affine-circulant maps: "
              f"{affine_circulant_commutant(tab)}/65536\n")

    # Anashin criteria and the rotation ablation ------------------------
    print("=== Anashin: single-orbit checks and the rotation ablation ===")
    ks = lambda x: x + (x * x | 5)
    print(f"  x + 0x{KA:2x}: single 2^k-cycle for all k <= "
          f"{single_cycle_lengths(lambda x: x + KA, 16)} (checked to 16)")
    print(f"  x + (x^2|5): single 2^k-cycle for all k <= "
          f"{single_cycle_lengths(ks, 16)} (checked to 16)")
    k = 16
    n, m = 1 << k, (1 << k) - 1
    stats = []
    for r in range(1, k):
        tab = [rotr(ks(x) & m, r, k) for x in range(n)]
        stats.append(cycle_stats(tab))
    ncyc = [s[0] for s in stats]
    lfrac = [s[1] for s in stats]
    h_n = sum(1.0 / i for i in range(1, n + 1))
    print(f"  ablation ROTR^r(x + (x^2|5)) mod 2^16, r = 1..15:")
    print(f"    cycles: mean {sum(ncyc)/len(ncyc):.2f} "
          f"range [{min(ncyc)}, {max(ncyc)}]  "
          f"(random-perm expectation H_65536 = {h_n:.2f})")
    print(f"    longest-cycle fraction: mean {sum(lfrac)/len(lfrac):.3f} "
          f"(Golomb-Dickman lambda = 0.624)")
    rp = []
    for _ in range(100):
        p = list(range(n))
        random.shuffle(p)
        rp.append(cycle_stats(p)[0])
    print(f"    (100 random perms of 2^16: cycles mean "
          f"{sum(rp)/len(rp):.2f}, range [{min(rp)}, {max(rp)}])")
    t_spread = truncation_spread(lambda x: ks(x) & m, k, 8)
    r_spread = truncation_spread(lambda x: rotr(ks(x) & m, 7, k), k, 8)
    print(f"  truncation: distinct low bytes per residue class mod 2^8:")
    print(f"    T-function x+(x^2|5): {t_spread:.1f}   "
          f"with ROTR^7: {r_spread:.1f}   "
          f"(random expectation 256(1-1/e) = {256*(1-math.exp(-1)):.1f})")


if __name__ == "__main__":
    main()
