#!/usr/bin/env python3
"""Regroup a .bbl file into thematic clusters.

Run automatically after every bibtex invocation (see .latexmkrc),
which passes the .bbl path of whichever variant is being built
(main.bbl, main-algebra.bbl, ...); defaults to main.bbl.
Reads cluster membership from the bib/s*.bib files, then rewrites
the .bbl as one thebibliography block per cluster, each preceded by a
\\refcluster{title}{note} header (macro defined in backmatter.tex).
Idempotent: a marker comment prevents double processing.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BBL = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "main.bbl"
MARKER = "% clustered-by-clusterbbl"

# (bibfile, title, one-line thematic note) in presentation order
CLUSTERS = [
    ("s1-solvable",
     "Patterns, exact solutions, and the algebra of iteration",
     r"Where seemingly procedural discrete dynamics collapsed to closed "
     r"forms---the genre of \cref{sec:josephus,sec:rivest,sec:hope}."),
    ("s2-prehistory",
     "Hashing before cryptography, and the birth of one-wayness",
     r"Scatter storage in early computing, and the adversary's gradual "
     r"arrival (\cref{sec:falling})."),
    ("s3-family",
     "The SHA family: standards, design theory, breaks, hardware",
     r"The object itself: specifications, the Merkle--Damg{\aa}rd "
     r"theorem, the fossil record of attacks, and silicon."),
    ("s4-uses",
     "Uses: fingerprints, trees, proofs, money",
     r"What the world does with a function it cannot analyze "
     r"(\cref{sec:uses})."),
    ("s5-mappings",
     "Random mappings and analytic combinatorics",
     r"Flajolet's telescope (\cref{sec:randommaps}): the exact "
     r"statistical silhouette of a structureless map."),
    ("s6-prng",
     "Pseudorandom numbers, finite fields, and mixing",
     r"The finished classical theories (\cref{sec:prngs}): linear "
     r"generators over finite fields, and random walks on groups."),
    ("s7-foundations",
     "Foundations: definitions, idealized models, limits of proof",
     r"What ``pseudorandom'' and ``secure'' officially mean, and why "
     r"unconditional theorems are out of reach (\cref{sec:frameworks})."),
    ("s10-provable",
     "Hashes with theorems: the provable road",
     r"Hash functions whose collision resistance is a theorem---and what "
     r"that trade costs (\cref{sec:provable})."),
    ("s9-differential",
     "The adversary's calculus: differential and dedicated cryptanalysis",
     r"The method that felled MD5 and SHA-1, and the state of its siege "
     r"of \SHA{} (\cref{sec:differential})."),
    ("s8-computability",
     "Computability, complexity, and barriers",
     r"Where breaking \SHA{} sits in the complexity landscape "
     r"(\cref{sec:computability}): \textsc{np}-hardness, \textsc{sat} "
     r"solvers, quantum limits, and proof barriers."),
]


def cluster_of_key():
    mapping = {}
    for fname, _, _ in CLUSTERS:
        text = (ROOT / "bib" / f"{fname}.bib").read_text()
        for m in re.finditer(r"^@\w+\{([^,\s]+)\s*,", text, re.M):
            mapping[m.group(1)] = fname
    return mapping


def main():
    if not BBL.exists():
        return 0
    src = BBL.read_text()
    if MARKER in src:
        return 0

    head, sep, rest = src.partition("\\begin{thebibliography}")
    if not sep:
        return 0
    m = re.match(r"\{[^}]*\}", rest)
    widest = m.group(0) if m else "{99}"
    body = rest[m.end():] if m else rest
    body = body.replace("\\end{thebibliography}", "")

    # split into preamble (\providecommand block) and \bibitem chunks
    first_item = body.find("\\bibitem")
    preamble, items_txt = body[:first_item], body[first_item:]
    chunks = re.split(r"(?=\\bibitem)", items_txt)
    chunks = [c for c in chunks if c.strip()]

    keymap = cluster_of_key()
    grouped = {fname: [] for fname, _, _ in CLUSTERS}
    for c in chunks:
        km = re.search(r"\\bibitem\[[^\]]*\]\{([^}]+)\}", c, re.S) or \
             re.search(r"\\bibitem\{([^}]+)\}", c)
        key = km.group(1) if km else None
        grouped.setdefault(keymap.get(key, CLUSTERS[0][0]), []).append(c)

    # ONE thebibliography environment (multiple ones break hyperref's
    # backref stream); cluster headers are label-less items inside it.
    out = [MARKER + "\n" + head.rstrip() + "\n"]
    out.append(f"\\begin{{thebibliography}}{widest}\n")
    out.append(preamble.strip() + "\n\n")
    for fname, title, note in CLUSTERS:
        items = grouped.get(fname)
        if not items:
            continue
        out.append(f"\\refclusteritem{{{title}}}{{{note}}}\n\n")
        out.extend(item.rstrip() + "\n\n" for item in items)
    out.append("\\end{thebibliography}\n")
    BBL.write_text("".join(out))
    print(f"clusterbbl: regrouped {len(chunks)} entries into clusters")
    return 0


if __name__ == "__main__":
    sys.exit(main())
