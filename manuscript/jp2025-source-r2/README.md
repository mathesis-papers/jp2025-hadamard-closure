# Hadamard-closed polygons and cones

The note answers both parts of Johnson and Paparella's Question 9.1, arXiv:2409.07682v4, Section 9, p. 27, negatively. The rational family has `m+4` facets in three coordinates; every integer `m >= 1` gives a counterexample, beginning with an explicit pentagon. All its spectra admit a common symmetric stochastic realization. Products with cubes give full-dimensional counterexamples in every coordinate size `n >= 3`, with a common stochastic eigenbasis.

## Contents and reproduction

- `note.tex`: complete mathematical statement and proof, bibliography included.
- `note.pdf`: typeset note, provided alongside the source archive.
- `verify.py`: standard-library Python 3 exact rational certificates.
- `verification.json`: expected output of `python3 verify.py`.
- `metadata.json`: scientific description, authorship, license and verification link.

Build with Tectonic 0.17.0:

```sh
tectonic note.tex
python3 verify.py > verification-new.json
cmp verification.json verification-new.json
```

Alternatively, run `pdflatex note.tex` twice. The bibliography is embedded, so BibTeX and external bibliography files are unnecessary.

## Mathematical coverage

The note gives written proofs of the family for every positive integer `m`, its exact facet count, Hadamard closure, common stochastic realization, and both nonrepresentation conclusions for every `m >= 1`. The obstruction uses the diagonal entries for facets omitting the identity, together with the trace identity. The product extension has `m+2n-2` facets in coordinate size `n`.

The declaration `JohnsonPaparellaHadamardPolyhedraRefutation.result` in the [pinned Lean source](https://github.com/the-omega-institute/trureturing/blob/e277c80d045570daf43bd03b4e313f7ce457b4e2/D5/S0/FiniteGeometry/Polytopes/JohnsonPaparellaHadamardPolyhedraRefutation.lean) proves the negation of the joint universal converse using the `m=10` cone. Its fixed-cone scope is distinct from the parameter-family and polytope proofs and all three corollaries in the note. `verify.py` checks the finite examples `m=1,...,16,32,64`, including all supporting-face incidences, ordered generator products, explicit diamond convex coefficients, and stochastic matrix/eigenvector identities. The rectangular generator matrices are row Hadamard conic of rank three, yet their row cones have no same-coordinate spectracone representation.

The final corollary concerns the source's Conjecture 9.2, p. 27: in the complex affine space with first coordinate 1, `(1,0,...,0)` is a stochastic spectrum on the boundary but is not radially extremal. The matrix with all entries `1/n` realizes it; a perturbation with nonreal trace certifies boundary membership. The same conclusion holds for the boundary in complex Euclidean space. This does not assert a result for a different topology relative to conjugation-symmetric spectral lists or a real parameterization.

The finite checks also cover the normalized pentagonal trace discrepancy `19/9 != 1`, the cube extensions for `m=1,2` and `n=3,...,8`, and the stochastic projector spectra for `n=2,...,10`. The general identity-facet and ambient-boundary arguments are written proofs.

## Authors and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

Copyright 2026 Haobo Ma and Wenlin Zhang. Licensed under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/) (`cc-by-4.0`).
