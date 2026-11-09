# Split Primes and the Exponent of Imaginary Quadratic Class Groups

Deep Bhattacharjee. DOI (all versions): [10.5281/zenodo.24093927](https://doi.org/10.5281/zenodo.24093927). Paper in `paper/` (LaTeX source, PDF, figures); programs and data in `anc/`.

Under a zero-free half-plane for Dirichlet L-functions (Re s > 11/12, or Re s > 7/8), every imaginary quadratic field has a split prime below (2 log|d|)^12 (resp. ^8), the class group exponent grows like log|d|/log log|d|, and Euler's 65 idoneal numbers are all of them. The half-planes are claimed in unrefereed OpenAI preprints; every result states the hypothesis it uses.

Check: `python3 -I anc/certify.py` (needs python-flint) ends with `all assertions passed`.
