# The modular evolution of chitinases is governed by coevolution of auxiliary and catalytic domains

Analysis code and derived data for Pulsford *et al.* (2026), *Nature Communications*.


Correspondence: sbp92@cornell.edu · colin.jackson@anu.edu.au

---

## What this is

We mapped the sequence space of prokaryotic GH18 chitinases (39,582 UniProt sequences), reconciled their InterProScan domain annotations, quantified auxiliary-domain promiscuity with a GH18-specific adaptation of the Domain Versatility Index (φ), and tested how domain architecture, secretion and taxonomy partition across a 773-sequence maximum-likelihood phylogeny.

This repository holds the **custom analysis code** and the **derived tables small enough to version**. Large primary matrices, full IQ-TREE run directories and the SSN node/edge tables live at Zenodo; raw sequences come from public repositories. See [Where everything lives](#where-everything-lives).

## Where everything lives

| What | Where |
|---|---|
| Analysis code, derived tables < 25 MB, replicate trees | **This repository** |
| Master annotation matrix, full IQ-TREE outputs, SSN node/edge tables, per-domain φ table, AlphaFold3 models, InDel HMM library | **Zenodo — [10.5281/zenodo.21828480](https://doi.org/10.5281/zenodo.21828480)** |
| Supplementary Data Files 1–3, Supplementary Tables 1–24, Supplementary Figures 1–11 | Published Supplementary Information |
| Prokaryotic + archaeal GH18 sequences (n = 39,582) | UniProt, InterPro entry [IPR001223](https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001223/), retrieved 7 Oct 2025 |
| Eukaryotic GH18 sequences (n = 46,524) | UniProt, IPR001223, retrieved 2 Jul 2026 |
| Long-branch stabilising sequences | NCBI nr clustered database |
| Characterised GH18 entries | [CAZy GH18](https://www.cazy.org/GH18.html), accessed 7 Jul 2026 |
| Crystal structures | PDB [9BUF](https://www.rcsb.org/structure/9BUF), [9BUG](https://www.rcsb.org/structure/9BUG), [9OJL](https://www.rcsb.org/structure/9OJL) |

`docs/zenodo_manifest.md` lists exactly which files were too large to version here and belong in the Zenodo deposit.

## Running the analyses

**Python** (notebooks 01–07), Python 3.13:

```bash
conda env create -f environment.yml && conda activate gh18
jupyter lab
```

or `pip install -r requirements.txt`.

**R** (`code/phylo_traits/`), R 4.6.0:

```r
install.packages(c("ape", "caper", "phylolm", "phytools", "phangorn",
                   "corHMM", "nlme", "dplyr", "tidyr"))
```

Exact versions used for the published results are in `R_sessionInfo.txt`. Run the R scripts from `data/phylogeny/` (they read `tree.nwk` and the tip-trait table from the working directory) — see `docs/known_gaps.md` for the path and column-name fixes still outstanding.

**External tools** (not run from this repo): InterProScan 6.0.0, MMseqs2, MAFFT, trimAl, IQ-TREE 3.1.3, SignalP 6.0, AlphaFold3, and [ProteinClusterTools](https://github.com/johnchen93/ProteinClusterTools) for SSN construction, cluster assignment and representative selection. Exact IQ-TREE invocations are in `code/phylogenetics/iqtree_commands.sh`.


## Licence

**Not yet set** — see `LICENSE`. This must be resolved before the repository is made public; *Nature Communications* requires custom code to be available under a recognised licence.

## Contact

Sacha B. Pulsford — sbp92@cornell.edu
