# Data sources

Every external dataset used, with retrieval date and terms.

## Sequence data

| Source | Query | Retrieved | n |
|---|---|---|---|
| UniProt | InterPro entry [IPR001223](https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001223/), bacteria + archaea | 7 October 2025 | 39,582 sequences |
| UniProt | IPR001223, eukaryotes | 2 July 2026 | 46,524 proteins / 51,593 GH18 domains |
| NCBI nr (clustered) | Targeted BLAST to stabilise long branches | — | included in the 773-sequence tree; no SSN cluster assignment |
| CAZy | [GH18 family](https://www.cazy.org/GH18.html) | 7 July 2026 | 258 characterised bacterial/archaeal entries after filtering |

UniProt data are released under CC-BY-4.0. InterPro is freely available. CAZy content is the property of its authors — cite Drula *et al.* and check the CAZy terms of use before redistributing `data/cazy/`.

## Structures

| PDB | Origin |
|---|---|
| 9BUF, 9BUG, 9OJL | Determined in this work (MX2 beamline, Australian Synchrotron) |
| 1KFW, 4W5U | Previously determined; shown in Fig. 3b |
| AlphaFold3 models, InDel classes F and G | Predicted here; deposited at Zenodo |

## Software

| Tool | Version | Use |
|---|---|---|
| InterProScan | 6.0.0 | Domain annotation, default thresholds |
| MMseqs2 | — *(record the version)* | `easy-search` for the SSN; `easy-cluster` for dereplication (`--min-seq-id 0.30–0.90 -c 0.8 --cov-mode 0 -s 7.5`) |
| MAFFT | — *(record the version)* | Alignment, `--dash --localpair` |
| trimAl | — *(record the version)* | Alignment trimming |
| FastTree | — *(record the version)* | Initial tree visualisation during curation |
| IQ-TREE | 3.1.3 | ML inference, Q.PFAM+F+R10; UFBoot, SH-aLRT, sCF, AU root tests |
| SignalP | 6.0 | Signal peptide prediction (installed locally; **academic use only**) |
| AlphaFold3 | — | Structure prediction where no experimental structure existed |
| ProteinClusterTools | [johnchen93/ProteinClusterTools](https://github.com/johnchen93/ProteinClusterTools) | SSN construction, cluster assignment, HMM-based representative selection |
| Python | 3.13.5 | NumPy, SciPy, pandas, statsmodels, matplotlib, seaborn, Biopython |
| R | 4.6.0 | ape 5.8-1, caper 1.0.4, phylolm 2.6.5, phytools 2.5-2, phangorn 2.12.1, corHMM 2.8, dplyr 1.2.1 — full list in `R_sessionInfo.txt` |
| XDS, AIMLESS, Phaser (CCP4), phenix.refine, Coot 1.1.17, MolProbity | — | Crystallographic data processing and refinement |

Versions marked *(record the version)* should be filled in before release — reviewers of a *Nature Communications* Code Availability statement routinely ask for them.

## Redistribution note

`data/` contains tables **derived** from UniProt, InterPro, CAZy and SignalP output. Whichever licence is chosen for this repository (see `LICENSE`) applies to the code and to our own derived analysis, not to the upstream annotations, which retain their own terms.
