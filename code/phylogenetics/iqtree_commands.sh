#!/usr/bin/env bash
# =============================================================================
# Phylogenetic inference — GH18 catalytic domains
# IQ-TREE 3.1.3. Five independent replicates were run for each model regime;
# all five reversible replicates converged on Q.PFAM+F+R10.
#
# Input alignment: data/phylogeny/gh18_msa_trimmed.fasta
#   (773 sequences: 623 HMM-based SSN representatives across bit scores
#    400-800 in steps of 25, plus UniProt "Reviewed" entries, plus targeted
#    NCBI nr sequences added to stabilise long branches.
#    MAFFT --dash --localpair, manually curated, trimAl.)
#
# The alignment is referred to below by its original run-time name,
# addedpathog_trimal.fasta. Outputs kept in this repository are the trees and
# support tables only; full run directories (~1.7 GB of .state, .ufboot,
# .ckp.gz, .trees and .mldist files) are deposited at Zenodo
# (10.5281/zenodo.21828480).
#
# Runs were performed on NCI Australia. Paths below are as executed.
# =============================================================================

set -euo pipefail
IQTREE=/programs/iqtree-3.1.3-Linux/bin/iqtree3
IQTREE_INTEL=/programs/iqtree-3.1.3-Linux/bin/iqtree3_intel
ALN=addedpathog_trimal.fasta

# -----------------------------------------------------------------------------
# 1. Main ML search, reversible models — ModelFinder Plus.
#    Run independently five times (replicates 1-5 -> data/phylogeny/replicates_mfp/).
#    All replicates selected Q.PFAM+F+R10.
# -----------------------------------------------------------------------------
"$IQTREE" \
    -s "$ALN" \
    -m MFP \
    -nt AUTO -ntmax 16 \
    -bb 1000 -bnni \
    -pers 0.3 -nstop 200 -ninit 10 -nbest 5 \
    -wbtl -asr \
    -safe -T AUTO --threads-max 8

# -----------------------------------------------------------------------------
# 2. Main ML search, non-reversible model (rooting).
#    Run independently five times -> data/phylogeny/replicates_nq/.
# -----------------------------------------------------------------------------
"$IQTREE" \
    -s "$ALN" \
    -mset NQ.pfam+F+R10 \
    -nt AUTO -ntmax 16 \
    -bb 1000 -bnni \
    -pers 0.3 -nstop 200 -ninit 10 -nbest 5 \
    -wbtl -asr \
    -safe -T AUTO --threads-max 8

# -----------------------------------------------------------------------------
# 3. SH-aLRT branch support (1000 replicates) on the fixed ML topology.
#    -> replicates_mfp/repN/sh_alrt.treefile
# -----------------------------------------------------------------------------
"$IQTREE_INTEL" \
    -s "$ALN" \
    -te addedpathog_trimalout.fasta.treefile \
    -m Q.PFAM+F+R10 \
    --alrt 1000 \
    -safe -T AUTO --threads-max 8 \
    --prefix run_alrt

# -----------------------------------------------------------------------------
# 4. Site concordance factors (100 quartets).
#    -> replicates_mfp/repN/site_concordance.cf.tree
# -----------------------------------------------------------------------------
"$IQTREE_INTEL" \
    -s "$ALN" \
    -te "$ALN.treefile" \
    -m Q.PFAM+F+R10 \
    --scfl 100 \
    -safe -T AUTO --threads-max 8 \
    --prefix scfl_

# -----------------------------------------------------------------------------
# 5. AU topology test over the combined replicate topologies.
#    Input list: data/phylogeny/root_tests/combined_topologies.treels
# -----------------------------------------------------------------------------
"$IQTREE" -s "$ALN" -z combinedtrees.treels \
    -m Q.PFAM+F+R10 -zb 10000 -zw -au

# -----------------------------------------------------------------------------
# 6. Root placement. Non-reversible NQ.PFAM+F+R10 root test.
#    -> data/phylogeny/root_tests/repN_roottest.csv (Supplementary Table 17)
# -----------------------------------------------------------------------------
"$IQTREE_INTEL" \
    -s "$ALN" -te addedpathog_trimalout.fasta.treefile \
    -m NQ.PFAM+F+R10 --root-test -au -zb 10000 --prefix roottest

"$IQTREE_INTEL" \
    -s "$ALN" -te addedpathog_trimalout.fasta.treefile \
    -m NQ.PFAM+F+R10 --root-test -au -zb 10000 \
    -nt 8 --prefix roottest8

# -----------------------------------------------------------------------------
# 7. AU test over seven candidate root positions (Supplementary Table 17).
#    Trees differ ONLY in root placement.
#    Input list: data/phylogeny/root_tests/candidate_roots.treels
#    The displayed root (InDel classes A+B vs C-G) was not rejected
#    (p-AU = 0.71) and is a presentational choice, not an inference.
# -----------------------------------------------------------------------------
"$IQTREE_INTEL" \
    -s "$ALN" -z candidate_roots.treels \
    -m NQ.PFAM+F+R10 -n 0 -zb 10000 -zw -au --prefix rootAU -nt AUTO

# -----------------------------------------------------------------------------
# 8. Minimal ancestral deviation rooting was applied separately
#    -> replicates_mfp/repN/mad_rooted.treefile
#    Internal node labels are stripped before rooting with
#    code/phylogenetics/strip_internal_labels.py
# -----------------------------------------------------------------------------
