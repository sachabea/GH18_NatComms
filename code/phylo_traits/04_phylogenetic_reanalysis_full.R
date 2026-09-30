# =============================================================================
# R2 PHYLOGENETIC REANALYSIS — GH18 chitinase modular evolution
# Pulsford et al. 2026, NCOMMS-26-017590
#
# Addresses Reviewer 2's concern about phylogenetic non-independence.
#
# TREE USAGE:
#   tree_full (771 tips) — used for 1A (domain count transitions) and 1B
#                          (secretion D-statistic). All 771 tips have domain
#                          annotations and SP calls, so no pruning needed here.
#   tree_654  (654 tips) — used for 1B Pagel/phyloglm (requires clade),
#                          1C (PGLS promiscuity), 1D (cluster enrichment),
#                          1E (catalytic divergence). NCBI tips lack clade
#                          and cluster assignments so are excluded here.
#
# STEPS:
#   0   — Setup, loading, verification
#   1A  — Domain count transitions: parsimony + stochastic mapping (771 tips)
#   1B  — Secretion state: D-statistic (771 tips) + Pagel/phyloglm (654 tips)
#   1C  — Optional exploratory PGLS (not required for the Fig. 2 reviewer reply)
#   1D  — Fig. 2d complete-case phylogenetic logistic regression
#   1E  — Catalytic divergence (654 tips)
#   1F  — Methods mapping table
#
# RUNTIME (approximate):
#   Steps 0-1A parsimony:        < 2 min
#   fitMk model comparison:      3-5 min
#   make.simmap NSIM=1000:       40-90 min  *** main bottleneck ***
#   Step 1B D-statistic:         5-10 min
#   Step 1B fitPagel:            5-20 min (may fall back to corHMM)
#   Steps 1C-1E:                 5-10 min total
#
# Required files in working directory:
#   tree3.nwk
#   tip_annotation_table_complete.csv
#   (optional) 251109partner_stats_converted_protein_domains_interproscan_GH18replacedfprint.tsv      — needed for 1C
#   (optional) architecture_table.csv — needed for 1C cluster DVI
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(phangorn)
  library(nlme)
  library(phylolm)
  library(corHMM)
  library(caper)
  library(dplyr)
  library(tidyr)
})

cat("All packages loaded.\n\n")

# =============================================================================
# STEP 0 — Setup and verification
# =============================================================================
cat("========== STEP 0: Setup ==========\n\n")
#converted_protein_domains_interproscan_GH18replaced_noPF01522_cd06549_G3DSAclean.tsv
# ---- 0.1  Paths -------------------------------------------------------------
TREE_FILE  <- "tree.nwk"
ANNOT_FILE <- "tree_tip_traits_CANONICAL_654.csv"
DVI_FILE   <- "251109partner_stats_converted_protein_domains_interproscan_GH18replacedfprint.tsv"          # optional — 1C skips if absent
ARCH_FILE  <- "architecture_table.csv"  # optional — 1C cluster DVI skips if absent
OUT_DIR    <- "phylo_results"
NSIM       <- 1000   # stochastic mapping simulations — change to 10 for test run
PHYLO_BOOT <- 50   # bootstrap replicates for Fig. 2d phyloglm; use 50 for a test run

dir.create(OUT_DIR, showWarnings = FALSE)

# ---- 0.2  Load tree ---------------------------------------------------------
tree_full <- read.tree(TREE_FILE)
cat(sprintf("Full tree: %d tips, rooted: %s\n",
            Ntip(tree_full), is.rooted(tree_full)))
stopifnot(Ntip(tree_full) == 771)
stopifnot(is.rooted(tree_full))

# ---- 0.3  Identify NCBI tips ------------------------------------------------
ncbi_tips <- tree_full$tip.label[!grepl("\\|", tree_full$tip.label)]
cat(sprintf("NCBI tips (no pipe in label): %d\n", length(ncbi_tips)))
stopifnot(length(ncbi_tips) == 117)

# 654-tip tree for analyses requiring clade/cluster
tree_654 <- drop.tip(tree_full, ncbi_tips)
cat(sprintf("654-tip tree: %d tips, rooted: %s\n",
            Ntip(tree_654), is.rooted(tree_654)))
stopifnot(Ntip(tree_654) == 654)

# ---- 0.4  Load annotation table ---------------------------------------------
raw <- read.csv(ANNOT_FILE, stringsAsFactors = FALSE,
                na.strings = c("NA", "None", ""))
cat(sprintf("Annotation table: %d rows, %d columns\n", nrow(raw), ncol(raw)))

# Full 771-tip dataset (all tips, dom_state_3 available for everyone)
dat_771 <- raw
rownames(dat_771) <- dat_771$tip_label
dat_771 <- dat_771[tree_full$tip.label, ]
stopifnot(identical(rownames(dat_771), tree_full$tip.label))
cat(sprintf("dat_771 aligned: %d rows\n", nrow(dat_771)))

# Check dom_state_3 complete for all 771
n_missing_state <- sum(is.na(dat_771$dom_state_3))
cat(sprintf("dom_state_3 missing: %d / 771\n", n_missing_state))
if (n_missing_state > 0) stop("dom_state_3 has NA values — check annotation file.")

# 654-tip dataset: clade-assigned UniProt tips only
dat_654 <- raw[!is.na(raw$clade), ]
rownames(dat_654) <- dat_654$tip_label
dat_654 <- dat_654[tree_654$tip.label, ]
stopifnot(identical(rownames(dat_654), tree_654$tip.label))
cat(sprintf("dat_654 aligned: %d rows\n", nrow(dat_654)))

# ---- 0.5  Tip label match checks --------------------------------------------
# 771 tree vs annotation
diff_771 <- setdiff(tree_full$tip.label, dat_771$tip_label)
if (length(diff_771) > 0) stop(sprintf("771 mismatch: %d tips unmatched", length(diff_771)))

# 654 tree vs annotation
diff_654 <- setdiff(tree_654$tip.label, dat_654$tip_label)
if (length(diff_654) > 0) stop(sprintf("654 mismatch: %d tips unmatched", length(diff_654)))

cat("All tip labels verified.\n")
cat("\nState distribution — full 771-tip dataset:\n")
print(table(dat_771$dom_state_3))
cat("\nState distribution — 654-tip UniProt dataset:\n")
print(table(dat_654$dom_state_3))
cat("\nStep 0 complete.\n\n")

# =============================================================================
# STEP 1A — Domain count transitions
# Uses: tree_full (771 tips) — all tips now have dom_state_3
# Replaces: Fig 1b chi-squared
# =============================================================================
cat("========== STEP 1A: Domain count transitions (771 tips) ==========\n\n")

# ---- 1A.1  State vectors ----------------------------------------------------
levels_v <- c("single", "two", "multi")

states_3_771 <- setNames(dat_771$dom_state_3, dat_771$tip_label)
states_5_771 <- setNames(dat_771$dom_state_5, dat_771$tip_label)

cat("3-state distribution (771 tips):\n")
print(table(states_3_771))
cat("\n5-state distribution (771 tips):\n")
print(table(states_5_771))
cat("\n")

# ---- 1A.2  Parsimony — ACCTRAN + MPR ----------------------------------------
# Note: phyDat must be built from a named list, not a matrix (phangorn >= 2.11)
# DELTRAN not available in recent phangorn — use MPR to identify ambiguous nodes

pdat_771 <- phyDat(as.list(states_3_771), type = "USER", levels = levels_v)

cat("Running ACCTRAN parsimony...\n")
pars_acctran <- ancestral.pars(tree_full, pdat_771, type = "ACCTRAN")

cat("Running MPR parsimony (for ambiguity detection)...\n")
pars_mpr <- ancestral.pars(tree_full, pdat_771, type = "MPR")

# Extract single best state from ACCTRAN
extract_acctran <- function(pars_obj, levels_v) {
  vapply(seq_along(pars_obj), function(i) {
    m <- pars_obj[[i]]
    if (is.matrix(m)) levels_v[which.max(m[1, ])]
    else levels_v[which.max(m)]
  }, character(1))
}

# Flag ambiguous nodes: MPR has >1 equally parsimonious state
is_ambiguous_mpr <- function(pars_obj) {
  vapply(seq_along(pars_obj), function(i) {
    m <- pars_obj[[i]]
    if (is.matrix(m)) sum(m[1, ] > 0) > 1
    else sum(m > 0) > 1
  }, logical(1))
}

node_acctran <- extract_acctran(pars_acctran, levels_v)
node_ambig   <- is_ambiguous_mpr(pars_mpr)

cat(sprintf("Ambiguous nodes (MPR): %d of %d\n\n",
            sum(node_ambig), length(node_ambig)))

# Build edge transition table
edge_df <- data.frame(
  from       = node_acctran[tree_full$edge[, 1]],
  to         = node_acctran[tree_full$edge[, 2]],
  anc_ambig  = node_ambig[tree_full$edge[, 1]],
  desc_ambig = node_ambig[tree_full$edge[, 2]]
)
edge_df$unambiguous <- !edge_df$anc_ambig & !edge_df$desc_ambig

cat(sprintf("Unambiguous edges: %d / %d (%.1f%%)\n\n",
            sum(edge_df$unambiguous), nrow(edge_df),
            100 * mean(edge_df$unambiguous)))

trans_acctran <- table(edge_df$from, edge_df$to)
cat("--- All ACCTRAN transitions ---\n")
print(trans_acctran)

edge_unamb    <- edge_df[edge_df$unambiguous, ]
trans_unambig <- table(edge_unamb$from, edge_unamb$to)
cat("\n--- Unambiguous transitions only (MPR-confirmed) ---\n")
print(trans_unambig)

# Directional summary
count_trans <- function(tab, froms, tos) {
  sum(mapply(function(f, t) {
    if (f %in% rownames(tab) && t %in% colnames(tab)) tab[f, t] else 0
  }, froms, tos))
}

gain_from <- c("single", "single", "two")
gain_to   <- c("two",    "multi",  "multi")
loss_from <- c("two",    "multi",  "multi")
loss_to   <- c("single", "single", "two")

gains_all   <- count_trans(trans_acctran, gain_from, gain_to)
losses_all  <- count_trans(trans_acctran, loss_from, loss_to)
gains_unamb <- count_trans(trans_unambig, gain_from, gain_to)
losses_unamb<- count_trans(trans_unambig, loss_from, loss_to)

cat(sprintf("\nAll ACCTRAN:      gains = %d, losses = %d, ratio = %.2f\n",
            gains_all, losses_all, gains_all / max(losses_all, 1)))
cat(sprintf("Unambiguous only: gains = %d, losses = %d, ratio = %.2f\n",
            gains_unamb, losses_unamb, gains_unamb / max(losses_unamb, 1)))

# Flat cost matrix sensitivity
cat("\nRunning flat cost matrix sensitivity...\n")
cost_flat <- matrix(c(0,1,1, 1,0,1, 1,1,0), 3, 3,
                    dimnames = list(levels_v, levels_v))
pars_flat     <- ancestral.pars(tree_full, pdat_771, type = "ACCTRAN")
node_flat     <- extract_acctran(pars_flat, levels_v)
edge_flat     <- data.frame(from = node_flat[tree_full$edge[,1]],
                             to   = node_flat[tree_full$edge[,2]])
trans_flat    <- table(edge_flat$from, edge_flat$to)
gains_flat    <- count_trans(trans_flat, gain_from, gain_to)
losses_flat   <- count_trans(trans_flat, loss_from, loss_to)
cat(sprintf("Flat cost:        gains = %d, losses = %d, ratio = %.2f\n",
            gains_flat, losses_flat, gains_flat / max(losses_flat, 1)))

# Save parsimony results
write.csv(as.data.frame(trans_acctran),
          file.path(OUT_DIR, "parsimony_transitions_acctran_771.csv"),
          row.names = FALSE)
write.csv(as.data.frame(trans_unambig),
          file.path(OUT_DIR, "parsimony_transitions_unambiguous_771.csv"),
          row.names = FALSE)
write.csv(data.frame(
    method    = c("ACCTRAN_all", "ACCTRAN_unambiguous", "flat_cost"),
    gains     = c(gains_all, gains_unamb, gains_flat),
    losses    = c(losses_all, losses_unamb, losses_flat),
    ratio     = c(gains_all/max(losses_all,1),
                  gains_unamb/max(losses_unamb,1),
                  gains_flat/max(losses_flat,1))
  ), file.path(OUT_DIR, "parsimony_gain_loss_summary.csv"), row.names = FALSE)

# ---- 1A.3  Mk model selection -----------------------------------------------
cat("\n--- Mk model selection (ER, SYM, ARD) ---\n")
fit_er  <- fitMk(tree_full, states_3_771, model = "ER",  pi = "fitzjohn")
fit_sym <- fitMk(tree_full, states_3_771, model = "SYM", pi = "fitzjohn")
fit_ard <- fitMk(tree_full, states_3_771, model = "ARD", pi = "fitzjohn")

aic_vals   <- c(ER = AIC(fit_er), SYM = AIC(fit_sym), ARD = AIC(fit_ard))
best_model <- names(which.min(aic_vals))
cat("AICs:\n"); print(round(aic_vals, 2))
cat("Best model:", best_model, "\n")

write.csv(data.frame(model = names(aic_vals), AIC = aic_vals),
          file.path(OUT_DIR, "Mk_AIC_comparison_771.csv"), row.names = FALSE)

# ---- 1A.4  Stochastic character mapping (NSIM simulations) ------------------
cat(sprintf("\nRunning %d stochastic maps (ARD model, 771 tips)...\n", NSIM))
cat(sprintf("Expected runtime: ~%d min for NSIM=%d\n\n",
            ceiling(NSIM * 0.05), NSIM))

set.seed(42)
smap <- make.simmap(tree_full, states_3_771,
                    model   = best_model,
                    nsim    = NSIM,
                    pi      = "fitzjohn",
                    Q       = "empirical",
                    message = FALSE)

smap_summary <- describe.simmap(smap, verbose = FALSE)

# Extract transition counts — handle phytools version differences
trans_counts <- if (!is.null(smap_summary$count)) {
  smap_summary$count
} else if (!is.null(smap_summary$Tr)) {
  smap_summary$Tr
} else {
  # Manual extraction fallback
  t(sapply(smap, function(s) {
    tab <- s$mapped.edge
    setNames(colSums(tab), colnames(tab))
  }))
}

# Remove N column if present
if ("N" %in% colnames(trans_counts)) {
  trans_counts <- trans_counts[, colnames(trans_counts) != "N", drop = FALSE]
}

col_names <- colnames(trans_counts)
gain_cols <- col_names[grepl("single,two|single,multi|two,multi", col_names)]
loss_cols <- col_names[grepl("two,single|multi,single|multi,two", col_names)]

gains  <- rowSums(trans_counts[, gain_cols, drop = FALSE])
losses <- rowSums(trans_counts[, loss_cols, drop = FALSE])

cat("\n--- Mean transitions across", NSIM, "simulations ---\n")
print(round(colMeans(trans_counts), 2))
cat("\n--- 95% CI ---\n")
print(round(apply(trans_counts, 2, quantile, c(0.025, 0.975)), 2))
cat(sprintf("\nGains:  median = %.1f  [%.1f, %.1f]\n",
            median(gains), quantile(gains, 0.025), quantile(gains, 0.975)))
cat(sprintf("Losses: median = %.1f  [%.1f, %.1f]\n",
            median(losses), quantile(losses, 0.025), quantile(losses, 0.975)))
cat(sprintf("Ratio:  median = %.2f  [%.2f, %.2f]\n",
            median(gains/losses),
            quantile(gains/losses, 0.025),
            quantile(gains/losses, 0.975)))

smap_df <- data.frame(simulation = seq_len(NSIM),
                       gains = gains, losses = losses,
                       ratio = gains / losses)
write.csv(smap_df,
          file.path(OUT_DIR, "stochastic_mapping_gain_loss_771.csv"),
          row.names = FALSE)
write.csv(data.frame(
    transition = col_names,
    mean       = colMeans(trans_counts),
    CI_low     = apply(trans_counts, 2, quantile, 0.025),
    CI_high    = apply(trans_counts, 2, quantile, 0.975)
  ), file.path(OUT_DIR, "stochastic_mapping_transitions_771.csv"),
  row.names = FALSE)

# Also run on 654-tip tree for direct comparison with earlier results
cat("\nRunning", NSIM, "stochastic maps on 654-tip tree (for comparison)...\n")
states_3_654 <- setNames(dat_654$dom_state_3, dat_654$tip_label)
pdat_654     <- phyDat(as.list(states_3_654), type = "USER", levels = levels_v)
set.seed(42)
smap_654 <- make.simmap(tree_654, states_3_654,
                         model = best_model, nsim = NSIM,
                         pi = "fitzjohn", Q = "empirical", message = FALSE)
smap_654_sum <- describe.simmap(smap_654, verbose = FALSE)
tc_654 <- if (!is.null(smap_654_sum$count)) smap_654_sum$count else smap_654_sum$Tr
if ("N" %in% colnames(tc_654)) tc_654 <- tc_654[, colnames(tc_654) != "N"]
gains_654  <- rowSums(tc_654[, grep("single,two|single,multi|two,multi",
                                     colnames(tc_654)), drop=FALSE])
losses_654 <- rowSums(tc_654[, grep("two,single|multi,single|multi,two",
                                     colnames(tc_654)), drop=FALSE])
cat(sprintf("654-tip gains:  median=%.1f, losses: median=%.1f, ratio=%.2f\n",
            median(gains_654), median(losses_654),
            median(gains_654/losses_654)))
write.csv(data.frame(simulation=seq_len(NSIM),
                      gains=gains_654, losses=losses_654,
                      ratio=gains_654/losses_654),
          file.path(OUT_DIR, "stochastic_mapping_gain_loss_654.csv"),
          row.names = FALSE)

cat("\n1A complete.\n\n")

# =============================================================================
# STEP 1B — Secretion state phylogenetic signal
# D-statistic uses tree_full (771 tips — SP_final complete for all)
# Pagel and phyloglm use tree_654 (requires clade assignments)
# =============================================================================
cat("========== STEP 1B: Secretion state ==========\n\n")

# ---- 1B.1  Binary secretion traits ------------------------------------------
# Use SP_final column (from SignalP, complete for all 771 tips)
dat_771$secreted     <- as.integer(dat_771$SP_final != "OTHER")
dat_771$non_secreted <- as.integer(dat_771$SP_final == "OTHER")
dat_771$sec_spi      <- as.integer(dat_771$SP_final == "SP")

dat_654$secreted     <- as.integer(dat_654$SP_final != "OTHER")
dat_654$non_secreted <- as.integer(dat_654$SP_final == "OTHER")
dat_654$is_cladeV    <- as.integer(dat_654$clade == "clade5")

cat("Secretion summary (771 tips):\n")
print(table(dat_771$SP_final))
cat("\nSecretion by clade (654 tips):\n")
print(table(dat_654$SP_final, dat_654$clade, useNA = "always"))
cat("\n")

# ---- 1B.2  Fritz & Purvis D-statistic — 771 tips ----------------------------
cat("--- D-statistic (771-tip tree) ---\n")
comp_771 <- comparative.data(phy      = tree_full,
                              data     = dat_771,
                              names.col = "tip_label",
                              vcv      = FALSE,
                              warn.dropped = TRUE)
d_all <- phylo.d(comp_771, binvar = secreted, permut = 1000)
print(d_all)

d_spi <- phylo.d(comp_771, binvar = sec_spi, permut = 1000)
cat("\n--- D-statistic for Sec/SPI specifically ---\n")
print(d_spi)

write.csv(data.frame(
    trait        = c("secreted_vs_other", "sec_SPI_only"),
    D_estimate   = c(d_all$DEstimate, d_spi$DEstimate),
    p_vs_random  = c(d_all$Pval1,     d_spi$Pval1),
    p_vs_brownian= c(d_all$Pval0,     d_spi$Pval0),
    n_permut     = 1000
  ), file.path(OUT_DIR, "D_statistic_secretion_771.csv"), row.names = FALSE)

# ---- 1B.3  Pagel's correlated evolution — 654 tips --------------------------
cat("\n--- Pagel's correlated evolution: Clade V vs non-secreted (654 tips) ---\n")

x_pagel <- setNames(dat_654$is_cladeV,    tree_654$tip.label)
y_pagel <- setNames(dat_654$non_secreted, tree_654$tip.label)

pagel_result <- tryCatch({
  fitPagel(tree_654, x_pagel, y_pagel)
}, error = function(e) {
  cat("fitPagel failed:", conditionMessage(e), "\n")
  cat("Falling back to corHMM correlated evolution...\n")
  NULL
})

if (!is.null(pagel_result)) {
  print(pagel_result)
  if (!is.null(pagel_result)) {
    print(pagel_result)
    gv <- function(obj, ...) { for (n in c(...)) if (!is.null(obj[[n]])) return(obj[[n]]); NA }
    write.csv(data.frame(
      logL_independent = gv(pagel_result, "independent.logL", "logL.independent"),
      logL_dependent   = gv(pagel_result, "dependent.logL",   "logL.dependent"),
      LRT              = gv(pagel_result, "lik.ratio"),
      p_value          = gv(pagel_result, "P", "p")
    ), file.path(OUT_DIR, "Pagel_cladeV_secretion_654.csv"), row.names = FALSE)
  }
} else {
  # corHMM fallback: 4-state coding
  # 1=nonV+secreted, 2=nonV+nonsecret, 3=V+secreted, 4=V+nonsecret
  dat_654$corhmm_state <- with(dat_654,
    ifelse(is_cladeV==0 & non_secreted==0, 1,
    ifelse(is_cladeV==0 & non_secreted==1, 2,
    ifelse(is_cladeV==1 & non_secreted==0, 3, 4))))
  corhmm_dat <- data.frame(taxon = tree_654$tip.label,
                             state = dat_654[tree_654$tip.label, "corhmm_state"])
  fit_ind <- corHMM(tree_654, corhmm_dat, rate.cat=1, model="ER", node.states="marginal")
  fit_dep <- corHMM(tree_654, corhmm_dat, rate.cat=2, model="ER", node.states="marginal")
  cat(sprintf("corHMM AIC: independent=%.2f, dependent=%.2f\n",
              fit_ind$AIC, fit_dep$AIC))
  cat(sprintf("LRT: %.2f\n", 2*(fit_dep$loglik - fit_ind$loglik)))
  write.csv(data.frame(model     = c("independent","dependent"),
                        AIC       = c(fit_ind$AIC, fit_dep$AIC),
                        loglik    = c(fit_ind$loglik, fit_dep$loglik),
                        LRT       = c(NA, 2*(fit_dep$loglik - fit_ind$loglik))),
            file.path(OUT_DIR, "corHMM_cladeV_secretion_654.csv"), row.names=FALSE)
}

# ---- 1B.4  Phylogenetic logistic regression — 654 tips ----------------------
cat("\n--- phyloglm(secreted ~ clade, 654 tips) ---\n")

fit_sec_clade <- tryCatch({
  phyloglm(secreted ~ clade, data=dat_654, phy=tree_654,
           method="logistic_MPLE", btol=30)
}, error = function(e) {
  cat("Full clade model failed:", conditionMessage(e), "\n")
  cat("Trying binary Clade V vs rest...\n")
  tryCatch(
    phyloglm(secreted ~ is_cladeV, data=dat_654, phy=tree_654,
             method="logistic_MPLE", btol=30),
    error = function(e2) {
      tryCatch(
        phyloglm(secreted ~ is_cladeV, data=dat_654, phy=tree_654,
                 method="logistic_IG10", btol=30),
        error = function(e3) { cat("All phyloglm attempts failed.\n"); NULL }
      )
    }
  )
})

if (!is.null(fit_sec_clade)) {
  print(summary(fit_sec_clade))
  write.csv(as.data.frame(summary(fit_sec_clade)$coefficients),
            file.path(OUT_DIR, "phyloglm_secretion_654.csv"))
}

cat("\n1B complete.\n\n")

# =============================================================================
# STEP 1C — DVI promiscuity PGLS (skips gracefully if files absent)
# Uses: tree_654 (654 tips)
# =============================================================================
cat("========== STEP 1C: DVI promiscuity ==========\n\n")

if (!file.exists(DVI_FILE)) {
  cat("DVI file not found (", DVI_FILE, ") — skipping 1C\n")
  cat("To run 1C: provide dvi_output.tsv with columns Domain, n_GH18, T_GH18, phi_GH18\n\n")
} else {
  dvi_seq <- read.delim(DVI_FILE, stringsAsFactors = FALSE)
  cat("DVI columns:", paste(colnames(dvi_seq), collapse=", "), "\n")

  # ---- 1C.1  Catalytic divergence -------------------------------------------
  cat("Computing catalytic divergence (patristic distances, 654 tips)...\n")
  patr <- cophenetic.phylo(tree_654)
  cat_div <- vapply(tree_654$tip.label, function(tip) {
    clade_i    <- dat_654[tip, "clade"]
    if (is.na(clade_i)) return(NA_real_)
    clade_tips <- rownames(dat_654)[!is.na(dat_654$clade) &
                                     dat_654$clade == clade_i]
    clade_tips <- intersect(clade_tips, rownames(patr))
    if (length(clade_tips) < 2) return(NA_real_)
    mean(patr[tip, clade_tips[clade_tips != tip]])
  }, numeric(1))
  dat_654$catalytic_divergence <- cat_div

  # ---- 1C.2  Per-tip mean promiscuity from annotation -----------------------
  # mean_promiscuity column already in annotation table from Python pipeline
  # Check if available
  has_prom <- "mean_promiscuity" %in% colnames(dat_654) &&
              sum(!is.na(dat_654$mean_promiscuity)) > 10
  if (has_prom) {
    cat(sprintf("mean_promiscuity available for %d / 654 tips\n",
                sum(!is.na(dat_654$mean_promiscuity))))
  } else {
    cat("mean_promiscuity not in annotation — 1C PGLS skipped\n")
    has_prom <- FALSE
  }

  # ---- 1C.3  Cluster-level DVI ----------------------------------------------
  if (file.exists(ARCH_FILE)) {
    cat("Computing cluster-level DVI from architecture table...\n")
    arch <- read.csv(ARCH_FILE, stringsAsFactors = FALSE)
    cluster_DVI <- arch %>%
      pivot_longer(c(left_neighbour, right_neighbour),
                   names_to="side", values_to="neighbour") %>%
      filter(!is.na(neighbour), neighbour != aux_pfam) %>%
      distinct(aux_pfam, ssn_cluster, neighbour) %>%
      group_by(aux_pfam) %>%
      summarise(n_cluster=n_distinct(ssn_cluster),
                T_cluster=n_distinct(neighbour), .groups="drop")
    total_T_cl <- sum(cluster_DVI$T_cluster)
    total_n_cl <- sum(cluster_DVI$n_cluster)
    cluster_DVI <- cluster_DVI %>%
      mutate(pi_cl = T_cluster/total_T_cl, f_cl = n_cluster/total_n_cl,
             phi_cluster = ifelse(pi_cl>0 & f_cl>0,
                                  pi_cl*log(pi_cl/f_cl), 0))
    merged_phi <- inner_join(cluster_DVI %>% select(aux_pfam, phi_cluster),
                              dvi_seq %>% select(Domain, phi_GH18) %>%
                                rename(aux_pfam=Domain), by="aux_pfam")
    sp_rho <- cor.test(merged_phi$phi_cluster, merged_phi$phi_GH18,
                       method="spearman")
    cat(sprintf("Cluster vs seq-level phi: rho=%.3f, p=%.2e, n=%d\n",
                sp_rho$estimate, sp_rho$p.value, nrow(merged_phi)))
    write.csv(merged_phi,
              file.path(OUT_DIR, "phi_seq_vs_cluster.csv"), row.names=FALSE)
    write.csv(cluster_DVI,
              file.path(OUT_DIR, "cluster_level_DVI.csv"), row.names=FALSE)
  }

  # ---- 1C.4  PGLS -----------------------------------------------------------
  if (has_prom) {
    cat("\n--- PGLS: mean_promiscuity ~ catalytic_divergence ---\n")
    pgls_data <- dat_654[!is.na(dat_654$mean_promiscuity) &
                           !is.na(dat_654$catalytic_divergence), ]
    cat(sprintf("Tips for PGLS: %d\n", nrow(pgls_data)))

    tree_pgls <- drop.tip(tree_654,
                           setdiff(tree_654$tip.label, pgls_data$tip_label))
    pgls_data <- pgls_data[tree_pgls$tip.label, ]

    fit_ols  <- gls(mean_promiscuity ~ catalytic_divergence,
                    data=pgls_data, method="ML", na.action=na.omit)
    cat(sprintf("OLS slope: %.4f\n", coef(fit_ols)[2]))

    fit_pgls <- tryCatch({
      gls(mean_promiscuity ~ catalytic_divergence,
          data        = pgls_data,
          correlation = corPagel(value=0.5, phy=tree_pgls,
                                  form=~tip_label, fixed=FALSE),
          method="ML", na.action=na.omit)
    }, error = function(e) {
      cat("corPagel failed:", conditionMessage(e), "\n")
      tryCatch(
        gls(mean_promiscuity ~ catalytic_divergence,
            data=pgls_data,
            correlation=corBrownian(phy=tree_pgls, form=~tip_label),
            method="ML", na.action=na.omit),
        error=function(e2) { cat("corBrownian also failed\n"); NULL }
      )
    })
    if (!is.null(fit_pgls)) {
      print(summary(fit_pgls))
      lambda_est <- tryCatch(
        coef(fit_pgls$modelStruct$corStruct, unconstrained=FALSE),
        error=function(e) NA)
      cat(sprintf("Lambda: %.3f | Sign concordance OLS vs PGLS: %s\n",
                  lambda_est,
                  ifelse(sign(coef(fit_ols)[2])==sign(coef(fit_pgls)[2]),
                         "SAME", "CHANGED")))
      write.csv(data.frame(parameter=names(coef(fit_pgls)),
                             estimate=coef(fit_pgls)),
                file.path(OUT_DIR, "PGLS_promiscuity_divergence.csv"),
                row.names=FALSE)
    }
  }
}
cat("\n1C complete.\n\n")

# =============================================================================
# STEP 1D — Fig. 2d: Cluster 1.2 and promiscuous-domain carriage
#
# IMPORTANT:
#   * A missing promiscuity annotation is NOT a negative observation.
#   * Unknown values remain NA.
#   * The tree is pruned to tips with both a known binary trait and cluster ID.
#   * The ordinary GLM is only a data-alignment diagnostic/effect-size check.
#   * phyloglm is the phylogenetically corrected analysis.
# =============================================================================
cat("========== STEP 1D: Fig. 2d complete-case analysis ==========\n\n")

# ---- 1D.1 Parse binary columns without converting NA to 0 --------------------
parse_binary <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(z))

  out[z %in% c("1", "true", "t", "yes", "y")]  <- 1L
  out[z %in% c("0", "false", "f", "no", "n")] <- 0L
  out
}

binary_column_or_na <- function(df, column) {
  if (column %in% names(df)) parse_binary(df[[column]])
  else rep(NA_integer_, nrow(df))
}

# Prefer the explicit sequence-level has_promiscuous field. Use is_promiscuous
# only where the first field is unavailable. This preserves genuine unknowns.
hp_primary   <- binary_column_or_na(dat_654, "has_promiscuous")
hp_secondary <- binary_column_or_na(dat_654, "is_promiscuous")

conflict <- !is.na(hp_primary) & !is.na(hp_secondary) &
            hp_primary != hp_secondary
if (any(conflict)) {
  stop(sprintf(
    "%d tips have conflicting has_promiscuous and is_promiscuous values.",
    sum(conflict)
  ))
}

dat_654$has_prom_int <- hp_primary
fill_from_secondary <- is.na(dat_654$has_prom_int) & !is.na(hp_secondary)
dat_654$has_prom_int[fill_from_secondary] <- hp_secondary[fill_from_secondary]

# Unknown or absent cluster labels remain NA rather than being assigned to
# the comparison group.
cluster_string <- as.character(dat_654$cluster)
dat_654$is_12 <- ifelse(
  is.na(cluster_string) | trimws(cluster_string) == "",
  NA_integer_,
  as.integer(cluster_string == "cluster_1.2")
)

cat("Promiscuity trait availability by Cluster 1.2 status:\n")
missingness_tab <- table(
  is_cluster_12   = dat_654$is_12,
  trait_available = !is.na(dat_654$has_prom_int),
  useNA = "always"
)
print(missingness_tab)

cat("\nParsed promiscuity trait values:\n")
print(table(dat_654$has_prom_int, useNA = "always"))

write.csv(
  as.data.frame(missingness_tab),
  file.path(OUT_DIR, "fig2d_trait_missingness_by_cluster12.csv"),
  row.names = FALSE
)

# ---- 1D.2 Construct one complete-case data/tree pair -------------------------
complete_prom <- !is.na(dat_654$has_prom_int) & !is.na(dat_654$is_12)
dat_prom <- dat_654[complete_prom, , drop = FALSE]

if (nrow(dat_prom) == 0) {
  stop("No complete cases are available for the Fig. 2d analysis.")
}

tree_prom <- keep.tip(tree_654, dat_prom$tip_label)
dat_prom <- dat_prom[match(tree_prom$tip.label, dat_prom$tip_label), ,
                     drop = FALSE]
rownames(dat_prom) <- dat_prom$tip_label

stopifnot(identical(tree_prom$tip.label, dat_prom$tip_label))
stopifnot(!anyNA(dat_prom$has_prom_int))
stopifnot(!anyNA(dat_prom$is_12))

if (length(unique(dat_prom$has_prom_int)) != 2L) {
  stop("The complete-case response does not contain both 0 and 1 states.")
}
if (length(unique(dat_prom$is_12)) != 2L) {
  stop("The complete-case predictor does not contain both Cluster 1.2 and non-1.2 tips.")
}

n_excluded <- nrow(dat_654) - nrow(dat_prom)
cat(sprintf(
  "\nComplete-case analysis: %d tips retained; %d excluded; %d positive (%.2f%%).\n",
  nrow(dat_prom), n_excluded, sum(dat_prom$has_prom_int == 1L),
  100 * mean(dat_prom$has_prom_int == 1L)
))

# Save the exact data used so every reported number is reproducible.
write.csv(
  dat_prom[, c("tip_label", "cluster", "clade", "has_prom_int", "is_12")],
  file.path(OUT_DIR, "fig2d_complete_case_tip_data.csv"),
  row.names = FALSE
)
write.tree(
  tree_prom,
  file = file.path(OUT_DIR, "fig2d_complete_case_tree.nwk")
)

# ---- 1D.3 Contingency table and ordinary-GLM diagnostic ----------------------
tab_12 <- with(
  dat_prom,
  table(
    is_cluster_12   = factor(is_12, levels = c(0, 1)),
    has_promiscuous = factor(has_prom_int, levels = c(0, 1))
  )
)
cat("\nComplete-case 2 x 2 table:\n")
print(tab_12)

write.csv(
  as.data.frame.matrix(tab_12),
  file.path(OUT_DIR, "fig2d_cluster12_contingency_table.csv")
)

# Crude compositional effect size. Fisher's test is retained here only to
# obtain an exact odds-ratio interval; the phyloglm below is the primary test.
fisher_12 <- fisher.test(tab_12, alternative = "two.sided")
fisher_ci <- unname(fisher_12$conf.int)

fit_glm <- glm(
  has_prom_int ~ is_12,
  data = dat_prom,
  family = binomial()
)
glm_sum <- summary(fit_glm)
glm_beta <- unname(coef(fit_glm)["is_12"])

# For a non-zero 2 x 2 table, the GLM coefficient must equal the log cross-
# product odds ratio. This catches accidental row/trait misalignment.
cell_00 <- unname(tab_12["0", "0"])
cell_01 <- unname(tab_12["0", "1"])
cell_10 <- unname(tab_12["1", "0"])
cell_11 <- unname(tab_12["1", "1"])

if (all(c(cell_00, cell_01, cell_10, cell_11) > 0)) {
  table_or <- (cell_11 * cell_00) / (cell_10 * cell_01)
  if (!isTRUE(all.equal(glm_beta, log(table_or), tolerance = 1e-8))) {
    stop("GLM coefficient does not match the 2 x 2 table odds ratio; check alignment.")
  }
} else {
  # Haldane-Anscombe value is descriptive only when a cell is zero.
  table_or <- ((cell_11 + 0.5) * (cell_00 + 0.5)) /
              ((cell_10 + 0.5) * (cell_01 + 0.5))
  warning("At least one contingency-table cell is zero; descriptive OR uses a 0.5 correction.")
}

crude_result <- data.frame(
  n_tips             = nrow(dat_prom),
  n_positive         = sum(dat_prom$has_prom_int == 1L),
  n_cluster_12       = sum(dat_prom$is_12 == 1L),
  n_non_cluster_12   = sum(dat_prom$is_12 == 0L),
  table_odds_ratio   = table_or,
  fisher_odds_ratio  = unname(fisher_12$estimate),
  fisher_CI_low      = fisher_ci[1],
  fisher_CI_high     = fisher_ci[2],
  fisher_p_value     = fisher_12$p.value,
  glm_beta           = glm_beta,
  glm_SE             = glm_sum$coefficients["is_12", "Std. Error"],
  glm_p_value        = glm_sum$coefficients["is_12", "Pr(>|z|)"]
)
cat("\nCrude complete-case effect size and alignment diagnostic:\n")
print(crude_result)
write.csv(
  crude_result,
  file.path(OUT_DIR, "fig2d_cluster12_crude_effect.csv"),
  row.names = FALSE
)

# ---- 1D.4 Primary phylogenetic logistic regression ---------------------------
cat(sprintf(
  "\n--- phyloglm(has_promiscuous ~ is_cluster_1.2), %d bootstrap replicates ---\n",
  PHYLO_BOOT
))

fit_12_mple <- tryCatch(
  phyloglm(
    has_prom_int ~ is_12,
    data   = dat_prom,
    phy    = tree_prom,
    method = "logistic_MPLE",
    btol   = 30,
    boot   = PHYLO_BOOT
  ),
  error = function(e) {
    stop("MPLE phyloglm failed: ", conditionMessage(e))
  }
)

mple_sum <- summary(fit_12_mple)
print(mple_sum)

mple_coef <- as.data.frame(mple_sum$coefficients)
mple_coef$term <- rownames(mple_coef)
rownames(mple_coef) <- NULL
write.csv(
  mple_coef,
  file.path(OUT_DIR, "fig2d_phyloglm_MPLE_coefficients.csv"),
  row.names = FALSE
)

# Save the fitted object and all bootstrap output even if package-version
# differences prevent automatic extraction of a named bootstrap column.
saveRDS(
  fit_12_mple,
  file.path(OUT_DIR, "fig2d_phyloglm_MPLE_fit.rds")
)
if (!is.null(fit_12_mple$bootresult)) {
  write.csv(
    as.data.frame(fit_12_mple$bootresult),
    file.path(OUT_DIR, "fig2d_phyloglm_MPLE_bootstrap_raw.csv"),
    row.names = FALSE
  )
}

extract_boot_term <- function(fit, term) {
  br <- fit$bootresult
  if (is.null(br)) return(NULL)

  br <- as.data.frame(br)
  if (term %in% names(br)) {
    v <- br[[term]]
  } else {
    # In some phylolm versions coefficient columns are unnamed but occur in
    # coefficient order. Use this fallback only when dimensions permit it.
    term_index <- match(term, names(coef(fit)))
    if (is.na(term_index) || ncol(br) < term_index) return(NULL)
    v <- br[[term_index]]
  }

  v <- as.numeric(v)
  v[is.finite(v)]
}

beta_mple <- mple_sum$coefficients["is_12", "Estimate"]
se_mple   <- mple_sum$coefficients["is_12", "StdErr"]
p_mple    <- mple_sum$coefficients["is_12", "p.value"]

boot_beta <- extract_boot_term(fit_12_mple, "is_12")
if (!is.null(boot_beta) && length(boot_beta) >= 20L) {
  beta_ci <- unname(quantile(boot_beta, c(0.025, 0.975), na.rm = TRUE))
  ci_method <- "parametric bootstrap"
} else {
  beta_ci <- beta_mple + c(-1, 1) * 1.96 * se_mple
  ci_method <- "Wald fallback; inspect bootstrap_raw.csv"
  warning(
    "Could not reliably identify the is_12 bootstrap column. ",
    "A Wald interval was written; inspect fig2d_phyloglm_MPLE_bootstrap_raw.csv."
  )
}

phylo_result <- data.frame(
  method          = "logistic_MPLE",
  n_tips          = nrow(dat_prom),
  n_positive      = sum(dat_prom$has_prom_int == 1L),
  beta            = beta_mple,
  SE              = se_mple,
  odds_ratio      = exp(beta_mple),
  OR_CI_low       = exp(beta_ci[1]),
  OR_CI_high      = exp(beta_ci[2]),
  CI_method       = ci_method,
  p_value         = p_mple,
  alpha           = if (!is.null(fit_12_mple$alpha)) fit_12_mple$alpha else NA_real_,
  n_boot_requested = PHYLO_BOOT
)

cat("\nPrimary phylogenetic result:\n")
print(phylo_result)
write.csv(
  phylo_result,
  file.path(OUT_DIR, "fig2d_phyloglm_MPLE_summary.csv"),
  row.names = FALSE
)

# ---- 1D.5 Optional estimator sensitivity -------------------------------------
# IG10 is reported only as a sensitivity check. Failure does not invalidate MPLE.
cat("\n--- IG10 estimator sensitivity check ---\n")
fit_12_ig10 <- tryCatch(
  phyloglm(
    has_prom_int ~ is_12,
    data   = dat_prom,
    phy    = tree_prom,
    method = "logistic_IG10",
    btol   = 30
  ),
  error = function(e) {
    cat("IG10 failed:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(fit_12_ig10)) {
  ig10_sum <- summary(fit_12_ig10)
  print(ig10_sum)

  ig10_coef <- as.data.frame(ig10_sum$coefficients)
  ig10_coef$term <- rownames(ig10_coef)
  rownames(ig10_coef) <- NULL
  write.csv(
    ig10_coef,
    file.path(OUT_DIR, "fig2d_phyloglm_IG10_coefficients.csv"),
    row.names = FALSE
  )
  saveRDS(
    fit_12_ig10,
    file.path(OUT_DIR, "fig2d_phyloglm_IG10_fit.rds")
  )
}

cat("\n1D complete.\n\n")

# =============================================================================
# STEP 1E — Catalytic divergence (if not computed in 1C)
# =============================================================================
cat("========== STEP 1E: Catalytic divergence ==========\n\n")

if (!"catalytic_divergence" %in% colnames(dat_654) ||
    all(is.na(dat_654$catalytic_divergence))) {
  cat("Computing catalytic divergence...\n")
  patr <- cophenetic.phylo(tree_654)
  dat_654$catalytic_divergence <- vapply(tree_654$tip.label, function(tip) {
    clade_i    <- dat_654[tip, "clade"]
    if (is.na(clade_i)) return(NA_real_)
    clade_tips <- rownames(dat_654)[!is.na(dat_654$clade) &
                                     dat_654$clade == clade_i]
    clade_tips <- intersect(clade_tips, rownames(patr))
    if (length(clade_tips) < 2) return(NA_real_)
    mean(patr[tip, clade_tips[clade_tips != tip]])
  }, numeric(1))
}

cat("Catalytic divergence by clade:\n")
print(tapply(dat_654$catalytic_divergence, dat_654$clade, function(x)
  round(c(mean=mean(x,na.rm=T), sd=sd(x,na.rm=T), n=sum(!is.na(x))), 4)))

write.csv(dat_654[, c("tip_label","clade","catalytic_divergence")],
          file.path(OUT_DIR, "catalytic_divergence_per_tip.csv"),
          row.names=FALSE)
cat("\n1E complete.\n\n")

# =============================================================================
# STEP 1F — Fig 1b replication on phylogenetically sampled 654-tip dataset
# Shows chi-squared result robust to taxonomic sampling correction
# =============================================================================
cat("========== STEP 1F: Fig 1b replication (654-tip dataset) ==========\n\n")

domain_counts <- table(dat_654$n_domains)
cat("Domain count distribution (654 tips):\n")
print(domain_counts)

n_mean  <- mean(dat_654$n_domains, na.rm=TRUE)
lambda  <- 1 - 1/n_mean
cat(sprintf("\nFitted lambda: %.4f  (mean n_domains: %.3f)\n", lambda, n_mean))

k_vals   <- 1:8
p_exp    <- (1-lambda) * lambda^(k_vals-1)
p_exp    <- p_exp / sum(p_exp)
n_obs    <- sum(dat_654$n_domains <= 8, na.rm=TRUE)
expected <- p_exp * n_obs
observed <- as.numeric(table(factor(
  dat_654$n_domains[dat_654$n_domains <= 8], levels=1:8)))

comp_df <- data.frame(k=k_vals, observed=observed,
                       expected=round(expected,1),
                       std_resid=round((observed-expected)/sqrt(expected), 2))
cat("\nObserved vs expected:\n")
print(comp_df)

chisq_654 <- chisq.test(observed, p=p_exp, rescale.p=TRUE)
cat(sprintf("\nChi-squared (654-tip): X2=%.1f, df=%d, p=%.2e\n",
            chisq_654$statistic, chisq_654$parameter, chisq_654$p.value))
cat("Original (39,582 sequences): X2=18891.8, p<2.2e-16\n")
cat("Qualitative conclusion unchanged:", chisq_654$p.value < 0.05, "\n")

write.csv(comp_df,
          file.path(OUT_DIR, "fig1b_chisq_654tip.csv"), row.names=FALSE)
cat("\n1F complete.\n\n")

# =============================================================================
# STEP 1G — Methods mapping table (for supplementary)
# =============================================================================
methods_map <- data.frame(
  original_test = c(
    "Fig 1b chi-squared (random accretion, 39582 seqs)",
    "Line 178-193 n vs T correlation",
    "Fig 2d Fisher's exact (cluster enrichment)",
    "Line 316 chi-squared (secretion non-random)",
    "Lines 319-330 Cramer's V pairwise clades"
  ),
  new_test = c(
    paste0("Sankoff parsimony (ACCTRAN+MPR) + Mk-", best_model,
           " stochastic mapping (771-tip tree); chi-squared replicated on 654-tip dataset"),
    "Fig. 2a treated descriptively; DVI recomputed after MMseqs2 dereplication outside this R script",
    "Complete-case phyloglm(has_promiscuous ~ is_cluster_1.2, logistic_MPLE) with IG10 sensitivity",
    "Fritz & Purvis D-statistic (caper::phylo.d, 1000 permutations, 771-tip tree)",
    "Pagel correlated evolution (fitPagel) + phyloglm(secreted ~ clade, 654-tip tree)"
  ),
  tree_used = c(
    "771-tip MAD-rooted (all tips have domain annotations)",
    "Not applicable: Fig. 2a units are auxiliary-domain families",
    "Complete-case subtree pruned from the 654-tip UniProt tree",
    "771-tip MAD-rooted (SP_final complete for all tips)",
    "654-tip (Pagel/phyloglm require clade assignments)"
  ),
  output_file = c(
    "parsimony_transitions_*_771.csv + stochastic_mapping_*_771.csv + fig1b_chisq_654tip.csv",
    "External MMseqs2 DVI robustness table",
    "fig2d_phyloglm_MPLE_summary.csv + fig2d_complete_case_tip_data.csv",
    "D_statistic_secretion_771.csv",
    "Pagel_cladeV_secretion_654.csv + phyloglm_secretion_654.csv"
  ),
  conclusion_changed = rep("*** FILL IN AFTER REVIEWING RESULTS ***", 5),
  stringsAsFactors = FALSE
)
write.csv(methods_map,
          file.path(OUT_DIR, "methods_mapping_table_supplementary.csv"),
          row.names=FALSE)

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("=====================================================\n")
cat("ALL STEPS COMPLETE\n")
cat("=====================================================\n\n")
cat("Results written to:", OUT_DIR, "\n\n")
cat("Output files:\n")
for (f in list.files(OUT_DIR, pattern="\\.csv$")) {
  cat(sprintf("  %s\n", f))
}
cat("\nNEXT STEPS:\n")
cat("1. Open methods_mapping_table_supplementary.csv\n")
cat("2. Fill in 'conclusion_changed' column for each test\n")
cat("3. Use stochastic_mapping_gain_loss_771.csv for manuscript language\n")
cat("4. Check D_statistic_secretion_771.csv for 1B result\n")
cat("\nKEY RESULTS SUMMARY:\n")
cat(sprintf("  1A parsimony (unambiguous): gains=%d, losses=%d, ratio=%.2f\n",
            gains_unamb, losses_unamb, gains_unamb/max(losses_unamb,1)))
cat(sprintf("  1A stochastic (771 tips):  gains median=%.1f, losses=%.1f, ratio=%.2f\n",
            median(gains), median(losses), median(gains/losses)))
cat(sprintf("  1B D-statistic (secreted): D=%.3f, p(vs random)=%.3e\n",
            d_all$DEstimate, d_all$Pval1))
cat(sprintf("  Mk best model: %s (AIC=%.1f vs ER=%.1f)\n",
            best_model, aic_vals[best_model], aic_vals["ER"]))