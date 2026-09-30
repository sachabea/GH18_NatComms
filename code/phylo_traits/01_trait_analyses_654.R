## =============================================================================
## All trait analyses on the 654-tip canonical set
##
## 654 = all UniProt tree tips with a clade/cluster assignment (771 - 117 NCBI).
## Trait sources on this set:
##   secretion, promiscuity, clade/cluster : from tip annotation (654/654)
##   auxiliary-domain category carriage    : DERIVED from the tip 'domains'
##       column via signature_broad_categories_2.csv (100% concordant with the
##       master-annotation categories)
##   phylum                                : 630/654 (24 tips unassigned)
##
## PROMISCUITY COLUMN: uses 'promiscuous' (= the DVI-pipeline is_promiscuous flag
## that drives Figure 2; 28 tips). The stricter 'promiscuous_stricter' column
## (21 tips) is kept in the table for reference only and is NOT analysed.
##
## => secretion, promiscuity, function-category and catalytic-divergence
##    analyses run at n = 654; phylum one-vs-rest traits run at n = 630
##    (phylum-less tips are NA and dropped, NOT coded as absent).
##
## InDel class is NOT used by any statistic here (D is tree-based; phyloglm uses
## clade), so the 80 tips lacking an InDel assignment are retained.
##
## Kept on their own N by design:
##   domain-count parsimony/stochastic mapping -> 771 (observed for all tips)
##   continuous mean_promiscuity ~ divergence PGLS -> its own N (score sparse)
## =============================================================================

suppressPackageStartupMessages({
  library(ape); library(caper); library(phylolm); library(phytools); library(dplyr)
})

TREE_FILE  <- "tree.nwk"
TRAIT_FILE <- "tree_tip_traits_CANONICAL_654.csv"
PROM_COL   <- "promiscuous"        # authoritative promiscuity flag (= is_promiscuous)
OUT_DIR    <- "canonical654_results"; dir.create(OUT_DIR, showWarnings = FALSE)
NPERM      <- 1000
MIN_TIPS   <- 15
set.seed(1)

tree <- read.tree(TREE_FILE)
d    <- read.csv(TRAIT_FILE, stringsAsFactors = FALSE, na.strings = c("NA","None",""))
common <- intersect(tree$tip.label, d$tip_label)
tree <- keep.tip(tree, common)
d <- d[match(tree$tip.label, d$tip_label), ]; rownames(d) <- d$tip_label
stopifnot(identical(tree$tip.label, d$tip_label))
cat(sprintf("Canonical tree: %d tips\n", length(common)))
cat("Clade sizes:\n"); print(table(d$clade))
cat(sprintf("Phylum present: %d / %d | promiscuous=1: %d\n",
            sum(!is.na(d$phylum) & d$phylum != ""), nrow(d), sum(d[[PROM_COL]] == 1, na.rm = TRUE)))

## one-vs-rest phylum columns: NA where phylum is missing (so those tips drop out)
d$phylum[d$phylum == ""] <- NA
for (p in names(which(table(d$phylum) >= MIN_TIPS))) {
  col <- paste0("phy_", make.names(p))
  d[[col]] <- ifelse(is.na(d$phylum), NA_integer_, as.integer(d$phylum == p))
}
cat_cols <- grep("^cat_", names(d), value = TRUE)
phy_cols <- grep("^phy_", names(d), value = TRUE)
d$prom   <- d[[PROM_COL]]
traits   <- c("secreted", "prom", cat_cols, phy_cols)

## ---- A + B : D-statistic and phyloglm(trait ~ clade) per trait ---------------
run_trait <- function(trait) {
  y <- d[[trait]]
  keep <- !is.na(y)                                   # drop NA tips (phylum-less)
  if (min(sum(y == 1, na.rm = TRUE), sum(y == 0, na.rm = TRUE)) < MIN_TIPS) return(NULL)
  dsub <- d[keep, c("tip_label", trait)]; names(dsub)[2] <- "y"
  tr   <- keep.tip(tree, dsub$tip_label)
  cp <- comparative.data(tr, dsub, names.col = "tip_label")
  dd <- phylo.d(cp, binvar = y, permut = NPERM)
  pg <- NA_real_
  fit <- tryCatch(phyloglm(as.formula(paste(trait, "~ clade")), data = d[keep, ],
                           phy = keep.tip(tree, d$tip_label[keep]),
                           method = "logistic_MPLE", btol = 30), error = function(e) NULL)
  if (!is.null(fit)) {
    co <- summary(fit)$coefficients; rr <- grep("^clade", rownames(co))
    if (length(rr)) pg <- min(co[rr, "p.value"], na.rm = TRUE)
  }
  data.frame(trait = trait, n = sum(keep), n_pos = sum(y == 1, na.rm = TRUE),
             D = unname(dd$DEstimate), p_vs_random = unname(dd$Pval1),
             p_vs_brownian = unname(dd$Pval0), phyloglm_min_clade_p = pg)
}
sig <- do.call(rbind, lapply(traits, run_trait)) %>% arrange(D)
cat("\n=== A/B  Phylogenetic signal + clade association ===\n")
print(sig, row.names = FALSE, digits = 3)
write.csv(sig, file.path(OUT_DIR, "signal_and_clade_654.csv"), row.names = FALSE)

## ---- C : promiscuity phyloglm (n = 654) + parametric bootstrap CI -----------
cat("\n=== C  Promiscuity (n = 654) ===\n")
cat("is_12 counts:\n"); print(table(d$is_12))
NBOOT <- 1000
fit_12 <- tryCatch(phyloglm(prom ~ is_12, data = d, phy = tree,
                            method = "logistic_MPLE", btol = 30, boot = NBOOT),
                   error = function(e) NULL)
if (!is.null(fit_12)) {
  s <- summary(fit_12); print(s)
  co   <- s$coefficients
  beta <- co["is_12", "Estimate"]; p <- co["is_12", "p.value"]
  # version-robust extraction of the is_12 bootstrap interval
  if (all(c("lowerbootCI","upperbootCI") %in% colnames(co))) {
    b_lo <- co["is_12","lowerbootCI"]; b_hi <- co["is_12","upperbootCI"]
  } else if (!is.null(fit_12$bootstrap)) {
    bs  <- as.data.frame(fit_12$bootstrap)
    col <- if ("is_12" %in% names(bs)) bs[["is_12"]] else bs[[ which(names(coef(fit_12)) == "is_12") ]]
    qq  <- quantile(col, c(0.025, 0.975), na.rm = TRUE); b_lo <- qq[1]; b_hi <- qq[2]
  } else { b_lo <- b_hi <- NA_real_ }
  cat(sprintf("\nis_12: beta=%.4f  OR=%.3f  OR 95%% CI [%.3f, %.3f]  P=%.4f  alpha=%.3f\n",
              beta, exp(beta), exp(b_lo), exp(b_hi), p, fit_12$alpha))
  write.csv(data.frame(model="promiscuous ~ is_12", n=nrow(d), beta=beta,
                       OR=exp(beta), OR_CI_low=exp(b_lo), OR_CI_high=exp(b_hi),
                       p=p, alpha=fit_12$alpha, n_boot=NBOOT),
            file.path(OUT_DIR, "promiscuity_phyloglm_654.csv"), row.names = FALSE)
}
## ---- D : secretion Pagel + phyloglm (n = 654) --------------------------------
cat("\n=== D  Secretion (n = 654) ===\n")
x_cV <- setNames(as.integer(d$clade == "clade5"), tree$tip.label)
y_ns <- setNames(as.integer(d$secreted == 0),     tree$tip.label)
pag <- tryCatch(fitPagel(tree, x_cV, y_ns), error = function(e) NULL)
if (!is.null(pag)) { print(pag)
  write.csv(data.frame(LRT = pag$lik.ratio, p = pag$P),
            file.path(OUT_DIR, "secretion_pagel_654.csv"), row.names = FALSE) }
fit_sec <- tryCatch(phyloglm(secreted ~ clade, data = d, phy = tree,
                             method = "logistic_MPLE", btol = 30), error = function(e) NULL)
if (!is.null(fit_sec)) { print(summary(fit_sec))
  write.csv(as.data.frame(summary(fit_sec)$coefficients),
            file.path(OUT_DIR, "secretion_phyloglm_654.csv")) }

## ---- E : catalytic divergence per clade (n = 654) ----------------------------
cat("\n=== E  Catalytic divergence per clade ===\n")
patr <- cophenetic.phylo(tree)
d$cat_div <- vapply(tree$tip.label, function(tp) {
  cl <- d[tp, "clade"]; ct <- setdiff(rownames(d)[d$clade == cl], tp)
  if (length(ct) < 2) return(NA_real_); mean(patr[tp, ct])
}, numeric(1))
print(do.call(rbind, tapply(d$cat_div, d$clade, function(x)
  c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), n = sum(!is.na(x))))))
write.csv(data.frame(tip_label = d$tip_label, clade = d$clade, cat_div = d$cat_div),
          file.path(OUT_DIR, "catalytic_divergence_654.csv"), row.names = FALSE)

cat("\nAll analyses complete (n = 654; phylum traits n = 630). Results in", OUT_DIR, "\n")