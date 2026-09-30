## =============================================================================
## Figure 4 — phylogenetically-aware "how traits partition across the classes"
##
## Replaces the non-phylogenetic enrichment (chi2 / pairwise Fisher / one-vs-rest
## Fisher+FDR) behind Supp Tables 15-16 and Fig 9 with a per-tip, tree-aware
## treatment. This is the SAME template you already used for secretion:
##   headline  = Fritz & Purvis D (is the trait phylogenetically clustered?)
##   secondary = phyloglm(trait ~ clade)  (residual association after ancestry)
##
## Input: tree_tip_traits_for_phylo.csv  (574 tips; built by joining the master
##        annotation fastaID -> tree tip_label). One binary column per trait:
##          secreted, cat_* (auxiliary-function carriage), and phylum handled
##          as one-vs-rest below.
##        tree3.nwk  (your 771-tip tree; pruned to the 574 joined tips)
##
## WHY THIS ANSWERS REVIEWER 2:
##   * The unit is the SEQUENCE (a tip), not exploded domain-category rows, so
##     multi-category sequences are no longer counted multiple times.
##   * D quantifies phylogenetic clustering directly; a strongly clustered trait
##     means the class-wise "enrichment" is largely shared ancestry, and a naive
##     per-sequence Fisher p-value overstates the evidence.
## =============================================================================

suppressPackageStartupMessages({ library(ape); library(caper); library(phylolm); library(dplyr) })

TREE_FILE  <- "tree.nwk"
TRAIT_FILE <- "tree_tip_traits_CANONICAL_654.csv"
OUT        <- "fig4_trait_partition_phylo_results.csv"
NPERM      <- 1000
MIN_TIPS   <- 15      # minimum tips in the minority state to bother testing
set.seed(1)

tree <- read.tree(TREE_FILE)
d    <- read.csv(TRAIT_FILE, stringsAsFactors = FALSE)
common <- intersect(tree$tip.label, d$tip_label)
tree <- keep.tip(tree, common)
d <- d[match(tree$tip.label, d$tip_label), ]
rownames(d) <- d$tip_label
cat(sprintf("Tips analysed: %d\n", length(common)))

## ---- build the trait set ----------------------------------------------------
## (a) secretion + each auxiliary-function category (already binary)
cat_cols <- grep("^cat_", names(d), value = TRUE)
bin_traits <- c("secreted", cat_cols)

## (b) taxa: one-vs-rest binary for each sufficiently sampled phylum
top_phyla <- names(which(table(d$phylum) >= MIN_TIPS))
top_phyla <- setdiff(top_phyla, c("", NA))
for (p in top_phyla) {
  d[[paste0("phy_", make.names(p))]] <- as.integer(d$phylum == p)
}
phy_traits <- grep("^phy_", names(d), value = TRUE)

all_traits <- c(bin_traits, phy_traits)

## ---- D-statistic + phyloglm per trait ---------------------------------------
run_trait <- function(trait) {
  y <- d[[trait]]
  n1 <- sum(y == 1, na.rm = TRUE); n0 <- sum(y == 0, na.rm = TRUE)
  if (min(n1, n0) < MIN_TIPS) return(NULL)

  # caper::phylo.d uses non-standard evaluation for binvar, so give the trait a
  # fixed column name ("y") and pass it as a bare symbol.
  dsub <- d[, c("tip_label", trait)]
  names(dsub)[2] <- "y"
  cp <- comparative.data(tree, dsub, names.col = "tip_label")
  dd <- phylo.d(cp, binvar = y, permut = NPERM)

  # secondary: residual association with clade after ancestry
  pg_p <- NA_real_
  fit <- tryCatch(
    phyloglm(as.formula(paste(trait, "~ clade")), data = d, phy = tree,
             method = "logistic_MPLE", btol = 30),
    error = function(e) NULL)
  if (!is.null(fit)) {
    co <- summary(fit)$coefficients
    clade_rows <- grep("^clade", rownames(co))
    if (length(clade_rows)) pg_p <- min(co[clade_rows, "p.value"], na.rm = TRUE)
  }

  data.frame(trait = trait, n_pos = n1, n_neg = n0,
             D = unname(dd$DEstimate),
             p_vs_random   = unname(dd$Pval1),
             p_vs_brownian = unname(dd$Pval0),
             phyloglm_min_clade_p = pg_p)
}

res <- do.call(rbind, lapply(all_traits, run_trait))
res <- res %>% arrange(D)
cat("\n--- Trait partitioning: phylogenetic-signal + clade association ---\n")
print(res, row.names = FALSE, digits = 3)
write.csv(res, OUT, row.names = FALSE)
cat("\nWritten to:", OUT, "\n")

## =============================================================================
## HOW TO READ / REPORT
##   D << 1 (near 0 or negative), p_vs_random ~ 0, p_vs_brownian large
##       -> trait is strongly phylogenetically clustered (Brownian-like).
##          Report the class/clade pattern DESCRIPTIVELY (proportions in Fig 4)
##          and state it is phylogenetically structured (cite D). Do NOT attach
##          a per-sequence Fisher p-value; it overstates evidence.
##   D ~ 1, p_vs_random large
##       -> trait is phylogenetically random; a class-wise difference, if any,
##          would be closer to independent, but this is unlikely here.
##   phyloglm_min_clade_p:  small  -> some clade association survives ancestry
##                          large  -> the partition is collinear with the tree
##          (expected for traits that define clades; mirrors your is_12 p=0.9994).
##
## Bottom line for the manuscript: present the Fig 4 taxa + function heatmaps as
## proportions, cite D per trait as the phylogenetic-signal statistic, and drop
## the chi2 / pairwise-Fisher / one-vs-rest-Fisher+FDR enrichment tables from any
## inferential claim (retain only as descriptive effect sizes if at all).
## =============================================================================