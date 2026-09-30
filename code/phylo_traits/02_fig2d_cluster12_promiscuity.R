## ============================================================================
## Fig. 2d — COMPLETE-CASE PHYLOGENETIC LOGISTIC REGRESSION
##
## Tests whether carriage of >=1 promiscuous auxiliary domain is associated
## with Cluster 1.2 after accounting for shared ancestry.
##
## Missing trait values remain NA and are excluded. They are never converted
## to biological absences. The tree is pruned to the exact complete-case set.
## ============================================================================
suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
})

## ---- EDIT ------------------------------------------------------------------
TREE_FILE  <- "tree.nwk"
TRAIT_FILE <- "tree_tip_traits_CANONICAL_654.csv"

TIP_COL    <- "tip"
TRAIT_COL  <- "has_prom"
CL12_COL   <- "is_12"

OUT_DIR    <- "fig2d_phylo_results"
NBOOT      <- 1000  # use 50 for a test run
RUN_IG10   <- FALSE # optional; MPLE is the primary analysis
## ---------------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE)

parse_binary <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(z))

  out[z %in% c("1", "true", "t", "yes", "y")]  <- 1L
  out[z %in% c("0", "false", "f", "no", "n")] <- 0L
  out
}

tr_raw <- read.tree(TREE_FILE)
dat_raw <- read.csv(
  TRAIT_FILE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "None", "")
)

required <- c(TIP_COL, TRAIT_COL, CL12_COL)
missing_columns <- setdiff(required, names(dat_raw))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

dat_raw$has_prom <- parse_binary(dat_raw[[TRAIT_COL]])
dat_raw$is_12    <- parse_binary(dat_raw[[CL12_COL]])

cat(sprintf("Raw tree: %d tips\n", Ntip(tr_raw)))
cat(sprintf("Trait table before collapsing: %d rows\n", nrow(dat_raw)))

# The input may contain more than one row per tree tip, for example one row
# per auxiliary-domain record. For the sequence-level response used in Fig. 2d:
#   * any known positive row makes the tip positive;
#   * otherwise any known negative row makes the tip negative;
#   * a tip remains NA only when every row is unknown.
# Cluster 1.2 status must be consistent across duplicate rows.
collapse_any_positive <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_integer_)
  if (any(x == 1L)) return(1L)
  0L
}

collapse_consistent_binary <- function(x, tip) {
  x <- unique(x[!is.na(x)])
  if (length(x) == 0L) return(NA_integer_)
  if (length(x) > 1L) {
    stop(
      "Conflicting Cluster 1.2 assignments among duplicate rows for tip: ",
      tip
    )
  }
  as.integer(x)
}

tip_groups <- split(seq_len(nrow(dat_raw)), dat_raw[[TIP_COL]])

dat_collapsed <- do.call(
  rbind,
  lapply(names(tip_groups), function(tip) {
    idx <- tip_groups[[tip]]
    data.frame(
      tip = tip,
      has_prom = collapse_any_positive(dat_raw$has_prom[idx]),
      is_12 = collapse_consistent_binary(dat_raw$is_12[idx], tip),
      n_input_rows = length(idx),
      n_known_positive_rows = sum(dat_raw$has_prom[idx] == 1L, na.rm = TRUE),
      n_known_negative_rows = sum(dat_raw$has_prom[idx] == 0L, na.rm = TRUE),
      n_unknown_rows = sum(is.na(dat_raw$has_prom[idx])),
      stringsAsFactors = FALSE
    )
  })
)
rownames(dat_collapsed) <- NULL
names(dat_collapsed)[names(dat_collapsed) == "tip"] <- TIP_COL

n_duplicate_tips <- sum(dat_collapsed$n_input_rows > 1L)
cat(sprintf(
  "Collapsed to %d unique tips; %d tips had duplicate input rows.\n",
  nrow(dat_collapsed), n_duplicate_tips
))

write.csv(
  dat_collapsed,
  file.path(OUT_DIR, "collapsed_trait_table.csv"),
  row.names = FALSE
)
write.csv(
  dat_collapsed[dat_collapsed$n_input_rows > 1L, , drop = FALSE],
  file.path(OUT_DIR, "duplicate_tip_aggregation_summary.csv"),
  row.names = FALSE
)

dat_raw <- dat_collapsed

excluded_tips <- setdiff(tr_raw$tip.label, dat_raw[[TIP_COL]])
cat(sprintf(
  "Tree tips absent from collapsed trait table: %d\n",
  length(excluded_tips)
))

excluded_df <- data.frame(
  tip = excluded_tips,
  no_pipe_label = !grepl("\\|", excluded_tips),
  stringsAsFactors = FALSE
)
write.csv(
  excluded_df,
  file.path(OUT_DIR, "tips_absent_from_trait_table.csv"),
  row.names = FALSE
)
cat(sprintf(
  "Of the absent tips, %d have no-pipe labels (the expected NCBI-style group).\n",
  sum(excluded_df$no_pipe_label)
))

# Retain only trait rows whose labels occur in the tree.
dat_tree <- dat_raw[dat_raw[[TIP_COL]] %in% tr_raw$tip.label, , drop = FALSE]

cat("\nTrait missingness by Cluster 1.2 status before pruning:\n")
missingness_tab <- table(
  is_cluster_12   = dat_tree$is_12,
  trait_available = !is.na(dat_tree$has_prom),
  useNA = "always"
)
print(missingness_tab)
write.csv(
  as.data.frame(missingness_tab),
  file.path(OUT_DIR, "trait_missingness_by_cluster12.csv"),
  row.names = FALSE
)

# The analysis requires both the response and predictor.
dat <- dat_tree[
  !is.na(dat_tree$has_prom) &
  !is.na(dat_tree$is_12),
  ,
  drop = FALSE
]

if (nrow(dat) == 0L) stop("No complete cases remain.")

tr <- keep.tip(tr_raw, dat[[TIP_COL]])
dat <- dat[match(tr$tip.label, dat[[TIP_COL]]), , drop = FALSE]
rownames(dat) <- dat[[TIP_COL]]

stopifnot(identical(tr$tip.label, dat[[TIP_COL]]))
stopifnot(!anyNA(dat$has_prom))
stopifnot(!anyNA(dat$is_12))

if (length(unique(dat$has_prom)) != 2L) {
  stop("The response does not contain both 0 and 1 among complete cases.")
}
if (length(unique(dat$is_12)) != 2L) {
  stop("The predictor does not contain both 0 and 1 among complete cases.")
}

cat(sprintf(
  "\nComplete-case tree: %d tips | has_prom=1: %d (%.2f%%) | excluded: %d\n",
  Ntip(tr), sum(dat$has_prom == 1L), 100 * mean(dat$has_prom == 1L),
  Ntip(tr_raw) - Ntip(tr)
))

write.csv(
  dat[, c(TIP_COL, "has_prom", "is_12")],
  file.path(OUT_DIR, "complete_case_tip_data.csv"),
  row.names = FALSE
)
write.tree(tr, file = file.path(OUT_DIR, "complete_case_tree.nwk"))

## ---- 2 x 2 table and non-phylogenetic alignment diagnostic -----------------
tab <- with(
  dat,
  table(
    is_cluster_12   = factor(is_12, levels = c(0, 1)),
    has_promiscuous = factor(has_prom, levels = c(0, 1))
  )
)
cat("\nComplete-case 2 x 2 table:\n")
print(tab)
write.csv(
  as.data.frame.matrix(tab),
  file.path(OUT_DIR, "cluster12_contingency_table.csv")
)

fit_glm <- glm(has_prom ~ is_12, data = dat, family = binomial())
glm_sum <- summary(fit_glm)
glm_beta <- unname(coef(fit_glm)["is_12"])

n00 <- unname(tab["0", "0"])
n01 <- unname(tab["0", "1"])
n10 <- unname(tab["1", "0"])
n11 <- unname(tab["1", "1"])

if (all(c(n00, n01, n10, n11) > 0)) {
  table_or <- (n11 * n00) / (n10 * n01)
  if (!isTRUE(all.equal(glm_beta, log(table_or), tolerance = 1e-8))) {
    stop("GLM coefficient and table odds ratio disagree; check tip alignment.")
  }
} else {
  table_or <- ((n11 + 0.5) * (n00 + 0.5)) /
              ((n10 + 0.5) * (n01 + 0.5))
  warning("A table cell is zero; descriptive OR uses a 0.5 correction.")
}

fisher_result <- fisher.test(tab, alternative = "two.sided")
fisher_ci <- unname(fisher_result$conf.int)

crude <- data.frame(
  n_tips            = nrow(dat),
  n_positive        = sum(dat$has_prom == 1L),
  n_cluster_12      = sum(dat$is_12 == 1L),
  n_non_cluster_12  = sum(dat$is_12 == 0L),
  table_odds_ratio  = table_or,
  fisher_odds_ratio = unname(fisher_result$estimate),
  fisher_CI_low     = fisher_ci[1],
  fisher_CI_high    = fisher_ci[2],
  fisher_p_value    = fisher_result$p.value,
  glm_beta          = glm_beta,
  glm_SE            = glm_sum$coefficients["is_12", "Std. Error"],
  glm_p_value       = glm_sum$coefficients["is_12", "Pr(>|z|)"]
)
cat("\nCrude compositional effect/alignment diagnostic:\n")
print(crude)
write.csv(crude, file.path(OUT_DIR, "crude_effect.csv"), row.names = FALSE)

## ---- Primary phylogenetic logistic regression -------------------------------
cat(sprintf(
  "\nRunning MPLE phyloglm with %d bootstrap replicates...\n",
  NBOOT
))

fit_mple <- phyloglm(
  has_prom ~ is_12,
  data   = dat,
  phy    = tr,
  method = "logistic_MPLE",
  btol   = 30,
  boot   = NBOOT
)

mple_sum <- summary(fit_mple)
print(mple_sum)
saveRDS(fit_mple, file.path(OUT_DIR, "phyloglm_MPLE_fit.rds"))

mple_coef <- as.data.frame(mple_sum$coefficients)
mple_coef$term <- rownames(mple_coef)
rownames(mple_coef) <- NULL
write.csv(
  mple_coef,
  file.path(OUT_DIR, "phyloglm_MPLE_coefficients.csv"),
  row.names = FALSE
)

if (!is.null(fit_mple$bootresult)) {
  write.csv(
    as.data.frame(fit_mple$bootresult),
    file.path(OUT_DIR, "phyloglm_MPLE_bootstrap_raw.csv"),
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
    term_index <- match(term, names(coef(fit)))
    if (is.na(term_index) || ncol(br) < term_index) return(NULL)
    v <- br[[term_index]]
  }

  v <- as.numeric(v)
  v[is.finite(v)]
}

beta <- mple_sum$coefficients["is_12", "Estimate"]
se   <- mple_sum$coefficients["is_12", "StdErr"]
pval <- mple_sum$coefficients["is_12", "p.value"]

# Recent phylolm versions place coefficient bootstrap intervals directly in
# summary(fit)$coefficients. Prefer those named columns when available.
coef_names <- colnames(mple_sum$coefficients)
if (all(c("lowerbootCI", "upperbootCI") %in% coef_names)) {
  beta_ci <- c(
    mple_sum$coefficients["is_12", "lowerbootCI"],
    mple_sum$coefficients["is_12", "upperbootCI"]
  )
  ci_method <- "parametric bootstrap from summary"
} else {
  boot_beta <- extract_boot_term(fit_mple, "is_12")
  if (!is.null(boot_beta) && length(boot_beta) >= 20L) {
    beta_ci <- unname(quantile(boot_beta, c(0.025, 0.975), na.rm = TRUE))
    ci_method <- "parametric bootstrap from raw output"
  } else {
    beta_ci <- beta + c(-1, 1) * 1.96 * se
    ci_method <- "Wald fallback; inspect bootstrap_raw.csv"
    warning(
      "Could not identify a bootstrap interval for is_12. ",
      "A Wald interval was written; inspect phyloglm_MPLE_bootstrap_raw.csv."
    )
  }
}

primary_result <- data.frame(
  method           = "logistic_MPLE",
  n_tips           = nrow(dat),
  n_positive       = sum(dat$has_prom == 1L),
  beta             = beta,
  SE               = se,
  odds_ratio       = exp(beta),
  OR_CI_low        = exp(beta_ci[1]),
  OR_CI_high       = exp(beta_ci[2]),
  CI_method        = ci_method,
  p_value          = pval,
  alpha            = if (!is.null(fit_mple$alpha)) fit_mple$alpha else NA_real_,
  n_boot_requested = NBOOT
)
cat("\nPrimary phylogenetic result:\n")
print(primary_result)
write.csv(
  primary_result,
  file.path(OUT_DIR, "phyloglm_MPLE_summary.csv"),
  row.names = FALSE
)

## ---- Optional IG10 estimator sensitivity -----------------------------------
if (isTRUE(RUN_IG10)) {
  cat("\nRunning IG10 estimator sensitivity check...\n")
  fit_ig10 <- tryCatch(
    suppressWarnings(
      phyloglm(
        has_prom ~ is_12,
        data   = dat,
        phy    = tr,
        method = "logistic_IG10",
        btol   = 30
      )
    ),
    error = function(e) {
      cat("IG10 failed:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(fit_ig10)) {
    ig10_sum <- summary(fit_ig10)
    print(ig10_sum)

    # Treat non-convergence as a failed sensitivity analysis, not as a result.
    converged_flag <- if (!is.null(fit_ig10$convergence)) {
      identical(fit_ig10$convergence, 0L)
    } else {
      NA
    }

    saveRDS(fit_ig10, file.path(OUT_DIR, "phyloglm_IG10_fit.rds"))

    ig10_coef <- as.data.frame(ig10_sum$coefficients)
    ig10_coef$term <- rownames(ig10_coef)
    ig10_coef$converged <- converged_flag
    rownames(ig10_coef) <- NULL
    write.csv(
      ig10_coef,
      file.path(OUT_DIR, "phyloglm_IG10_coefficients.csv"),
      row.names = FALSE
    )
  }
} else {
  cat("\nIG10 sensitivity check skipped (RUN_IG10 = FALSE).\n")
}

writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo.txt"))

cat("\nAnalysis complete. Results written to:", OUT_DIR, "\n")
