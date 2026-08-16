###############################################################################
##########  CARMA on T2D and fasting insulin, signed and harmonised LD  #######
###############################################################################
#
# Fine maps the C5orf67 locus for two traits, T2D and fasting insulin, over a shared signed
# LD matrix, and reports the credible sets and any variant with high posterior inclusion
# probability in both.
#
# Inputs
# ------
# All four come from 2a_colocalization/SharePro_signed.ipynb and share one variant order:
#
#   T2D_bse_signed.txt              SNP, BETA, SE, N
#   fasting_insulin_bse_signed.txt  SNP, BETA, SE, N
#   rs3843467_signed.ld             247 x 247 signed correlation matrix, no labels
#   variant_order.tsv               the matrix row index, and the allele every sign refers to
#
# CARMA models association statistics as a function of the signed correlation between
# variants, so the matrix must hold r rather than R2, and every beta must already be flipped
# onto the matrix's reference allele. Section 2 refuses to run if either is untrue: an
# unsigned matrix and unharmonised betas both produce credible sets that look like many
# independent signals at one locus.
#
# What to expect
# --------------
# The matrix is rank-deficient, 149 independent dimensions across 247 variants, with 347
# pairs at |r| = 1. Variants inside a collinear group cannot be ranked on association
# evidence alone, so read the credible set as a whole and not its ordering. See the note at
# section 4 on why the seed is fixed.
#
# Run:  Rscript CARMA_T2D_signed.R
###############################################################################

# install.packages(c("data.table", "magrittr", "dplyr", "remotes", "R.utils"))
# remotes::install_github("ZikunY/CARMA", ref = "master")

library(CARMA)
library(data.table)
library(magrittr)
library(dplyr)

###############################################################################
## 0. Layout
##
##   2_colocalization_and_finemapping/
##     dat/                inputs, read only. Written by SharePro_signed.ipynb.
##     2b_finemapping/
##       interim/          CARMA's own per-locus dumps, safe to delete
##       results/          the CSVs this script produces
##
## Runs from either 2_colocalization_and_finemapping/ or 2b_finemapping/.
###############################################################################

if (dir.exists("dat")) {
  ROOT <- "."; HERE <- "2b_finemapping"
} else if (dir.exists("../dat")) {
  ROOT <- ".."; HERE <- "."
} else {
  stop("cannot find the dat folder. Run from 2_colocalization_and_finemapping/ ",
       "or from 2b_finemapping/.")
}

DAT     <- file.path(ROOT, "dat")
INTERIM <- file.path(HERE, "interim")
RESULTS <- file.path(HERE, "results")
dir.create(INTERIM, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)

cat("reading from", normalizePath(DAT), "\n")
cat("interim    ", normalizePath(INTERIM), "\n")
cat("results    ", normalizePath(RESULTS), "\n")

###############################################################################
## 1. Inputs
###############################################################################

sumstats1 <- read.table(file.path(DAT, "T2D_bse_signed.txt"), header = TRUE)
sumstats2 <- read.table(file.path(DAT, "fasting_insulin_bse_signed.txt"), header = TRUE)
ld_matrix <- as.matrix(read.table(file.path(DAT, "rs3843467_signed.ld"), header = FALSE))
variants  <- read.table(file.path(DAT, "variant_order.tsv"), header = TRUE, sep = "\t")

###############################################################################
## 2. Guards
##
## Bad inputs here produce plausible-looking output rather than an error, so each assumption
## the model depends on is checked and made noisy.
###############################################################################

stopifnot(
  nrow(ld_matrix) == ncol(ld_matrix),
  nrow(ld_matrix) == nrow(sumstats1),
  nrow(ld_matrix) == nrow(sumstats2),
  nrow(ld_matrix) == nrow(variants)
)

# The matrix, the variant table and both traits must be in one order.
stopifnot(identical(as.character(sumstats1$SNP), as.character(variants$SNP)))
stopifnot(identical(as.character(sumstats2$SNP), as.character(variants$SNP)))

# A correlation matrix has a unit diagonal. An R2 matrix does too, so this is necessary but
# not sufficient, which is what the next two checks are for.
stopifnot(isTRUE(all.equal(diag(ld_matrix), rep(1, nrow(ld_matrix)), tolerance = 1e-6)))

# An R2 matrix has no negative entries, which is impossible for a real correlation matrix over
# this many variants. This is what catches an unsigned matrix reaching ld.list.
off <- ld_matrix[upper.tri(ld_matrix)]
cat(sprintf("LD matrix: %d x %d, mean |off-diagonal| = %.4f, range [%.3f, %.3f], %d negative\n",
            nrow(ld_matrix), ncol(ld_matrix), mean(abs(off)), min(off), max(off), sum(off < 0)))
if (sum(off < 0) == 0) {
  stop("the matrix has no negative entries, so it is still unsigned. ",
       "Rerun 2a_colocalization/SharePro_signed.ipynb, which rewrites dat/rs3843467_signed.ld.")
}
# Squaring also shrinks every off-diagonal, so an R2 file looks unusually flat.
if (mean(abs(off)) < 0.05) {
  warning("the off-diagonal entries look unusually small. Check that this is r and not R2.")
}

# The betas must be harmonised. If nothing was flipped, the harmonisation step did not run.
if (!"T2D_beta_flipped" %in% names(variants)) {
  stop("variant_order.tsv has no T2D_beta_flipped column: rerun the preprocessing notebook")
}
cat(sprintf("harmonisation: %d of %d T2D betas were flipped onto the reference allele\n",
            sum(variants$T2D_beta_flipped %in% c(TRUE, "True", "TRUE")), nrow(variants)))

###############################################################################
## 3. Z scores
###############################################################################

sumstats1$Z <- sumstats1$BETA / sumstats1$SE
sumstats2$Z <- sumstats2$BETA / sumstats2$SE

# Under a single causal variant the Z of a correlated variant tracks r times the lead Z, so at
# high |r| the product Z_i * Z_j should carry the sign of r_ij. This compares the matrix
# against the summary statistics, so it is the one check that an unharmonised set of betas
# still fails even when the matrix itself is correctly signed. A small minority of
# contradicting pairs is expected; a large fraction means the alleles are not aligned.
lead <- which.max(abs(sumstats1$Z))
strong <- which(abs(ld_matrix[lead, ]) >= 0.8)
strong <- setdiff(strong, lead)
if (length(strong) > 0) {
  contradicting <- sum(sign(sumstats1$Z[strong] * sumstats1$Z[lead]) !=
                       sign(ld_matrix[lead, strong]))
  cat(sprintf("LD mismatch check against %s: %d of %d variants at |r| >= 0.8 contradict the LD sign\n",
              sumstats1$SNP[lead], contradicting, length(strong)))
  if (contradicting > 0.2 * length(strong)) {
    stop("more than a fifth of the strongly linked variants still contradict the LD sign. ",
         "The alleles are not harmonised. Do not interpret the output.")
  }
}

###############################################################################
## 4. CARMA
###############################################################################

# CARMA explores models by shotgun stochastic search, calling sample() throughout, and sets no
# seed of its own. Without a fixed seed, two runs on identical inputs return different credible
# sets. On this locus the difference is large: across runs the T2D lead moved between rs9686661
# and rs9687846, and rs256903 between PIP 0.27 and 0.61. The variants that swap sit in groups at
# |r| = 1, where the likelihood cannot separate them, so the sampler is choosing a
# representative rather than finding a signal.
#
# The seed makes a run reproducible. It does not make the ordering within a collinear group
# meaningful, and any claim about a single variant here needs evidence from outside the
# summary statistics.
#
# CARMA_seed_stability.R runs this analysis across 20 seeds and reports mean PIP, its spread,
# and how often each variant lands in a credible set. Those are the numbers to quote: a single
# run's PIP is one draw. Credible-set membership is far more stable than the PIP ordering
# within a set, so prefer it when summarising.
set.seed(1)

z.list      <- list(sumstats1$Z, sumstats2$Z)
ld.list     <- list(ld_matrix, ld_matrix)
lambda.list <- list(1, 1)

# outlier.switch is CARMA's own test for summary statistics that disagree with the supplied
# LD. Variants it flags are informative about the inputs rather than a failure.
#
# output.labels defaults to ".", which drops post_locus_* files wherever the script happens to
# be run from. Pointed at interim/ instead. label.list names them after the traits rather than
# locus_1 and locus_2, so the two runs cannot be confused.
CARMA.results <- CARMA(z.list, ld.list,
                       lambda.list    = lambda.list,
                       outlier.switch = TRUE,
                       output.labels  = INTERIM,
                       label.list     = list("T2D", "fasting_insulin"))

###############################################################################
## 5. Results
###############################################################################

annotate <- function(sumstats, res) {
  out <- sumstats %>% mutate(PIP = res$PIPs, CS = 0)
  cs <- res$`Credible set`[[2]]
  if (length(cs) != 0) {
    for (l in seq_along(cs)) out$CS[cs[[l]]] <- l
  }
  out
}

res1 <- annotate(sumstats1, CARMA.results[[1]])
res2 <- annotate(sumstats2, CARMA.results[[2]])

summarise_run <- function(res, label) {
  sets <- sort(unique(res$CS[res$CS > 0]))
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("credible sets: %d\n", length(sets)))
  cat(sprintf("variants at PIP > 0.95: %d\n", sum(res$PIP > 0.95)))
  for (l in sets) {
    members <- res[res$CS == l, ] %>% arrange(desc(PIP))
    cat(sprintf("  set %d: %d variant(s), total PIP %.3f -> %s\n",
                l, nrow(members), sum(members$PIP),
                paste(sprintf("%s(%.2f)", members$SNP, members$PIP),
                      collapse = ", ")))
  }
  # Many single-variant sets at one locus is the signature of an unsigned matrix or
  # unharmonised betas, not of many independent signals.
  singletons <- sum(sapply(sets, function(l) sum(res$CS == l) == 1))
  if (singletons >= 5) {
    cat(sprintf("  NOTE: %d single-variant credible sets at one locus. Check the LD mismatch\n",
                singletons))
    cat("        line above before interpreting these as separate signals.\n")
  }
  invisible(res)
}

summarise_run(res1, "T2D")
summarise_run(res2, "fasting insulin")

cat("\nTop 10 by PIP, T2D:\n")
print(res1 %>% arrange(desc(PIP)) %>%
        select(SNP, BETA, SE, Z, PIP, CS) %>% head(10))

###############################################################################
## 6. Save
###############################################################################

write.csv(res1, file.path(RESULTS, "carma_results_T2D_signed.csv"), row.names = FALSE)
write.csv(res2, file.path(RESULTS, "carma_results_fasting_insulin_signed.csv"), row.names = FALSE)

write.csv(res1 %>% filter(CS > 0) %>% arrange(CS, desc(PIP)),
          file.path(RESULTS, "credible_set_variants_T2D_signed.csv"), row.names = FALSE)
write.csv(res2 %>% filter(CS > 0) %>% arrange(CS, desc(PIP)),
          file.path(RESULTS, "credible_set_variants_fasting_insulin_signed.csv"), row.names = FALSE)

shared <- inner_join(
  res1 %>% filter(PIP > 0.5) %>% select(SNP, PIP_T2D = PIP),
  res2 %>% filter(PIP > 0.5) %>% select(SNP, PIP_FI  = PIP),
  by = "SNP")
write.csv(shared, file.path(RESULTS, "shared_high_pip_variants_signed.csv"), row.names = FALSE)

cat("\nwritten to", normalizePath(RESULTS), "\n")
for (f in sort(list.files(RESULTS, full.names = TRUE))) {
  cat(sprintf("   %-46s %7.1f KB\n", basename(f), file.size(f) / 1024))
}

###############################################################################
## Reading the output
##
## Each trait resolves to a small number of credible sets, one of which covers the 3' end
## haplotype. Read that set as a unit. rs459193, rs256903, rs173964 and rs256904 sit at
## pairwise R2 above 0.99, so no summary-statistic method can rank them on association
## evidence alone, and the member CARMA places first changes with the seed.
##
## shared_high_pip_variants_signed.csv lists variants above PIP 0.5 in both traits. It is
## empty whenever the two traits' posteriors settle on different members of the same
## collinear group, which happens here, so an empty file is not evidence against a shared
## signal. Compare the credible sets themselves, and the SharePro colocalization result in
## 2a_colocalization/results/, before concluding anything about sharing.
##
## Prioritising a single variant within the haplotype needs evidence the LD matrix cannot
## provide, such as the sequence models.
###############################################################################
