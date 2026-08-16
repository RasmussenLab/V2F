###############################################################################
##########  How much of the CARMA result survives the sampler  ################
###############################################################################
#
# CARMA explores models by shotgun stochastic search and sets no seed of its own, so a single
# run reports one draw from the posterior rather than the posterior. This script runs the same
# analysis across N_SEEDS seeds and summarises what is stable across them.
#
# CARMA_T2D_signed.R is the analysis and owns the input guards. This script re-reads the same
# four files and assumes they have already passed those guards. Run it after that script.
#
# Reads   dat/                    the same inputs
# Writes  results/carma_seed_stability_<trait>.csv    per-variant summary across seeds
#         results/carma_seed_credible_sets.csv        every credible set from every seed
#
# Nothing is written to interim/: output.labels is NULL, which turns off CARMA's per-run dumps.
#
# Run:  Rscript CARMA_seed_stability.R          (about 4 minutes per seed)
###############################################################################

suppressMessages({
  library(CARMA); library(data.table); library(magrittr); library(dplyr)
})

N_SEEDS <- as.integer(Sys.getenv("N_SEEDS", "20"))

if (dir.exists("dat")) {
  ROOT <- "."; HERE <- "2b_finemapping"
} else if (dir.exists("../dat")) {
  ROOT <- ".."; HERE <- "."
} else {
  stop("cannot find the dat folder. Run from 2_colocalization_and_finemapping/ ",
       "or from 2b_finemapping/.")
}
DAT     <- file.path(ROOT, "dat")
RESULTS <- file.path(HERE, "results")
dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)

###############################################################################
## 1. Inputs
###############################################################################

sumstats1 <- read.table(file.path(DAT, "T2D_bse_signed.txt"), header = TRUE)
sumstats2 <- read.table(file.path(DAT, "fasting_insulin_bse_signed.txt"), header = TRUE)
ld_matrix <- as.matrix(read.table(file.path(DAT, "rs3843467_signed.ld"), header = FALSE))

sumstats1$Z <- sumstats1$BETA / sumstats1$SE
sumstats2$Z <- sumstats2$BETA / sumstats2$SE

z.list      <- list(sumstats1$Z, sumstats2$Z)
ld.list     <- list(ld_matrix, ld_matrix)
lambda.list <- list(1, 1)

traits <- c("T2D", "fasting_insulin")
snps   <- sumstats1$SNP

cat(sprintf("%d variants, %d seeds\n\n", length(snps), N_SEEDS))

###############################################################################
## 2. Run every seed
###############################################################################

# pip[[trait]] is seeds x variants; incs counts credible-set membership.
pip  <- lapply(traits, function(x) matrix(NA_real_, N_SEEDS, length(snps),
                                          dimnames = list(NULL, snps)))
incs <- lapply(traits, function(x) setNames(integer(length(snps)), snps))
names(pip) <- names(incs) <- traits
cs_rows <- list()

started <- Sys.time()
for (s in seq_len(N_SEEDS)) {
  set.seed(s)
  # CARMA prints per-locus progress; keep the log readable.
  invisible(capture.output(
    res <- CARMA(z.list, ld.list, lambda.list = lambda.list,
                 outlier.switch = TRUE, output.labels = NULL)
  ))

  for (ti in seq_along(traits)) {
    tr <- traits[ti]
    pip[[tr]][s, ] <- res[[ti]]$PIPs
    cs <- res[[ti]]$`Credible set`[[2]]
    if (length(cs) != 0) {
      for (l in seq_along(cs)) {
        members <- cs[[l]]
        incs[[tr]][members] <- incs[[tr]][members] + 1L
        cs_rows[[length(cs_rows) + 1L]] <- data.frame(
          seed = s, trait = tr, set = l, size = length(members),
          total_PIP = sum(res[[ti]]$PIPs[members]),
          members = paste(sprintf("%s(%.2f)", snps[members],
                                  res[[ti]]$PIPs[members]), collapse = "; "))
      }
    }
  }
  el <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  cat(sprintf("seed %2d/%d done  (%.1f min elapsed, ~%.0f min left)\n",
              s, N_SEEDS, el, el / s * (N_SEEDS - s)))
}

###############################################################################
## 3. Per-variant summary
###############################################################################

summarise_trait <- function(tr) {
  m <- pip[[tr]]
  data.frame(
    SNP        = snps,
    mean_PIP   = colMeans(m),
    sd_PIP     = apply(m, 2, sd),
    min_PIP    = apply(m, 2, min),
    max_PIP    = apply(m, 2, max),
    # how often the variant would have been called, one seed at a time
    n_PIP_gt50 = colSums(m > 0.5),
    n_in_CS    = as.integer(incs[[tr]]),
    n_seeds    = N_SEEDS,
    row.names  = NULL
  ) %>% arrange(desc(mean_PIP))
}

for (tr in traits) {
  out <- summarise_trait(tr)
  write.csv(out, file.path(RESULTS, sprintf("carma_seed_stability_%s.csv", tr)),
            row.names = FALSE)

  cat(sprintf("\n--- %s, across %d seeds ---\n", tr, N_SEEDS))
  cat("variant       meanPIP   sdPIP   minPIP  maxPIP   PIP>0.5   in a CS\n")
  for (i in seq_len(min(10, nrow(out)))) {
    r <- out[i, ]
    cat(sprintf("%-12s  %7.3f %7.3f  %7.3f %7.3f   %2d/%-2d      %2d/%-2d\n",
                r$SNP, r$mean_PIP, r$sd_PIP, r$min_PIP, r$max_PIP,
                r$n_PIP_gt50, N_SEEDS, r$n_in_CS, N_SEEDS))
  }
}

cs_all <- do.call(rbind, cs_rows)
write.csv(cs_all, file.path(RESULTS, "carma_seed_credible_sets.csv"), row.names = FALSE)

###############################################################################
## 4. Is there a shared signal?
##
## The single-run output applies PIP > 0.5 to one draw. Counting how often a variant clears
## that bar in both traits separates a stable shared signal from a threshold artefact.
###############################################################################

shared <- data.frame(
  SNP      = snps,
  n_both   = colSums(pip[[1]] > 0.5 & pip[[2]] > 0.5),
  n_T2D    = colSums(pip[[1]] > 0.5),
  n_FI     = colSums(pip[[2]] > 0.5),
  meanPIP_T2D = colMeans(pip[[1]]),
  meanPIP_FI  = colMeans(pip[[2]]),
  row.names = NULL
) %>% arrange(desc(n_both), desc(meanPIP_T2D))

write.csv(shared, file.path(RESULTS, "carma_seed_shared_high_pip.csv"), row.names = FALSE)

cat(sprintf("\n--- PIP > 0.5 in BOTH traits, out of %d seeds ---\n", N_SEEDS))
top <- shared %>% filter(n_T2D > 0 | n_FI > 0) %>% head(10)
if (nrow(top) == 0) {
  cat("  no variant cleared PIP > 0.5 in either trait under any seed\n")
} else {
  cat("variant       both   T2D    FI    meanPIP_T2D  meanPIP_FI\n")
  for (i in seq_len(nrow(top))) {
    r <- top[i, ]
    cat(sprintf("%-12s  %2d/%-2d  %2d/%-2d %2d/%-2d      %6.3f      %6.3f\n",
                r$SNP, r$n_both, N_SEEDS, r$n_T2D, N_SEEDS, r$n_FI, N_SEEDS,
                r$meanPIP_T2D, r$meanPIP_FI))
  }
}

cat("\nwritten to", normalizePath(RESULTS), "\n")

###############################################################################
## Reading the output
##
## mean_PIP over seeds is a better summary than any single run's PIP, and sd_PIP says how much
## the sampler moved. A variant with high mean_PIP and low sd_PIP is a stable call. A group of
## variants that share the posterior between them, each with large sd_PIP, is one signal the
## data cannot resolve, not several.
##
## n_in_CS is the more robust statistic: credible-set membership is stable even where the
## ranking inside the set is not.
###############################################################################
