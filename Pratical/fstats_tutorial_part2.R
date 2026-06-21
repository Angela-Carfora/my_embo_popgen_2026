################################################################################
# F-STATISTICS TUTORIAL PART 2
# EMBO Population Genomics
#
# This tutorial covers:
#   Part 0: Setup and data provenance
#   Part 1: Building admixture graphs with find_graphs()
#   Part 2: Robustness of the admixture graph results
#   Part 3: qpWave and qpAdm


################################################################################


################################################################################
# PART 0: SETUP AND DATA PROVENANCE
################################################################################

# load the packages
library(admixtools)
library(ggplot2)

# Set a working directory within your home directory
setwd("~/test/")

# File paths used throughout this tutorial
info_file  <- "info_embo_1240k.txt"                 # metadata (one row per individual)
prefix     <- "v62.0_1240k_public_embosubset"       # EIGENSTRAT prefix (.geno/.snp/.ind)
f2_dir     <- "f2data_1240k"                        # will store f2 blocks


################################################################################
# PART 1: Building admixture graphs with find_graphs()
################################################################################

# find_graphs() automates the search for admixture graphs by exploring many
# random topologies and optimising branch lengths and admixture weights for each.
# It returns a ranked list of graphs by score.

# Here we will investigate the relationships between ancient samples from Africa
# that predate migrations associated with food production in the continent.

# Select the testing groups
pops = c("Chimp", "South_Africa_2000BP", "Ethiopia_4500BP", 
         "Cameroon_3100-7800BP", "Malawi_2500-16000BP", "Morocco_14000BP")

# Find more context about the selected groups in the info file
info <- read.table(info_file, sep="\t", header=TRUE, stringsAsFactors=FALSE,
                   quote="", comment.char="")
info

# Compute f2 for the subset needed
f2_blocks_strict <- f2_from_geno(prefix, pops = pops)

# Let's find graphs for 0 to 4 admixture events
results <- lapply(0:4, function(n) {

  # Run find_graphs
  g      <- find_graphs(f2_blocks_strict, numadmix = n, outpop = "Chimp")
  
  # Get winner (lowest score) for each numadmix
  winner <- g[which.min(g$score), ]
  
  # Run the original qpGraph function to obtain the worst residual
  qpg    <- qpgraph(f2_blocks_strict, edges_to_igraph(winner$edges[[1]]), return_fstats=TRUE)
  
  # Plot winning graph
  p <- plot_graph(winner$edges[[1]], textsize = 4, fix = TRUE) +
    ggtitle(paste("numadmix =", n, "| score =", round(winner$score, 4),
                  "| worst residual =", round(qpg$worst_residual, 2)))
  print(p)
  list(winner = winner, qpg = qpg)
})

# Questions: What is the minimum number of admixture events required to obtain a 
# graph with the worst residuals (Z scores) below 3 in your results?
# Does the score improve substantially for higher number of admixture events?


################################################################################
# PART 2: ROBUSTNESS OF THE ADMIXTURE GRAPH RESULTS
################################################################################
#
# The graph score varies across SNP subsets, so small score differences
# between models may not be meaningful. ADMIXTOOLS 2 provides tools to test
# whether one graph fits significantly better than another by comparing
# score distributions across bootstrap-resampled SNP blocks.
#
#   1. Run find_graphs several times with numadmix=3 — each run may find a
#      different topology (local optimum)
#   2. Evaluate each candidate on bootstrap-resampled SNP blocks with
#      qpgraph_resample_snps()
#   3. Compare score distributions pairwise with compare_fits()


# --- Step 1: Find multiple candidate graphs (numadmix = 3) ---

n_runs <- 5

candidates <- lapply(1:n_runs, function(i) {
  
  # Run find_graphs and extract winner
  g      <- find_graphs(f2_blocks_strict, numadmix = 3, outpop = "Chimp")
  winner <- g[which.min(g$score), ]
  
  # Print score for this run
  cat("Run", i, "| score =", round(winner$score, 4), "\n")
  
  # Plot the winning graph
  p <- plot_graph(winner$edges[[1]], textsize = 4, fix = TRUE) +
    ggtitle(paste("Candidate", i, "| score =", round(winner$score, 4)))
  print(p)
  
  # Return edges, igraph and score for later use
  list(edges  = winner$edges[[1]],
       igraph = edges_to_igraph(winner$edges[[1]]),
       score  = winner$score)
})

# Questions: Are the scores of the winner graphs similar? Are the graphs 
# converging to similar solutions?
# Note that all graphs are build with the same SNP set.


# --- Step 2: Evaluate candidates on bootstrap-resampled SNP blocks ---
#
# qpgraph_resample_snps() fits each graph on many bootstrap-resampled sets of
# SNP blocks. The variability in scores across bootstrap samples reflects
# uncertainty due to the finite number of SNPs.

fits <- lapply(candidates, function(cand) {
  qpgraph_resample_snps(f2_blocks_strict, graph = cand$igraph, boot = 100)
})
#this shows the score variability of a fixed graph across bootstrap SNP blocks
fits[[1]]$score


# --- Step 3: Pairwise model comparison ---
#
# compare_fits() takes the bootstrap score distribution of two graphs on the same
# populations and tests whether the scores of one graph are significantly higher
# or lower than the scores of the other graph.
# The key output is p_emp: a two-sided bootstrap p-value.
# A non-significant result (p_emp > 0.05) means that we cannot distinguish
# the two models — they are equally good fits.
# ci_low: The 2.5% quantile of distribution of score differences
# ci_high: The 97.5% quantile of distribution of score differences

cat("\nPairwise model comparisons:\n")
for (i in 1:(n_runs - 1)) {
  for (j in (i + 1):n_runs) {
    cf <- compare_fits(fits[[i]]$score, fits[[j]]$score)
    cat("  Run", i, "vs Run", j,
        "| p_emp =", round(cf$p_emp, 3),
        "| 95% CI of score diff: [", round(cf$ci_low, 4), ",",
        round(cf$ci_high, 4), "]\n")
  }
}

# Question: Are any two graphs significantly different in fit?
# A non-significant p_emp for most pairs illustrates a key limitation of
# admixture graph modelling: many different topologies can explain the the data
# equally well.


# --- Step 4: Summary results across multiple graphs ---
#
# the number of admixture events can be summarized with summarize_numadmix()

nadmix_list <- lapply(seq_along(candidates), function(i) {
  df <- summarize_numadmix(candidates[[i]]$igraph)
  setNames(df$nadmix, df$pop)
})

nadmix_table <- do.call(cbind, lapply(nadmix_list, function(x) x[pops]))
rownames(nadmix_table) <- pops
colnames(nadmix_table) <- paste0("run", seq_along(candidates))
nadmix_table

# the order of admixture events can be summarized with summarize_eventorder()
lapply(candidates, function(cand) summarize_eventorder(cand$igraph))


################################################################################
# PART 3: qpWave and qpAdm
################################################################################

# qpWave tests the RANK of the matrix of f4 statistics between left and right
# populations. The rank indicates the minimum number of independent waves
# between left and right populations needed to explain the data.
#
#   rank = 0  →  all left pops share a single ancestral source
#   rank = 1  →  two independent streams are needed
#   rank = k  →  k+1 independent streams are needed
#
# Strategy:
#   1. One-wave test: pair target with each candidate source on the left.
#      If rank 0 is not rejected, that source alone explains the target.
#   2. Two-wave test: pair target with each combination of two sources.
#      If rank 0 is rejected but rank 1 is not, two streams are needed.
#
# We will use the same test populations as before, except for Chimp, and add
# Mbuti (a deeply-divergent African group not well represented by any of the 
# ancient groups in the previous test), French (for cases of Eurasian admixture),
# and the Denisova, which should not contribute gene flow to any of the modern 
# human groups included. The present day Ju'|hoan North (Khoisan group from 
# Southern Africa) will be used as target population.

all_pops <- c("Denisova", "Morocco_14000BP", "Ethiopia_4500BP", "Cameroon_3100-7800BP",
              "Malawi_2500-16000BP", "South_Africa_2000BP", "French", "Mbuti", "Ju_hoan_North")

f2_blocks_strict   <- f2_from_geno(prefix, pops = all_pops)
target             <- "Ju_hoan_North"
right_all          <- setdiff(all_pops, target)
non_denisova_right <- setdiff(right_all, "Denisova")   # candidate sources


# ── One-wave test ─────────────────────────────────────────────────────────────
# For each candidate source, we place it on the left together with the target.
# left  = c(target, source)       → 2 populations being modelled
# right = Denisova + remaining    → reference populations
#
# Remember: a high p-value for rank 0 means target and source are consistent with
# sharing a single ancestral stream → one wave suffices.

one_wave <- lapply(non_denisova_right, function(src) {
  left   <- c(target, src)
  right  <- c("Denisova", setdiff(non_denisova_right, src))
  result <- qpwave(f2_blocks_strict, left = left, right = right)
  list(source = src, left = left, right = right, result = result)
})

# Summary table — all ranks and p-values, sorted by p_rank0
one_wave_table <- do.call(rbind, lapply(one_wave, function(x) {
  rd  <- x$result$rankdrop
  row <- data.frame(source = x$source)
  for (r in rd$f4rank) row[[paste0("p_rank", r)]] <- rd$p[rd$f4rank == r]
  row
}))
one_wave_table[order(-one_wave_table$p_rank0), ]

# Question: Can we describe the Ju'|hoan North as deriving from one of the 
# populations tested alone, given the chosen right populations?


# ── Two-wave test ─────────────────────────────────────────────────────────────
# For each pair of candidate sources, we place both on the left with the target.
# left  = c(target, src1, src2)   → 3 populations being modelled
# right = Denisova + remaining    → reference populations
#
# Interpretation:
#   p_rank0 < 0.05 AND p_rank1 > 0.05  →  two streams needed
#   p_rank1 < 0.05                      →  even two streams are insufficient
#
# The table is sorted by p_rank1 so the best-fitting two-source models
# appear at the top.

two_wave <- lapply(combn(non_denisova_right, 2, simplify = FALSE), function(pair) {
  left   <- c(target, pair)
  right  <- c("Denisova", setdiff(non_denisova_right, pair))
  result <- qpwave(f2_blocks_strict, left = left, right = right)
  list(source1 = pair[1], source2 = pair[2], left = left, right = right, result = result)
})

# Summary table — all ranks and p-values, sorted by p_rank1
two_wave_table <- do.call(rbind, lapply(two_wave, function(x) {
  rd  <- x$result$rankdrop
  row <- data.frame(source1 = x$source1, source2 = x$source2)
  for (r in rd$f4rank) row[[paste0("p_rank", r)]] <- rd$p[rd$f4rank == r]
  row
}))
two_wave_table[order(-two_wave_table$p_rank1), ]

#show only models that are not rejected
good_pairs <- two_wave_table[two_wave_table$p_rank1>0.05, ]

# Full results accessible via e.g.:
# one_wave[[1]]$result      — full qpwave output for first source
# two_wave[[1]]$result      — full qpwave output for first source pair

# Question: Can we describe the Ju'|hoan North as deriving from two of the 
# populations tested alone, given the chosen right populations? Is there a
# single good-fiting model?


# ── qpAdm: estimate admixture proportions for non-rejected two-source models ─
# qpAdm uses the same left/right structure as qpWave but additionally estimates
# the proportion of ancestry from each source population.

qpadm_results <- lapply(seq_len(nrow(good_pairs)), function(i) {
  src1  <- good_pairs$source1[i]
  src2  <- good_pairs$source2[i]
  left  <- c(src1, src2)
  right <- c("Denisova", setdiff(non_denisova_right, c(src1, src2)))
  
  res <- qpadm(f2_blocks_strict, target = target, left = left, right = right)
  
  cat("Sources:", src1, "+", src2, "\n")
  print(res$weights)   # admixture proportions + standard errors
  cat("\n")
  
  list(source1 = src1, source2 = src2, result = res)
})

# Question: In previous studies, the Ju'|hoan North was shown to have some East
# African-related ancestry associated with pastoralist groups that arrived in 
# southern Africa < 2000 ya. How do you interpret the proportions estimated 
# for these particular sources based on all you learned about them today?


# OPTIONAL CHALLENGE: Select a new set of populations from our full dataset, 
# find suitable graphs, and use what you learned so far to evaluate the 
# robustness of the graphs. Alternatively, select a modern population and 
# try to describe it as a 1-wave or multiple-wave mixture of other populations.


################################################################################
# END OF TUTORIAL
#
# For further reading:
#   - Patterson et al. (2012) Genetics — original f-statistics paper
#   - Maier et al. (2023) eLife - On the limits of fitting complex models of 
#     population history to f-statistics
#   - Flegontov et al. (2023) PLOS Genetics - Modeling of African population 
#     history using f-statistics is biased when applying all previously proposed 
#     SNP ascertainment schemes
#   - ADMIXTOOLS v2 documentation
################################################################################

