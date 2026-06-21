################################################################################
# F-STATISTICS TUTORIAL PART 1
# EMBO Population Genomics
#
# This tutorial covers:
#   Part 0: Setup and data provenance
#   Part 1: Exploring the dataset
#   Part 2: f2 statistic and FST — genetic divergence
#   Part 3: Outgroup f3 statistics — shared drift
#   Part 4: Admixture f3 statistics — testing admixture
#   Part 5: f4 statistics — detecting excess allele sharing
#   Part 6: f4-ratio — quantifying admixture proportions

################################################################################


################################################################################
# PART 0: SETUP AND DATA PROVENANCE
################################################################################

# --- 0.1  Setup ---
#
# For this practical, we will use ADMIXTOOLS v2 
# (https://uqrmaie1.github.io/admixtools/index.html), an R package that 
# re-implements and expands the methods originally developed in the ADMIXTOOLS 
# program (https://github.com/DReichLab/AdmixTools), which required a 
# combination of bash scripting and manual editing of input files.
#
# All necessary packages for this tutorial have been pre-installed:
# devtools::install_github("uqrmaie1/admixtools")
# install.packages("ggplot2") # for plotting
# install.packages("maps")    # for plotting on a world map
#
# load the packages
library(admixtools)
library(ggplot2)
library(maps)
#
# Set a working directory within your home directory
setwd("~/test/")


# --- 0.2  Data provenance ---
#
# The data comes from the Allen Ancient DNA Resource (AADR), a curated repository
# of ancient and modern human genotype data maintained by the Reich Lab
# (https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW).
# The AADR includes modern DNA genotyped on the Affymetrix Human Origins SNP
# Array (~600K SNPs) as well as target-captured modern and ancient data for an
# extended SNP panel (1240K SNPs). Additionally, some publicly available full
# genomes (shotgun data) are included for the same 1240K SNP subset.
#
# The AADR stores data in EIGENSTRAT format (.geno / .snp / .ind). This format
# is also used in ADMIXTOOLS v2. More information on format can be found
# here: https://reich.hms.harvard.edu/software/InputFileFormats
#
# The subset of samples selected for this tutorial consists of various African
# populations as well as some non-African representatives from the Simon Genome
# diversity Project (SGDP) and the 1000 Genomes Project. It also contains one
# chimpanzee genome (Pan troglodytes) to be used as a distant outgroup
# and four high-coverage archaic genomes:
#   - Altai Neanderthal
#   - Chagyrskaya Neanderthal
#   - Vindija Neanderthal
#   - Denisova

# --- 0.3  Data types: diploid vs pseudo-haploid ---
#
# The AADR contains two main types of genotype data:
# - diploid -> both alleles are called at each SNP
# - pseudo-haploid -> a single allele is drawn at random at each site
# At typical ancient DNA coverages, diploid calling would systematically inflate
# homozygosity — sites where only one allele is observed due to low coverage
# get called homozygous even when the individual is truly heterozygous. 
# Randomly drawing one allele per site instead yields unbiased allele frequency
# estimates at the population level.
# In our subset, modern samples are diploid, while ancient samples can be either
# diploid or pseudo-haploid.
#
# Taking into account the type of data is important for accurate f-statistic
# estimation. In the original ADMIXTOOLS, this was controlled by the inbreed flag.
# ADMIXTOOLS v2 detects ploidy automatically, making mixed-ploidy analyses more
# straightforward for the user.

# File paths used throughout this tutorial
info_file  <- "info_embo_1240k.txt"                 # metadata (one row per individual)
prefix     <- "v62.0_1240k_public_embosubset"       # EIGENSTRAT prefix (.geno/.snp/.ind)
f2_dir     <- "f2data_1240k"                        # will store f2 blocks


################################################################################
# PART 1: EXPLORING THE DATASET
################################################################################

# --- 1.1  Load the metadata ---

info <- read.table(info_file, sep="\t", header=TRUE, stringsAsFactors=FALSE,
                   quote="", comment.char="")
head(info)

# How many individuals, populations, and language groups are in the dataset?

cat("Individuals:", nrow(info), "\n")
cat("Populations:", length(unique(info$Group)), "\n")
cat("Language groups:", length(unique(info$Language[!is.na(info$Language)])), "\n")

# --- 1.2  Sample sizes per population ---

pop_counts <- sort(table(info$Group))
pop_counts

# --- 1.3  Which populations belong to which language family? ---
# Get one representative row per population (to avoid duplicates)

pop_meta <- info[!duplicated(info$Group),
                 c("Group","Continent","Language","Lat","Long","Date_mean_BP")]
pop_meta

# --- 1.4  Define population lists ---

# Archaic groups and outgroup
neanderthals <- c("Altai_Neanderthal", "Vindija_Neanderthal", "Chagyrskaya_Neanderthal")
denisovan    <- "Denisova"
outgroup     <- "Chimp"
archaic_outgroup <- c(neanderthals, denisovan, outgroup)

# Modern and ancient populations derived from the info file.
# Order follows the info file; Date_mean_BP == 0 flags modern samples.
modern         <- pop_meta$Group[pop_meta$Date_mean_BP == 0 &
                                   !pop_meta$Group %in% archaic_outgroup]
ancient <- pop_meta$Group[pop_meta$Date_mean_BP >  0 &
                                   !pop_meta$Group %in% archaic_outgroup]

# Colour palette. For modern populations color by language family
lang_cols <- c("Afro-Asiatic"  = "#e41a1c",
               "Indo-European" = "#377eb8",
               "Khoisan"       = "#4daf4a",
               "Niger-Congo"   = "#984ea3",
               "Nilo-Saharan"  = "#ff7f00",
               "Papuan"        = "#a65628",
               "Sino-Tibetan"  = "#f781bf")
label_cols <- lang_cols[pop_meta$Language[pop_meta$Group %in% modern]]
names(label_cols) <- modern
ancient_col    <- "black"   # colour used for ancient populations in all plots
label_cols_all <- c(label_cols, setNames(rep(ancient_col, length(ancient)), ancient))


# --- 1.5  World map of sampling locations ---

# Map data: use modern and ancient populations with valid coordinates 
modern_meta <- pop_meta[pop_meta$Group %in% modern & !is.na(pop_meta$Lat), ]
ancient_meta <- pop_meta[pop_meta$Group %in% ancient & !is.na(pop_meta$Lat), ]

world <- map_data("world")

ggplot() +
  geom_polygon(data=world, aes(x=long, y=lat, group=group),
               fill="grey85", colour="white", linewidth=0.15) +
  geom_point(data=modern_meta, aes(x=Long, y=Lat, fill=Language),
             shape=21, size=2.5, stroke=0.5) +
  geom_point(data=ancient_meta, aes(x=Long, y=Lat, colour="Ancient samples"),
             shape=18, size=2.5) +
  scale_fill_manual(values=lang_cols, name="Modern samples by language family") +
  scale_colour_manual(values=c("Ancient samples" = ancient_col), name=NULL) +
  coord_quickmap(xlim=range(modern_meta$Long) + c(-5, 5),
                 ylim=range(modern_meta$Lat)  + c(-5, 5)) +
  labs(title="Sampling locations", x=NULL, y=NULL) +
  theme_minimal()


# --- 1.6  Number of SNPs and distribution of private SNPs ---

n_snps <- length(readLines(paste0(prefix, ".snp")))
cat("SNPs:", n_snps, "\n")

afdat  <- eigenstrat_to_afs(prefix, pops = modern)
afs    <- afdat$afs

n_ind        <- table(info$Group[info$Group %in% modern])
total_called <- colSums(!is.na(afs[, modern]))

private_counts <- sapply(modern, function(pop) {
  others      <- setdiff(modern, pop)
  has_alt     <- !is.na(afs[, pop]) & afs[, pop] > 0
  others_zero <- rowSums(!is.na(afs[, others]) & afs[, others] > 0) == 0
  sum(has_alt & others_zero)
})

# Table 1: private SNPs per population
private_df <- data.frame(
  pop         = modern,
  n_ind       = as.integer(n_ind[modern]),
  non_missing = total_called,
  private     = private_counts
)
private_df

# Table 2: private SNPs per language family
lang_map <- setNames(pop_meta$Language[pop_meta$Group %in% modern], modern)

private_lang <- sapply(unique(lang_map), function(lang) {
  fam_pops   <- names(lang_map)[lang_map == lang]
  other_pops <- names(lang_map)[lang_map != lang]
  has_alt     <- rowSums(!is.na(afs[, fam_pops,   drop=FALSE]) & afs[, fam_pops,   drop=FALSE] > 0) > 0
  others_zero <- rowSums(!is.na(afs[, other_pops,  drop=FALSE]) & afs[, other_pops,  drop=FALSE] > 0) == 0
  sum(has_alt & others_zero)
})

lang_df <- data.frame(
  language    = unique(lang_map),
  n_ind       = sapply(unique(lang_map), function(l) {
                  sum(n_ind[names(lang_map)[lang_map == l]])
                }),
  non_missing = sapply(unique(lang_map), function(l) {
                  round(mean(total_called[names(lang_map)[lang_map == l]]))
                }),
  private     = private_lang
)
lang_df

################################################################################
# PART 2: F2 STATISTIC AND FST — GENETIC DIVERGENCE
################################################################################
#
# f2(A, B) estimates the genetic drift that has accumulated between populations
# A and B since they diverged from their common ancestor. Briefly, it is the
# expected squared allele frequency difference. The actual estimator implemented
# in ADMIXTOOLS v2 additionally accounts for low sample counts, missing data, 
# and differences in ploidy. Larger f2 means more accumulated drift (greater
# genetic distance).

# All f-statistics in ADMIXTOOLS v2 are derived from f2.
# The package uses a two-step workflow:
#   1. extract_f2() — reads genotype data, computes f2 for all population
#      pairs, and writes them to disk. Slow but only run once.
#   2. f2_from_precomp() — reads the stored f2 blocks. Fast.

# --- 2.1  Compute pairwise f2 ---

# By default maxmix is set to 0. Since we have ancient samples with substantial 
# missingness, in order not to lose too many SNPs, we allow 10% of missingness
extract_f2(prefix, f2_dir, maxmiss = 0.1)

# read in f2 blocks
# use afprod=TRUE to compute allele frequency products in addition to f2.
# This option will handle negative f2-statistic across blocks (useful if
# there are too many missing or rare SNPs in populations with low sample size)
f2_blocks <- f2_from_precomp(f2_dir, afprod = T)

dim(f2_blocks)
# The resulting object is a 3D array with f2-statistics for each population pair:
#   dim 1 × dim 2 = population pairs
#   dim 3         = SNP blocks
# The purpose of having separate estimates for each SNP block is to compute
# jackknife or bootstrap standard errors

f2_blocks[,,1]   # f2 matrix for the first SNP block (populations × populations)

# --- 2.2  Pairwise f2 heatmap (modern populations) ---

pairwise_f2 <- f2(f2_blocks, pop1=modern, pop2=modern)

# OPTIONAL: As a sanity check you can also compute f2 from genotypes for the 
# modern set only. The strict missingness filter still returns more than 1 
# million SNPs for this subset. In this case afprod=FALSE.
# f2_blocks_strict <- f2_from_geno(prefix, pops = c(modern, outgroup))
# pairwise_f2 <- f2(f2_blocks_strict, pop1=modern, pop2=modern)

# Exclude within-population pairs: f2(A,A)
f2_plot      <- pairwise_f2[pairwise_f2$pop1 != pairwise_f2$pop2, ]
f2_plot$pop1 <- factor(f2_plot$pop1, levels = modern)
f2_plot$pop2 <- factor(f2_plot$pop2, levels = modern)

ggplot(f2_plot, aes(x=pop1, y=pop2, fill=est)) +
  geom_tile() +
  scale_fill_gradientn(colours=hcl.colors(64, "Blues", rev=TRUE), name="f2") +
  labs(x=NULL, y=NULL, title="Pairwise f2 among modern human populations") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7,
                                  colour=label_cols_all),
        axis.text.y=element_text(size=7, colour=label_cols_all))

# Questions: How does the f2 pattern relate to geography and language family?
# Which populations are most diverged from each other? Which are closest?
# (OPTIONAL) Are there differences in the f2 results computed with different 
# SNP sets?


# --- 2.3  Pairwise f2 heatmap (modern + ancient populations) ---

all_pops <- c(modern, ancient)

pairwise_f2_all <- f2(f2_blocks, pop1=all_pops, pop2=all_pops)

f2_all_plot      <- pairwise_f2_all[pairwise_f2_all$pop1 != pairwise_f2_all$pop2, ]
f2_all_plot$pop1 <- factor(f2_all_plot$pop1, levels=all_pops)
f2_all_plot$pop2 <- factor(f2_all_plot$pop2, levels=all_pops)

ggplot(f2_all_plot, aes(x=pop1, y=pop2, fill=est)) +
  geom_tile() +
  scale_fill_gradientn(colours=hcl.colors(64, "Blues", rev=TRUE), name="f2") +
  labs(x=NULL, y=NULL, title="Pairwise f2 — modern and ancient populations") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7,
                                  colour=label_cols_all),
        axis.text.y=element_text(size=7, colour=label_cols_all))


# Questions: How do the ancient samples relate to each other?
# How do they related to the modern samples?


# --- 2.4  Pairwise FST among modern populations ---
#
# FST is the most widely used summary of population differentiation. Unlike f2,
# it is normalised by total diversity. FST per-block estimates were stored
# alongside f2 when extract_f2() was run (fst=TRUE is the default).
# ADMIXTOOLS v2 uses the Hudson FST estimator (Bhatia et al. 2013),
# which is unbiased for unequal sample sizes:
#
#   numerator   = (p1 - p2)^2 - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1)
#   denominator = p1 + p2 - 2·p1·p2
#   FST         = mean(numerator) / mean(denominator)   [averaged over SNPs]


# read FST from the f2 directory.
pairwise_fst <- fst(f2_dir, pop1=modern, pop2=modern)

# Exclude within-population pairs
fst_plot      <- pairwise_fst[pairwise_fst$pop1 != pairwise_fst$pop2, ]
fst_plot$pop1 <- factor(fst_plot$pop1, levels = modern)
fst_plot$pop2 <- factor(fst_plot$pop2, levels = modern)

ggplot(fst_plot, aes(x=pop1, y=pop2, fill=est)) +
  geom_tile() +
  scale_fill_gradientn(colours=hcl.colors(64, "Blues", rev=TRUE), name="FST") +
  labs(x=NULL, y=NULL, title="Pairwise FST among modern human populations") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7,
                                  colour=label_cols_all),
        axis.text.y=element_text(size=7, colour=label_cols_all))

# Question: Compare the FST heatmap to the f2 heatmap.
# Do the same pairs stand out? Are there any where FST seems high but f2 is not?


# --- 2.5  Comparing f2 and FST ---
#
# f2 and FST both quantify population divergence, but they capture different things:
#
#   f2(A,B)  = total accumulated drift between A and B (unbounded)
#   FST(A,B) = divergence normalised by within-population diversity (0–1)
#
# They should correlate positively: more divergence → higher FST.
# Deviations arise when populations differ in within-pop diversity (e.g.
# bottlenecked populations have low heterozygosity, inflating their FST
# relative to f2).

stats_df <- merge(
  pairwise_f2[pairwise_f2$pop1 < pairwise_f2$pop2,   c("pop1","pop2","est")],
  pairwise_fst[pairwise_fst$pop1 < pairwise_fst$pop2, c("pop1","pop2","est")],
  by = c("pop1","pop2"))
names(stats_df)[3:4] <- c("f2","fst")

# Identify outliers: FST disproportionately high relative to f2
# First, fit a simple linear regression
fit <- lm(fst ~ f2, data=stats_df)
# Second, find outliers based on the residuals (more than 2 standard deviations
# away from zero). For each population pair, the residual is the difference 
# between the observed FST and the FST predicted by the linear model.
outliers <- stats_df[abs(residuals(fit)) > 2 * sd(residuals(fit)), ]
outliers$pair <- paste0(outliers$pop1, " / ", outliers$pop2)

ggplot(stats_df, aes(x=f2, y=fst)) +
  geom_point(alpha=0.4, size=1.5, colour="steelblue") +
  geom_smooth(method="lm", se=FALSE, colour="grey40", linetype="dashed",
              linewidth=0.6) +
  geom_point(data=outliers, colour="#e41a1c", size=2) +
  geom_text(data=outliers, aes(label=pair), size=2.3, vjust=-0.7,
            colour="#e41a1c") +
  labs(x="f2", y="FST",
       title="f2 vs FST") +
  theme_minimal()

# Question: Are there populations that recurrently appear in the outlier pairs?
# Did these populations go through a recent bottleneck?


################################################################################
# PART 3: OUTGROUP F3 STATISTICS — SHARED DRIFT
################################################################################
#
# f3(O; A, B) measures the branch length from O to the common ancestor of A
# and B — the more drift A and B have accumulated together relative to O,
# the higher the f3 value. This is a measure of genetic similarity that is not 
# inflated by drift specific to A or B individually (e.g. population bottlenecks).
#
# Larger f3(O; A, B)  →  A and B share more drift (closer relatives)
# Smaller f3(O; A, B) →  A and B diverged earlier


# --- 3.1  Pairwise outgroup f3 heatmap (modern + ancient populations) ---

outgroup_f3_all <- f3(f2_blocks, pop1 = outgroup, pop2 = all_pops, pop3 = all_pops)
# We use Chimp as outgroup for all modern and archaic humans.
# Note that pop1 should always be the outgroup. This may not be the case outside
# ADMIXTOOLS v2 (e.g. admixr package)

# Exclude diagonal (same pop vs. same pop)
f3_all_plot      <- outgroup_f3_all[outgroup_f3_all$pop2 != outgroup_f3_all$pop3, ]
f3_all_plot$pop2 <- factor(f3_all_plot$pop2, levels=all_pops)
f3_all_plot$pop3 <- factor(f3_all_plot$pop3, levels=all_pops)

ggplot(f3_all_plot, aes(x=pop2, y=pop3, fill=est)) +
  geom_tile() +
  scale_fill_gradientn(colours=hcl.colors(64, "Blues", rev=TRUE), name="f3") +
  labs(x=NULL, y=NULL, title="Outgroup f3 — modern and ancient populations") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7,
                                  colour=label_cols_all),
        axis.text.y=element_text(size=7, colour=label_cols_all))

# Questions: Which populations share the most drift? Does this match what you
# know about human population history?


# --- 3.2  Targeted outgroup f3: who shares most drift with South_Africa_2000BP? ---
#
# There is a lot of information in the heatmap. We can also visualize results 
# for a specific ancient individual. Let's use South_Africa_2000BP — two 
# ancient individuals dated to ~2,000 years BP — as the target, and check which 
# populations share the most drift with it.

# Extract the relevant rows from the already-computed outgroup_f3_all object
sa_f3 <- outgroup_f3_all[outgroup_f3_all$pop2 == "South_Africa_2000BP" &
                           outgroup_f3_all$pop3 != "South_Africa_2000BP", ]

# Merge with coordinates for all populations
sa_map <- merge(sa_f3, pop_meta[, c("Group", "Lat", "Long")],
                by.x="pop3", by.y="Group")
sa_map$is_ancient <- sa_map$pop3 %in% ancient

# Location of the target ancient sample
target_loc <- pop_meta[pop_meta$Group == "South_Africa_2000BP", ]

# Visualize results in a map (shared drift projected onto modern sampling
# locations) or as a ranked dot plot with standard error bars.

ggplot() +
  geom_polygon(data=world, aes(x=long, y=lat, group=group),
               fill="grey85", colour="white", linewidth=0.15) +
  geom_point(data=sa_map, aes(x=Long, y=Lat, fill=est, shape=is_ancient),
             size=4, stroke=0.4) +
  geom_point(data=target_loc,
             aes(x=Long, y=Lat), shape=23, size=4, fill=NA, colour="red", stroke=1) +
  scale_shape_manual(values=c("FALSE"=21, "TRUE"=23),
                     labels=c("FALSE"="modern", "TRUE"="ancient"), name=NULL) +
  scale_fill_gradientn(colours=hcl.colors(64, "Blues", rev=TRUE),
                       name="f3") +
  coord_quickmap(xlim=range(modern_meta$Long) + c(-5, 5),
                 ylim=range(modern_meta$Lat)  + c(-5, 5)) +
  labs(title="Shared drift with South_Africa_2000BP — f3(Chimp; South_Africa_2000BP, X)",
       x=NULL, y=NULL) +
  theme_minimal()

# Ranked dot plot with SE bars
sa_f3 <- sa_f3[order(sa_f3$est), ]
sa_f3$pop3       <- factor(sa_f3$pop3, levels=sa_f3$pop3)
sa_f3$is_ancient <- sa_f3$pop3 %in% ancient
sa_cols <- label_cols_all[levels(sa_f3$pop3)]

ggplot(sa_f3, aes(x=est, y=pop3, shape=is_ancient)) +
  geom_point(size=2) +
  geom_errorbarh(aes(xmin=est - 2*se, xmax=est + 2*se), height=0.3) +
  scale_shape_manual(values=c("FALSE"=19, "TRUE"=18),
                     labels=c("FALSE"="modern", "TRUE"="ancient"),
                     name=NULL) +
  labs(x="f3(Chimp; South_Africa_2000BP, X)", y=NULL,
       title="Shared drift with South_Africa_2000BP") +
  theme_minimal() +
  theme(axis.text.y=element_text(colour=sa_cols))

# Question: Which populations share the most genetic history with South_Africa_2000BP?
# What does the ranking tell you about population continuity across space and 
# time in Africa?


################################################################################
# PART 4: ADMIXTURE F3 STATISTICS — TESTING ADMIXTURE
################################################################################
#
# A significantly negative f3(A; B, C) indicates that population A is admixed
# between sources related to B and C.

# --- 4.1  f3 admixture test for Neanderthal ancestry in non-Africans ---
#
# f3(non-African; Neanderthal, African)
# A negative value indicates the non-African carries ancestry from both a
# Neanderthal-related and a modern human source (here represented by an African
# group)
#
# neand_f3 <- f3(f2_blocks, pop1 = "French", pop2 = "Altai_Neanderthal", pop3 = "Mbuti")
# We could use pre-computed blocks, but since the test involves few pops, we
# can perform it on a stricter set of non-missing SNPs across these specific pops
f2_blocks_strict <- f2_from_geno(prefix, maxmiss = 0,
                                pops = c("French", "Altai_Neanderthal", "Mbuti"))
neand_f3 <- f3(f2_blocks_strict,
               pop1 = "French", pop2 = "Altai_Neanderthal", pop3 = "Mbuti")
neand_f3

# A non-significant result does NOT rule out admixture — the f3 test simply
# lacks power at this scale. Common reasons for lack of power include a low 
# admixture proportion and genetic drift since admixture. This test is most
# powerful for recent admixture events and high admixture levels (e.g., 20–50%).


################################################################################
# PART 5: F4 STATISTICS — DETECTING EXCESS ALLELE SHARING
################################################################################
#
# f4(A, B; C, D) measures the average product of allele frequency 
# differences between the pair (A vs. B) and (C vs. D). The f4 statistic is the
# unnormalised equivalent of the D-statistic (ABBA-BABA test). It directly
# measures excess allele sharing due to gene flow without being diluted by drift
# specific to any individual population. Under a tree with no gene flow, f4 = 0.
# A significantly non-zero value signals that the tree topology is inadequate
# (i.e. gene flow has occurred).
#
# f4(A, B; C, D) > 0  →  A and C share more alleles than B and C
#                         equivalently: B and D share more alleles than A and D
# f4(A, B; C, D) < 0  →  B and C share more alleles than A and C
#                         equivalently: A and D share more alleles than B and D
#
# Symmetry: f4(A, B; C, D) = −f4(B, A; C, D) = −f4(A, B; D, C)
#
# Neanderthal and Denisovan introgression are well-documented cases of
# archaic gene flow into modern humans:
#   - All non-African populations carry ~1–3% Neanderthal ancestry, acquired
#     after the Out-of-Africa dispersal.
#   - Papuans carry an additional ~3–5% Denisovan ancestry.
#
# Let's find evidence for archaic introgression in our dataset!
# Here we use Yoruba (African population) as the unadmixed modern human
# and Vindija as the Neanderthal source


# --- 5.1  Neanderthal and Denisovan introgression tests ---

non_africans <- pop_meta$Group[pop_meta$Continent != "Africa" & pop_meta$Group %in% modern]

# Compute f2 blocks for all relevant populations in one pass
f2_blocks_strict <- f2_from_geno(prefix, maxmiss = 0,
                                pops = c(outgroup, "Vindija_Neanderthal",
                                "Altai_Neanderthal", denisovan,
                                "Yoruba", non_africans))

neand_f4 <- f4(f2_blocks_strict, maxmiss = 0,
               pop1 = outgroup,
               pop2 = "Vindija_Neanderthal",
               pop3 = "Yoruba",
               pop4 = non_africans)
neand_f4$archaic <- "Neanderthal"

denisova_f4 <- f4(f2_blocks_strict, maxmiss = 0,
                  pop1 = outgroup,
                  pop2 = denisovan,
                  pop3 = "Yoruba",
                  pop4 = non_africans)
denisova_f4$archaic <- "Denisova"

archaic_f4        <- rbind(neand_f4, denisova_f4)
archaic_f4$sig    <- abs(archaic_f4$z) > 3
archaic_f4$label  <- factor(archaic_f4$pop4,
                             levels = neand_f4$pop4[order(neand_f4$est)])

ggplot(archaic_f4, aes(x=label, y=est, colour=archaic, shape=sig)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey40") +
  geom_point(position=position_dodge(width=0.5), size=2.5) +
  geom_errorbar(aes(ymin=est - 2*se, ymax=est + 2*se),
                position=position_dodge(width=0.5), width=0.25) +
  scale_colour_manual(values=c("Neanderthal"="steelblue", "Denisova"="#ff7f00"),
                      name="Archaic") +
  scale_shape_manual(values=c("FALSE"=1, "TRUE"=16),
                     labels=c("FALSE"="not significant", "TRUE"="|z| > 3"),
                     name="") +
  labs(x=NULL, y="f4 (Chimp, Archaic; Yoruba, X)",
       title="Archaic introgression across non-African populations") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# Question: Are there differences in the Neanderthal signal?
# What could cause these differences?



################################################################################
# PART 6: F4-RATIO — QUANTIFYING ADMIXTURE PROPORTIONS
################################################################################
#
# The f4-ratio (qpf4ratio in ADMIXTOOLS v2) estimates the admixture
# proportion alpha from an f4 statistic divided by another f4 statistic.
#
# The function qpf4ratio() computes:
#   alpha = f4(pop1, pop2; pop3, pop4) / f4(pop1, pop2; pop5, pop4)
#
# In the Neanderthal introgression context (Petr et al. 2019):
#   pop1 = Altai Neanderthal   (sister group of introgressing Neanderthal)
#   pop2 = Chimp               (outgroup)
#   pop3 = test population     (e.g., French)
#   pop4 = an African without Neanderthal ancestry (e.g., Yoruba)
#   pop5 = Vindija Neanderthal (introgressing Neanderthal)


# --- 6.1  Estimate Neanderthal admixture for non-African populations and plot ---

# qpf4ratio() accepts a matrix with 5 columns (one row per test population)
pops_matrix <- cbind("Altai_Neanderthal", outgroup,
                     non_africans, "Yoruba", "Vindija_Neanderthal")

f4ratio_res <- qpf4ratio(f2_blocks_strict, pops=pops_matrix)
f4ratio_res$label <- non_africans

# Plot
par(mar=c(6, 5, 4, 2))
bp3 <- barplot(f4ratio_res$alpha * 100,
               names.arg = f4ratio_res$label,
               las=2, cex.names=0.8,
               col="steelblue",
               ylab="Estimated Neanderthal ancestry (%)",
               main="Neanderthal admixture proportion\n(f4-ratio, Vindija as proxy)",
               ylim=c(0, max(f4ratio_res$alpha + 2*f4ratio_res$se) * 100 * 1.1))
arrows(bp3,
       (f4ratio_res$alpha - 2*f4ratio_res$se) * 100,
       bp3,
       (f4ratio_res$alpha + 2*f4ratio_res$se) * 100,
       angle=90, code=3, length=0.05, lwd=1)
par(mar=c(5, 4, 4, 2))


################################################################################
# END OF TUTORIAL
#
# For further reading:
#   - Patterson et al. (2012) Genetics
#   - Peter (2016) Genetics
#   - Petr et al. (2019) PNAS
#   - ADMIXTOOLS v2 documentation
################################################################################
