# EMBO Practical Course: Genomic Diversity & Natural Selection Scan (Worksheet)

This worksheet guides you through the analysis of genomic data in humans and canines to identify genomic regions under natural selection. The tutorial is divided into two main parts:
  
1. **Part 1: Human Genomic Diversity and Natural Selection**: Investigating selection signatures in the candidate gene ***EDAR*** (associated with ectodermal traits in East Asians and Native Americans) using population differentiation ($F_{ST}$, PBS) and haplotype-based metrics (EHH, iHS, XP-EHH).
2. **Part 2: Genomic Selection Scan in Canines**: Identifying the selective sweep at the ***IGF1*** body-size locus by comparing small vs. large dog breeds using PCA, PCAdapt (outlier scan), and haplotype homozygosity methods (XP-nSL and Rsb).

---
  
  # Part 1: Human Genomic Diversity and Natural Selection
  
  ## 1. Background and Dataset
  
  ### Goal
  Our goal is to explore approaches and methods which seek to identify regions of the genome with signatures of natural selection. We will use real genomic data and two classes of tests: one based on population differentiation ($F_{ST}$ / PBS) and another based on extended haplotype homozygosity (EHH / iHS / XP-EHH).

### Dataset
Whole-genome sequencing data from the 1000 Genomes Project Phase III. The full database can be accessed via:
  <ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>
  
  ### Data Pre-processing
  We will analyze a pre-processed dataset for chromosome 2 corresponding to individuals sampled from the African (AFR: 504 individuals), European (EUR: 503 individuals), and East Asian (EAS: 504 individuals) populations. In this dataset, INDELs, singletons, and SNPs with MAF < 0.05 have been removed. The pairwise $F_{ST}$ was then estimated using `vcftools`.

All data files are located in the `input/` directory:
# - `input/Part_1_HumanDiversity/AFR_EAS.weir.fst` (Fst between Africans and East Asians)
- `input/Part_1_HumanDiversity/AFR_EUR.weir.fst` (Fst between Africans and Europeans)
- `input/Part_1_HumanDiversity/EAS_EUR.weir.fst` (Fst between East Asians and Europeans)
- `input/Part_1_HumanDiversity/Chr2_EDAR_LWK_500K.recode.vcf` (Phased African haplotypes around *EDAR*)
- `input/Part_1_HumanDiversity/Chr2_EDAR_CHS_500K.recode.vcf` (Phased East Asian haplotypes around *EDAR*)
- `input/Part_1_HumanDiversity/Chr2_NAM_EAS.weir.fst` (Fst between Native Americans and East Asians in candidate region)
- `input/Part_1_HumanDiversity/Chr2_NAM_EUR.weir.fst` (Fst between Native Americans and Europeans in candidate region)
- `input/Part_1_HumanDiversity/Chr2_EUR_EAS.weir.fst` (Fst between Europeans and East Asians in candidate region)

---
  # ==============================================================================
# R Code Exercise: Pairwise FST Calculation (Soluzione Completa)
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# PUNTO 1: Read the pairwise FST files from input/
# ------------------------------------------------------------------------------
# Leggiamo i file specificando la cartella 'input/' come richiesto
fst_afr_eas <- read_table("AFR_EAS.weir.fst", col_names = FALSE, skip = 1)
fst_afr_eur <- read_table("AFR_EUR.weir.fst", col_names = FALSE, skip = 1)
fst_eas_eur <- read_table("EAS_EUR.weir.fst", col_names = FALSE, skip = 1)

# Assegniamo i nomi corretti alle colonne per poter lavorare
colnames(fst_afr_eas) <- c("CHROM", "POS", "FST_AFR_EAS")
colnames(fst_afr_eur) <- c("CHROM", "POS", "FST_AFR_EUR")
colnames(fst_eas_eur) <- c("CHROM", "POS", "FST_EAS_EUR")

clean_afr_eas <- fst_afr_eas %>% filter(!is.na(FST_AFR_EAS)) %>% distinct(POS, .keep_all = TRUE) %>% mutate(FST_AFR_EAS = ifelse(FST_AFR_EAS < 0, 0, FST_AFR_EAS))
clean_afr_eur <- fst_afr_eur %>% filter(!is.na(FST_AFR_EUR)) %>% distinct(POS, .keep_all = TRUE) %>% mutate(FST_AFR_EUR = ifelse(FST_AFR_EUR < 0, 0, FST_AFR_EUR))
clean_eas_eur <- fst_eas_eur %>% filter(!is.na(FST_EAS_EUR)) %>% distinct(POS, .keep_all = TRUE) %>% mutate(FST_EAS_EUR = ifelse(FST_EAS_EUR < 0, 0, FST_EAS_EUR))

# ------------------------------------------------------------------------------
# 3. Allineamento dei dataset (Inner Join)
# ------------------------------------------------------------------------------
aligned_data <- clean_afr_eas %>%
  inner_join(clean_afr_eur, by = c("CHROM", "POS")) %>%
  inner_join(clean_eas_eur, by = c("CHROM", "POS"))

# ------------------------------------------------------------------------------
# 4. Verifica posizione target e calcolo Quantili (Outlier)
# ------------------------------------------------------------------------------
target_pos <- 109513601
target_snp <- aligned_data %>% filter(POS == target_pos)

cat("\n--- VALORI DI FST ALLA POSIZIONE TARGET ---\n")
print(target_snp)

# Quantili complessivi per capire se è un outlier nel confronto principale (AFR_EAS)
quantiles_afr_eas <- quantile(aligned_data$FST_AFR_EAS, probs = c(0.5, 0.9, 0.95, 0.99, 0.999))
cat("\n--- QUANTILI DELLA DISTRIBUZIONE (AFR_EAS) ---\n")
print(quantiles_afr_eas)

# ------------------------------------------------------------------------------
# 5. Preparazione dati per il Regional Manhattan Plot (Tutti i confronti insieme)
# ------------------------------------------------------------------------------
# Definiamo i confini della finestra di 10kb (±5000 bp dal target)
window_size <- 5000 
start_pos <- target_pos - window_size
end_pos <- target_pos + window_size

# Filtriamo i dati nella finestra e li trasformiamo in formato LONG per ggplot
window_data_long <- aligned_data %>%
  filter(POS >= start_pos & POS <= end_pos) %>%
  pivot_longer(cols = starts_with("FST_"), names_to = "Confronto", values_to = "FST") %>%
  mutate(Confronto = str_replace(Confronto, "FST_", "")) # Pulisce i nomi in AFR_EAS, AFR_EUR, ecc.

# Facciamo lo stesso per il punto target da evidenziare
target_snp_long <- target_snp %>%
  pivot_longer(cols = starts_with("FST_"), names_to = "Confronto", values_to = "FST") %>%
  mutate(Confronto = str_replace(Confronto, "FST_", ""))

# ------------------------------------------------------------------------------
# 6. GENERAZIONE REGIONAL MANHATTAN PLOT
# ------------------------------------------------------------------------------
p <- ggplot(window_data_long, aes(x = POS, y = FST, color = Confronto)) +
  # Linee di andamento per ciascun confronto di popolazione
  geom_line(alpha = 0.4, size = 0.6) +
  
  # Punti stile Manhattan divisi per colore del confronto di popolazione
  geom_point(alpha = 0.8, size = 2.5) +
  
  # Evidenzia lo SNP candidato rs3827760 con un rombo con bordo nero su tutti i confronti
  geom_point(data = target_snp_long, aes(x = POS, y = FST), color = "black", size = 4.5, shape = 23, fill = "red") +
  
  # Mette il testo del gene solo sul punto più alto (AFR vs EAS) per non sovrapporsi
  geom_text(data = target_snp_long %>% filter(Confronto == "AFR_EAS"), 
            aes(x = POS, y = FST, label = "rs3827760 (EDAR)"), 
            vjust = -1.5, color = "black", fontface = "bold", size = 3.8) +
  
  # Impostazioni grafiche e assi
  scale_x_continuous(labels = scales::comma) +
  scale_color_manual(values = c("AFR_EAS" = "firebrick2", "AFR_EUR" = "dodgerblue3", "EAS_EUR" = "forestgreen")) +
  ylim(0, 1.05) +
  labs(
    title = "Regional Manhattan Plot (Locus EDAR)",
    subtitle = paste("Finestra Genomica di 10kb su Cromosoma", target_snp$CHROM[1]),
    x = "Posizione Genomica (bp)",
    y = expression("Pairwise " ~ italic(F)[ST]),
    color = "Confronto"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# Mostra il grafico a schermo
print(p)



## 2. Genetic Differentiation ($F_{ST}$ and PBS)
  
  ### Investigating the Candidate Gene *EDAR*
  The human Ectodysplasin A receptor gene, or ***EDAR***, is part of the EDA signaling pathway which specifies prenatally the location, size, and shape of ectodermal appendages (such as hair follicles, teeth, and glands). *EDAR* is a textbook example of positive selection in East Asians. A specific non-synonymous variant, **rs3827760** (chr2:109,513,601 A>G), results in a Val370Ala substitution and is strongly associated with thicker hair shafts and shovel-shaped incisors. Another hypothesis states that *EDAR* acted along with *FADS* and *VDR* in the Beringia Standstill, allowing Native American ancestors to survive in extreme arctic environments.

### Questions for Students

1. **The estimate of $F_{ST}$ by the Weir and Cockerham metric can sometimes generate negative values and "NA". What does that mean? How can this interfere with the results?**
  * *Answer*:
  2. **The $F_{ST}$ values observed between pairs of populations for the SNP rs3827760 (position 109,513,601) fall within which distribution quantiles of $F_{ST}$ values for the studied chromosome? Can they be considered outliers?**
  * *Answer*:
  3. **From the observed $F_{ST}$ values between population pairs and the significance estimates, what can we say about the rs3827760 SNP differentiation between populations?**
  * *Answer*:
  4. **Discuss how these results justify performing another type of analysis based on PBS (Population Branch Statistics).**
  * *Answer*:
  5. **What does the PBS analysis reveal? What is the difference between PBS and $F_{ST}$ analysis?**
  * *Answer*:
  
  ---
  
  ### R Code Exercise: Pairwise $F_{ST}$ Calculation
  Write the R code necessary to perform the following:
  1. Read the pairwise $F_{ST}$ files from `input/`.
2. Filter duplicate SNP positions and exclude NA values.
3. Align the datasets by overlapping positions.
4. Set negative $F_{ST}$ values to zero.
5. Check Fst values at position `109513601`.
6. Calculate distribution quantiles to determine if rs3827760 is an outlier.
7. Plot pairwise Fst around `109513601` in a 10kb window, highlighting the candidate SNP.

**Write your R code here:**
  ```R
# 


```

---
  
  ### R Code Exercise: Population Branch Statistics (PBS)
  Write the R code necessary to:
  1. Estimate the Population Branch Statistic for East Asians ($PBS_{EAS}$) using the AFR, EAS, and EUR populations.
2. Convert negative branch lengths to zero.
3. Check the PBS value for the candidate SNP rs3827760.
4. Calculate distribution quantiles to determine if it is an outlier.
5. Plot PBS values around the candidate SNP in a 10kb window.

**Write your R code here:**
  ```R
pbs_data <- aligned_data %>%
  mutate(
    # Trasformazione in distanze evolutive (T)
    T_AFR_EAS = -log(1 - ifelse(FST_AFR_EAS >= 1, 0.9999, FST_AFR_EAS)),
    T_AFR_EUR = -log(1 - ifelse(FST_AFR_EUR >= 1, 0.9999, FST_AFR_EUR)),
    T_EAS_EUR = -log(1 - ifelse(FST_EAS_EUR >= 1, 0.9999, FST_EAS_EUR)),
    
    # Formula matematica della PBS per il ramo East Asian (EAS)
    PBS_EAS = (T_AFR_EAS + T_EAS_EUR - T_AFR_EUR) / 2
  ) %>%
  # Converti le lunghezze dei rami negative a zero (come richiesto dalla traccia)
  mutate(PBS_EAS = ifelse(PBS_EAS < 0, 0, PBS_EAS))

target_pos <- 109513601

target_pbs_snp <- pbs_data %>% 
  filter(POS == target_pos)

cat("\n--- VALORE DI PBS_EAS ALLA POSIZIONE 109513601 ---\n")
print(target_pbs_snp %>% select(CHROM, POS, PBS_EAS))

# ------------------------------------------------------------------------------
# 3. Calcolo dei quantili della distribuzione per determinare l'outlier
# ------------------------------------------------------------------------------
quantiles_pbs <- quantile(pbs_data$PBS_EAS, probs = c(0.5, 0.9, 0.95, 0.99, 0.999))

cat("\n--- QUANTILI DELLA DISTRIBUZIONE GENERALE DI PBS_EAS ---\n")
print(quantiles_pbs)

# Calcola l'esatto rango percentile genomico per rs3827760
percentile_rank_pbs <- mean(pbs_data$PBS_EAS < target_pbs_snp$PBS_EAS[1]) * 100
cat("\nLo SNP rs3827760 si trova nel", round(percentile_rank_pbs, 4), "° percentile genomico per la PBS_EAS.\n")

# ------------------------------------------------------------------------------
# 4. Filtro per la finestra di 10kb intorno al candidato
# ------------------------------------------------------------------------------
window_size <- 5000 
start_pos <- target_pos - window_size
end_pos <- target_pos + window_size

pbs_window_data <- pbs_data %>%
  filter(POS >= start_pos & POS <= end_pos)

# Estraiamo il valore del 99° percentile da usare come linea di soglia nel grafico
soglia_pbs_99 <- quantiles_pbs["99%"]

# ------------------------------------------------------------------------------
# 5. Grafico dei valori di PBS (Regional Manhattan Plot)
# ------------------------------------------------------------------------------
p_pbs <- ggplot() +
  # Linea di andamento grigia per unire i valori di PBS nella finestra
  geom_line(data = pbs_window_data, aes(x = POS, y = PBS_EAS), color = "grey75", size = 0.5) +
  
  # Punti stile Manhattan colorati in gradiente in base all'intensità della PBS
  geom_point(data = pbs_window_data, aes(x = POS, y = PBS_EAS, color = PBS_EAS), alpha = 0.8, size = 2.5) +
  scale_color_gradientn(colors = c("cadetblue3", "darkorange", "firebrick2")) +
  
  # Linea tratteggiata orizzontale che indica il Top 1% del genoma (99° percentile)
  geom_hline(yintercept = soglia_pbs_99, linetype = "dashed", color = "darkgrey", size = 0.7) +
  annotate("text", x = start_pos + 1500, y = soglia_pbs_99 + (max(pbs_window_data$PBS_EAS)*0.03), 
           label = "Soglia Top 1% Genoma", color = "grey40", size = 3.5, fontface = "italic") +
  
  # Evidenzia graficamente lo SNP candidato rs3827760 (EDAR) con un rombo rosso
  geom_point(data = target_pbs_snp, aes(x = POS, y = PBS_EAS), color = "black", size = 5, shape = 23, fill = "red") +
  
  # Etichetta di testo sopra lo SNP candidato
  geom_text(data = target_pbs_snp, aes(x = POS, y = PBS_EAS, label = "rs3827760 (EDAR)"), 
            vjust = -1.5, color = "black", fontface = "bold", size = 4) +
  
  # Ottimizzazione degli assi e formattazione con i punti delle migliaia
  scale_x_continuous(labels = scales::comma) +
  # Adatta dinamicamente l'altezza dell'asse Y per far spazio alle scritte
  ylim(0, max(pbs_window_data$PBS_EAS) * 1.15) + 
  
  # Etichette del grafico
  labs(
    title = "Regional PBS Plot - Locus EDAR (East Asian Selection)",
    subtitle = paste("Finestra Genomica di 10kb su Cromosoma", target_pbs_snp$CHROM[1]),
    x = "Posizione Genomica (bp)",
    y = expression("Population Branch Statistic (" * italic(PBS)[EAS] * ")"),
    color = expression(italic(PBS)[EAS])
  ) +
  
  # Tema pulito per report scientifici
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# Mostra il grafico della PBS a schermo
print(p_pbs)
```

---
  
  ## 3. Extended Haplotype Homozygosity (EHH)
  
  ### Extended Haplotype Homozygosity (EHH) and Haplotype Sweeps
  Different approaches can detect genomic signatures of selection at different timescales. More recent selection signals can be detected from haplotype-based tests. Positive selection causes a rapid rise in the frequency of the selected allele, such that recombination does not have enough time to break down the haplotype on which the mutation arose. This creates a signature of **Extended Haplotype Homozygosity (EHH)** extending over a long physical distance.

### Questions for Students

1. **How is the haplotype profile of genetic variants under recent positive selection?**
  * *Answer*:
  2. **What is the profile of ancestral and derived haplotypes of the rs3827760 SNP in AFR and EAS?**
  * *Answer*:
  3. **The iHS score observed for the SNP rs3827760 falls within which distribution quantiles of iHS values for the studied chromosome? Can it be considered an outlier? How can we make this analysis more robust?**
  * *Answer*:
  4. **What information does the XP-EHH analysis add about natural selection in the candidate SNP?**
  * *Answer*:
  
  ---
  
  ### R Code Exercise: EHH & Furcation Trees
  Write the R code necessary to:
  1. Convert the VCF databases to `rehh` format using `data2haplohh()`.
2. Estimate the EHH decay for rs3827760 in both populations.
3. Plot the EHH decay and furcation trees for both AFR and EAS.

**Write your R code here:**
  ```R
# 


```

---
  
  ### R CodeExercise: iHS & XP-EHH (Window-based)
  Write the R code necessary to:
  1. Perform a genome-wide scan of homozygosity using `scan_hh()` for AFR and EAS.
2. Calculate iHS scores for both populations using `ihh2ihs()`.
3. Check the iHS score at rs3827760 and generate a single-site iHS plot in EAS.
4. Create a function to estimate the average absolute iHS in sliding windows (50 SNPs/40 step) and plot the results.
5. Estimate cross-population XP-EHH between EAS and AFR using `ies2xpehh()`, calculate window-based averages, and plot them.

**Write your R code here:**
  ```R
# 


```

---
  
  ## 4. Native American Selection Analysis
  
  ### Background
  Hlusko et al. (2018), using morphological data, found a strong selection signal in the *EDAR* gene in Native Americans. Using the additional database from the 1000 Genomes Project (Peruvian samples with over 95% Native American Ancestry, represented as **NAM**), we evaluate genomic signatures of selection at the functional variant rs3827760.

### Questions for Students

1. **Is the functional allele in East Asian at high frequency in other human populations (e.g. Native Americans)?**
  * *Answer*:
  2. **Can we identify signatures of natural selection on EDAR in Native Americans using PBS?**
  * *Answer*:
  3. **Is selection targeting the same functional variant?**
  * *Answer*:
  4. **What is your conclusion based on the data generated?**
  * *Answer*:
  
  ---
  
  ### R Code Exercise: PBS in Native Americans (NAM)
  Write the R code necessary to:
  1. Read the pairwise $F_{ST}$ files involving Native Americans (`input/Part_1_HumanDiversity/Chr2_NAM_EAS.weir.fst` and `input/Part_1_HumanDiversity/Chr2_NAM_EUR.weir.fst`) and Europeans-East Asians (`input/Part_1_HumanDiversity/Chr2_EUR_EAS.weir.fst`).
2. Filter duplicates, exclude NA values, and align positions.
3. Convert negative $F_{ST}$ values to zero.
4. Estimate $PBS_{NAM}$ using NAM, EAS, and EUR.
5. Check PBS value at rs3827760, check quantiles, and plot the PBS scan.

**Write your R code here:**
  ```R
# 


```

---
  ---
# Part 2: Genomic Selection Sweep Scan in Canines
  
  ## 1. Background and Dataset
  
  The dataset is sourced from the **Dog10K** consortium ([Download Link](https://dog10k.kiz.ac.cn/Home/Download)). The original genomic dataset is a high-coverage phased BCF file containing 1,929 individuals and over 29 million SNPs:
  - Original file: `AutoAndXPAR.Dog10K.phased.bcf`
- Metadata table: `dog10K-alignment-sample-table.2022-02-23-v7.txt`

### Sample Selection
For this analysis, we will use a subset of **130 individuals** representing body size extremes:
  
  | Group | Dog Breed (Breed.Type) | Number of Samples |
  | :--- | :--- | :---: |
  | **Small** (61) | Dachshund | 17 |
  | | Toy Fox Terrier | 10 |
  | | Pomeranian | 8 |
  | | Brussels Griffon | 7 |
  | | Yorkshire Terrier | 5 |
  | | Shih Tzu | 5 |
  | | Maltese | 4 |
  | | Pekingese | 3 |
  | | Chihuahua | 2 |
  | **Large** (69) | Saint Bernard | 13 |
  | | Leonberger | 11 |
  | | Bernese Mountain Dog | 10 |
  | | Greater Swiss Mountain Dog | 10 |
  | | Great Pyrenees | 7 |
  | | Bullmastiff | 6 |
  | | Mastiff | 6 |
  | | Newfoundland | 6 |
  
  ---
  

bcftools view \
-S subset_dogs.txt \
-r chr15 \
-q 0.05:minor \
-O b \
-o subset_chr15.bcf \
AutoAndXPAR.Dog10K.phased.bcf
  ## 2. Preprocessing and Filtering
  
  We filter the massive BCF file to include only our 130 samples and chromosome 15, while removing low-frequency SNPs (MAF < 0.05) that are not informative for breed/size differentiation. The coordinates correspond to the CanFam3 reference genome, where the ***IGF1*** gene is located approximately in the **41.2 Mb to 44.5 Mb** region.

### Bash Code Exercise 1: Extraction & Format Conversion
```bash
# 1. Extract chr15 and filter for samples and MAF >= 0.05
bcftools view \
-S input/Part_2_CanidDiversity/subset_dogs.txt \
-r chr15 \
-q 0.05:minor \
-O b \
-o input/Part_2_CanidDiversity/subset_chr15.bcf \
input/Part_2_CanidDiversity/AutoAndXPAR.Dog10K.phased.bcf

# 2. Index the subset BCF
bcftools index input/Part_2_CanidDiversity/subset_chr15.bcf

# 3. Convert BCF to PLINK binary format (.bed/.bim/.fam)
plink1.9 \
--bcf input/Part_2_CanidDiversity/subset_chr15.bcf \
--dog \
--keep-allele-order \
--make-bed \
--out output/subset_chr15
```

> **Premise**: The filtered BCF file contains **177,953 SNPs** on chromosome 15 across the 130 samples. The coordinates correspond to the CanFam3 reference genome, where the ***IGF1*** gene is located approximately in the **41.2 Mb to 44.5 Mb** region.

---
  
  ## 3. Population Structure Analysis (PCA)
  
  Before checking for selection outliers, we must examine the genetic structure of our subset. We will run a PCA using **PLINK 1.9** and visualize it in R.

### R Code Exercise 2: PCA Visualization in R
Write the R code necessary to:
  1. Load the eigenvectors and eigenvalues generated by PLINK.
2. Merge them with the sample metadata (`input/Part_2_CanidDiversity/sample_info.txt`).
3. Calculate the percentage of variance explained by PC1 and PC2.
4. Generate a scatter plot of PC1 vs. PC2, coloring the points by breed and shaping them by size group (small vs. large).

**Write your R code here:**
  ```R
# 
eigenvec <- read.table("plink_pca.eigenvec", header = FALSE)
eigenval <- scan("plink_pca.eigenval")
categories <- read.table("sample_info.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
num_pcs <- ncol(eigenvec) - 2
colnames(eigenvec) <- c("FID", "IID", paste0("PC", 1:num_pcs))
colnames(categories)[colnames(categories) == "sampleName"] <- "IID"

# Unisci i dati per IID
df <- merge(eigenvec, categories, by = "IID")

# Leggi i valori propri
var_explained <- eigenval / sum(eigenval) * 100

library(ggplot2)
PCA <- ggplot(df, aes(x = PC1, y = PC2)) +
  geom_point(data = subset(df, !is.na(IID)), aes(color = group, shape = breed), alpha = 0.9, size = 8) +
  geom_point(data = subset(df, is.na(IID)), aes(shape = breed), color = "grey50", alpha = 0.9, size = 8) +
  scale_color_gradient(low = "green", high = "red") +  
  labs(
    x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_minimal() +
  theme(
    # axis.text: colore nero e font bold
    axis.text = element_text(size = 18, face = "bold", color = "black"), 
    
    # axis.title: colore nero, font bold e aggiunta di margin
    # t = top, r = right, b = bottom, l = left
    axis.title.x = element_text(size = 20, face = "bold", color = "black", 
                                margin = margin(t = 20)), # Spazio sopra il titolo X
    axis.title.y = element_text(size = 20, face = "bold", color = "black", 
                                margin = margin(r = 20)), # Spazio a destra del titolo Y
    
    legend.title = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 20, color = "black"),
    
    # Opzionale: aggiunge spazio bianco intorno a tutto il plot
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20) 
  )
PCA
```

### Questions for Students
1. **What pattern do you observe along the first principal component (PC1)?**
  * *Answer*:
  2. **Why does PC1 capture body size differences in this particular dataset?**
  * *Answer*:
  
  ---
  
## 4. Genomic Outlier Detection using PCAdapt
  
  **PCAdapt** is a method designed to find SNPs that are exceptionally related to population structure (PCs) rather than neutral drift.

### The Role of Linkage Disequilibrium (LD) and Clumping
> [!IMPORTANT]
> **Key Concept**: PCA is highly sensitive to Linkage Disequilibrium (LD). If a region contains many highly correlated markers (due to a selective sweep or low recombination), that single region will dominate the principal components, biasing the PCA and masking other genomic signals (such as the *IGF1* gene sweep).
> 
  > To resolve this, we must enable **LD Clumping** in PCAdapt. Thinning out redundant SNPs in strong LD allows the global genomic structure to be correctly computed and helps locate narrow selection sweeps.

### R Code Exercise 3: PCAdapt with LD Clumping
Write the R script necessary to:
  1. Load the genotypes into PCAdapt.
2. Execute `pcadapt()` using $K = 2$ components and enable **LD clumping** (for example, with a window size of 500 SNPs and an $r^2$ threshold of 0.1).
3. Merge the resulting p-values with the physical genomic positions from the `.bim` file.
4. Generate a Manhattan Plot of the results, plotting physical position (in Mb) on the X-axis and $-\log_{10}(\text{p-value})$ on the Y-axis.

**Write your R code here:**
  ```R
# 
# Load required libraries
library(pcadapt)
library(ggplot2)

# ==========================================
# 1. Load Genotypes and Enable LD Clumping
# ==========================================

# Define paths to your PLINK files
bed <- "subset_chr15.bed"
bim <- "subset_chr15.bim"

# Load the genotype data into pcadapt
genotypes <- read.pcadapt(bed, type = "bed")

# Run pcadapt with K = 2 and enable LD Clumping
# window.size = 500 SNPs, r2 threshold = 0.1
res_pcadapt <- pcadapt(genotypes, 
                       K = 2, 
                       LD.clumping = list(size = 500, thr = 0.1))

# ==========================================
# 2. Merge P-values with .bim physical data
# ==========================================

# Read the PLINK .bim file (contains Chromosome and Physical Positions)
# PLINK .bim format has no header: 
# V1 = Chromosome, V2 = SNP ID, V3 = Morgans, V4 = Base Pair Position, V5 = Allele 1, V6 = Allele 2
bim_data <- read.table(bim, header = FALSE, stringsAsFactors = FALSE)
colnames(bim_data) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

# Create a data frame combining the SNP positions with pcadapt p-values
# Note: pcadapt returns p-values for ALL SNPs; those clumped out will have NA p-values
results_df <- data.frame(
  CHR = bim_data$CHR,
  BP = bim_data$BP,
  P_VAL = res_pcadapt$pvalues
)

# Remove NA values (the SNPs that were thinned out by LD clumping)
results_df <- subset(results_df, !is.na(P_VAL))

# Calculate -log10(p-value) and convert physical position BP to Megabases (Mb)
results_df$logP <- -log10(results_df$P_VAL)
results_df$Position_Mb <- results_df$BP / 1e6

# ==========================================
# 3. Generate Manhattan Plot
# ==========================================

# Ensure Chromosome is treated as a factor for clean coloring
results_df$CHR <- as.factor(results_df$CHR)

manhattan_plot <- ggplot(results_df, aes(x = Position_Mb, y = logP, color = CHR)) +
  geom_point(alpha = 0.7, size = 1.5) +
  
  # Alternating colors for chromosomes (optional but makes it clean)
  scale_color_manual(values = rep(c("slateblue4", "orange2"), length.out = length(unique(results_df$CHR)))) +
  
  # Group points by chromosome if plotting genome-wide
  facet_grid(. ~ CHR, scales = "free_x", space = "free_x") +
  
  # Labels
  labs(
    title = "PCAdapt Manhattan Plot (with LD Clumping)",
    x = "Physical Position (Mb)",
    y = expression(-log[10](p-value))
  ) +
  
  # Styling matching your previous plots
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none", # Hide legend since facets show chromosomes
    panel.spacing = unit(0.2, "lines"), # Reduce spacing between chromosome columns
    strip.text = element_text(size = 12, face = "bold") # Chromosome labels on top
  )

# Display the plot
print(manhattan_plot)

```

### Questions for Students
1. **If you ran PCAdapt without enabling LD clumping, a single massive peak at ~61 Mb would dominate the plot, hiding other regions. What is the effect of LD clumping on outlier detection and why is it necessary here?**
  * *Answer*:
  2. **Did you detect the outlier peak at the *IGF1* locus? What is the approximate coordinate of the peak and what is its biological significance?**
  * *Answer*:
  
  ---
  
  ## 5. Cross-Population Selection Scan using XP-nSL
  
  To confirm that the outlier peak on chromosome 15 is indeed driven by a selective sweep, we will perform a haplotype-based selection scan. Specifically, we will run **XP-nSL (Cross-Population Number of Segregating Sites by Length)** to compare the haplotype homozygosity decay between small dogs and large dogs.

### Key Concepts: nSL and XP-nSL
- **nSL**: A within-population selection scan metric similar to iHS. However, instead of measuring haplotype decay in terms of genetic distance (which requires a genetic map), nSL measures distance by counting the number of segregating sites (segregating site count by length). This makes it highly robust to recombination rate variation and suitable for genomes without well-defined genetic maps.
- **XP-nSL**: A cross-population statistic that compares nSL profiles between a target population and a reference population. A high positive score indicates a selective sweep specific to the target population (longer haplotypes around the derived allele).
- **Phased Mode**: Since our input Dog10K BCF file is already phased (containing haplotype data formatted as `0|0`, `1|0`, etc.), we will perform a phased XP-nSL scan. This utilizes the precise haplotype sequences, which provides a significantly stronger selection signal compared to unphased analyses.

### Outgroup Allele Polarization
Haplotype selection scans require knowing which allele is **ancestral** (original) and which is **derived** (new mutant). `selscan` expects a VCF file where `0` is the ancestral allele and `1` is the derived allele.
To polarize our dataset, we use the **gray wolves** in the Dog10K metadata as an outgroup:
  - Gray wolves are the evolutionary ancestor of domestic dogs.
- For each SNP, the most common (major) allele in the gray wolf population is designated as the ancestral allele.
- If the ALT allele in the original VCF is the major allele in wolves, we must physically swap the REF/ALT alleles and swap genotypes (`0` becomes `1`, and `1` becomes `0`) for all individuals.

---

  # Definizione dei percorsi dei file VCF (Modifica i percorsi se sono diversi)
vcf_afr_path <- "Chr2_EDAR_LWK_500K.recode.vcf"
vcf_eas_path <- "Chr2_EDAR_CHS_500K.recode.vcf"  
# Popolazione Africana (AFR)

hap_afr <- data2haplohh(
  hap_file = vcf_afr_path,
  chr.name = 15,            # Cambia il numero se rs3827760 si trova su un altro CHR (es. 2)
  polarize_vcf = FALSE,     # Evita errori se manca l'info sull'allele ancestrale "AA"
)

# Popolazione Est-Asiatica (EAS)
hap_eas <- data2haplohh(
  hap_file = vcf_eas_path,
  chr.name = 15, 
  polarize_vcf = FALSE,
)

ehh_afr <- calc_ehh(hap_afr, mrk = target_marker)
ehh_eas <- calc_ehh(hap_eas, mrk = target_marker)

target_marker <- "rs3827760"

# Calcolo EHH pulito
ehh_afr <- calc_ehh(hap_afr, mrk = target_marker)
ehh_eas <- calc_ehh(hap_eas, mrk = target_marker)

# Generazione dei Grafici 2x2
par(mfrow = c(2, 2))

# 1. Grafico EHH Decay AFR
plot(ehh_afr, main = paste("EHH Decay AFR (LWK) -", target_marker), col = c("blue", "lightblue"), lwd = 2)

# 2. Grafico EHH Decay EAS
plot(ehh_eas, main = paste("EHH Decay EAS (CHS) -", target_marker), col = c("red", "pink"), lwd = 2)

# 3. Furcation Tree AFR
bif_afr <- calc_furcation(hap_afr, mrk = target_marker)
plot(bif_afr, main = paste("Furcation Tree AFR (LWK) -", target_marker), col = c("blue", "lightblue"))

# 4. Furcation Tree EAS
bif_eas <- calc_furcation(hap_eas, mrk = target_marker)
plot(bif_eas, main = paste("Furcation Tree EAS (CHS) -", target_marker), col = c("red", "pink"))

# Reset layout grafico
par(mfrow = c(1, 1))

library(rehh)
library(ggplot2)

target_marker <- "rs3827760"

# ==============================================================================
# 1. Genome-wide Scan con scan_hh()
# ==============================================================================
# Calcola l'omozigosi integrata per tutti i marker di ciascuna popolazione
scan_afr <- scan_hh(hap_afr)
scan_eas <- scan_hh(hap_eas)

# ==============================================================================
# 2. Calcolo dei punteggi iHS con ihh2ihs()
# ==============================================================================
# Calcola e standardizza i punteggi iHS per AFR ed EAS
ihs_afr <- ihh2ihs(scan_afr, freqbin = 0.025)
ihs_eas <- ihh2ihs(scan_eas, freqbin = 0.025)

# Estraiamo i data frame dei risultati standardizzati
ihs_afr_df <- ihs_afr$ihs
ihs_eas_df <- ihs_eas$ihs

# ==============================================================================
# 3. Controllo iHS a rs3827760 e Single-Site Plot in EAS
# ==============================================================================
# Stampa a console il valore di iHS per il marker target in EAS
target_ihs_eas <- ihs_eas_df[rownames(ihs_eas_df) == target_marker, ]
print("iHS score per rs3827760 in EAS:")
print(target_ihs_eas)

# Generazione del Single-Site Plot per EAS (Manhattan plot di iHS)
# rehh ha una funzione nativa per visualizzare i p-value di iHS: -log10(2 * (1 - Phi(|iHS|)))
plot(ihs_eas, main = paste("iHS Genome Scan in EAS - Target:", target_marker))


# ==============================================================================
# 4. Funzione Medie iHS in Sliding Windows (50 SNPs / step 40)
# ==============================================================================
calculate_window_ihs <- function(ihs_data, window_size = 50, step_size = 40) {
  # Rimuoviamo i valori NA per sicurezza
  clean_data <- na.omit(ihs_data[, c("POSITION", "IHS")])
  clean_data <- clean_data[order(clean_data$POSITION), ]
  
  n_snps <- nrow(clean_data)
  starts <- seq(1, n_snps - window_size + 1, by = step_size)
  
  window_positions <- numeric(length(starts))
  window_avg_abs_ihs <- numeric(length(starts))
  
  for (i in seq_along(starts)) {
    idx <- starts[i]:(starts[i] + window_size - 1)
    window_positions[i] <- mean(clean_data$POSITION[idx])
    # Calcolo della media del valore ASSOLUTO di iHS
    window_avg_abs_ihs[i] <- mean(abs(clean_data$IHS[idx]))
  }
  
  return(data.frame(Position_Mb = window_positions / 1e6, Avg_Abs_iHS = window_avg_abs_ihs))
}

# Calcolo delle finestre per EAS
eas_windows_ihs <- calculate_window_ihs(ihs_eas_df, window_size = 50, step_size = 40)

# Plot dei risultati iHS a finestre in EAS con ggplot2
ggplot(eas_windows_ihs, aes(x = Position_Mb, y = Avg_Abs_iHS)) +
  geom_line(color = "darkred", size = 1) +
  geom_point(color = "red", alpha = 0.5) +
  labs(
    title = "EAS Sliding Window Average |iHS| (50 SNPs / 40 Step)",
    x = "Position (Mb)",
    y = "Mean |iHS|"
  ) +
  theme_minimal()


# ==============================================================================
# 5. Calcolo di XP-EHH (EAS vs AFR) e relative sliding windows
# ==============================================================================
# Calcolo del cross-population EHH usando i risultati di scan_hh()
xpehh_res <- ies2xpehh(scan_pop1 = scan_eas, scan_pop2 = scan_afr, 
                       popname1 = "EAS", popname2 = "AFR")

# Estraiamo il data frame
xpehh_df <- xpehh_res$xpehh

# Creazione finestre mobili per XP-EHH (riutilizzando la stessa logica)
calculate_window_xpehh <- function(xpehh_data, window_size = 50, step_size = 40) {
  clean_data <- na.omit(xpehh_data[, c("POSITION", "XPEHH")])
  clean_data <- clean_data[order(clean_data$POSITION), ]
  
  n_snps <- nrow(clean_data)
  starts <- seq(1, n_snps - window_size + 1, by = step_size)
  
  window_positions <- numeric(length(starts))
  window_avg_xpehh <- numeric(length(starts))
  
  for (i in seq_along(starts)) {
    idx <- starts[i]:(starts[i] + window_size - 1)
    window_positions[i] <- mean(clean_data$POSITION[idx])
    # Per XP-EHH usiamo il valore reale (i picchi positivi indicano selezione in EAS)
    window_avg_xpehh[i] <- mean(clean_data$XPEHH[idx])
  }
  
  return(data.frame(Position_Mb = window_positions / 1e6, Avg_XPEHH = window_avg_xpehh))
}


# 1. Definizione della funzione per le finestre mobili di iHS
calculate_window_ihs <- function(ihs_data, window_size = 50, step_size = 40) {
  clean_data <- na.omit(ihs_data[, c("POSITION", "IHS")])
  clean_data <- clean_data[order(clean_data$POSITION), ]
  
  n_snps <- nrow(clean_data)
  starts <- seq(1, n_snps - window_size + 1, by = step_size)
  
  window_positions <- numeric(length(starts))
  window_avg_abs_ihs <- numeric(length(starts))
  
  for (i in seq_along(starts)) {
    idx <- starts[i]:(starts[i] + window_size - 1)
    window_positions[i] <- mean(clean_data$POSITION[idx])
    window_avg_abs_ihs[i] <- mean(abs(clean_data$IHS[idx]))
  }
  
  return(data.frame(Position_Mb = window_positions / 1e6, Avg_Abs_iHS = window_avg_abs_ihs))
}

# 2. Applicazione della funzione ai dati EAS (estratti in precedenza da ihs_eas$ihs)
eas_windows_ihs <- calculate_window_ihs(ihs_eas_df, window_size = 50, step_size = 40)

# 3. Plot nativo di R
plot(eas_windows_ihs$Position_Mb, eas_windows_ihs$Avg_Abs_iHS, 
     pch = 16,             # Punto pieno
     col = "red",          # Colore rosso per i punti
     xlab = "Position (Mb)", 
     ylab = "Mean |iHS|",
     main = "EAS Sliding Window Average |iHS| (50 SNPs / 40 Step)")
eas_afr_windows_xpehh <- calculate_window_xpehh(xpehh_df, window_size = 50, step_size = 40)
# Sostituisci il vecchio plot() con la funzione nativa corretta di rehh
manhattanplot(ihs_eas, main = paste("iHS Genome Scan in EAS - Target:", target_marker))


### Bash Code Exercise 4: Extracting and Polarizing Alleles
  ```bash
# 1. Extract wolf samples and combine with dog samples
bcftools query -l input/Part_2_CanidDiversity/AutoAndXPAR.Dog10K.phased.bcf | grep -E '^CLUP' > input/Part_2_CanidDiversity/wolves.txt
cat input/Part_2_CanidDiversity/subset_dogs.txt input/Part_2_CanidDiversity/wolves.txt > input/Part_2_CanidDiversity/dogs_and_wolves.txt

# 2. Extract polymorphic sites for both dogs and wolves
bcftools query -f '%CHROM\t%POS\n' input/Part_2_CanidDiversity/subset_chr15.bcf > input/Part_2_CanidDiversity/subset_chr15_positions.txt
bcftools view \
-S input/Part_2_CanidDiversity/dogs_and_wolves.txt \
-T input/Part_2_CanidDiversity/subset_chr15_positions.txt \
-O z \
-o input/Part_2_CanidDiversity/subset_chr15_with_wolves.vcf.gz \
input/Part_2_CanidDiversity/AutoAndXPAR.Dog10K.phased.bcf

# 3. Polarize alleles using Python (major allele in wolves = 0)
python3 scripts/polarize_by_wolves.py

# 4. Re-compress to block gzip format and index polarized VCF
bcftools view input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz -O z -o output/subset_chr15_polarized_bgzf.vcf.gz
mv output/subset_chr15_polarized_bgzf.vcf.gz input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz
bcftools index input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz

# 5. Extract polarized VCFs for small and large dogs separately
bcftools view -S input/Part_2_CanidDiversity/small_dogs.txt -O z -o input/Part_2_CanidDiversity/small_dogs_polarized.vcf.gz input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz
bcftools index input/Part_2_CanidDiversity/small_dogs_polarized.vcf.gz

bcftools view -S input/Part_2_CanidDiversity/large_dogs.txt -O z -o input/Part_2_CanidDiversity/large_dogs_polarized.vcf.gz input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz
bcftools index input/Part_2_CanidDiversity/large_dogs_polarized.vcf.gz
```

---
  
  ### Selection Scan Execution and Normalization
  We will execute the selection scan using the compiled `selscan` binary, comparing small dogs (target) against large dogs (reference) in phased mode, and normalize the raw scores.

### Bash Code Exercise 5: Running XP-nSL and Normalizing
Write the bash commands necessary to:
  1. Execute the `selscan` program in phased XP-nSL mode using the polarized target (small dogs) and reference (large dogs) VCFs. Use 4 threads and output results to `output/xpnsl_phased`.
2. Normalize the raw XP-nSL scores using the `selscan norm` utility.
3. Calculate window-based statistics (fraction of outliers) in 100 Kb non-overlapping windows.

**Write your Bash code here:**
  ```bash
# 

fst_nam_eas <- read_table("Chr2_NAM_EAS.weir.fst", col_names = FALSE, skip = 1)
fst_nam_eur <- read_table("Chr2_NAM_EUR.weir.fst", col_names = FALSE, skip = 1)
fst_eur_eas <- read_table("Chr2_EUR_EAS.weir.fst", col_names = FALSE, skip = 1)

# Assegniamo i nomi corretti alle colonne per poter lavorare
colnames(fst_nam_eas) <- c("CHROM", "POS", "FST_NAM_EAS")
colnames(fst_nam_eur) <- c("CHROM", "POS", "FST_NAM_EUR")
colnames(fst_eur_eas) <- c("CHROM", "POS", "FST_EUR_EAS")
 overlap_NAM <- fst_nam_eas[fst_nam_eas$POS %in% fst_nam_eur$POS, ]
 overlap_NAM_all <- overlap_NAM[overlap_NAM$POS %in% fst_eur_eas$POS, ]

 

```

---
  
  ### R Code Exercise 6: XP-nSL Haplotype Manhattan Plots
  Write the R code necessary to:
  1. Load the normalized SNP-level XP-nSL scores from `input/Part_2_CanidDiversity/xpnsl_phased.xpnsl.out.norm`.
2. Load the 100 Kb window-level data from `input/Part_2_CanidDiversity/xpnsl_phased.xpnsl.out.norm.100kb.windows`.
3. Highlight the *IGF1* region (between 41 Mb and 45.5 Mb).
4. Generate two Manhattan plots: one for the raw SNP-level scores (calculating $-\log_{10}(p\text{-value})$ from the normalized Z-scores), and another for the window-based fraction of extreme positive SNPs (`frac_top`).

**Write your R code here:**
  ```R
# 


```

### Questions for Students
1. **Why do we need to polarize alleles using an outgroup like the gray wolf? What does the `0` vs `1` coding represent in `selscan`?**
  * *Answer*:
  2. **Why does the raw SNP-level XP-nSL scan look like a noisy cloud at individual sites? What is the effect of calculating window-based scores (e.g. 100 Kb)?**
  * *Answer*:
  
  ---
  
  ## 6. Alternative Haplotype Selection Scan: Rsb using `rehh`
  
  As a complementary approach to XP-nSL, we will run **Rsb**, another widely used cross-population EHH-based statistic.

### Key Concepts & Comparison
- **Rsb**: Compares the integrated Extended Haplotype Homozygosity (iHH) between two populations. It is calculated as $\ln(iES_{pop1} / iES_{pop2})$, where $iES$ is the integrated EHH over physical distance (bp). A high positive score indicates selection in population 1 (small dogs), while a negative score indicates selection in population 2 (large dogs).
- **Difference from XP-nSL**: 
  - **XP-nSL** integrates the nSL metric over the number of segregating sites (SNP count). This makes it highly robust to local recombination rate variation.
- **Rsb** integrates EHH over physical distance (bp). In species with very strong selective sweeps and long-range Linkage Disequilibrium (like domestic dogs), Rsb can produce exceptionally high, clear peaks at sweep loci like *IGF1*.

### Phase Preservation during Outgroup Polarization
A common question in haplotype selection scans is: **Does swapping the REF and ALT alleles (and flipping 0 to 1 and 1 to 0) to align with the outgroup corrupt or destroy the phasing information?** The answer is **No**. In a phased VCF, genotypes are represented as `0|1` or `1|0` to denote which allele lies on which homologous chromosome (haplotype). Swapping REF and ALT alleles and swapping `0` and `1` (converting `0|1` to `1|0` and vice-versa) is a mathematically symmetric operation. It maintains the exact same physical haplotype alignment across all sites on each chromosome, merely updating the label of which allele is ancestral and which is derived. Thus, outgroup polarization is fully compatible with phased haplotype scans and no data is lost.

---
  
  ### R Code Exercise 7: Running Rsb in R using `rehh`
  Write the R code necessary to:
  1. Load the polarized Small Dogs and Large Dogs VCFs into `rehh` using `data2haplohh()`.
2. Compute the haplotype homozygosity scan for both populations using `scan_hh()`.
3. Compute the cross-population Rsb statistic using `ines2rsb()`.
4. Plot the Rsb Manhattan plot using `ggplot2`, highlighting the *IGF1* region between 41 Mb and 45.5 Mb in red.

**Write your R code here:**
  ```R
# 
input_dir <- "/disk/home/carfora/my_embo_popgen_2026/Tabita_Hunemeier/Practical_Session_Selection_EMBO/input/Part_2_CanidDiversity"

norm_files <- file.path(input_dir, "xpnsl_phased.xpnsl.out.norm")
xpnls_data <- read.table(norm_files, header = TRUE, sep ="\t")
xpnls_data <- na.omit(xpnls_data)
              
```

---
  
  ### Questions for Students
  1. **Explain the physical and mathematical difference between Rsb and XP-nSL. Why does Rsb show a much higher, less noisy peak at the *IGF1* locus in dogs compared to XP-nSL?**
  * *Answer*:
  2. **Does outgroup polarization of a phased VCF file corrupt or destroy the phasing information? Why or why not?**
  * *Answer*:
  