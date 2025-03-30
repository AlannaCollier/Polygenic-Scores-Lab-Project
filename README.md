# Polygenic Scores Lab Project

Analyses performed for 4th year lab research project "
Genetic Contributions to ADHD: Implications for Sleep and Circadian Rhythms in The 1000 Genomes Project"

## Objectives
The main aims of this project are to: 

- Compute polygenic scores (PGS) for ADHD, chronotype, and insomnia.

- Examine PGS distributions across diverse populations using the 1000 Genomes 
  Project.

- Assess the impact of population structure on PGS interpretation.


## Data Preparation

### Download and Convert 1000 Genomes Dataset

The 1000 Genomes dataset must be formatted for compatability SBayesRC in order to calculate PGS.

<pre> ```bash 
  # Download The 1000 Genomes data wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/*.vcf.gz 
  # Convert Target data from VCF to PLINK format
  plink --vcf ALL.chr*.vcf.gz --make-bed --out 1000G
  ``` </pre>

  ### Quality Control (QC)

  Quality control was performed to remove low-quality variants and ensure high-quality data for further analyses.

  ```bash
# Filter SNPs and individuals based on missingness and MAF
plink --bfile 1000G --geno 0.05 --mind 0.05 --maf 0.01 --make-bed --out 1000G_QC

#Perform LD pruning to reduce SNP correlation
plink --bfile 1000G_QC --indep-pairwise 50 5 0.2 --out 1000G_pruned
plink --bfile 1000G_QC --extract 1000G_pruned.prune.in --make-bed --out 1000G_LDpruned
```

### Principal Component Analysis (PCA)

Principal Components were extracted to capture population structure, which has an influence on interpretation of PGS

```bash
plink --bfile 1000G_LDpruned --pca 10 --out 1000G_PCA
```

## Convert GWAS Summary Statistics to COJO format

This is done to ensure SNP IDs, alleles and effect sizes match The 1000 Genomes Project dataset. The example below is for chronotype summary statistics, but similar editing was performed for all GWAS summary statistics.

```R
library(data.table)

# load the data into the environment in a table format
chronotype_data <-fread("chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt")

library(dplyr)
library(magrittr)


#Reformat the data to create a table in cojo format
chronotype_cojo <- chronotype_data %>%
  mutate (
    b = BETA, #rename BETA to b
    freq = A1FREQ, #average fequency of A1 (effect allele) across cases and controls
    N = 449734, #sample size (cases and controls)
    p = P_BOLT_LMM , #rename p value variable
    A1 = ALLELE1, #rename allele 1 variable
    A2 = ALLELE0, #rename allele 2 variable
    se = SE #change se from lowercase to uppercase
  ) %>%
  select(SNP, A1, A2, freq, b, SE, p,N) #selects appropriate columns to be used in SBayesRC

#save COJO table as a .txt file that is in tab-delimited format for running in SBayesRC
write.table(adhd_cojo, "adhd_cojo_tab.txt", sep = "\t", row.names = FALSE, quote = FALSE ) 
```

### Filter GWAS to 1000 Genomes SNPs

This was done to ensure that only SNPs preset in the 1000 Genomes data were used

```bash
plink --bfile 1000G_QC --extract <GWAS_SNP_list.txt> --make-bed --out GWAS_SNPs_filtered
```

## Running SBayesRC

The code below is an example of how SBayesRC was run using the chronotype data in COJO format in order to produce weighted effects estimates. This was carried out in the same way for both insomnia and chronotype.

``` bash
# Tidy
Rscript -e "SBayesRC::tidy(mafile=‘chronotype, LDdir='../ukbEUR_Imputed', output=‘chronotype_tidy.ma', log2file=TRUE)"

# Impute
Rscript -e "SBayesRC::impute(mafile=‘chronotype_tidy.ma', LDdir='../ukbEUR_Imputed', output=‘chronotype_imp.ma', log2file=TRUE)"

# Run model

Rscript -e "SBayesRC::sbayesrc(mafile=‘chronotype_imp.ma', LDdir='../ukbEUR_Imputed', outPrefix=‘chronotype_tidy_sbrc', annot='../annot_baseline2.2.txt', log2file=TRUE)"
```

## Computing PGS in 1000 Genomes 

Individual PGS were calculated based on GWAS effect sizes, this was carried out for each trait.

```bash
plink --bfile 1000G_QC \
--score Chronotype_PGS_weights.txt 1 2 3 header \
--out 1000G_Chronotype_PGS
```

## Performance metrics and Data Visualisation

To analyse PGS distributions across populations, the PGS scores were combined with population labels and PCA components to make a single dataset for each trait. Below is an example of how this was done for ADHD

```R
#load in adhd data and rename columns
adhd_pgs_data <- read.table("adhd_1kg_PRS.sscore", header=FALSE)

head(adhd_pgs_data)

colnames(adhd_pgs_data)[1] <- "FID"
colnames(adhd_pgs_data)[2] <- "IID"
colnames(adhd_pgs_data)[3] <- "ALLELE_CT"
colnames(adhd_pgs_data)[4] <- "NAMED_ALLELE_DOSAGE_SUM"
colnames(adhd_pgs_data)[5] <- "SCORE1_AVG"

#load in population label data
pop_data <- read.table("1000G_Merged_Population_Corrected.txt", header=TRUE)

#load in pca data and rename columns 
pca_data <- read.table("1kg_pca.eigenvec", header=FALSE)

colnames(pca_data) <- c("FID", "IID", paste0("PC", 1:(ncol(pca_data)-2)))

pca_data <- pca_data[, -1]

#merge  ADHD PRS, PCA and population labels

merged_data <- adhd_pgs_data %>%
  inner_join(pop_data, by="IID") %>%
  inner_join(pca_data, by="IID") %>%
  drop_na
```

The merged datasets for each trait were then combined together into one large dataset 

```R
#combine all merged datasets using the mutate function 
combined_data <- bind_rows(
  chronotype_merged %>% mutate(Trait = "Chronotype"),
  insomnia_merged %>% mutate(Trait = "Insomnia"),
  merged_data %>% mutate(Trait = "ADHD")
)
```

### PCA Plot 

A PCA scatterplot was generated to confirm the genetic diversity within The 1000 Genomes Project dataset.

```R
ggplot(insomnia_merged, aes(x=PC1, y=PC2, color=super_pop)) +
  geom_point(size=3, alpha=0.7) +
  labs(title="PCA of 1000 Genomes with Insomnia PGS", x="PC1", y="PC2") +
  theme_minimal()
```

### Boxplots of PGS distribution for each trait 

To visualise how PGS distributions vary across traits (ADHD, insomnia and chronotype).

```R
ggplot(combined_data, aes(x=Trait, y=SCORE1_AVG, fill=Trait)) +
  geom_boxplot() +
  labs(title="PGS Across Phenotype", x="Phenotype", y="Polygenic Score") +
  theme_minimal()
```

### Boxplots of PGS grouped by ancestry 

To viasualise how PGS distributions vary across traits that are grouped by superpopulation.

```R
#side by side boxplots grouped by ancestry (population)
ggplot(combined_data, aes(x=super_pop, y=SCORE1_AVG, fill=Trait)) +
  geom_boxplot() +
  labs(title="PGS Across Populations", x="Population", y="Polygenic Score") +
  theme_minimal()
```

### Side by side scatterplots of PGS distribution across superpopulations with overlying boxplots

Side by side scatterplots were created to reinforce and confirm the findings from the boxplots using a different visualisation type.

```R
# Determine global min and max PGS values across all datasets
y_min <- min(merged_data$SCORE1_AVG, insomnia_merged$SCORE1_AVG, chronotype_merged$SCORE1_AVG, na.rm = TRUE)
y_max <- max(merged_data$SCORE1_AVG, insomnia_merged$SCORE1_AVG, chronotype_merged$SCORE1_AVG, na.rm = TRUE)

# Define function to create scatter plots with fixed y-axis range
create_scatter_plot <- function(data, trait_name) {
  ggplot(data, aes(x = super_pop, y = SCORE1_AVG, color = super_pop)) +
    geom_jitter(width = 0.2, alpha = 0.6) +  # Jitter to avoid overlap
    geom_boxplot(alpha = 0.2, outlier.shape = NA) +  # Boxplot for distribution
    stat_smooth(method = "lm", se = FALSE, aes(group = 1), color = "black") +  # Regression line
    theme_minimal() +
    labs(title = paste( trait_name, "PGS"),
         x = "Superpopulation",
         y = "Polygenic Score (PGS)") +
    ylim(y_min, y_max) +  # Set consistent y-axis range
    theme(legend.position = "none")  # Hide legend since it's redundant
}

# Generate plots for each trait
plot_adhd <- create_scatter_plot(merged_data, "ADHD")
plot_insomnia <- create_scatter_plot(insomnia_merged, "Insomnia")
plot_chronotype <- create_scatter_plot(chronotype_merged, "Chronotype")

# Arrange plots side by side
(plot_adhd | plot_insomnia | plot_chronotype)

```

### Linear Regressions

Individual linear regressions were carried out for each trait, the results were then summarised into a table. 2 types of linear regressions were carried out on each trait and this resulted in 2 summary tables.

An example of the linear regressions carried out insomnia can be seen below, similar code was used for each trait.

```R
#linear regression to assess population structure bias in insomnia 

insomnia_lm <- lm(SCORE1_AVG ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=insomnia_merged)

summary(insomnia_lm)


insomnia_lm_pop <- lm(SCORE1_AVG ~ super_pop + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=insomnia_merged)

summary(insomnia_lm_pop)
```

### Scree plot 

A scree plot was produced to see which PCs captured the largest proportion of variance in the data.

```R
# Extract eigenvalues (variance explained by each PC)
eigenvalues <- pca_result$sdev^2

# Create a dataframe for plotting
scree_df <- data.frame(
  PC = paste0("PC", 1:length(eigenvalues)),
  Variance_Explained = eigenvalues / sum(eigenvalues)  # Proportion of variance
)

# Print eigenvalues and variance explained for each PC
print(scree_df)



plot(eigenvalues, type = "b",
     xlab = "Principal Component",
     ylab = "Eigenvalue")
```
