# Script header ----
# Programmed by: Kairi Tanaka
# Programmed on: 2023-04-04
# Programmed to: Create a template to run dada2 pipeline
# Last modified by: 
# Last modified on: 
# Last modified to: 

# Install packages ----
# Uncomment below if not installed yet
# bioc_packs <- c(
#   "dada2",
#   "phyloseq"
# )

# Load packages ----
packages <- c(
  "tidyverse",
  "dada2"
)

lapply(
  packages,
  library,
  character.only = TRUE
)

# Switch working directory to the project folder ----
setwd('./sewage_health_equity/') # CHANGE this to the project directory

# Run DADA2 pipeline ----
## Specify path ----
path <- "intermediate_data/16s_clean_fastq/" # CHANGE this to the directory containing the fastq files
list.files(path)

## Read in the names of fastq files ----
# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq.gz and SAMPLENAME_2.fastq.gz
fnFs <- 
  sort(
    list.files(
      path,
      # pattern = "_L001_R1_001.fastq.gz", # CHANGE to the appropriate filename pattern
      pattern = "_clean_1.fq.gz", # CHANGE to the appropriate filename pattern
      full.names = TRUE
    )
  )
fnFs

fnRs <- 
  sort(
    list.files(
      path, 
      # pattern = "_L001_R2_001.fastq.gz", # CHANGE to the appropriate filename pattern
      pattern = "_clean_2.fq.gz", # CHANGE to the appropriate filename pattern
      full.names = TRUE
    )
  )
fnRs

# Extract sample names
# ADJUST below to obtain unique sample names based on filenames
sample.names <- 
  sapply(
    strsplit(basename(fnFs), "_"), `[`, 3
  )
sample.names

## Plot quality profiles ----
# Create a list to loop through
qual_vec <- 
  list(
    c(1:9), c(10:18), c(19:27), c(28:36), c(37:43)
  )

# Plot the quality scores of the sequencing read files via for loop
for (i in 1:length(qual_vec)) {
  # Forward read quality pre QC
  jpeg(
    file = paste0("intermediate_data/dada2/F_qual_",i,".jpg"),
    units = "in", width = 8, height = 8, res = 300
  )
  plotQualityProfile(fnFs[qual_vec[[i]]])
  dev.off()
  
  # Reverse read quality pre QC
  jpeg(
    file = paste0("intermediate_data/dada2/R_qual_",i,".jpg"),
    units = "in", width = 8, height = 8, res = 300
  )
  plotQualityProfile(fnRs[qual_vec[[i]]])
  dev.off()
}

# Place filtered files in filtered/ subdirectory
filtFs <- 
  file.path(
    # path, 
    "./intermediate_data/clean_filtered", 
    paste0(
      sample.names, 
      "_filt_1.fq.gz"
    )
  )

filtRs <- 
  file.path(
    # path, 
    "./intermediate_data/clean_filtered", 
    paste0(
      sample.names, 
      "_filt_2.fq.gz"
    )
  )

# Assign sample names
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## Filter and trim ----
# CHANGE BELOW with the appropriate primer sequences
# Or comment out if QC program (e.g. fastp) was run beforehand
# 16S rRNA gene V3-V4 region 
# 341F and 785R
# n_fwd <- nchar("CCTACGGGNGGCWGCAG") 
# n_fwd
# 
# n_rev <- nchar("GACTACVHGGGTATCTAATCC")
# n_rev

out <- 
  filterAndTrim(
    fnFs, 
    filtFs, 
    fnRs, 
    filtRs, 
    # trimLeft = c(
    #   # n_fwd+3, # CHANGE to the length of fwd primer used
    #   n_fwd, # CHANGE to the length of fwd primer used
    #   n_rev # CHANGE to the length of rev primer used
    # ),
    truncLen = c(
      240, # CHANGE based on the quality of fwd reads
      200 # CHANGE based on the quality of rev reads
    ),
    maxN = 0,
    maxEE = c(2, 2), 
    truncQ = 2, 
    rm.phix = TRUE,
    compress = TRUE, 
    multithread = TRUE
  )
head(out)
# View(out)

# Save output as rds
# saveRDS(
#   out,
#   file = "data/out.rds"
# )

## Plot quality profiles of filtered reads ----
# Plot the quality scores of the sequencing read files via for loop
for (i in 1:length(qual_vec)) {
  # Forward read quality pre QC
  jpeg(
    file = paste0("intermediate_data/dada2/filtF_qual_",i,".jpg"),
    units = "in", width = 8, height = 8, res = 300
  )
  plotQualityProfile(filtFs[qual_vec[[i]]])
  dev.off()
  
  # Reverse read quality pre QC
  jpeg(
    file = paste0("intermediate_data/dada2/filtR_qual_",i,".jpg"),
    units = "in", width = 8, height = 8, res = 300
  )
  plotQualityProfile(filtRs[qual_vec[[i]]])
  dev.off()
}

## Learn error rates ----
# Forward reads
errF <- 
  learnErrors(
    filtFs,
    multithread = TRUE
  )

# Reverse reads
errR <- 
  learnErrors(
    filtRs,
    multithread = TRUE
  )

# Save as rds
saveRDS(
  errF,
  file = "intermediate_data/dada2/errF.rds"
)

saveRDS(
  errR,
  file = "intermediate_data/dada2/errR.rds"
)

# Visualize the error rates
jpeg(file = "intermediate_data/dada2/F_error_plot.jpg",
     units = "in", width = 9, height = 9, res = 300)
plotErrors(
  errF, 
  nominalQ = TRUE
)
dev.off()

jpeg(file = "intermediate_data/dada2/R_error_plot.jpg",
     units = "in", width = 9, height = 9, res = 300)
plotErrors(
  errR, 
  nominalQ = TRUE
)
dev.off()

## Infer ASVs ----
# Sample inference for forward reads
dadaFs <- 
  dada(
    filtFs, 
    err = errF, 
    multithread = TRUE
  )

# Sample inference for reverse reads
dadaRs <- 
  dada(
    filtRs, 
    err = errR,
    multithread = TRUE
  )

# Save as rds
saveRDS(
  dadaFs,
  file = "intermediate_data/dada2/dadaFs.rds"
)

saveRDS(
  dadaRs,
  file = "intermediate_data/dada2/dadaRs.rds"
)

## Merge reads ----
mergers <- 
  mergePairs(
    dadaFs,
    filtFs,             
    dadaRs,
    filtRs, 
    verbose = TRUE
  )
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Save data
saveRDS(
  mergers,
  file = "intermediate_data/dada2/mergers.rds"
)

## Construct sequence table ----
# Create a sequence table
seqtab <- 
  makeSequenceTable(mergers) 

# Save as rds
saveRDS(
  seqtab,
  file = "intermediate_data/dada2/seqtab.rds"
)

# Check dimension
dim(seqtab)

## Inspect distribution of sequence lengths 
table(nchar(getSequences(seqtab)))

## Remove chimeras ----
seqtab.nochim <- 
  removeBimeraDenovo(
    seqtab,
    method = "consensus",
    multithread = TRUE,
    verbose = TRUE
  )

# Check dimension
dim(seqtab.nochim)
# Chimeric sequences make up about ??% of the merged sequence variants

## Find out frequency of chimeric sequences 
sum(seqtab.nochim/sum(seqtab))
# ??% was chimeric sequences 

# Save the data
saveRDS(
  seqtab.nochim,
  file = "intermediate_data/dada2/seqtab_nochim.rds"
)

## Track reads ----
getN <- function(x) sum(getUniques(x))

track <- 
  cbind(
    out,
    sapply(
      dadaFs,
      getN
    ),
    rowSums(seqtab.nochim)
  )
# If processing a single sample, remove the sapply calls:
# e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c(
  "input",
  "filtered",
  "denoisedF",
  "nonchim"
)

rownames(track) <- sample.names

head(track)
# View(track)

# Save data
saveRDS(
  track,
  file = 'intermediate_data/dada2/track.rds'
)

# Save n of reads tracking table as an xlsx file
write.xlsx(
  track, 
  file = 'processed_data/summary_tables//dada2_read_tracking.xlsx'
)

## Assign taxonomy ----
# Silva 138.2 database
taxa <- 
  assignTaxonomy(
    seqtab.nochim,
    "databases/dada2/silva_nr99_v138.2_toGenus_trainset.fa.gz", # Silva 138.2
    # Alternatively, GreenGenes2 (comment out above and comment in below)
    # "databases/dada2/gg2_2024_09_toGenus_trainset.fa.gz", 
    tryRC = TRUE,
    multithread = TRUE,
    verbose = TRUE
  )

# Save as rds
saveRDS(
  taxa,
  file = "intermediate_data/dada2/taxa.rds"
)

# Inspect the taxonomic assignments
# Copy the dataset for display purpose only
taxa.print <- taxa

# Removing sequence rownames for display only
rownames(taxa.print) <- NULL

# Show the first few
head(taxa.print)
View(taxa.print)

# Save the outputs ----
# seqtab.nochim_untrim <- seqtab.nochim
# taxa_untrim <- taxa
# track_untrim <- track

# seqtab.nochim_cutadapt <- seqtab.nochim
# taxa_cutadapt <- taxa
# track_cutadapt <- track


# Create metadata ----
metadata <- 
  data.frame(
    SampleNames = sample.names
  ) |> 
  mutate(
    Chamber = str_split_i(string = SampleNames, pattern = "-", i = 1),
    Medium = str_split_i(string = SampleNames, pattern = "-", i = 2),
    Location = str_split_i(string = SampleNames, pattern = "-", i = 3),
    SampleID = str_split_i(string = SampleNames, pattern = "-", i = 4),
    Wash = str_split_i(string = SampleNames, pattern = "-", i = 5),
    Rownames = SampleNames
  ) |> 
  column_to_rownames(var = "Rownames")
metadata

# saveRDS(
#   metadata,
#   file = 'data/metadata.rds'
# )

# Create a dataframe containing the information above ----
df <- 
  metadata |> 
  rownames_to_column(var = "SampleName") |> 
  merge(
    seqtab.nochim |> 
      data.frame() |> 
      rownames_to_column(var = "SampleName"),
    by = "SampleName"
  ) |> 
  pivot_longer(
    cols = where(is.integer),
    names_to = "Sequence",
    values_to = "Abundance"
  ) |> 
  merge(
    taxa |> 
      data.frame() |> 
      rownames_to_column(var = "Sequence"),
    by = "Sequence"
  ) |> 
  mutate(
    GenusSpecies = paste0("*",Genus," ",Species,"*"),
    Chamber = factor(
      Chamber,
      levels = c("A", "D", "B", "C"),
      labels = c("A", "A1", "B", "C")
    ),
    O2_ppm = factor(
      Chamber,
      levels = c("A", "A1", "B", "C"),
      labels = c("50000", "10000", "10-50", "0")
    ),
    Wash = factor(
      Wash,
      levels = c("W0", "W1"),
      labels = c("Undiluted", "Diluted")
    )
  ) |> 
  # filter out A1 chamber for presentation
  filter(Chamber != "A1") |> 
  # calculate the total abundance (# of reads) per sample
  mutate(
    TotalAbund = sum(Abundance, na.rm = TRUE),
    .by = c("SampleName"),
    .keep = "all"
  ) |> 
  # calculate relative abundance per sample
  mutate(
    RelAbund = Abundance / TotalAbund
  ) |> 
  # calculate prevalence per sample
  mutate(
    Prevalence = Abundance/Abundance
  ) |> 
  # calculate prevalence per chamber and medium type, across samples
  mutate(
    PerChamPrev = sum(Prevalence, na.rm = TRUE),
    .by = c("Sequence", "Chamber", "Medium")
  )
df |> View()

# saveRDS(
#   df,
#   file = "data/SewageProj/amplicon_data/df.rds"
# )

df <- readRDS(
  file = "data/SewageProj/amplicon_data/df.rds"
)
# df |> View()


## Take a look at AW02 only ----
df |> 
  filter(SampleID == 'AW02') |>
  select(Sequence, Kingdom:Genus, SampleName, RelAbund) |>
  pivot_wider(
    names_from = SampleName,
    values_from = RelAbund
  ) |> 
  write_csv(
    file = 'data/SewageProj/bmk_wgs/AW02_16s_taxa_table.csv'
  )

saveRDS(
  df |> filter(SampleID == 'AW02'),
  file = 'data/AW02_amplicon_tax_table.rds'
)

### Subset for analysis ----
df <- 
  df |> 
  filter(
    Wash == "Diluted"
  )

## Looking at archaea ----
df |> 
  filter(Kingdom =="Archaea") |> 
  ggplot(
    mapping = aes(
      x = SampleID,
      y = RelAbund*100,
      fill = Chamber
    )
  ) +
  geom_col(
    mapping = aes(
      fill = Chamber
    ),
    color = "black",
    position = position_dodge(
      width = 0.75
    )
  ) +
  facet_grid(
    rows = vars(Medium, Wash),
    cols = vars(GenusSpecies)
  ) +
  labs(
    title = "Relative abundance of archaeal amplicon sequence variants",
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(
      angle = 30, vjust = .75
    )
  )

# Save plot
# ggsave(
#   filename = "figures/methanogen_asv_relabund.png",
#   units = "in", width = 8, height = 6, dpi = 300
# )

## Looking at bacteria ----
# Set a grouping parameter
group_level <- "Phylum"

n_group <- df$Phylum |> unique() |> length()

# Expand color palette
library(RColorBrewer)
my_col <- colorRampPalette(brewer.pal(8, "Dark2"))

df |> 
  summarize(
    GroupRelAbund = sum(RelAbund, na.rm = TRUE),
    .by = c(
      "SampleID", "Chamber", "Medium", "Wash", "Phylum"
    )
  ) |> 
  mutate(
    Phylum = if_else(
      GroupRelAbund <= 0.01,
      "≤ 1%",
      Phylum
    )
  ) |> 
  ggplot(
    mapping = aes(
      x = Chamber,
      y = GroupRelAbund
    )
  ) +
  geom_col(
    mapping = aes(
      fill = reorder(Phylum, -GroupRelAbund),
      group = Chamber
    ),
    color = "black",
    position = position_fill(reverse = TRUE)
  ) +
  facet_grid(
    cols = vars(SampleID),
    rows = vars(Medium),
    switch = "x"
  ) +
  labs(
    # title = "Stacked bar chart showing the relative abundance of\nbacterial and archaeal amplicon sequence variants"
    x = "Sample ID & Chamber",
    y = "Relative abundance (%)",
    fill = "Phylum"
  ) +
  scale_y_continuous(labels = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0,0))) +
  scale_fill_manual(values = my_col(n_group)) + 
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.ticks = element_blank()
    # axis.text.x = element_text(
    #   angle = 30, vjust = .75
    # )
  )

# Save plot
# ggsave(
#   filename = "figures/asv_relabund.png",
#   units = "in", width = 8, height = 5, dpi = 1200
# )

# Load data ----
seqtab.nochim <- readRDS(file = 'data/seqtab_nochim.rds')
metadata <- readRDS(file = 'data/metadata.rds')
taxa <- readRDS(file = 'data/taxa.rds')
fitGTR <- readRDS(file = 'data/fitGTR.rds')

# Root the phylogenetic tree ----
rooted_tree <- phangorn::midpoint(fitGTR$tree)
ape::is.rooted(rooted_tree) # Output should be TRUE

# Create a phyloseq object ----
ps <- 
  phyloseq(
    otu_table(
      seqtab.nochim,
      taxa_are_rows = FALSE 
    ),
    sample_data(metadata),
    tax_table(taxa),
    phy_tree(rooted_tree)
  )
ps

# Save as rds
# saveRDS(
#   ps,
#   file = "data/SewageProj/amplicon_data/ps.rds"
# )

ps <- 
  readRDS(
    file = "data/SewageProj/amplicon_data/ps.rds"
  )
ps

# Remove Chamber D samples for presentation
ps <- ps |> 
  phyloseq::subset_samples(
    Chamber != "D" & Wash == "W1"
  )
ps |> sample_data()

# Decontaminate ----
# Load decontam
# BiocManager::install('decontam')
# install.packages('decontam')
library(decontam)

# Load the DNA concentration data ----
conc <- 
  read_csv(
    file = "data/conc_table.csv"
  ) |> 
  filter(
    SampleName %in% c(df$SampleName |> unique())
  )
# conc |> View()

# Run decontam
contam_df <- 
  decontam::isContaminant(
    seqtab = ps,
    conc = conc$Concentration_ng_uL,
    method = "frequency"
  )
# contam_df

ps_dec <- prune_taxa(
  !contam_df$contaminant,
  ps
)
ps_dec
# Now with 3838 taxa (3849 previously)

# save as rds
# saveRDS(
#   ps_dec,
#   file = "data/SewageProj/amplicon_data/ps_dec.rds"
# )

# load rds
ps_dec <-
  readRDS(
    file = 'data/SewageProj/amplicon_data/ps_dec.rds'
  )

## Update df ----
df <- 
  metadata |> 
  rownames_to_column(var = "SampleName") |> 
  merge(
    otu_table(ps_dec) |> 
      data.frame() |> 
      rownames_to_column(var = "SampleName"),
    by = "SampleName"
  ) |> 
  pivot_longer(
    cols = where(is.integer),
    names_to = "Sequence",
    values_to = "Abundance"
  ) |> 
  merge(
    taxa |> 
      data.frame() |> 
      rownames_to_column(var = "Sequence"),
    by = "Sequence"
  ) |> 
  mutate(
    GenusSpecies = paste0("*",Genus," ",Species,"*"),
    Chamber = factor(
      Chamber,
      levels = c("A", "D", "B", "C"),
      labels = c("A", "A1", "B", "C")
    ),
    O2_ppm = factor(
      Chamber,
      levels = c("A", "A1", "B", "C"),
      labels = c("50000", "10000", "10-50", "0")
    ),
    Wash = factor(
      Wash,
      levels = c("W0", "W1"),
      labels = c("Undiluted", "Diluted")
    )
  ) |> 
  # filter out A1 chamber for presentation
  filter(Chamber != "A1") |> 
  # calculate the total abundance (# of reads) per sample
  mutate(
    TotalAbund = sum(Abundance, na.rm = TRUE),
    .by = c("SampleName"),
    .keep = "all"
  ) |> 
  # calculate relative abundance per sample
  mutate(
    RelAbund = Abundance / TotalAbund
  ) |> 
  # calculate prevalence per sample
  mutate(
    Prevalence = Abundance/Abundance
  ) |> 
  # calculate prevalence per chamber and medium type, across samples
  mutate(
    PerChamPrev = sum(Prevalence, na.rm = TRUE),
    .by = c("Sequence", "Chamber", "Medium")
  )

# Create a rarefaction curve ----
library(vegan)
# jpeg(
#   filename = "figures/16s_rarecurve.jpg",
#   units = "in", width = 6, height = 6, res = 300
# )
rarecurve(
  otu_table(ps_dec) |> matrix() |> t(),
  step = 10,
  label = FALSE
)
# dev.off()

# Calculate Good's coverage ----
df_goods <- 
  df |> 
  summarize(
    n_seq = sum(Abundance),
    n_sing = sum(Abundance == 1),
    goods = 100* (1 - n_sing / n_seq),
    .by = c("Chamber", "SampleID", "Medium", "Wash")
  )

## Visualize Good's coverage ----
df_goods |> 
  ggplot(
    mapping = aes(
      x = n_seq,
      y = goods
    )
  ) +
  geom_point(
    mapping = aes(
      color = Chamber
    ),
    size = 5
  ) +
  labs(
    y = "Good's Coverage (%)",
    x = "Number of ASVs"
  ) +
  theme_classic()
# Probably no need to cut off samples

# Alpha diversity ----
## Calculate Shannon diversity ----
shannon_df <- 
  df |> 
  summarize(
    PerSpecies = RelAbund * log(RelAbund),
    .by = c(
      "Sequence", "Chamber", "SampleID", "Medium", "Wash"
    )
  ) |> 
  summarize(
    Shannon = -sum(PerSpecies, na.rm = TRUE),
    .by = c(
      "Chamber", "SampleID", "Medium", "Wash"
    )
  )

## Calculate Simpson diversity ----
simpson_df <- 
  df |> 
  summarize(
    n = Abundance * (Abundance - 1),
    .by = c(
      "Sequence", "Chamber", "SampleID", "Medium", "Wash", "TotalAbund"
    )
  ) |> 
  summarize(
    n_sum = sum(n, na.rm = TRUE),
    .by = c(
      "Chamber", "SampleID", "Medium", "Wash", "TotalAbund"
    ),
  ) |> 
  mutate(
    PreSimpson = n_sum / (TotalAbund * (TotalAbund - 1)),
    Simpson = 1 - PreSimpson,
    `Inverse Simpson` = 1 / PreSimpson,
    .by = c("Chamber", "SampleID", "Medium", "Wash")
  ) |> 
  select(-TotalAbund, -n_sum, -PreSimpson)
simpson_df

## Calculate abundance-based coverage estimator (ACE) ----
ace_df <- 
  df |> 
  summarize(
    Observed = sum(Prevalence, na.rm = TRUE),
    n_sing = sum(Abundance == 1),
    n_doub = sum(Abundance == 2),
    ACE = Observed + (n_sing^2) / (2 * (n_doub-1)),
    .by = c("Chamber", "SampleID", "Medium", "Wash")
  ) |> 
  select(-n_sing, -n_doub)
ace_df

## Merge alpha diveristy results
alpha_div_df <- 
  shannon_df |> 
  merge(
    simpson_df,
    by = c(
      "Chamber", "SampleID", "Medium", "Wash"
    )
  ) |> 
  merge(
    ace_df,
    by = c(
      "Chamber", "SampleID", "Medium", "Wash"
    )
  )
alpha_div_df

# Calculate alpha diversity using phyloseq ----
alpha_div <- 
  estimate_richness(
    ps_dec,
    split = TRUE,
    measures = c(
      "Observed",
      "Shannon"
    )
  )
alpha_div |> View()

## Statistically analyze ----
### Check for normality ----
sapply(
  alpha_div,
  function(x){shapiro.test(x)}
)
# p = 0.054

### Check for homogeneity ----
sapply(
  alpha_div,
  function(x){
    bartlett.test(x ~ metadata$Chamber[metadata$Chamber != "D" & metadata$Wash != "W0"])
  }
)
# p = 0.645

### Observed ----
obs_mm <- 
  RRPP::lm.rrpp(
    # alpha_div$Observed ~ Chamber * Medium + SampleID*Wash,
    alpha_div$Observed ~ Chamber * Medium + SampleID,
    # alpha_div$Observed ~ Chamber * Medium,
    data = metadata |> filter(Chamber != "D" & Wash != "W0"),
    SS.type = "III",
    iter = 9999,
    Parallel = TRUE
  )

# PERMANOVA 
obs_anova <- 
  anova(
    obs_mm,
    effect.type = "F",
    error = c(
      "Medium",
      "Chamber",
      "Residuals",
      "Residuals"
    )
  )
obs_anova
# SampleID is significantly different

# Null model for pairwise comparison
obs_null <- 
  RRPP::lm.rrpp(
    # alpha_div$Observed ~ SampleID,
    alpha_div$Observed ~ SampleID,
    data = metadata |> filter(Chamber != "D" & Wash != "W0"),
    SS.type = "III",
    iter = 9999,
    Parallel = TRUE
  )

# Pairwise comparison
obs_pw <- RRPP::pairwise(
  fit = obs_mm,
  fit.null = obs_null,
  groups = metadata |> filter(Chamber != "D" & Wash != "W0") |> _$Chamber
) |> 
  summary()
obs_pw

p.adjust(
  obs_pw$summary.table$`Pr > d`,
  method = "fdr"
)
# A and C are significantly different
# rest tended to differ

### Shannon ----
shannon_mm <- 
  RRPP::lm.rrpp(
    alpha_div$Shannon ~ Chamber * Medium + SampleID,
    data = metadata |> filter(Chamber != "D" & Wash != "W0"),
    SS.type = "III",
    iter = 9999,
    Parallel = TRUE
  )

# PERMANOVA 
shannon_anova <- 
  anova(
    shannon_mm,
    effect.type = "F",
    error = c(
      "Medium",
      "Chamber",
      "Residuals",
      "Residuals"
    )
  )
shannon_anova
# Chamber tended to differ
# SampleID is significant

# Null model for pairwise comparison
shan_null <- 
  RRPP::lm.rrpp(
    alpha_div$Shannon ~ SampleID,
    data = metadata |> filter(Chamber != "D" & Wash != "W0"),
    SS.type = "III",
    iter = 9999,
    Parallel = TRUE
  )

# Pairwise comparision
shannon_pw <- RRPP::pairwise(
  fit = shannon_mm,
  fit.null = shan_null,
  groups = metadata |> filter(Chamber != "D" & Wash != "W0") |> _$Chamber
) |> 
  summary()
shannon_pw

p.adjust(
  shannon_pw$summary.table$`Pr > d`,
  method = "fdr"
)
# A is different from both B and C
# B and C are not different

## Visualize alpha diversity ----
### Set parameters ----
anno_colors <- 
  list(
    # "Sample" = c("AW03"="#FFC20A", "AO02"="#0C7BDC"),
    # "Chamber" = c("A5"="grey88", "A1"="grey66", "B"="grey25","C"="grey11"),
    "Chamber" = c("All"="white", "A"="grey88", "B"="grey55","C"="grey11"),
    "Medium" = c("Y1"="grey55","YC"="grey11"),
    "Wash" = c("Raw" = "#E66100", "Diluted" = "#5D3A9B"),
    "O2_ppm" = c("50000"="grey88", "10000"="grey66", "10-50"="grey25","0"="grey11")
    # "O2_ppm" = c("50000"="grey99", "10-50"="grey55","0"="grey0")
  )

### Visualize ----
alpha_gg <- 
  alpha_div_df |>
  # mutate(
  #   O2_ppm = factor(
  #     Chamber,
  #     levels = c("A5", "A1", "B", "C"),
  #     labels = c("50000", "10000", "10-50", "0")
  #   )
  # ) |> 
  pivot_longer(
    cols = c(
      "Observed", "ACE", 
      "Shannon", 
      "Simpson", "Inverse Simpson"
    ),
    names_to = "Metric",
    values_to = "Estimate"
  ) |>
  filter(Metric %in% c("Observed", "Shannon")) |> 
  ggplot(
    mapping = aes(
      x = Chamber,
      y = Estimate,
      fill = Chamber
    )
  ) +
  geom_point(
    shape = 21,
    alpha = 0.75,
    size = 3.5,
    position = position_jitter(
      width = 0.3
    ),
    show.legend = FALSE
  ) +
  geom_boxplot(
    alpha = 0.75,
    outliers = FALSE
  ) +
  scale_fill_manual(
    values = anno_colors$Chamber
  ) +
  facet_grid(
    cols = vars(Medium),
    rows = vars(Metric),
    scales = "free"
  ) +
  labs(
    y = "Diversity estimate",
    x = "Chamber",
    fill = "Chamber"
  ) +
  # theme_classic() +
  theme_bw() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
alpha_gg

# Save plot
# ggsave(
#   filename = "figures/sewageProj/ampl_alpha_div_gg.png",
#   width = 6,
#   height = 6,
#   units = "in",
#   dpi = 300
# )

### Make a df for significance ----
sig_df <- 
  data.frame(
    padj = c(
      p.adjust(obs_pw$summary.table$`Pr > d`, method = "fdr"),
      p.adjust(shannon_pw$summary.table$`Pr > d`, method = "fdr")
    ),
    Chamber1 = rep(c("A","A","B"), times = 2),
    Chamber2 = rep(c("B", "C", "C"), times = 2),
    y = c(
      max(alpha_div_df$Observed) * 1.15,
      max(alpha_div_df$Observed) * 1.05,
      max(alpha_div_df$Observed) * 1.15,
      max(alpha_div_df$Shannon) * 1.155,
      max(alpha_div_df$Shannon) * 1.05,
      max(alpha_div_df$Shannon) * 1.155
    ),
    Medium = rep("NA", times = 6),
    Metric = rep(c("Observed", "Shannon"), each = 3)
  ) |> 
  mutate(
    padj = round(padj, digits = 2),
    padj = if_else(padj==0, 0.01, padj)
  )
sig_df

### Viz 2 ----
alpha_div_df |>
  pivot_longer(
    cols = c(
      "Observed", "ACE", 
      "Shannon", 
      "Simpson", "Inverse Simpson"
    ),
    names_to = "Metric",
    values_to = "Estimate"
  ) |>
  filter(Metric %in% c("Observed", "Shannon")) |> 
  rbind(
    alpha_div_df |>
      pivot_longer(
        cols = c(
          "Observed", "ACE", 
          "Shannon", 
          "Simpson", "Inverse Simpson"
        ),
        names_to = "Metric",
        values_to = "Estimate"
      ) |>
      filter(Metric %in% c("Observed", "Shannon")) |> 
      mutate(Chamber = "All")
  ) |> 
  mutate(
    Chamber = factor(
      Chamber,
      levels = c("All", "A", "B", "C")
    )
  ) |> 
  ggplot(
    mapping = aes(
      x = Chamber,
      y = Estimate,
      fill = Medium
    )
  ) +
  geom_point(
    mapping = aes(
      shape = Medium
    ),
    size = 3,
    alpha = 0.5,
    position = position_jitterdodge(
      dodge.width = 0.75,
      jitter.width = 0.25
    )
  ) +
  geom_point(
    mapping = aes(
      y = Estimate * 1.2
    ),
    alpha = 0
  ) +
  geom_boxplot(
    mapping = aes(
      fill = Medium
    ),
    outliers = FALSE,
    alpha = 0.5 
  ) +
  ggpubr::geom_bracket(
    data = sig_df,
    mapping = aes(
      xmin = Chamber1,
      xmax = Chamber2,
      y.position = y,
      label = paste0("p ≤ ", padj)
    )
  ) +
  facet_grid(
    rows = vars(Metric),
    scales = "free_y"
  ) +
  scale_fill_manual(values = anno_colors$Medium) +
  scale_shape_manual(values = c(21, 24)) +
  theme_bw() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

# Save plot
# ggsave(
#   filename = "figures/sewageProj/ampl_alpha_div_gg_kt.png",
#   width = 8,
#   height = 6,
#   units = "in",
#   dpi = 1200
# )

# Transform data ----
ps_ra <- 
  microbiome::transform(
    ps_dec,
    transform = 'compositional'
  )

# Beta diversity ----
# Calculate BC distance ----
bray_dist <- phyloseq::distance(
  ps_ra, 
  method = "bray"
)

# PERMANOVA 
bray_perm <- 
  RRPP::lm.rrpp(
    bray_dist ~ Chamber * Medium + SampleID,
    data = metadata |> filter(Chamber != "D" & Wash != "W0"),
    SS.type = "III",
    seed = 957,
    iter = 9999,
    Parallel = TRUE
  )

bray_permanova <- 
  anova(
    bray_perm,
    error = c(
      "Medium",
      "Chamber",
      "Residuals",
      "Residuals"
    ),
    effect.type = "F"
  )
bray_permanova
# Chamber and SampleID differ (p < 0.037 and 0.0001)

# Construct a null model
bray_null <- 
  RRPP::lm.rrpp(
    bray_dist ~ SampleID,
    data = metadata |> filter(Chamber != "D" & Wash != "W0"),
    SS.type = "III",
    seed = 957,
    iter = 9999,
    Parallel = TRUE
  )
bray_null

# Pairwise comparison
RRPP::pairwise(
  fit = bray_perm,
  fit.null = bray_null,
  groups = metadata |> filter(Chamber != "D" & Wash != "W0") |> _$Chamber
) |> summary()
# A differed from both B and C

# Ordinate via PCoA ----
pcoa_bray <- 
  phyloseq::ordinate(
    ps_ra,
    method = "PCoA",
    distance = "bray"
  )

# Plot scree plot
plot_scree(pcoa_bray)

# Obtain the first three principal components and make data frame
pcoa_bray_df <- 
  pcoa_bray$vectors[, 1:3] |> 
  data.frame() |>
  merge(
    metadata,
    by = "row.names"
  )
pcoa_bray_df

## Visualize ----
pcoa_bray_df |> 
  ggplot(
    mapping = aes(
      x = Axis.1,
      y = Axis.2
    )
  ) +
  stat_ellipse(
    mapping = aes(color = Chamber),
    alpha = 0.3
  ) +
  geom_point(
    mapping = aes(
      fill = Chamber,
      # shape = Medium
    ),
    shape = 21,
    size = 5,
    alpha = 0.7
  ) +
  # geom_text(mapping = aes(label = run)) +
  # scale_shape_manual(values = c(21,23)) +
  # scale_fill_manual(values = c("grey33", "red")) +
  # scale_color_manual(values = c("grey33", "red")) +
  labs(
    color = "Age",
    fill = "Age",
    shape = "Sex",
    # lty = "Run",
    x = paste0(
      "PCo1 [Variance explained: ",
      round(
        pcoa_bray$values$Relative_eig[1] * 100,
        digits = 2
      ),
      "%]"
    ),
    y = paste0(
      "PCo2 [Variance explained: ",
      round(
        pcoa_bray$values$Relative_eig[2] * 100,
        digits = 2
      ),
      "%]"
    )
  ) +
  theme(panel.grid = element_blank())
# guides(
#   fill = guide_legend(
#     override.aes = list(
#       fill = c("grey33", "red"),
#       color = c("grey33", "red"),
#       shape = 21
#     )
#   )
# )

# ggsave(
#   filename = "../figures/bray_pcoa.png",
#   units = "in",
#   width = 10,
#   height = 8,
#   dpi = 320
# )

# Ordinate via NMDS ----
bray_dist <- 
  df |> 
  select(SampleName, Sequence, RelAbund) |> 
  pivot_wider(
    names_from = "Sequence",
    values_from = "RelAbund"
  ) |> 
  column_to_rownames(var = "SampleName") |> 
  vegan::vegdist(method = "bray")

set.seed(957)
bray_nmds <- 
  vegan::metaMDS(bray_dist) |> 
  vegan::scores() |> 
  as_tibble(rownames = "SampleName")

### Visualize ----
bray_nmds_df <- 
  bray_nmds |> 
  separate(
    col = SampleName,
    into = c("Chamber", "Medium", "Location", "SampleID", "Wash", "five"),
    sep = "-"
  ) |> 
  select(-five) |> 
  mutate(
    O2_ppm = factor(
      Chamber,
      levels = c("A", "D", "B", "C"),
      labels = c("50000", "10000", "10-50", "0")
    )
  )

bray_nmds_gg <- 
  bray_nmds_df |> 
  ggplot(
    mapping = aes(
      x = NMDS1,
      y = NMDS2
    )
  ) +
  geom_point(
    mapping = aes(
      fill = O2_ppm
    ),
    shape = 21,
    size = 5
  ) +
  stat_ellipse(
    mapping = aes(
      color = O2_ppm
    )
  ) +
  scale_fill_manual(values = anno_colors$O2_ppm) +
  scale_color_manual(values = anno_colors$O2_ppm) +
  labs(
    fill = "O~2~ (ppm)",
    color = "O~2~ (ppm)"
  ) +
  theme_classic() +
  theme(
    legend.title = ggtext::element_markdown(),
    legend.position = "bottom"
  )
bray_nmds_gg

# Save plot
# ggsave(
#   filename = "figures/sewageProj/ampl_beta_bray_nmds_gg.png",
#   width = 6,
#   height = 6,
#   units = "in",
#   dpi = 300
# )

# Create a density curve plot per nmds axis
library(cowplot)
nmds1_curve_gg <- 
  bray_nmds_df |> 
  ggplot(
    mapping = aes(
      x = NMDS1,
      fill = O2_ppm
    )
  ) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(values = anno_colors$O2_ppm) +
  cowplot::theme_nothing()
nmds1_curve_gg

nmds2_curve_gg <- 
  bray_nmds_df |> 
  ggplot(
    mapping = aes(
      y = NMDS2,
      fill = O2_ppm
    )
  ) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(values = anno_colors$O2_ppm) +
  cowplot::theme_nothing()
nmds2_curve_gg

v <- cowplot::plot_grid(
  nmds1_curve_gg,
  bray_nmds_gg,
  ncol = 1,
  align = "v",
  rel_heights = c(1,3)
)
v

plot_grid(
  ggdraw(v_aligned[[1]]),
  NULL,
  ggdraw(v_aligned[[2]]),
  ggdraw(h_aligned[[2]]),
  align = "hv",
  ncol = 2,
  rel_heights = c(1,3),
  rel_widths = c(3,1)
)

v_aligned <- cowplot::align_plots(
  nmds1_curve_gg,
  bray_nmds_gg,
  align = "v")

h_aligned <- cowplot::align_plots(
  ggdraw(v_aligned[[2]]),
  nmds2_curve_gg,
  align = "h")

ggdraw(h_aligned[[2]])

# Calculate weighted UniFrac distances ----
wunif_dist <- 
  phyloseq::distance(
    ps_ra,
    method = 'wunifrac'
  )

# PERMANOVA 
wunif_perm <- 
  RRPP::lm.rrpp(
    wunif_dist ~ Chamber * Medium + SampleID,
    data = metadata |> filter(Chamber != "D" & Wash != "W0"),
    SS.type = "III",
    seed = 957,
    iter = 9999,
    Parallel = TRUE
  )

wunif_permanova <- 
  anova(
    wunif_perm,
    error = c(
      "Medium",
      "Chamber",
      "Residuals",
      "Residuals"
    ),
    effect.type = "F"
  )
wunif_permanova
# SampleID differ (p < 0.0001)
# Interaction between Chamber and Medium was tendency (p < 0.09)

# Construct a null model
wunif_null <- 
  RRPP::lm.rrpp(
    wunif_dist ~ SampleID,
    data = metadata |> filter(Chamber != "D" & Wash != "W0"),
    SS.type = "III",
    seed = 957,
    iter = 9999,
    Parallel = TRUE
  )
wunif_null

# Pairwise comparison
RRPP::pairwise(
  fit = wunif_perm,
  fit.null = wunif_null,
  groups = metadata |> filter(Chamber != "D" & Wash != "W0") |> _$Chamber
) |> summary()
# A differed from both B and C

# Ordinate via PCoA ----
pcoa_wunif <- 
  phyloseq::ordinate(
    ps_ra,
    method = "PCoA",
    distance = "wunifrac"
  )

# Plot scree plot
plot_scree(pcoa_wunif)

# Obtain the first three principal components and make data frame
pcoa_wunif_df <- 
  pcoa_wunif$vectors[, 1:3] |> 
  data.frame() |>
  merge(
    metadata,
    by = "row.names"
  )
pcoa_wunif_df

## Visualize ----
pcoa_wunif_df |> 
  ggplot(
    mapping = aes(
      x = Axis.1,
      y = Axis.2
    )
  ) +
  stat_ellipse(
    mapping = aes(color = Chamber),
    alpha = 0.3
  ) +
  geom_point(
    mapping = aes(
      fill = Chamber,
      # shape = Medium
    ),
    shape = 21,
    size = 5,
    alpha = 0.7
  ) +
  # geom_text(mapping = aes(label = run)) +
  # scale_shape_manual(values = c(21,23)) +
  # scale_fill_manual(values = c("grey33", "red")) +
  # scale_color_manual(values = c("grey33", "red")) +
  labs(
    color = "Chamber",
    fill = "Chamber",
    # lty = "Run",
    x = paste0(
      "PCo1 [Variance explained: ",
      round(
        pcoa_bray$values$Relative_eig[1] * 100,
        digits = 2
      ),
      "%]"
    ),
    y = paste0(
      "PCo2 [Variance explained: ",
      round(
        pcoa_bray$values$Relative_eig[2] * 100,
        digits = 2
      ),
      "%]"
    )
  ) +
  theme(panel.grid = element_blank())
# guides(
#   fill = guide_legend(
#     override.aes = list(
#       fill = c("grey33", "red"),
#       color = c("grey33", "red"),
#       shape = 21
#     )
#   )
# )

# ggsave(
#   filename = "../figures/wunif_pcoa.png",
#   units = "in",
#   width = 10,
#   height = 8,
#   dpi = 320
# )


# Differential abundance analysis ----
library(ANCOMBC)

# Create another ps object but with chamber B as the first rank of the chamber factor
ps_ancom <- 
  phyloseq(
    otu_table(
      otu_table(ps_ra)
    ),
    sample_data(
      metadata |> 
        filter(
          Chamber != "D" & Wash != "W0"
        ) |> 
        mutate(
          Chamber = factor(
            Chamber,
            levels = c("B", "A", "C")
          )
        )
    ),
    tax_table(
      tax_table(ps_ra)
    )
  )

ancombc_sp <- 
  ANCOMBC::ancombc2(
    data = ps_ancom,
    tax_level = "Phylum",
    fix_formula = "Chamber",
    # rand_formula = "(1|SampleID)",
    p_adj_method = "fdr",
    pseudo_sens = TRUE,
    # prv_cut = 0.10,
    group = "Chamber",
    struc_zero = TRUE,
    n_cl = 10,
    verbose = TRUE,
    global = TRUE,
    pairwise = TRUE,
    dunnet = TRUE
  )
ancombc_sp

## Show results of ancombc ----
ancombc_sp$res |> View()

# Prepare input for LEfSe ----
df_ra<- 
  sample_data(ps_ra) |> 
  data.frame() |> 
  merge(
    otu_table(ps_ra) |>
      data.frame() |> 
      rownames_to_column(var = "SampleNames") |>
      pivot_longer(
        cols = where(is.numeric),
        names_to = "Sequence",
        values_to = "RelativeAbundance"
      ),
    by = "SampleNames"
  ) |> 
  merge(
    tax_table(ps_ra) |> 
      data.frame() |> 
      rownames_to_column(var = "Sequence"),
    by = "Sequence"
  )

# Save as rds
# saveRDS(
#   df_ra,
#   file = "data/df_ra.rds"
# )

lefse_input <-
  df_ra |> 
  summarize(
    KingdomSum = sum(RelativeAbundance),
    .by = c("SampleID", "Chamber", "Medium", "Kingdom")
  ) |> 
  arrange(Chamber, Medium, SampleID) |> 
  pivot_wider(
    names_from = Kingdom,
    values_from = KingdomSum
  ) |>
  t() |>
  rbind(
    df_ra |> 
      mutate(
        KtoPhylum = if_else(is.na(Phylum)==FALSE, paste0(Kingdom,"|",Phylum), NA)
      ) |> 
      na.omit() |>
      summarize(
        PhylumSum = sum(RelativeAbundance),
        .by = c("SampleID", "Chamber", "Medium", "KtoPhylum")
      ) |> 
      arrange(Chamber, Medium, SampleID) |> 
      pivot_wider(
        names_from = KtoPhylum,
        values_from = PhylumSum
      ) |>
      t() |> 
      tail(-3)
  ) |> 
  rbind(
    df_ra |> 
      mutate(
        KtoClass = if_else(is.na(Class)==FALSE, paste0(Kingdom,"|",Phylum,"|",Class), NA)
      ) |> 
      na.omit() |>
      summarize(
        ClassSum = sum(RelativeAbundance),
        .by = c("SampleID", "Chamber", "Medium", "KtoClass")
      ) |> 
      arrange(Chamber, Medium, SampleID) |> 
      pivot_wider(
        names_from = KtoClass,
        values_from = ClassSum
      ) |>
      t() |> 
      tail(-3)
  ) |> 
  rbind(
    df_ra |> 
      mutate(
        KtoOrder = if_else(is.na(Order)==FALSE, paste0(Kingdom,"|",Phylum,"|",Class,"|",Order), NA)
      ) |> 
      na.omit() |>
      summarize(
        OrderSum = sum(RelativeAbundance),
        .by = c("SampleID", "Chamber", "Medium", "KtoOrder")
      ) |> 
      arrange(Chamber, Medium, SampleID) |> 
      pivot_wider(
        names_from = KtoOrder,
        values_from = OrderSum
      ) |>
      t() |> 
      tail(-3)
  ) |> 
  rbind(
    df_ra |> 
      mutate(
        KtoFamily = if_else(is.na(Family)==FALSE, paste0(Kingdom,"|",Phylum,"|",Class,"|",Order,"|",Family), NA)
      ) |> 
      na.omit() |>
      summarize(
        FamilySum = sum(RelativeAbundance),
        .by = c("SampleID", "Chamber", "Medium", "KtoFamily")
      ) |> 
      arrange(Chamber, Medium, SampleID) |> 
      pivot_wider(
        names_from = KtoFamily,
        values_from = FamilySum
      ) |>
      t() |> 
      tail(-3)
  ) |> 
  rbind(
    df_ra |> 
      mutate(
        KtoGenus = if_else(is.na(Genus)==FALSE, paste0(Kingdom,"|",Phylum,"|",Class,"|",Order,"|",Family,"|",Genus), NA)
      ) |> 
      na.omit() |>
      summarize(
        GenusSum = sum(RelativeAbundance),
        .by = c("SampleID", "Chamber", "Medium", "KtoGenus")
      ) |> 
      arrange(Chamber, Medium, SampleID) |> 
      pivot_wider(
        names_from = KtoGenus,
        values_from = GenusSum
      ) |>
      t() |> 
      tail(-3)
  )
lefse_input |> head(15)

# Output as TSV
# write.table(
#   lefse_input,
#   file = "data/lefse_input.tsv",
#   sep = "\t"
# )

# Make a heatmap ----
filt_seq <- 
  df |> 
  filter(Chamber != "D" & Wash != "W0") |> 
  mutate(
    PrevPerCham = sum(Prevalence, na.rm = TRUE),
    .by = c("Sequence", "Chamber")
  ) |> 
  filter(
    PrevPerCham >= 3 & RelAbund >= 0.01
  ) |> 
  _$Sequence |> 
  unique()

## Filter data frame ----
df_filt <- 
  df |> 
  filter(
    Sequence %in% filt_seq
  )
df_filt

## Set parameters ----
anno_colors <- 
  list(
    # "Sample" = c("AW03"="#FFC20A", "AO02"="#0C7BDC"),
    # "Chamber" = c("A5"="grey88", "A1"="grey66", "B"="grey25","C"="grey11"),
    "Chamber" = c("A"="grey88", "B"="grey55","C"="grey11"),
    "Medium" = c("Y1"="#e1be6a","YC"="#40b0a6"),
    # "Wash" = c("Undiluted" = "#E66100", "Diluted" = "#5D3A9B"),
    "O2_ppm" = c("50000"="grey88", "10000"="grey66", "10-50"="grey25","0"="grey11")
    # "O2_ppm" = c("50000"="grey99", "10-50"="grey55","0"="grey0")
  )

## Visualize ----
df_filt |> 
  filter(is.na(Genus) == FALSE & is.na(Species) == FALSE) |> 
  select(GenusSpecies, SampleName, RelAbund) |> 
  summarize(
    RelAbund = sum(RelAbund, na.rm = TRUE),
    .by = c("GenusSpecies", "SampleName")
  ) |> 
  pivot_wider(
    names_from = "SampleName",
    values_from = "RelAbund"
  ) |> 
  column_to_rownames(
    var = "GenusSpecies"
  ) |> 
  as.matrix() |> 
  pheatmap::pheatmap(
    color = viridis::inferno(100),
    scale = "row",
    display_numbers = TRUE,
    number_color = 'dodgerblue',
    border_color = FALSE,
    # labels_row = as.expression(ital_all_chamber_tax),
    annotation_col = metadata |> filter(Chamber != "D" & Wash != "W0"),
    annotation_colors = anno_colors,
    show_colnames = FALSE
  )
