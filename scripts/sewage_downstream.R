# Script header ----
# Programmed by: Kairi Tanaka
# Programmed on: 04-16-2025
# Programmed to: Compile downstream analysis into one script
# Last modified by: Kairi Tanaka
# Last modified on: 08-04-2025
# Last modified to: Test out other options for DA

# Load packages ----
sewage_packs <- c(
  # For GM Repo
  "httr", "jsonlite", "xml2",
  # For general data handling
  "tidyverse",
  # For microbiome data handling
  "BiocManager", "dada2", "phyloseq", "microbiome", "vegan",
  # For statistical analysis
  "RRPP",
  # For data viz
  "ggpubr", "ggtext", "cowplot", "patchwork", "viridis", "ggpattern"
)

# Install packages if not already installed ----
# lapply(
#   sewage_packs,
#   function(x) {
#     if (!x %in% installed.packages()) {
#       install.packages(x)
#     }
#   }
# )

# lapply(
#   sewage_packs,
#   function(x) {
#     if (!x %in% installed.packages()) {
#       BiocManager::install(x)
#     }
#   }
# )

# Load packages ----
lapply(
  sewage_packs,
  library,
  character.only = TRUE
)

# Set working directory ----
setwd("~/Documents/LyuLab/sewage_health_equity/")
getwd() # Just to confirm

# Set theme for visualizations ----
theme_set(
  theme_cowplot() +
    theme(
      plot.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      legend.justification = "center"
    )
)

# Set options for plots ----
options(
  vsc.dev.args = list(
    width = 1500,
    height = 1500,
    pointsize = 12,
    res = 300
  )
)

# Set colors for covid conditions ----
covid_colors <- c(
  "Negative" = "grey55",
  "Positive" = "red"
)

# Set functions ----
## ps to df ----
ps2df <- function(ps, column.is.taxa = TRUE) {
  if (column.is.taxa == TRUE) {
    otu_table(ps) |>
      merge(sample_data(ps), by = "row.names") |>
      pivot_longer(
        cols = where(is.numeric),
        names_to = "ASV",
        values_to = "Abundance"
      ) |>
      merge(
        tax_table(ps) |>
          data.frame() |>
          rownames_to_column(var = "ASV"),
        by = "ASV"
      )
  } else {
    otu_table(ps) |>
      t() |>
      merge(sample_data(ps), by = "row.names") |>
      pivot_longer(
        cols = where(is.numeric),
        names_to = "ASV",
        values_to = "Abundance"
      ) |>
      merge(
        tax_table(ps) |>
          data.frame() |>
          rownames_to_column(var = "ASV"),
        by = "ASV"
      )
  }
}

# ps to df with taxonomic information ----
ps2df_tax <- function(ps_obj) {
  otu_table(ps_obj) |>
    data.frame() |>
    rownames_to_column(var = "Tax") |>
    pivot_longer(
      cols = where(is.numeric),
      names_to = "SampleID",
      values_to = "Abundance"
    ) |>
    merge(
      sample_data(ps_obj) |>
        data.frame() |>
        rownames_to_column(var = "SampleID"),
      by = "SampleID"
    )
}

# Function for adding logticks to one side only
# https://stackoverflow.com/questions/20128582/is-it-possible-to-have-annotation-logtics-appear-on-only-one-of-the-subplots-u
add_logticks <- function(
    base = 10, sides = "bl", scaled = TRUE,
    short = unit(0.1, "cm"), mid = unit(0.2, "cm"), long = unit(0.3, "cm"),
    colour = "black", size = 0.5, linetype = 1, alpha = 1, color = NULL,
    data = data.frame(x = NA), ...) {
  if (!is.null(color)) {
    colour <- color
  }
  layer(
    geom = "logticks",
    params = list(
      base = base, sides = sides, scaled = scaled, short = short,
      mid = mid, long = long, colour = colour, size = size,
      linetype = linetype, alpha = alpha, ...
    ),
    stat = "identity", data = data, mapping = NULL, inherit.aes = FALSE, position = "identity",
    show.legend = FALSE
  )
}

# Load data ----
ps <- readRDS(file = "intermediate_data/phyloseq/ps_noEuk.rds")
ps

# Filter data ----
filt_taxa <-
  ps |>
  microbiome::transform(transform = "compositional") |>
  phyloseq::filter_taxa(
    function(x) mean(x) > 0.001,
    prune = FALSE
  )
filt_taxa

# Apply to ps
ps_filt <- prune_taxa(filt_taxa, ps)
ps_filt

# Beta diversity ----
## Run statistics using PERMANOVA ----
### Covid positive communities ----
cov_pos_bray <-
  vegdist(
    otu_table(ps_filt |> subset_samples(Covid == "Positive")),
    method = "bray"
  )

cov_pos_meta <-
  sample_data(ps_filt |> subset_samples(Covid == "Positive")) |>
  data.frame()

cov_pos_bray_lm_rrpp <-
  lm.rrpp(
    cov_pos_bray ~ Community2 / Site,
    data = cov_pos_meta,
    iter = 999,
    seed = 1234,
    Parallel = TRUE
  )

cov_pos_bray_perm <-
  anova(
    cov_pos_bray_lm_rrpp,
    effect.type = "F",
    error = c("Residuals", "Community2:Site")
  )
cov_pos_bray_perm

cov_pos_bray_lm_null <-
  lm.rrpp(
    cov_pos_bray ~ 1 / Site,
    data = cov_pos_meta,
    iter = 999,
    seed = 1234,
    Parallel = TRUE
  )

cov_pos_pw <-
  pairwise(
    cov_pos_bray_lm_rrpp,
    cov_pos_bray_lm_null,
    groups = cov_pos_meta$Community2
  ) |>
  summary(show.vectors = TRUE) |>
  _$summary.table |>
  mutate(
    padj = p.adjust(`Pr > d`, method = "BH"),
    sig = case_when(
      padj <= 0.05 & padj > 0.01 ~ "*",
      padj <= 0.01 & padj > 0.005 ~ "**",
      padj <= 0.005 ~ "***",
      TRUE ~ NA
    )
  )
cov_pos_pw

### Covid negative communities ----
cov_neg_bray <-
  vegdist(
    otu_table(ps_filt |> subset_samples(Covid == "Negative")),
    method = "bray"
  )

cov_neg_meta <-
  sample_data(ps_filt |> subset_samples(Covid == "Negative")) |>
  data.frame()

cov_neg_bray_lm_rrpp <-
  lm.rrpp(
    cov_neg_bray ~ Community2 / Site,
    data = cov_neg_meta,
    iter = 999,
    seed = 1234,
    Parallel = TRUE
  )

cov_neg_bray_perm <-
  anova(
    cov_neg_bray_lm_rrpp,
    effect.type = "F",
    error = c("Residuals", "Community2:Site")
  )
cov_neg_bray_perm

cov_neg_bray_lm_null <-
  lm.rrpp(
    cov_neg_bray ~ 1 / Site,
    data = cov_neg_meta,
    iter = 999,
    seed = 1234,
    Parallel = TRUE
  )

cov_neg_pw <-
  pairwise(
    cov_neg_bray_lm_rrpp,
    cov_neg_bray_lm_null,
    groups = cov_neg_meta$Community2
  ) |>
  summary(show.vectors = TRUE) |>
  _$summary.table |>
  mutate(
    padj = p.adjust(`Pr > d`, method = "BH"),
    sig = case_when(
      padj <= 0.05 & padj > 0.01 ~ "*",
      padj <= 0.01 & padj > 0.005 ~ "**",
      padj <= 0.005 ~ "***",
      TRUE ~ NA
    )
  )
cov_neg_pw

### Elderly communities ----
com_eld_bray <-
  vegdist(
    otu_table(ps_filt |> subset_samples(Community2 == "Elderly" & Covid != "NA")),
    method = "bray"
  )

com_eld_meta <-
  sample_data(ps_filt |> subset_samples(Community2 == "Elderly" & Covid != "NA")) |>
  data.frame()

com_eld_bray_lm_rrpp <-
  lm.rrpp(
    com_eld_bray ~ Covid / Site,
    data = com_eld_meta,
    iter = 999,
    seed = 1234,
    Parallel = TRUE
  )

com_eld_bray_perm <-
  anova(
    com_eld_bray_lm_rrpp,
    effect.type = "F",
    error = c("Residuals", "Covid:Site")
  )
com_eld_bray_perm

com_eld_bray_lm_null <-
  lm.rrpp(
    com_eld_bray ~ 1 / Site,
    data = com_eld_meta,
    iter = 999,
    seed = 1234,
    Parallel = TRUE
  )

com_eld_pw <-
  pairwise(
    com_eld_bray_lm_rrpp,
    com_eld_bray_lm_null,
    groups = com_eld_meta$Covid
  ) |>
  summary(show.vectors = TRUE) |>
  _$summary.table |>
  mutate(
    padj = p.adjust(`Pr > d`, method = "BH"),
    sig = case_when(
      padj <= 0.05 & padj > 0.01 ~ "*",
      padj <= 0.01 & padj > 0.005 ~ "**",
      padj <= 0.005 ~ "***",
      TRUE ~ NA
    )
  )
com_eld_pw

### Marginalized communities ----
com_mar_bray <-
  vegdist(
    otu_table(ps_filt |> subset_samples(Community2 == "Marginalized")),
    method = "bray"
  )

com_mar_meta <-
  sample_data(ps_filt |> subset_samples(Community2 == "Marginalized")) |>
  data.frame()

com_mar_bray_lm_rrpp <-
  lm.rrpp(
    com_mar_bray ~ Covid / Site,
    data = com_mar_meta,
    iter = 999,
    seed = 1234,
    Parallel = TRUE
  )

com_mar_bray_perm <-
  anova(
    com_mar_bray_lm_rrpp,
    effect.type = "F",
    error = c("Residuals", "Covid:Site")
  )
com_mar_bray_perm

com_mar_bray_lm_null <-
  lm.rrpp(
    com_mar_bray ~ 1 / Site,
    data = com_mar_meta,
    iter = 999,
    seed = 1234,
    Parallel = TRUE
  )

com_mar_pw <-
  pairwise(
    com_mar_bray_lm_rrpp,
    com_mar_bray_lm_null,
    groups = com_mar_meta$Covid
  ) |>
  summary(show.vectors = TRUE) |>
  _$summary.table |>
  mutate(
    padj = p.adjust(`Pr > d`, method = "BH"),
    sig = case_when(
      padj <= 0.05 & padj > 0.01 ~ "*",
      padj <= 0.01 & padj > 0.005 ~ "**",
      padj <= 0.005 ~ "***",
      TRUE ~ NA
    )
  )
com_mar_pw

### Minority communities ----
com_min_bray <-
  vegdist(
    otu_table(ps_filt |> subset_samples(Community2 == "Minority")),
    method = "bray"
  )

com_min_meta <-
  sample_data(ps_filt |> subset_samples(Community2 == "Minority")) |>
  data.frame()

com_min_bray_lm_rrpp <-
  lm.rrpp(
    com_min_bray ~ Covid / Site,
    data = com_min_meta,
    iter = 999,
    seed = 1234,
    Parallel = TRUE
  )

com_min_bray_perm <-
  anova(
    com_min_bray_lm_rrpp,
    effect.type = "F",
    error = c("Residuals", "Covid:Site")
  )
com_min_bray_perm

com_min_bray_lm_null <-
  lm.rrpp(
    com_min_bray ~ 1 / Site,
    data = com_min_meta,
    iter = 999,
    seed = 1234,
    Parallel = TRUE
  )

com_min_pw <-
  pairwise(
    com_min_bray_lm_rrpp,
    com_min_bray_lm_null,
    groups = com_min_meta$Covid
  ) |>
  summary(show.vectors = TRUE) |>
  _$summary.table |>
  mutate(
    padj = p.adjust(`Pr > d`, method = "BH"),
    sig = case_when(
      padj <= 0.05 & padj > 0.01 ~ "*",
      padj <= 0.01 & padj > 0.005 ~ "**",
      padj <= 0.005 ~ "***",
      TRUE ~ NA
    )
  )
com_min_pw

### Min and Mar together ----
com_minmar_bray <-
  vegdist(
    otu_table(ps_filt |> subset_samples(Community2 != "Elderly")),
    method = "bray"
  )

com_minmar_meta <-
  sample_data(ps_filt |> subset_samples(Community2 != "Elderly")) |>
  data.frame()

com_minmar_bray_lm_rrpp <-
  lm.rrpp(
    com_minmar_bray ~ Covid / Site,
    data = com_minmar_meta,
    iter = 999,
    seed = 1234,
    Parallel = TRUE
  )

com_minmar_bray_perm <-
  anova(
    com_minmar_bray_lm_rrpp,
    effect.type = "F",
    error = c("Residuals", "Covid:Site")
  )
com_minmar_bray_perm

com_minmar_bray_lm_null <-
  lm.rrpp(
    com_minmar_bray ~ 1 / Site,
    data = com_minmar_meta,
    iter = 999,
    seed = 1234,
    Parallel = TRUE
  )

com_minmar_pw <-
  pairwise(
    com_minmar_bray_lm_rrpp,
    com_minmar_bray_lm_null,
    groups = com_minmar_meta$Covid
  ) |>
  summary(show.vectors = TRUE) |>
  _$summary.table |>
  mutate(
    padj = p.adjust(`Pr > d`, method = "BH"),
    sig = case_when(
      padj <= 0.05 & padj > 0.01 ~ "*",
      padj <= 0.01 & padj > 0.005 ~ "**",
      padj <= 0.005 ~ "***",
      TRUE ~ NA
    )
  )
com_minmar_pw

## Covid positive communities ----
### Ordinate using Bray-Curtis dissimilarity measures ----
cov_pos_bray_pcoa <-
  ps_filt |>
  # microbiome::transform(transform = 'rclr') |>
  subset_samples(Covid == "Positive") |>
  phyloseq::ordinate(
    method = "PCoA",
    distance = "bray"
  )

# Calculate centroids
cov_pos_bray_centroids <- aggregate(
  cbind(Axis.1, Axis.2) ~ Community2,
  data = cov_pos_bray_pcoa$vectors[, 1:4] |>
    merge(
      sample_data(ps_filt),
      by = "row.names"
    ),
  FUN = mean
)

# Add letters for significantly different centroids
cov_pos_bray_centroids$letter <- c("a", "b", "c")

### Visualize ----
cov_pos_bray_pcoa_gg <-
  cov_pos_bray_pcoa$vectors[, 1:4] |>
  merge(
    sample_data(ps_filt),
    by = "row.names"
  ) |>
  ggplot(mapping = aes(x = Axis.1, y = Axis.2)) +
  stat_stars(mapping = aes(color = Community2), geom = "segment") +
  geom_point(mapping = aes(fill = Community2, shape = Community2), size = 3) +
  stat_stars(mapping = aes(color = Community2), geom = "point", shape = 16, size = 9, show.legend = FALSE) +
  stat_stars(mapping = aes(color = Community2), fill = "grey33", geom = "point", shape = 21, size = 7, show.legend = FALSE) +
  geom_text(data = cov_pos_bray_centroids, mapping = aes(label = letter), color = "white", fontface = "bold", size = 5) +
  scale_y_continuous(sec.axis = dup_axis()) +
  scale_x_continuous(sec.axis = dup_axis()) +
  scale_shape_manual(values = c(24, 23, 22)) +
  # scale_fill_manual( values = viridis::turbo(n = 3, begin = 0.1, end = 0.9)) +
  # scale_color_manual(values = viridis::turbo(n = 3, begin = 0.1, end = 0.9)) +
  # scale_fill_manual( values = dichromat::dichromat(viridis::turbo(n = 3, begin = 0.1, end = 0.9))) +
  # scale_color_manual(values = dichromat::dichromat(viridis::turbo(n = 3, begin = 0.1, end = 0.9))) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  # scale_fill_manual( values = dichromat::dichromat(RColorBrewer::brewer.pal(3,'Dark2'))) +
  # scale_color_manual(values = dichromat::dichromat(RColorBrewer::brewer.pal(3,'Dark2'))) +
  labs(
    title = "SARS-CoV-2 Positive",
    x = paste0("PC1 (Variance explained: ", round(cov_pos_bray_pcoa$values$Relative_eig[1] * 100, 2), "%)"),
    y = paste0("PC2 (Variance explained: ", round(cov_pos_bray_pcoa$values$Relative_eig[2] * 100, 2), "%)"),
    fill = "Community",
    color = "Community",
    shape = "Community"
  ) +
  theme(
    axis.title.x.top = element_blank(),
    axis.title.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.text.y.right = element_blank(),
    legend.position = "bottom",
    legend.justification = "center"
  )
cov_pos_bray_pcoa_gg

## Covid negative communities ----
### Ordinate using Bray-Curtis dissimilarity measures ----
cov_neg_bray_pcoa <-
  ps_filt |>
  # microbiome::transform(transform = 'rclr') |>
  subset_samples(Covid == "Negative") |>
  phyloseq::ordinate(
    method = "PCoA",
    distance = "bray"
  )

# Calculate centroids
cov_neg_bray_centroids <- aggregate(
  cbind(Axis.1, Axis.2) ~ Community2,
  data = cov_neg_bray_pcoa$vectors[, 1:4] |>
    merge(
      sample_data(ps_filt),
      by = "row.names"
    ),
  FUN = mean
)
cov_neg_bray_centroids

# Add letters for significantly different centroids
cov_neg_bray_centroids$letter <- c("a", "b", "b")

### Visualize ----
cov_neg_bray_pcoa_gg <-
  cov_neg_bray_pcoa$vectors[, 1:4] |>
  merge(
    sample_data(ps_filt),
    by = "row.names"
  ) |>
  ggplot(mapping = aes(x = Axis.1, y = Axis.2)) +
  stat_stars(mapping = aes(color = Community2), geom = "segment") +
  geom_point(mapping = aes(fill = Community2, shape = Community2), size = 3) +
  stat_stars(mapping = aes(color = Community2), geom = "point", shape = 16, size = 9, show.legend = FALSE) +
  stat_stars(mapping = aes(color = Community2), fill = "grey33", geom = "point", shape = 21, size = 7, show.legend = FALSE) +
  geom_text(data = cov_neg_bray_centroids, mapping = aes(label = letter), color = "white", fontface = "bold", size = 5) +
  scale_y_continuous(sec.axis = dup_axis()) +
  scale_x_continuous(sec.axis = dup_axis()) +
  scale_shape_manual(values = c(24, 23, 22)) +
  # scale_fill_manual( values = viridis::turbo(n = 3, begin = 0.1, end = 0.9)) +
  # scale_color_manual(values = viridis::turbo(n = 3, begin = 0.1, end = 0.9)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  labs(
    title = "SARS-CoV-2 Negative",
    x = paste0("PC1 (Variance explained: ", round(cov_neg_bray_pcoa$values$Relative_eig[1] * 100, 2), "%)"),
    y = paste0("PC2 (Variance explained: ", round(cov_neg_bray_pcoa$values$Relative_eig[2] * 100, 2), "%)"),
    fill = "Community",
    color = "Community",
    shape = "Community"
  ) +
  theme(
    axis.title.x.top = element_blank(),
    axis.title.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.text.y.right = element_blank(),
    legend.position = "bottom",
    legend.justification = "center"
  )

cov_neg_bray_pcoa_gg

## Elderly community ----
### Ordinate using Bray-Curtis dissimilarity measures ----
eld_bray_pcoa <-
  ps_filt |>
  # microbiome::transform(transform = 'rclr') |>
  subset_samples(Community2 == "Elderly" & Covid != "NA") |>
  phyloseq::ordinate(
    method = "PCoA",
    distance = "bray"
  )

# Calculate centroids
com_eld_bray_centroids <- aggregate(
  cbind(Axis.1, Axis.2) ~ Covid,
  data = eld_bray_pcoa$vectors[, 1:4] |>
    merge(
      sample_data(ps_filt),
      by = "row.names"
    ),
  FUN = mean
)
com_eld_bray_centroids

# Add letters for significantly different centroids
com_eld_bray_centroids$letter <- c("a", "b")

### Visualize ----
eld_bray_pcoa_gg <-
  eld_bray_pcoa$vectors[, 1:4] |>
  merge(
    sample_data(ps_filt),
    by = "row.names"
  ) |>
  ggplot(mapping = aes(x = Axis.1, y = Axis.2)) +
  stat_stars(mapping = aes(color = Covid), geom = "segment") +
  geom_point(mapping = aes(fill = Covid), shape = 24, size = 3) +
  stat_stars(mapping = aes(color = Covid), geom = "point", shape = 16, size = 9, show.legend = FALSE) +
  stat_stars(mapping = aes(color = Covid), fill = "grey33", geom = "point", shape = 21, size = 7, show.legend = FALSE) +
  geom_text(data = com_eld_bray_centroids, mapping = aes(label = letter), color = "white", fontface = "bold", size = 5) +
  scale_y_continuous(sec.axis = dup_axis()) +
  scale_x_continuous(sec.axis = dup_axis()) +
  scale_fill_manual(values = covid_colors) +
  scale_color_manual(values = covid_colors) +
  # scale_fill_manual(values = viridis::inferno(n = 3, begin = 0.1, end = 0.9)) +
  # scale_color_manual(values = viridis::inferno(n = 3, begin = 0.1, end = 0.9)) +
  labs(
    title = "Elderly",
    x = paste0("PC1 (Variance explained: ", round(eld_bray_pcoa$values$Relative_eig[1] * 100, 2), "%)"),
    y = paste0("PC2 (Variance explained: ", round(eld_bray_pcoa$values$Relative_eig[2] * 100, 2), "%)"),
    fill = "SARS-CoV-2",
    color = "SARS-CoV-2"
  ) +
  theme(
    axis.title.x.top = element_blank(),
    axis.title.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.text.y.right = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.2, 0.85),
    legend.box.margin = margin(0, 0, 0, 0, "pt"),
    legend.margin = margin(5, 5, 5, 5, "pt"),
    legend.background = element_rect(color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
eld_bray_pcoa_gg

## Marginalized community ----
### Ordinate using Bray-Curtis dissimilarity measures ----
mar_bray_pcoa <-
  ps_filt |>
  # microbiome::transform(transform = 'rclr') |>
  subset_samples(Community2 == "Marginalized" & Covid != "NA") |>
  phyloseq::ordinate(
    method = "PCoA",
    distance = "bray"
  )

# Calculate centroids
com_mar_bray_centroids <- aggregate(
  cbind(Axis.1, Axis.2) ~ Covid,
  data = mar_bray_pcoa$vectors[, 1:4] |>
    merge(
      sample_data(ps_filt),
      by = "row.names"
    ),
  FUN = mean
)
com_mar_bray_centroids

# Add letters for significantly different centroids
com_mar_bray_centroids$letter <- c("a", "a")

### Visualize ----
mar_bray_pcoa_gg <-
  mar_bray_pcoa$vectors[, 1:4] |>
  merge(
    sample_data(ps_filt),
    by = "row.names"
  ) |>
  ggplot(mapping = aes(x = Axis.1, y = Axis.2)) +
  stat_stars(mapping = aes(color = Covid), geom = "segment") +
  geom_point(mapping = aes(fill = Covid), shape = 23, size = 3) +
  stat_stars(mapping = aes(color = Covid), geom = "point", shape = 16, size = 9, show.legend = FALSE) +
  stat_stars(mapping = aes(color = Covid), fill = "grey33", geom = "point", shape = 21, size = 7, show.legend = FALSE) +
  geom_text(data = com_mar_bray_centroids, mapping = aes(label = letter), color = "white", fontface = "bold", size = 5) +
  scale_y_continuous(sec.axis = dup_axis()) +
  scale_x_continuous(sec.axis = dup_axis()) +
  scale_fill_manual(values = covid_colors) +
  scale_color_manual(values = covid_colors) +
  # scale_fill_manual(values = viridis::inferno(n = 3, begin = 0.1, end = 0.9)) +
  # scale_color_manual(values = viridis::inferno(n = 3, begin = 0.1, end = 0.9)) +
  labs(
    title = "Marginalized",
    x = paste0("PC1 (Variance explained: ", round(mar_bray_pcoa$values$Relative_eig[1] * 100, 2), "%)"),
    y = paste0("PC2 (Variance explained: ", round(mar_bray_pcoa$values$Relative_eig[2] * 100, 2), "%)"),
    fill = "SARS-CoV-2",
    color = "SARS-CoV-2"
  ) +
  theme(
    axis.title.x.top = element_blank(),
    axis.title.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.text.y.right = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.55, 0.25),
    legend.box.margin = margin(0, 0, 0, 0, "pt"),
    legend.margin = margin(5, 5, 5, 5, "pt"),
    legend.background = element_rect(color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
mar_bray_pcoa_gg

## Minority community ----
### Ordinate using Bray-Curtis dissimilarity measures ----
min_bray_pcoa <-
  ps_filt |>
  # microbiome::transform(transform = 'rclr') |>
  subset_samples(Community2 == "Minority" & Covid != "NA") |>
  phyloseq::ordinate(
    method = "PCoA",
    distance = "bray"
  )

# Calculate centroids
com_min_bray_centroids <- aggregate(
  cbind(Axis.1, Axis.2) ~ Covid,
  data = min_bray_pcoa$vectors[, 1:4] |>
    merge(
      sample_data(ps_filt),
      by = "row.names"
    ),
  FUN = mean
)
com_min_bray_centroids

# Add letters for significantly different centroids
com_min_bray_centroids$letter <- c("a", "a")

### Visualize ----
min_bray_pcoa_gg <-
  min_bray_pcoa$vectors[, 1:4] |>
  merge(
    sample_data(ps_filt),
    by = "row.names"
  ) |>
  ggplot(mapping = aes(x = Axis.1, y = Axis.2)) +
  stat_stars(mapping = aes(color = Covid), geom = "segment") +
  geom_point(mapping = aes(fill = Covid), shape = 22, size = 3) +
  stat_stars(mapping = aes(color = Covid), geom = "point", shape = 16, size = 9, show.legend = FALSE) +
  stat_stars(mapping = aes(color = Covid), fill = "grey33", geom = "point", shape = 21, size = 7, show.legend = FALSE) +
  geom_text(data = com_min_bray_centroids, mapping = aes(label = letter), color = "white", fontface = "bold", size = 5) +
  scale_y_continuous(sec.axis = dup_axis()) +
  scale_x_continuous(sec.axis = dup_axis()) +
  scale_fill_manual(values = covid_colors) +
  scale_color_manual(values = covid_colors) +
  # scale_fill_manual(values = viridis::inferno(n = 3, begin = 0.1, end = 0.9)) +
  # scale_color_manual(values = viridis::inferno(n = 3, begin = 0.1, end = 0.9)) +
  labs(
    title = "Minority",
    x = paste0("PC1 (Variance explained: ", round(min_bray_pcoa$values$Relative_eig[1] * 100, 2), "%)"),
    y = paste0("PC2 (Variance explained: ", round(min_bray_pcoa$values$Relative_eig[2] * 100, 2), "%)"),
    fill = "SARS-CoV-2",
    color = "SARS-CoV-2"
  ) +
  theme(
    axis.title.x.top = element_blank(),
    axis.title.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.text.y.right = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.4, 0.15),
    legend.box.margin = margin(0, 0, 0, 0, "pt"),
    legend.margin = margin(5, 5, 5, 5, "pt"),
    legend.background = element_rect(color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
min_bray_pcoa_gg

## Marginalized and Minority combined ----
### Ordinate using Bray-Curtis dissimilarity measures ----
marmin_bray_pcoa <-
  ps_filt |>
  # microbiome::transform(transform = 'rclr') |>
  subset_samples(Community2 != "Elderly" & Covid != "NA") |>
  phyloseq::ordinate(
    method = "PCoA",
    distance = "bray"
  )

# Calculate centroids
com_marmin_bray_centroids <- aggregate(
  cbind(Axis.1, Axis.2) ~ Covid,
  data = marmin_bray_pcoa$vectors[, 1:4] |>
    merge(
      sample_data(ps_filt),
      by = "row.names"
    ),
  FUN = mean
)
com_marmin_bray_centroids

# Add letters for significantly different centroids
com_marmin_bray_centroids$letter <- c("a", "a")

covid_colors2 <- c(
  "Marginalized, Negative" = "grey55",
  "Marginalized, Positive" = "red",
  "Minority, Negative" = "grey55",
  "Minority, Positive" = "red"
)

### Visualize ----
marmin_bray_pcoa_gg <-
  marmin_bray_pcoa$vectors[, 1:4] |>
  merge(
    sample_data(ps_filt),
    by = "row.names"
  ) |>
  mutate(
    CovCom = factor(paste0(Community2, ", ", Covid), levels = c("Marginalized, Negative", "Marginalized, Positive", "Minority, Negative", "Minority, Positive"))
  ) |>
  ggplot(mapping = aes(x = Axis.1, y = Axis.2)) +
  stat_stars(mapping = aes(color = CovCom, group = Covid), geom = "segment", show.legend = FALSE) +
  geom_point(mapping = aes(fill = CovCom, shape = CovCom), size = 3, color = "black") +
  stat_stars(mapping = aes(color = Covid), geom = "point", shape = 16, size = 9, show.legend = FALSE) +
  stat_stars(mapping = aes(color = Covid), fill = "grey33", geom = "point", shape = 21, size = 7, show.legend = FALSE) +
  geom_text(data = com_marmin_bray_centroids, mapping = aes(x = Axis.1, y = Axis.2, label = letter), color = "white", fontface = "bold", size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)), sec.axis = dup_axis()) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1)), sec.axis = dup_axis()) +
  scale_fill_manual(values = rep(c("grey55", "red"), times = 2)) +
  scale_color_manual(values = rep(c("grey55", "red"), times = 3)) +
  # scale_fill_manual(
  #   values = rep(viridis::inferno(n = 3, begin = 0.1, end = 0.9)[1:2], times = 2)
  # ) +
  # scale_color_manual(
  #   values = rep(viridis::inferno(n = 3, begin = 0.1, end = 0.9)[1:2], times = 3)
  # ) +
  scale_shape_manual(
    values = c(23, 23, 22, 22)
  ) +
  labs(
    title = "Marginalized & Minority Combined",
    x = paste0("PC1 (Variance explained: ", round(marmin_bray_pcoa$values$Relative_eig[1] * 100, 2), "%)"),
    y = paste0("PC2 (Variance explained: ", round(marmin_bray_pcoa$values$Relative_eig[2] * 100, 2), "%)"),
    fill = "Community & SARS-CoV-2",
    color = "Community & SARS-CoV-2",
    shape = "Community & SARS-CoV-2"
  ) +
  theme(
    axis.title.x.top = element_blank(),
    axis.title.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.text.y.right = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.65, 0.85),
    legend.box.margin = margin(0, 0, 0, 0, "pt"),
    legend.margin = margin(5, 5, 5, 5, "pt"),
    legend.background = element_rect(color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
marmin_bray_pcoa_gg

## Combine figures ----
### Figure 1 - Community per covid ----
cov_neg_bray_pcoa_gg +
  cov_pos_bray_pcoa_gg +
  plot_annotation(tag_level = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.justification = "center")

# ggsave(
#   filename = "processed_data/figures/version4/figure1.png",
#   units = "in", width = 12, height = 6, dpi = 320
# )

### Figure 2 - Covid per community ----

layout <- "
  AABBCC
  AADDDD
"

eld_bray_pcoa_gg +
  mar_bray_pcoa_gg +
  min_bray_pcoa_gg +
  marmin_bray_pcoa_gg +
  plot_annotation(tag_level = "A") +
  plot_layout(design = layout, heights = c(0.6, 0.8, 0.6))

# ggsave(
#   filename = "processed_data/figures/version4/figure2.png",
#   units = "in", width = 12, height = 10, dpi = 320
# )

# Core microbiome ----
# ps_filt_ra <-
#   microbiome::transform(ps, transform = "compositional") |>
#   aggregate_taxa(level = "Genus")

# communities <- unique(sample_data(ps_filt_ra)$Community2)

# ## Parse out core microbiome using for-loop ----
# list_core <- c() # an empty object to store information

# for (n in communities) { # for each variable n in DiseaseState
#   # print(paste0("Identifying Core Taxa for ", n))

#   ps.sub <- subset_samples(ps_filt_ra, Community2 == n) # Choose sample from DiseaseState by n

#   core_m <- core_members(
#     ps.sub,
#     detection = 0,
#     prevalence = 0
#   )
#   print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
#   list_core[[n]] <- core_m # add to a list core taxa for each group.
#   # print(list_core)
# }

# comm_cols <- c(RColorBrewer::brewer.pal(3, "Dark2"))
# plot(
#   eulerr::venn(list_core),
#   fills = RColorBrewer::brewer.pal(3, "Dark2")
# )

# Differential abundance analysis via Lefse ----
library(microbiomeMarker)
## Community comparisons across the covid status ----
com_res <- run_lefse(
  ps_filt,
  group = "Community2",
  subgroup = "Covid",
  multigrp_strat = TRUE, # enable multi-class mode
  strict = "1",
  kw_cutoff = 0.05,
  wilcoxon_cutoff = 0.05,
  lda_cutoff = 4,
  norm = "CPM"
)
com_res

# Save the results ----
# saveRDS(
#   com_res,
#   file = "intermediate_data/lefse_results/com_lefse_results.rds"
# )

cov_all_lefse_res <-
  com_res@marker_table |>
  data.frame() |>
  mutate(
    feature = str_remove_all(string = feature, pattern = "\\w_{2,3}")
  ) |>
  tidyr::separate(
    col = feature,
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
    sep = "\\|"
  ) |>
  mutate_if(is.character, ~ na_if(., "")) |>
  rownames_to_column(var = "markers") |>
  mutate(
    # Genus = if_else(str_detect(Genus,'^\\W')==T,NA,Genus),
    # Genus = if_else(str_detect(Genus,'')==T,NA,Genus),
    # Genus = if_else(str_detect(Genus,'^_')==T,NA,Genus),
    # Family = if_else(str_detect(Genus,'^\\W')==T,NA,Family),
    feature = case_when(
      is.na(Genus) == FALSE ~ paste0("<i>", Genus, "</i> (Genus)"),
      is.na(Genus) & is.na(Family) == FALSE ~ paste0("<i>", Family, "</i> (Family)"),
      is.na(Genus) & is.na(Family) & is.na(Order) == FALSE ~ paste0("<i>", Order, "</i> (Order)"),
      is.na(Genus) & is.na(Family) & is.na(Order) & is.na(Class) == FALSE ~ paste0("<i>", Class, "</i> (Class)"),
      is.na(Genus) & is.na(Family) & is.na(Order) & is.na(Class) & is.na(Phylum) == FALSE ~ paste0("<i>", Phylum, "</i> (Phylum)"),
      TRUE ~ NA
    ),
    markerID = str_extract(markers, "\\d{1,2}$") |> as.numeric(),
    up_down = if_else(
      markerID %in% c(3, 4),
      "Depleted",
      "Enriched"
    ),
    grouping = paste0(up_down, " in ", enrich_group),
    grouping = factor(
      grouping,
      levels = c(
        "Depleted in Elderly",
        "Enriched in Elderly",
        "Depleted in Marginalized",
        "Enriched in Marginalized"
      )
    )
  )

# Visualize result
fig3 <-
  cov_all_lefse_res |>
  # filter(markerID != 23) |>
  ggplot(mapping = aes(x = ef_lda, y = reorder(feature, -markerID))) +
  geom_col(mapping = aes(fill = enrich_group), color = "black") +
  # ggpattern::geom_col_pattern(
  #   mapping = aes(pattern = grouping, fill = grouping),
  #   color = "black",
  #   pattern_density = 0.25,
  #   pattern_fill = "white",
  #   pattern_color = "white",
  #   pattern_key_scale_factor = 0.4,
  #   pattern_res = 320
  # ) +
  # facet_grid(rows = vars(enrich_group), scales = "free_y", space = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  # scale_fill_manual(values = c(RColorBrewer::brewer.pal(3, "Dark2")[1:2], RColorBrewer::brewer.pal(3, "Dark2")[2])) +
  # scale_pattern_manual(values = c("none", "stripe", "none")) +
  labs(
    x = "LDA score",
    fill = "Enrichment"
    # pattern = "Enrichment/Depletion"
  ) +
  theme(
    axis.text.y = ggtext::element_markdown(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(-0.22, "lines"),
    legend.position = "right",
  )
fig3
options(
  vsc.dev.args = list(
    width = 2500,
    height = 1500,
    pointsize = 12,
    res = 300
  )
)
# ggsave(
#   filename = 'processed_data/figures/version2/figure3.png',
#   units = 'in', width = 8, height = 6, dpi = 320
# )

# Plot the p-values
cov_all_lefse_res |>
  ggplot(mapping = aes(x = -log10(padj), y = reorder(feature, -markerID), fill = enrich_group)) +
  geom_point(mapping = aes(shape = enrich_group), size = 2, color = "white", show.legend = FALSE) +
  # facet_grid(rows = vars(enrich_group), scales = "free_y", space = "free_y") +
  scale_x_continuous(
    limits = c(0, 9),
    expand = expansion(mult = c(0, 0.1)),
    sec.axis = dup_axis()
  ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  scale_shape_manual(values = c(24, 23)) +
  labs(x = "-log<sub>10</sub>(p<sub>adj</sub>)", fill = "Community", title = "Overall") +
  theme(
    axis.text.y = ggtext::element_markdown(),
    axis.title.y = element_blank(),
    # axis.text.y.left = element_blank(),
    axis.title.x = ggtext::element_markdown(),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.line.x.bottom = element_blank(),
    axis.title.x.bottom = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(-0.22, "lines"),
    plot.background = element_blank(),
    legend.position = "right"
  ) -> fig3_p
fig3_p

# Combine figures
aligned_plots <- cowplot::align_plots(fig3, fig3_p, align = "hv", axis = "tblr")

com_ova_lefse_gg <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
com_ova_lefse_gg

## Plot relative abundance of significant taxa ----
### Stratify data by tax levels ----
tax_levels <- tax_table(ps_filt) |> colnames()
tax_levels

# Use for loop to create a list of dataframes
df_ps_filt_list <- list()

for (i in 1:length(tax_levels)) {
  tax_level <- tax_levels[i]
  df_ps_filt_list[[i]] <- ps_filt |>
    microbiome::transform(transform = "compositional") |>
    aggregate_taxa(level = tax_level) |>
    ps2df_tax() |>
    mutate(
      TaxLevel = tax_level,
      feature = paste0("<i>", Tax, "</i> (", tax_level, ")")
    )
  # summarize(
  #   TotalRA = sum(Abundance, na.rm = TRUE),
  #   .by = c('TaxLevelName', 'Community2')
  # )
}
# df_ps_filt_list

# Combine into one and visualize
cov_all_ra_box <-
  do.call(rbind, df_ps_filt_list) |>
  data.frame() |>
  # dplyr::filter(TaxLevelName %in% cov_all_lefse_res$feature) |>
  inner_join(cov_all_lefse_res, by = "feature") |>
  filter(Abundance > 0) |>
  ggplot(mapping = aes(x = Abundance * 100, y = reorder(feature, -markerID))) +
  geom_boxplot(mapping = aes(fill = Community2), position = position_dodge(width = 0.9), outliers = FALSE, alpha = 0.5) +
  geom_point(mapping = aes(fill = Community2, shape = Community2), position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2)) +
  scale_x_continuous(
    transform = "log10",
    sec.axis = dup_axis(),
    labels = scales::number_format(accuracy = 0.1),
    limits = c(0.01, 100)
  ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  scale_shape_manual(values = c(24, 23, 22)) +
  labs(color = "Community", fill = "Community", shape = "Community", title = "Relative abundance", x = "Relative abundance (%)") +
  annotation_logticks(sides = "tb") +
  theme(
    axis.title.y = element_blank(),
    # axis.text.y = ggtext::element_markdown()
    axis.text.y = element_blank(),
    legend.position = "bottom",
    plot.margin = unit(c(10, 8, 7, 0), "pt"),
  )
cov_all_ra_box

options(
  vsc.dev.args = list(
    width = 3000,
    height = 1500,
    pointsize = 12,
    res = 300
  )
)

# Put figures side by side
fig3a <- ggarrange(
  # fig3,
  com_ova_lefse_gg,
  cov_all_ra_box,
  # align = 'h',
  # axis = 'tb',
  widths = c(0.9, 1)
  # heights = c(0.5, 1)
  # legend = "bottom"
)
fig3a

# ggsave(filename = 'processed_data/figures/version4/cov_gg_c.png',
#        units = 'in', width = 15, height = 6, dpi = 320)

## Community comparisons covid negative ----
cov_neg_res <- run_lefse(
  ps_filt |> subset_samples(Covid == "Negative"),
  group = "Community2",
  multigrp_strat = TRUE, # enable multi-class mode
  strict = "1", # "1" = all‑against‑all (most stringent)
  kw_cutoff = 0.05,
  wilcoxon_cutoff = 0.05,
  lda_cutoff = 4
)
cov_neg_res

# Save the results
# saveRDS(
#   cov_neg_res,
#   file = "intermediate_data/lefse_results/cov_neg_lefse_results.rds"
# )

### Visualize ----
options(
  vsc.dev.args = list(
    width = 3000,
    height = 3000,
    pointsize = 12,
    res = 300
  )
)

cov_neg_lefse_res <-
  cov_neg_res@marker_table |>
  data.frame() |>
  mutate(
    feature = str_remove_all(string = feature, pattern = "\\w_{2,3}")
  ) |>
  tidyr::separate(
    col = feature,
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
    sep = "\\|"
  ) |>
  mutate_if(is.character, ~ na_if(., "")) |>
  rownames_to_column(var = "markers") |>
  mutate(
    # Genus = if_else(str_detect(Genus,'^\\W')==T,NA,Genus),
    # Genus = if_else(str_detect(Genus,'')==T,NA,Genus),
    # Genus = if_else(str_detect(Genus, "_$") == T, NA, Genus),
    Genus = if_else(str_detect(Genus, "_$") == T, "unknown_genera", Genus),
    # Family = if_else(str_detect(Genus,'^\\W')==T,NA,Family),
    feature = case_when(
      is.na(Genus) == FALSE & Genus != "unknown_genera" ~ paste0("<i>", Genus, "</i> (Genus)"),
      is.na(Genus) & is.na(Family) == FALSE ~ paste0("<i>", Family, "</i> (Family)"),
      is.na(Genus) & is.na(Family) & is.na(Order) == FALSE ~ paste0("<i>", Order, "</i> (Order)"),
      is.na(Genus) & is.na(Family) & is.na(Order) & is.na(Class) == FALSE ~ paste0("<i>", Class, "</i> (Class)"),
      is.na(Genus) & is.na(Family) & is.na(Order) & is.na(Class) & is.na(Phylum) == FALSE ~ paste0("<i>", Phylum, "</i> (Phylum)"),
      TRUE ~ NA
    ),
    markerID = str_extract(markers, "\\d{1,2}$") |> as.numeric(),
  ) |>
  filter(feature != "unknown_genera")
cov_neg_lefse_res

cov_neg_lefse_gg <-
  cov_neg_lefse_res |>
  # filter(markerID %notin% c(18,20)) |>
  ggplot(mapping = aes(x = ef_lda, y = reorder(feature, -markerID), fill = enrich_group)) +
  geom_col(color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  # facet_grid(rows = vars(enrich_group), scales = "free_y", space = "free_y") +
  labs(x = "LDA score", fill = "Enrichment") +
  theme(
    axis.text.y = ggtext::element_markdown(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(-0.19, "lines"),
    legend.position = "right"
    # legend.position = 'bottom'
  )
cov_neg_lefse_gg

cov_neg_lefse_ggp <-
  cov_neg_lefse_res |>
  ggplot(mapping = aes(x = -log10(padj), y = reorder(feature, -markerID), fill = enrich_group)) +
  geom_point(mapping = aes(shape = enrich_group), size = 2, color = "white", show.legend = FALSE) +
  # facet_grid(rows = vars(enrich_group), scales = "free_y", space = "free_y") +
  scale_x_continuous(
    limits = c(0, 9),
    expand = expansion(mult = c(0, 0.1)),
    sec.axis = dup_axis()
  ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  scale_shape_manual(values = c(24, 23, 22)) +
  labs(x = "-log<sub>10</sub>(p<sub>adj</sub>)", fill = "Community", title = "SARS-CoV-2 Negative") +
  theme(
    axis.text.y = ggtext::element_markdown(),
    axis.title.y = element_blank(),
    axis.text.y.left = element_blank(),
    axis.title.x = ggtext::element_markdown(),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.line.x.bottom = element_blank(),
    axis.title.x.bottom = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(-0.22, "lines"),
    plot.background = element_blank(),
    legend.position = "none"
  )
cov_neg_lefse_ggp

# cov_neg_lefse_ggp2 <- cov_neg_lefse_ggp + add_logticks(side = "t", data = data.frame(x = NA, enrich_group = "Elderly"))

cov_neg_aligned_plots <- cowplot::align_plots(cov_neg_lefse_gg, cov_neg_lefse_ggp, align = "hv", axis = "tblr")
cov_neg_aligned_plots2 <- ggdraw(cov_neg_aligned_plots[[1]]) + draw_plot(cov_neg_aligned_plots[[2]])
cov_neg_aligned_plots2

## Plot relative abundance of significant taxa ----
# Combine into one and visualize
cov_neg_ra_box <-
  do.call(rbind, df_ps_filt_list) |>
  data.frame() |>
  # dplyr::filter(TaxLevelName %in% cov_all_lefse_res$feature) |>
  inner_join(cov_neg_lefse_res, by = "feature") |>
  filter(Abundance > 0) |>
  ggplot(mapping = aes(x = Abundance * 100, y = reorder(feature, -markerID))) +
  geom_boxplot(mapping = aes(fill = Community2), position = position_dodge(width = 0.9), outliers = FALSE, alpha = 0.5) +
  geom_point(mapping = aes(fill = Community2, shape = Community2), position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2)) +
  scale_x_continuous(
    transform = "log10",
    sec.axis = dup_axis(),
    labels = scales::number_format(accuracy = 0.1),
    limits = c(0.01, 100)
  ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  scale_shape_manual(values = c(24, 23, 22)) +
  labs(color = "Community", fill = "Community", shape = "Community", title = "Relative abundance", x = "Relative abundance (%)") +
  annotation_logticks(sides = "tb") +
  theme(
    axis.title.y = element_blank(),
    # axis.text.y = ggtext::element_markdown()
    axis.text.y = element_blank(),
    legend.position = "bottom",
    plot.margin = unit(c(10, 8, 7, 0), "pt"),
  )
cov_neg_ra_box

# Put figures side by side
fig3b <- ggarrange(
  cov_neg_aligned_plots2,
  cov_neg_ra_box,
  # align = 'h',
  # axis = 'tb',
  widths = c(0.9, 1)
  # heights = c(0.5, 1)
  # legend = "bottom"
)
fig3b

# ggsave(filename = 'processed_data/figures/version4/cov_neg_gg.png',
#        units = 'in', width = 15, height = 15, dpi = 320)

## Community comparisons covid positive ----
cov_pos_res <- run_lefse(
  ps_filt |> subset_samples(Covid == "Positive"),
  group = "Community2",
  multigrp_strat = TRUE, # enable multi-class mode
  strict = "1", # "1" = all‑against‑all (most stringent)
  kw_cutoff = 0.05,
  wilcoxon_cutoff = 0.05,
  lda_cutoff = 4
)
cov_pos_res

# Save the results
# saveRDS(
#   cov_pos_res,
#   file = "intermediate_data/lefse_results/cov_pos_lefse_results.rds"
# )

### Visualize ----
cov_pos_lefse_res <-
  cov_pos_res@marker_table |>
  data.frame() |>
  mutate(
    feature = str_remove_all(string = feature, pattern = "\\w_{2,3}")
  ) |>
  tidyr::separate(
    col = feature,
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
    sep = "\\|"
  ) |>
  mutate_if(is.character, ~ na_if(., "")) |>
  rownames_to_column(var = "markers") |>
  mutate(
    # Genus = if_else(str_detect(Genus,'^\\W')==T,NA,Genus),
    # Genus = if_else(str_detect(Genus,'')==T,NA,Genus),
    Genus = if_else(str_detect(Genus, "_$") == T, NA, Genus),
    # Family = if_else(str_detect(Genus,'^\\W')==T,NA,Family),
    feature = case_when(
      is.na(Genus) == FALSE ~ paste0("<i>", Genus, "</i> (Genus)"),
      is.na(Genus) & is.na(Family) == FALSE ~ paste0("<i>", Family, "</i> (Family)"),
      is.na(Genus) & is.na(Family) & is.na(Order) == FALSE ~ paste0("<i>", Order, "</i> (Order)"),
      is.na(Genus) & is.na(Family) & is.na(Order) & is.na(Class) == FALSE ~ paste0("<i>", Class, "</i> (Class)"),
      is.na(Genus) & is.na(Family) & is.na(Order) & is.na(Class) & is.na(Phylum) == FALSE ~ paste0("<i>", Phylum, "</i> (Phylum)"),
      TRUE ~ NA
    ),
    markerID = str_extract(markers, "\\d{1,2}$") |> as.numeric()
  )

cov_pos_lefse_gg <-
  cov_pos_lefse_res |>
  filter(markerID != 19) |>
  ggplot(mapping = aes(x = ef_lda, y = reorder(feature, -markerID), fill = enrich_group)) +
  geom_col(color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")[2:3]) +
  labs(x = "LDA score", fill = "Enrichment") +
  theme(
    axis.text.y = ggtext::element_markdown(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(-0.19, "lines"),
    legend.position = "right"
  )
cov_pos_lefse_gg

cov_pos_lefse_ggp <-
  cov_pos_lefse_res |>
  ggplot(mapping = aes(x = -log10(padj), y = reorder(feature, -markerID), fill = enrich_group)) +
  geom_point(mapping = aes(shape = enrich_group), size = 2, color = "white", show.legend = FALSE) +
  scale_x_continuous(
    limits = c(0, 9),
    expand = expansion(mult = c(0, 0.1)),
    sec.axis = dup_axis()
  ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")[2:3]) +
  scale_shape_manual(values = c(23, 22)) +
  labs(x = "-log<sub>10</sub>(p<sub>adj</sub>)", fill = "Community", title = "SARS-CoV-2 Positive") +
  theme(
    axis.text.y = ggtext::element_markdown(),
    axis.title.y = element_blank(),
    axis.text.y.left = element_blank(),
    axis.title.x = ggtext::element_markdown(),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.line.x.bottom = element_blank(),
    axis.title.x.bottom = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(-0.22, "lines"),
    plot.background = element_blank(),
    legend.position = "none"
  )
cov_pos_lefse_ggp

# cov_pos_lefse_ggp2 <- cov_pos_lefse_ggp + add_logticks(side = "t", data = data.frame(x = NA, enrich_group = "Marginalized"))

cov_pos_aligned_plots <- cowplot::align_plots(cov_pos_lefse_gg, cov_pos_lefse_ggp, align = "hv", axis = "tblr")
cov_pos_aligned_plots2 <- ggdraw(cov_pos_aligned_plots[[1]]) + draw_plot(cov_pos_aligned_plots[[2]])
cov_pos_aligned_plots2

# Combine into one and visualize
cov_pos_ra_box <-
  do.call(rbind, df_ps_filt_list) |>
  data.frame() |>
  # dplyr::filter(TaxLevelName %in% cov_all_lefse_res$feature) |>
  inner_join(cov_pos_lefse_res, by = "feature") |>
  filter(Abundance > 0) |>
  ggplot(mapping = aes(x = Abundance * 100, y = reorder(feature, -markerID))) +
  geom_boxplot(mapping = aes(fill = Community2), position = position_dodge(width = 0.9), outliers = FALSE, alpha = 0.5) +
  geom_point(mapping = aes(fill = Community2, shape = Community2), position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2)) +
  scale_x_continuous(
    transform = "log10",
    sec.axis = dup_axis(),
    labels = scales::number_format(accuracy = 0.1),
    limits = c(0.01, 100)
  ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  scale_shape_manual(values = c(24, 23, 22)) +
  labs(color = "Community", fill = "Community", title = "Relative abundance", x = "Relative abundance (%)") +
  annotation_logticks(sides = "tb") +
  theme(
    axis.title.y = element_blank(),
    # axis.text.y = ggtext::element_markdown(),
    axis.text.y = element_blank(),
    legend.position = "bottom",
    plot.margin = unit(c(10, 8, 7, 0), "pt"),
  )
cov_pos_ra_box

options(
  vsc.dev.args = list(
    width = 3000,
    height = 1500,
    pointsize = 12,
    res = 300
  )
)

# Put figures side by side
fig3c <- ggarrange(
  cov_pos_aligned_plots2,
  cov_pos_ra_box,
  # align = 'h',
  # axis = 'tb',
  widths = c(0.9, 1)
  # heights = c(0.5, 1)
  # legend = "bottom"
)
fig3c

# ggsave(filename = 'processed_data/figures/version4/cov_pos_gg.png',
#        units = 'in', width = 15, height = 8, dpi = 320)

options(
  vsc.dev.args = list(
    width = 5000,
    height = 3500,
    pointsize = 12,
    res = 320
  )
)
# Combine figure 3 plots ----
fig3a +
  fig3b +
  fig3c +
  plot_annotation(tag_level = "A") +
  plot_layout(heights = c(0.2, 1, 0.4))

design <- "
  AABB
  CCBB
"

com_ova_lefse_gg +
  cov_neg_aligned_plots2 +
  cov_pos_aligned_plots2 +
  plot_annotation(tag_level = "A") +
  plot_layout(heights = c(0.2, 0.4, 1), design = design)

# ggsave(
#   filename = "processed_data/figures/version4/figure3.png",
#   units = "in", width = 12, height = 10, dpi = 320
# )

## Elderly covid comparisons ----
set.seed(1234)
eld_cov_res <- run_lefse(
  ps_filt |> subset_samples(Community2 == "Elderly" & Covid != "NA"),
  group = "Covid",
  subgroup = "SiteID",
  multigrp_strat = TRUE, # enable multi-class mode
  strict = "1", # "1" = all‑against‑all (most stringent)
  kw_cutoff = 0.05,
  wilcoxon_cutoff = 0.05,
  lda_cutoff = 4
)
eld_cov_res

# Save the results
# saveRDS(
#   eld_cov_res,
#   file = "intermediate_data/lefse_results/eld_cov_lefse_results.rds"
# )

### Visualize ----
eld_cov_lefse_res <- eld_cov_res@marker_table |>
  data.frame() |>
  mutate(
    feature = str_remove_all(string = feature, pattern = "\\w_{2,3}")
  ) |>
  tidyr::separate(
    col = feature,
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
    sep = "\\|"
  ) |>
  mutate_if(is.character, ~ na_if(., "")) |>
  rownames_to_column(var = "markers") |>
  mutate(
    # Genus = if_else(str_detect(Genus,'^\\W')==T,NA,Genus),
    # Genus = if_else(str_detect(Genus,'')==T,NA,Genus),
    Genus = if_else(str_detect(Genus, "_$") == T, NA, Genus),
    # Family = if_else(str_detect(Genus,'^\\W')==T,NA,Family),
    feature = case_when(
      is.na(Genus) == FALSE ~ paste0("<i>", Genus, "</i> (Genus)"),
      is.na(Genus) & is.na(Family) == FALSE ~ paste0("<i>", Family, "</i> (Family)"),
      is.na(Genus) & is.na(Family) & is.na(Order) == FALSE ~ paste0("<i>", Order, "</i> (Order)"),
      is.na(Genus) & is.na(Family) & is.na(Order) & is.na(Class) == FALSE ~ paste0("<i>", Class, "</i> (Class)"),
      is.na(Genus) & is.na(Family) & is.na(Order) & is.na(Class) & is.na(Phylum) == FALSE ~ paste0("<i>", Phylum, "</i> (Phylum)"),
      TRUE ~ NA
    ),
    markerID = str_extract(markers, "\\d{1,2}$") |> as.numeric(),
    up_down = if_else(
      markerID %in% c(35),
      "Depleted",
      "Enriched"
    ),
    grouping = paste0(up_down, " in ", enrich_group),
    grouping = factor(
      grouping,
      levels = c(
        "Depleted in Negative",
        "Enriched in Negative",
        "Depleted in Positive",
        "Enriched in Positive"
      )
    )
  ) |>
  filter(markerID != 43)

eld_cov_lefse_gg <-
  eld_cov_lefse_res |>
  ggplot(mapping = aes(x = ef_lda, y = reorder(feature, -markerID))) +
  # geom_col(color = "black") +
  geom_col(
    mapping = aes(fill = enrich_group),
    color = "black"
  ) +
  # facet_grid(rows = vars(enrich_group), scales = "free_y", space = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  # scale_fill_manual(values = c(viridis::inferno(n = 3, begin = 0.1, end = 0.9)[1], rep(viridis::inferno(n = 3, begin = 0.1, end = 0.9)[2], times = 2))) +
  scale_fill_manual(values = covid_colors) +
  labs(x = "LDA score", fill = "COVID-19") +
  theme(
    axis.text.y = ggtext::element_markdown(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(-0.19, "lines"),
    legend.position = "right"
  )
eld_cov_lefse_gg
# ggsave(
#   filename = 'processed_data/figures/figure4.png',
#   units = 'in', width = 8, height = 12, dpi = 320
# )

eld_cov_lefse_res |>
  ggplot(mapping = aes(x = -log10(padj), y = reorder(feature, -markerID), fill = enrich_group)) +
  geom_point(size = 2, color = "white", shape = 24, show.legend = FALSE) +
  # facet_grid(rows = vars(enrich_group), scales = "free_y", space = "free_y") +
  scale_x_continuous(
    limits = c(0, 7),
    expand = expansion(mult = c(0, 0.1)),
    sec.axis = dup_axis()
  ) +
  # scale_fill_manual(values = viridis::inferno(n = 3, begin = 0.1, end = 0.9)) +
  scale_fill_manual(values = covid_colors) +
  labs(x = "-log<sub>10</sub>(p<sub>adj</sub>)", fill = "SARS-CoV-2", title = "Elderly") +
  theme(
    axis.text.y = ggtext::element_markdown(),
    axis.title.y = element_blank(),
    # axis.text.y.left = element_blank(),
    axis.text.y.left = ggtext::element_markdown(),
    axis.title.x = ggtext::element_markdown(),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.line.x.bottom = element_blank(),
    axis.title.x.bottom = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(-0.22, "lines"),
    plot.background = element_blank(),
    legend.position = "none"
  ) -> fig4_p
fig4_p
# fig4_p2 <- fig4_p + add_logticks(side = "t", data = data.frame(x = NA, enrich_group = "Negative"))

aligned_plots <- cowplot::align_plots(eld_cov_lefse_gg, fig4_p, align = "hv", axis = "tblr")
eld_aligned_plots2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
eld_aligned_plots2

## Plot relative abundance of significant taxa ----
# Combine into one and visualize
eld_cov_ra_box <-
  do.call(rbind, df_ps_filt_list) |>
  data.frame() |>
  # dplyr::filter(TaxLevelName %in% cov_all_lefse_res$feature) |>
  inner_join(eld_cov_lefse_res, by = "feature") |>
  filter(Abundance > 0 & Community2 == "Elderly", Covid != "NA") |>
  ggplot(mapping = aes(x = Abundance * 100, y = reorder(feature, -markerID))) +
  geom_boxplot(mapping = aes(fill = Covid), position = position_dodge(width = 0.9, preserve = "single"), outliers = FALSE, alpha = 0.5) +
  geom_point(mapping = aes(fill = Covid), shape = 24, position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2)) +
  # stat_summary(aes(color = Covid), size = 5, geom = "point", position = position_dodge(width = 0.9)) +
  # stat_summary(aes(color = Covid), width = 0.5, geom = "errorbar", position = position_dodge(width = 0.9)) +
  scale_x_continuous(
    transform = "log10",
    sec.axis = dup_axis(),
    labels = scales::number_format(accuracy = 0.1),
    limits = c(0.01, 100)
  ) +
  scale_fill_manual(values = viridis::inferno(n = 3, begin = 0.1, end = 0.9)[1:2]) +
  scale_color_manual(values = viridis::inferno(n = 3, begin = 0.1, end = 0.9)[1:2]) +
  labs(color = "Community", fill = "Community", title = "Relative abundance", x = "Relative abundance (%)") +
  annotation_logticks(sides = "tb") +
  theme(
    axis.title.y = element_blank(),
    # axis.text.y = ggtext::element_markdown()
    axis.text.y = element_blank(),
    legend.position = "bottom",
    plot.margin = unit(c(10, 8, 7, 0), "pt"),
  )
eld_cov_ra_box

# Put figures side by side
fig4a <- ggarrange(
  eld_aligned_plots2,
  eld_cov_ra_box,
  # align = 'h',
  # axis = 'tb',
  widths = c(0.9, 1)
  # heights = c(0.5, 1)
  # legend = "bottom"
)
fig4a

# ggsave(
#   filename = 'processed_data/figures/version4/eld_cov_.png',
#   units = 'in', width = 15, height = 15, dpi = 320
# )

## Marginalized covid comparisons ----
mar_cov_res <- run_lefse(
  ps_filt |> subset_samples(Community2 == "Marginalized" & Covid != "NA"),
  group = "Covid",
  subgroup = "SiteID",
  multigrp_strat = TRUE, # enable multi-class mode
  strict = "2", # "1" = all‑against‑all (most stringent)
  kw_cutoff = 0.05,
  wilcoxon_cutoff = 0.05,
  lda_cutoff = 4
)
mar_cov_res

# Save the results
# saveRDS(
#   mar_cov_res,
#   file = "intermediate_data/lefse_results/mar_cov_lefse_results.rds"
# )

### Visualize ----
mar_lefse_res <-
  mar_cov_res@marker_table |>
  data.frame() |>
  mutate(
    feature = str_remove_all(string = feature, pattern = "\\w_{2,3}")
  ) |>
  tidyr::separate(
    col = feature,
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
    sep = "\\|"
  ) |>
  mutate_if(is.character, ~ na_if(., "")) |>
  rownames_to_column(var = "markers") |>
  mutate(
    # Genus = if_else(str_detect(Genus,'^\\W')==T,NA,Genus),
    # Genus = if_else(str_detect(Genus,'')==T,NA,Genus),
    Genus = if_else(str_detect(Genus, "_$") == T, NA, Genus),
    # Family = if_else(str_detect(Genus,'^\\W')==T,NA,Family),
    feature = case_when(
      is.na(Genus) == FALSE ~ paste0("<i>", Genus, "</i> (Genus)"),
      is.na(Genus) & is.na(Family) == FALSE ~ paste0("<i>", Family, "</i> (Family)"),
      is.na(Genus) & is.na(Family) & is.na(Order) == FALSE ~ paste0("<i>", Order, "</i> (Order)"),
      is.na(Genus) & is.na(Family) & is.na(Order) & is.na(Class) == FALSE ~ paste0("<i>", Class, "</i> (Class)"),
      is.na(Genus) & is.na(Family) & is.na(Order) & is.na(Class) & is.na(Phylum) == FALSE ~ paste0("<i>", Phylum, "</i> (Phylum)"),
      TRUE ~ NA
    ),
    markerID = str_extract(markers, "\\d{1,2}$") |> as.numeric()
  )

mar_lefse_gg <-
  mar_lefse_res |>
  ggplot(mapping = aes(x = ef_lda, y = reorder(feature, -markerID), fill = enrich_group)) +
  geom_col(color = "black") +
  # facet_grid(rows = vars(enrich_group), scales = "free_y", space = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  # scale_fill_manual(values = viridis::inferno(n = 3, begin = 0.1, end = 0.9)) +
  scale_fill_manual(values = covid_colors) +
  labs(x = "LDA score", fill = "COVID-19") +
  theme(
    axis.text.y = ggtext::element_markdown(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(-0.19, "lines"),
    legend.position = "right"
  )
mar_lefse_gg

mar_lefse_ggp <-
  mar_lefse_res |>
  ggplot(mapping = aes(x = -log10(padj), y = reorder(feature, -markerID), fill = enrich_group)) +
  geom_point(size = 2, color = "white", shape = 23, show.legend = FALSE) +
  # facet_grid(rows = vars(enrich_group), scales = "free_y", space = "free_y") +
  scale_x_continuous(
    limits = c(0, 7),
    expand = expansion(mult = c(0, 0.1)),
    sec.axis = dup_axis()
  ) +
  # scale_fill_manual(values = viridis::inferno(n = 3, begin = 0.1, end = 0.9)) +
  scale_fill_manual(values = covid_colors) +
  labs(x = "-log<sub>10</sub>(p<sub>adj</sub>)", fill = "SARS-CoV-2", title = "Marginalized") +
  theme(
    axis.text.y = ggtext::element_markdown(),
    axis.title.y = element_blank(),
    axis.text.y.left = element_blank(),
    axis.title.x = ggtext::element_markdown(),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.line.x.bottom = element_blank(),
    axis.title.x.bottom = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(-0.22, "lines"),
    plot.background = element_blank()
  )

mar_lefse_ggp2 <- mar_lefse_ggp + add_logticks(side = "t", data = data.frame(x = NA, enrich_group = "Negative"))

mar_aligned_plots <- cowplot::align_plots(mar_lefse_gg, mar_lefse_ggp, align = "hv", axis = "tblr")
mar_aligned_plots2 <- ggdraw(mar_aligned_plots[[1]]) + draw_plot(mar_aligned_plots[[2]])
mar_aligned_plots2

## Plot relative abundance of significant taxa ----
# Combine into one and visualize
mar_cov_ra_box <-
  do.call(rbind, df_ps_filt_list) |>
  data.frame() |>
  # dplyr::filter(TaxLevelName %in% cov_all_lefse_res$feature) |>
  inner_join(mar_lefse_res, by = "feature") |>
  filter(Abundance > 0 & Community2 == "Marginalized") |>
  ggplot(mapping = aes(x = Abundance * 100, y = reorder(feature, -markerID))) +
  geom_boxplot(mapping = aes(fill = Covid), position = position_dodge(width = 0.9, preserve = "single"), outliers = FALSE, alpha = 0.5) +
  geom_point(mapping = aes(fill = Covid), shape = 23, position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2)) +
  scale_x_continuous(
    transform = "log10",
    sec.axis = dup_axis(),
    labels = scales::number_format(accuracy = 0.1),
    limits = c(0.01, 100)
  ) +
  scale_fill_manual(values = viridis::inferno(n = 3, begin = 0.1, end = 0.9)[1:2]) +
  scale_color_manual(values = viridis::inferno(n = 3, begin = 0.1, end = 0.9)[1:2]) +
  labs(color = "Community", fill = "Community", title = "Relative abundance", x = "Relative abundance (%)") +
  annotation_logticks(sides = "tb") +
  theme(
    axis.title.y = element_blank(),
    # axis.text.y = ggtext::element_markdown()
    axis.text.y = element_blank(),
    legend.position = "bottom",
    plot.margin = unit(c(10, 8, 7, 0), "pt"),
  )
mar_cov_ra_box

# Put figures side by side
fig4b <- ggarrange(
  mar_aligned_plots2,
  mar_cov_ra_box,
  # align = 'h',
  # axis = 'tb',
  widths = c(0.9, 1)
  # heights = c(0.5, 1)
  # legend = "bottom"
)
fig4b

# ggsave(
#   filename = 'processed_data/figures/version4/mar_cov_gg.png',
#   units = 'in', width = 15, height = 6, dpi = 320
# )

## Minority covid comparisons ----
# min_cov_res <- run_lefse(
#   ps_filt |> subset_samples(Community2 == "Minority"),
#   group = "Covid",
#   subgroup = "SiteID",
#   multigrp_strat = TRUE, # enable multi-class mode
#   strict = "1", # "1" = all‑against‑all (most stringent)
#   kw_cutoff = 0.05,
#   wilcoxon_cutoff = 0.05,
#   lda_cutoff = 4
# )
# min_cov_res # NO SIGNIFICANT TAXA

## Combine figures ----
fig4a +
  fig4b +
  plot_annotation(tag_level = "A") +
  plot_layout(heights = c(1, 0.25), ncol = 1)

eld_aligned_plots2 +
  mar_aligned_plots2 +
  plot_annotation(tag_level = "A") +
  plot_layout(heights = c(1, 0.4), ncol = 1)

# ggsave(
#   filename = "processed_data/figures/version4/figure4.png",
#   units = "in", width = 8, height = 10, dpi = 320
# )

# Gut microbe filtering ----
## Load GM Repo taxa data ----
query <- POST("https://gmrepo.humangut.info/api/get_all_gut_microbes", body = NULL, encode = "json")
retrieved_contents <- content(query)
gm_repo_data <- fromJSON(xml_text(retrieved_contents))

## Extract genus data (present in ≥0.1% samples) ----
gm_gen_data <-
  gm_repo_data$all_genus |>
  filter(pct_of_all_samples >= 0.1)
gm_gen_data$name |>
  unique() |>
  length()
# 485 genera remaining

# gm_gen_data <-
#   gm_repo_data$all_genus

## Filter genera by 0.1% abundance filter ----
# Aggregate taxa at the genus level
ps_gen <-
  ps |>
  microbiome::aggregate_taxa(level = "Genus") |>
  microbiome::transform(transform = "compositional")

# Save data
# saveRDS(
#   ps_gen,
#   file = "intermediate_data/downstream/ps_gen.rds"
# )

# Filter by relative abundance
filt_gens <- filter_taxa(
  ps_gen,
  function(x) mean(x) > 0.001,
  prune = FALSE
)

ps_filt_gen <-
  prune_taxa(
    filt_gens,
    ps_gen
  )
ps_filt_gen

# Save data
# saveRDS(
#   ps_filt_gen,
#   file = "intermediate_data/downstream/ps_filt_gen.rds"
# )

## Fix taxa table ----
# Extract tax_table from ps object
tax_table_unfilt <- tax_table(ps_gen) |> data.frame()
tax_table_filt <- tax_table(ps_filt_gen) |> data.frame()

### Create a fucntion to adjust taxa names ----
extract_genus <- function(taxonomic_string) {
  # Handle NULL case
  if (is.null(taxonomic_string)) {
    return(NA_character_)
  }

  # Create a result vector of NA values with the same length as input
  result <- rep(NA_character_, length(taxonomic_string))

  # Process only non-NA values
  valid_indices <- which(!is.na(taxonomic_string))

  if (length(valid_indices) > 0) {
    valid_strings <- taxonomic_string[valid_indices]

    # Find which valid strings have brackets
    has_brackets <- str_detect(valid_strings, "\\[.*?\\]")

    # For strings with brackets
    if (any(has_brackets)) {
      bracket_indices <- valid_indices[has_brackets]
      result[bracket_indices] <- str_remove_all(
        str_extract(taxonomic_string[bracket_indices], "\\[(.*?)\\]"),
        "\\[|\\]"
      )
    }

    # For strings without brackets
    if (any(!has_brackets)) {
      no_bracket_indices <- valid_indices[!has_brackets]
      result[no_bracket_indices] <- str_extract(taxonomic_string[no_bracket_indices], "^[^ ]+")
    }
  }

  return(result)
}

### Revise taxa naming convention to match GM repo ----
taxrev_unfilt <-
  tax_table_unfilt |>
  mutate(
    Genus = case_when(
      str_detect(string = Genus, pattern = "group") & str_detect(string = Genus, pattern = "aceae") ~ paste0("unclassified ", extract_genus(Genus)),
      str_detect(string = Genus, pattern = "group") ~ extract_genus(Genus),
      # str_detect(string=Genus, pattern='aceae') ~ str_split_i(string=Genus, pattern=' ', i=1),
      str_detect(string = Genus, pattern = "sensu stricto") ~ str_split_i(string = Genus, pattern = " ", i = 1),
      str_detect(string = Genus, pattern = "aceae") ~ paste0("unclassified ", (str_split_i(string = Genus, pattern = " ", i = 1))),
      str_detect(string = Genus, pattern = "Prevotella_\\d") ~ "Prevotella",
      str_detect(string = Genus, pattern = "\\d+") ~ paste0("unclassified ", Family),
      is.na(Family) == FALSE & is.na(Genus) ~ paste0("unclassified ", Family),
      is.na(Order) == FALSE & is.na(Family) & is.na(Genus) ~ paste0("unclassified ", Order),
      TRUE ~ Genus
    )
  )
# taxrev_unfilt |> View()

taxrev_filt <-
  tax_table_filt |>
  mutate(
    Genus = case_when(
      str_detect(string = Genus, pattern = "group") & str_detect(string = Genus, pattern = "aceae") ~ paste0("unclassified ", extract_genus(Genus)),
      str_detect(string = Genus, pattern = "group") ~ extract_genus(Genus),
      # str_detect(string=Genus, pattern='aceae') ~ str_split_i(string=Genus, pattern=' ', i=1),
      str_detect(string = Genus, pattern = "sensu stricto") ~ str_split_i(string = Genus, pattern = " ", i = 1),
      str_detect(string = Genus, pattern = "aceae") ~ paste0("unclassified ", (str_split_i(string = Genus, pattern = " ", i = 1))),
      str_detect(string = Genus, pattern = "Prevotella_\\d") ~ "Prevotella",
      str_detect(string = Genus, pattern = "\\d+") ~ paste0("unclassified ", Family),
      is.na(Family) == FALSE & is.na(Genus) ~ paste0("unclassified ", Family),
      is.na(Order) == FALSE & is.na(Family) & is.na(Genus) ~ paste0("unclassified ", Order),
      TRUE ~ Genus
    )
  )

# Copy ps object
ps_taxrev <- ps_gen
ps_filt_taxrev <- ps_filt_gen

# Apply the new tax_table to ps object
tax_table(ps_taxrev) <- taxrev_unfilt |> as.matrix()
tax_table(ps_filt_taxrev) <- taxrev_filt |> as.matrix()

## Filter taxa by GM repo ----
ps_gm <- prune_taxa(
  gm_gen_data$name,
  ps_taxrev
)
ps_gm

# # Save data
# saveRDS(
#   ps_gm,
#   file = "intermediate_data/downstream/ps_gm.rds"
# )

ps_filt_gm <- prune_taxa(
  gm_gen_data$name,
  ps_filt_taxrev
)
ps_filt_gm

# Save data
# saveRDS(
#   ps_filt_gm,
#   file = "intermediate_data/downstream/ps_filt_gm.rds"
# )

filt_comp_df <-
  otu_table(ps_gm) |>
  t() |>
  rowSums() |>
  data.frame() |>
  dplyr::rename(total_ra = 1) |>
  mutate(abund_filt = "Unfiltered") |>
  rownames_to_column(var = "SampleID") |>
  full_join(
    otu_table(ps_filt_gm) |>
      t() |>
      rowSums() |>
      data.frame() |>
      dplyr::rename(total_ra = 1) |>
      mutate(abund_filt = "> 0.1%") |>
      rownames_to_column(var = "SampleID")
  ) |>
  mutate(
    abund_filt = factor(abund_filt, levels = c("Unfiltered", "> 0.1%"))
  ) |>
  dplyr::full_join(
    sample_data(ps_filt_taxrev) |>
      data.frame() |>
      rownames_to_column(var = "SampleID"),
    by = "SampleID"
  )
filt_comp_df

# Save data
# xlsx::write.xlsx(
#   filt_comp_df,
#   file = "processed_data/summary_tables/filt_comp_df.xlsx"
# )

fig5a <-
  filt_comp_df |>
  ggplot(mapping = aes(x = Community2, y = total_ra * 100, fill = abund_filt, group = abund_filt)) +
  stat_summary(geom = "col", color = "black", position = position_dodge(width = 0.9)) +
  stat_summary(geom = "errorbar", color = "black", width = 0.25, position = position_dodge(width = 0.9)) +
  stat_pwc(method = "wilcox_test", hide.ns = FALSE) +
  scale_fill_manual(values = viridis::rocket(n = 3, begin = 0.8, end = 0.45)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Community", y = "Total<br>relative abundance<br>(%)", fill = "Relative abundance filtering") +
  theme(
    legend.position = "bottom"
  )
fig5a

fig5b <-
  filt_comp_df |>
  filter(Covid != "NA") |>
  ggplot(mapping = aes(x = Covid, y = total_ra * 100, fill = abund_filt, group = abund_filt)) +
  stat_summary(geom = "col", color = "black", position = position_dodge(width = 0.9)) +
  stat_summary(geom = "errorbar", color = "black", width = 0.25, position = position_dodge(width = 0.9)) +
  stat_pwc(method = "wilcox_test", hide.ns = FALSE) +
  lemon::facet_rep_wrap(~Community2) +
  scale_fill_manual(values = viridis::rocket(n = 3, begin = 0.8, end = 0.45)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "SARS-CoV-2", y = "Total<br>relative abundance<br>(%)", fill = "Relative abundance filtering") +
  theme(
    axis.text.x = ggtext::element_markdown(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    panel.spacing = unit(-1.5, "lines")
  )
fig5b

fig5a +
  fig5b +
  plot_annotation(tag_level = "A") +
  plot_layout(guides = "collect", ncol = 1) &
  theme(
    axis.title.y = ggtext::element_markdown(size = 12),
    legend.position = "bottom",
    legend.justification = "center"
  )

# ggsave(
#   filename = "processed_data/figures/version4/figure5.png",
#   units = "in", width = 8, height = 6, dpi = 320
# )

## Create a histogram ----
fig5b2 <-
  filt_comp_df |>
  filter(Covid != "NA") |>
  ggplot(aes(x = total_ra * 100, fill = abund_filt)) +
  # geom_histogram(mapping=aes(y=after_stat(count/sum(count))), binwidth = 10, color = 'black') +
  # geom_histogram(binwidth = 10, color = 'black') +
  geom_histogram_pattern(mapping = aes(pattern = abund_filt), binwidth = 10, color = "black", pattern_fill = "white", pattern_key_scale_factor = 0.4, pattern_res = 320, pattern_colour = "white") +
  lemon::facet_rep_wrap(
    ~abund_filt,
    nrow = 1,
    repeat.tick.labels = TRUE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    sec.axis = sec_axis(transform = ~ .x / 81 * 100, name = "% of samples")
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  scale_fill_manual(values = viridis::rocket(n = 3, begin = 0.8, end = 0.45)) +
  scale_pattern_manual(values = c("none", "stripe")) +
  labs(
    x = "Total relative abundance (%)",
    y = "# of samples",
    fill = "Relative abundance filtering",
    pattern = "Relative abundance filtering"
    # fill = 'Relative<br>abundance<br>filtering',
    # pattern = 'Relative<br>abundance<br>filtering'
  ) +
  theme(
    axis.text.x = ggtext::element_markdown(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = ggtext::element_markdown(),
    panel.spacing = unit(1, "lines")
  )
fig5b2

## Combine previous 5A and 5B ----
fig5a2 <-
  filt_comp_df |>
  rbind(
    filt_comp_df |>
      mutate(Covid = "Overall")
  ) |>
  filter(Covid != "NA") |>
  mutate(Covid = factor(Covid, levels = c("Overall", "Negative", "Positive"))) |>
  ggplot(mapping = aes(x = Covid, y = total_ra * 100, fill = abund_filt, group = abund_filt)) +
  stat_summary(geom = "col_pattern", mapping = aes(pattern = abund_filt), color = "black", position = position_dodge(width = 0.9), pattern_fill = "white", pattern_key_scale_factor = 0.4, pattern_res = 320, pattern_colour = "white") +
  stat_summary(geom = "errorbar", color = "black", width = 0.25, position = position_dodge(width = 0.9)) +
  stat_pwc(method = "wilcox_test", hide.ns = FALSE) +
  lemon::facet_rep_wrap(~Community2) +
  scale_fill_manual(values = viridis::rocket(n = 3, begin = 0.8, end = 0.45)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_pattern_manual(values = c("none", "stripe")) +
  labs(
    x = "SARS-CoV-2",
    y = "Total relative<br>abundance (%)",
    fill = "Relative abundance filtering",
    pattern = "Relative abundance filtering"
  ) +
  theme(
    axis.text.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    panel.spacing = unit(-1.5, "lines")
  )
fig5a2

## Combine figures ----
# fig5a2 +
#   fig5b2 +
#   plot_annotation(tag_level = "A") +
#   plot_layout(guides = "collect", ncol = 1) &
#   theme(
#     axis.title.y = ggtext::element_markdown(size = 12),
#     legend.position = "bottom",
#     legend.justification = "center"
#   )

ggarrange(
  fig5a2,
  fig5b2 + theme(legend.position = "none"),
  ncol = 1,
  labels = "AUTO"
)

# Save figure
# ggsave(
#   filename = "processed_data/figures/version4/figure5.png",
#   units = "in", width = 8, height = 6, dpi = 320
# )

# Create summary tables for the LEfSe results ----
## All communities ----
# Transpose the ASV table
otu_table(ps_filt) <- otu_table(ps_filt) |> t()

# Normalize the ASV table using CPM (Counts Per Million)
ps_cpm <- microbiomeMarker::normalize(ps_filt, method = "CPM")

# Save data
# saveRDS(
#   ps_cpm,
#   file = "intermediate_data/downstream/ps_cpm.rds"
# )

# Summarize the CPM data
overall_cpm_summary_list <- list()
taxa_levels <- c("Order", "Family", "Genus")

for (i in 1:length(taxa_levels)) {
  overall_cpm_summary_list[[i]] <-
    ps2df(ps_cpm, column.is.taxa = FALSE) |>
    summarize(
      avg_cpm = mean(Abundance, na.rm = TRUE),
      se_cpm = sd(Abundance, na.rm = TRUE) / sqrt(n()),
      .by = c(taxa_levels[i], "Community2")
    ) |>
    dplyr::rename(taxon = taxa_levels[i]) |>
    mutate(
      feature = paste0("<i>", taxon, "</i> (", taxa_levels[i], ")"),
    )
}

# Combine the summaries into one data frame
overall_cpm_summary <-
  do.call(
    rbind,
    overall_cpm_summary_list
  )
overall_cpm_summary

# Reformat the summary df
overall_cpm_summary_2 <-
  overall_cpm_summary |>
  pivot_wider(
    names_from = Community2,
    values_from = c(avg_cpm, se_cpm),
    names_sep = "_"
  ) |>
  select(-taxon) |>
  merge(
    cov_all_lefse_res |>
      select(markerID, feature, ef_lda, padj),
    by = "feature"
  )

# Create a summary table using gt
library(gt)
overall_cpm_summary_gt <-
  overall_cpm_summary_2 |>
  arrange(markerID) |>
  gt() |>
  cols_hide(columns = "markerID") |>
  tab_spanner(
    label = "Elderly (CPM)",
    columns = c(avg_cpm_Elderly, se_cpm_Elderly)
  ) |>
  tab_spanner(
    label = "Marginalized (CPM)",
    columns = c(avg_cpm_Marginalized, se_cpm_Marginalized)
  ) |>
  tab_spanner(
    label = "Minority (CPM)",
    columns = c(avg_cpm_Minority, se_cpm_Minority)
  ) |>
  fmt_markdown(columns = feature) |>
  fmt_number(
    columns = c(avg_cpm_Elderly, avg_cpm_Marginalized, avg_cpm_Minority, se_cpm_Elderly, se_cpm_Marginalized, se_cpm_Minority),
    decimals = 2,
    scale_by = 100
  ) |>
  fmt_number(
    columns = c(ef_lda),
    decimals = 2
  ) |>
  fmt_scientific(
    columns = c(padj),
    decimals = 2
  ) |>
  cols_label(
    feature = "Taxon (level)",
    avg_cpm_Elderly = "Average",
    avg_cpm_Marginalized = "Average",
    avg_cpm_Minority = "Average",
    se_cpm_Elderly = "SE",
    se_cpm_Marginalized = "SE",
    se_cpm_Minority = "SE",
    ef_lda = "LDA score",
    padj = "Adjusted p-value"
  )
overall_cpm_summary_gt

# Save the summary table as an word file
# gt::gtsave(
#   overall_cpm_summary_gt,
#   filename = "processed_data/summary_tables/overall_cpm_summary_table.docx"
# )

## Covid negative all communities ----
# Summarize the CPM data
cov_neg_cpm_summary_list <- list()
taxa_levels <- c("Phylum", "Class", "Order", "Family", "Genus")

for (i in 1:length(taxa_levels)) {
  cov_neg_cpm_summary_list[[i]] <-
    ps2df(ps_cpm, column.is.taxa = FALSE) |>
    filter(Covid == "Negative") |>
    summarize(
      avg_cpm = mean(Abundance, na.rm = TRUE),
      se_cpm = sd(Abundance, na.rm = TRUE) / sqrt(n()),
      .by = c(taxa_levels[i], "Community2")
    ) |>
    dplyr::rename(taxon = taxa_levels[i]) |>
    mutate(
      feature = paste0("<i>", taxon, "</i> (", taxa_levels[i], ")"),
    )
}

# Combine the summaries into one data frame
cov_neg_cpm_summary <-
  do.call(
    rbind,
    cov_neg_cpm_summary_list
  )
cov_neg_cpm_summary

# Reformat the summary df
cov_neg_cpm_summary_2 <-
  cov_neg_cpm_summary |>
  pivot_wider(
    names_from = Community2,
    values_from = c(avg_cpm, se_cpm),
    names_sep = "_"
  ) |>
  select(-taxon) |>
  merge(
    cov_neg_lefse_res |>
      select(markerID, feature, ef_lda, padj),
    by = "feature"
  )

# Create a summary table using gt
cov_neg_cpm_summary_gt <-
  cov_neg_cpm_summary_2 |>
  arrange(markerID) |>
  gt() |>
  cols_hide(columns = "markerID") |>
  tab_spanner(
    label = "Elderly (CPM)",
    columns = c(avg_cpm_Elderly, se_cpm_Elderly)
  ) |>
  tab_spanner(
    label = "Marginalized (CPM)",
    columns = c(avg_cpm_Marginalized, se_cpm_Marginalized)
  ) |>
  tab_spanner(
    label = "Minority (CPM)",
    columns = c(avg_cpm_Minority, se_cpm_Minority)
  ) |>
  fmt_markdown(columns = feature) |>
  fmt_number(
    columns = c(avg_cpm_Elderly, avg_cpm_Marginalized, avg_cpm_Minority, se_cpm_Elderly, se_cpm_Marginalized, se_cpm_Minority),
    decimals = 2,
    scale_by = 100
  ) |>
  fmt_number(
    columns = c(ef_lda),
    decimals = 2
  ) |>
  fmt_scientific(
    columns = c(padj),
    decimals = 2
  ) |>
  cols_label(
    feature = "Taxon (level)",
    avg_cpm_Elderly = "Average",
    avg_cpm_Marginalized = "Average",
    avg_cpm_Minority = "Average",
    se_cpm_Elderly = "SE",
    se_cpm_Marginalized = "SE",
    se_cpm_Minority = "SE",
    ef_lda = "LDA score",
    padj = "Adjusted p-value"
  )
cov_neg_cpm_summary_gt

# Save the summary table as an word file
# gt::gtsave(
#   cov_neg_cpm_summary_gt,
#   filename = "processed_data/summary_tables/cov_neg_cpm_summary_table.docx"
# )

## Covid negative all communities ----
# Summarize the CPM data
cov_pos_cpm_summary_list <- list()
taxa_levels <- c("Phylum", "Class", "Order", "Family", "Genus")

for (i in 1:length(taxa_levels)) {
  cov_pos_cpm_summary_list[[i]] <-
    ps2df(ps_cpm, column.is.taxa = FALSE) |>
    filter(Covid == "Positive") |>
    summarize(
      avg_cpm = mean(Abundance, na.rm = TRUE),
      se_cpm = sd(Abundance, na.rm = TRUE) / sqrt(n()),
      .by = c(taxa_levels[i], "Community2")
    ) |>
    dplyr::rename(taxon = taxa_levels[i]) |>
    mutate(
      feature = paste0("<i>", taxon, "</i> (", taxa_levels[i], ")"),
    )
}

# Combine the summaries into one data frame
cov_pos_cpm_summary <-
  do.call(
    rbind,
    cov_pos_cpm_summary_list
  )
cov_pos_cpm_summary

# Reformat the summary df
cov_pos_cpm_summary_2 <-
  cov_pos_cpm_summary |>
  pivot_wider(
    names_from = Community2,
    values_from = c(avg_cpm, se_cpm),
    names_sep = "_"
  ) |>
  select(-taxon) |>
  merge(
    cov_pos_lefse_res |>
      select(markerID, feature, ef_lda, padj),
    by = "feature"
  )

# Create a summary table using gt
cov_pos_cpm_summary_gt <-
  cov_pos_cpm_summary_2 |>
  arrange(markerID) |>
  gt() |>
  cols_hide(columns = "markerID") |>
  tab_spanner(
    label = "Elderly (CPM)",
    columns = c(avg_cpm_Elderly, se_cpm_Elderly)
  ) |>
  tab_spanner(
    label = "Marginalized (CPM)",
    columns = c(avg_cpm_Marginalized, se_cpm_Marginalized)
  ) |>
  tab_spanner(
    label = "Minority (CPM)",
    columns = c(avg_cpm_Minority, se_cpm_Minority)
  ) |>
  fmt_markdown(columns = feature) |>
  fmt_number(
    columns = c(avg_cpm_Elderly, avg_cpm_Marginalized, avg_cpm_Minority, se_cpm_Elderly, se_cpm_Marginalized, se_cpm_Minority),
    decimals = 2,
    scale_by = 100
  ) |>
  fmt_number(
    columns = c(ef_lda),
    decimals = 2
  ) |>
  fmt_scientific(
    columns = c(padj),
    decimals = 2
  ) |>
  cols_label(
    feature = "Taxon (level)",
    avg_cpm_Elderly = "Average",
    avg_cpm_Marginalized = "Average",
    avg_cpm_Minority = "Average",
    se_cpm_Elderly = "SE",
    se_cpm_Marginalized = "SE",
    se_cpm_Minority = "SE",
    ef_lda = "LDA score",
    padj = "Adjusted p-value"
  )
cov_pos_cpm_summary_gt

# Save the summary table as an word file
# gt::gtsave(
#   cov_pos_cpm_summary_gt,
#   filename = "processed_data/summary_tables/cov_pos_cpm_summary_table.docx"
# )

## Elderly community + vs - ----
# Summarize the CPM data
eld_cov_cpm_summary_list <- list()
taxa_levels <- c("Phylum", "Class", "Order", "Family", "Genus")

for (i in 1:length(taxa_levels)) {
  eld_cov_cpm_summary_list[[i]] <-
    ps2df(ps_cpm, column.is.taxa = FALSE) |>
    filter(Community2 == "Elderly") |>
    filter(Covid != "NA") |>
    summarize(
      avg_cpm = mean(Abundance, na.rm = TRUE),
      se_cpm = sd(Abundance, na.rm = TRUE) / sqrt(n()),
      .by = c(taxa_levels[i], "Covid")
    ) |>
    dplyr::rename(taxon = taxa_levels[i]) |>
    mutate(
      feature = paste0("<i>", taxon, "</i> (", taxa_levels[i], ")"),
    )
}

# Combine the summaries into one data frame
eld_cov_cpm_summary <-
  do.call(
    rbind,
    eld_cov_cpm_summary_list
  )
eld_cov_cpm_summary

# Reformat the summary df
eld_cov_cpm_summary_2 <-
  eld_cov_cpm_summary |>
  pivot_wider(
    names_from = Covid,
    values_from = c(avg_cpm, se_cpm),
    names_sep = "_"
  ) |>
  select(-taxon) |>
  merge(
    eld_cov_lefse_res |>
      select(markerID, feature, ef_lda, padj),
    by = "feature"
  )

# Create a summary table using gt
eld_cov_cpm_summary_gt <-
  eld_cov_cpm_summary_2 |>
  arrange(markerID) |>
  gt() |>
  cols_hide(columns = "markerID") |>
  tab_spanner(
    label = "SARS-CoV-2 Positive (CPM)",
    columns = c(avg_cpm_Positive, se_cpm_Positive)
  ) |>
  tab_spanner(
    label = "SARS-CoV-2 Negative (CPM)",
    columns = c(avg_cpm_Negative, se_cpm_Negative)
  ) |>
  fmt_markdown(columns = feature) |>
  fmt_number(
    columns = c(avg_cpm_Positive, avg_cpm_Negative, se_cpm_Positive, se_cpm_Negative),
    decimals = 2,
    scale_by = 100
  ) |>
  fmt_number(
    columns = c(ef_lda),
    decimals = 2
  ) |>
  fmt_scientific(
    columns = c(padj),
    decimals = 2
  ) |>
  cols_label(
    feature = "Taxon (level)",
    avg_cpm_Positive = "Average",
    avg_cpm_Negative = "Average",
    se_cpm_Positive = "SE",
    se_cpm_Negative = "SE",
    ef_lda = "LDA score",
    padj = "Adjusted p-value"
  )
eld_cov_cpm_summary_gt

# Save the summary table as an word file
# gt::gtsave(
#   eld_cov_cpm_summary_gt,
#   filename = "processed_data/summary_tables/eld_cov_cpm_summary_table.docx"
# )

## Elderly community + vs - ----
# Summarize the CPM data
mar_cov_cpm_summary_list <- list()
taxa_levels <- c("Phylum", "Class", "Order", "Family", "Genus")

for (i in 1:length(taxa_levels)) {
  mar_cov_cpm_summary_list[[i]] <-
    ps2df(ps_cpm, column.is.taxa = FALSE) |>
    filter(Community2 == "Marginalized") |>
    filter(Covid != "NA") |>
    summarize(
      avg_cpm = mean(Abundance, na.rm = TRUE),
      se_cpm = sd(Abundance, na.rm = TRUE) / sqrt(n()),
      .by = c(taxa_levels[i], "Covid")
    ) |>
    dplyr::rename(taxon = taxa_levels[i]) |>
    mutate(
      feature = paste0("<i>", taxon, "</i> (", taxa_levels[i], ")"),
    )
}

# Combine the summaries into one data frame
mar_cov_cpm_summary <-
  do.call(
    rbind,
    mar_cov_cpm_summary_list
  )
mar_cov_cpm_summary

# Reformat the summary df
mar_cov_cpm_summary_2 <-
  mar_cov_cpm_summary |>
  pivot_wider(
    names_from = Covid,
    values_from = c(avg_cpm, se_cpm),
    names_sep = "_"
  ) |>
  select(-taxon) |>
  merge(
    mar_lefse_res |>
      select(markerID, feature, ef_lda, padj),
    by = "feature"
  )

# Create a summary table using gt
mar_cov_cpm_summary_gt <-
  mar_cov_cpm_summary_2 |>
  arrange(markerID) |>
  gt() |>
  cols_hide(columns = "markerID") |>
  tab_spanner(
    label = "SARS-CoV-2 Positive (CPM)",
    columns = c(avg_cpm_Positive, se_cpm_Positive)
  ) |>
  tab_spanner(
    label = "SARS-CoV-2 Negative (CPM)",
    columns = c(avg_cpm_Negative, se_cpm_Negative)
  ) |>
  fmt_markdown(columns = feature) |>
  fmt_number(
    columns = c(avg_cpm_Positive, avg_cpm_Negative, se_cpm_Positive, se_cpm_Negative),
    decimals = 2,
    scale_by = 100
  ) |>
  fmt_number(
    columns = c(ef_lda),
    decimals = 2
  ) |>
  fmt_scientific(
    columns = c(padj),
    decimals = 2
  ) |>
  cols_label(
    feature = "Taxon (level)",
    avg_cpm_Positive = "Average",
    avg_cpm_Negative = "Average",
    se_cpm_Positive = "SE",
    se_cpm_Negative = "SE",
    ef_lda = "LDA score",
    padj = "Adjusted p-value"
  )
mar_cov_cpm_summary_gt

# Save the summary table as an word file
# gt::gtsave(
#   mar_cov_cpm_summary_gt,
#   filename = "processed_data/summary_tables/mar_cov_cpm_summary_table.docx"
# )

# Create data sheets as follows ----
## Relative abundance of top 20 taxa from all communities ----
overall_df_top20 <- ps_filt |>
  microbiome::transform(transform = "compositional") |>
  ps2df() |>
  # filter(Abundance > 0) |>
  mutate(Abundance = Abundance * 100) |>
  summarize(
    avg_ra = mean(Abundance, na.rm = TRUE),
    med_ra = median(Abundance, na.rm = TRUE),
    min_ra = min(Abundance, na.rm = TRUE),
    max_ra = max(Abundance, na.rm = TRUE),
    sd_ra = sd(Abundance, na.rm = TRUE),
    se_ra = sd_ra / sqrt(n()),
    perc95_ci = qt(0.975, df = n() - 1) * (se_ra),
    .by = c("Genus")
  ) |>
  # mutate_if(is.numeric, function(x){x*100}) |>
  arrange(desc(avg_ra)) |>
  head(20)
# overall_df_top20 |> View()

# Save excel file ---
# write.csv(
#   overall_df_top20 |>
#   dplyr::rename(
#     Mean = "avg_ra",
#     Median = "med_ra",
#     Min = "min_ra",
#     Max = "max_ra",
#     SD = "sd_ra",
#     SE = "se_ra",
#     `95% CI` = "perc95_ci"
#   ),
#   file = "processed_data/summary_tables/overall_top20.csv",
#   row.names = FALSE
# )

### Boxplot ----
overall_top20_gg <-
  ps_filt |>
  microbiome::transform(transform = "compositional") |>
  ps2df() |>
  filter(Genus %in% overall_df_top20$Genus) |>
  # filter(Abundance > 0) |>
  mutate(
    Genus = factor(Genus, levels = overall_df_top20$Genus)
  ) |>
  mutate(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    .by = "Genus"
  ) |>
  ggplot(
    mapping = aes(
      x = Genus,
      y = Abundance * 100
    )
  ) +
  geom_point(alpha = 0.2, position = position_jitter(width = 0.25)) +
  geom_boxplot(alpha = 0.1, outliers = FALSE) +
  # stat_summary(fun = "mean", geom = "point", size = 3, color = "salmon", shape = 23) +
  geom_point(mapping = aes(y = mean_abundance * 100), size = 3, color = "salmon", shape = 23) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.25, color = "salmon", position = position_nudge(x = -0.2)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, color = "lightgreen") +
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.25, color = "dodgerblue", position = position_nudge(x = 0.2)) +
  labs(
    x = "Genus",
    y = "Relative abundance",
    title = "Top 20 taxa from all communities",
    caption = "Boxplot = five-number summary (min, Q1, median, Q3, max);<br>salmon diamond = mean; salmon error bar = mean ± 1 SD;<br>green error bar = mean ± SE; blue error bar = mean ± 95% CI"
  ) +
  scale_y_continuous(
    # transform = "log10",
    # transform = scales::pseudo_log_trans(base = 4),
    transform = "sqrt",
    # breaks = seq(0, 100, by = 4),
    breaks = c(0, 4, 16, 36, 64, 100),
    # limits = c(0.001, 100),
    limits = c(-10, 100),
    labels = scales::number_format(accuracy = 1, scale = 1, suffix = "%")
  ) +
  # annotation_logticks(sides = "l") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = ggtext::element_markdown()
  )
overall_top20_gg

# ggsave(
#   plot = overall_top20_gg,
#   filename = "processed_data/figures/version4/overall_top20_boxplot.png",
#   units = "px", width = 1500, height = 1500, dpi = 320
# )

## Relative abundance of top 20 taxa from nursing home samples only ----
overall_nh_top20 <- ps_filt |>
  microbiome::transform(transform = "compositional") |>
  ps2df() |>
  filter(Community2 == "Elderly") |>
  # filter(Abundance > 0) |>
  mutate(Abundance = Abundance * 100) |>
  summarize(
    avg_ra = mean(Abundance, na.rm = TRUE),
    med_ra = median(Abundance, na.rm = TRUE),
    min_ra = min(Abundance, na.rm = TRUE),
    max_ra = max(Abundance, na.rm = TRUE),
    sd_ra = sd(Abundance, na.rm = TRUE),
    se_ra = sd(Abundance, na.rm = TRUE) / sqrt(n()),
    perc95_ci = qt(0.975, df = n() - 1) * (se_ra),
    .by = c("Genus")
  ) |>
  arrange(desc(avg_ra)) |>
  # mutate_if(is.numeric, function(x){x*100}) |>
  head(20)
# overall_nh_top20 |> View()

# Save excel file ---
# write.csv(
#   overall_nh_top20 |>
#   dplyr::rename(
#     Mean = "avg_ra",
#     Median = "med_ra",
#     Min = "min_ra",
#     Max = "max_ra",
#     SD = "sd_ra",
#     SE = "se_ra",
#     `95% CI` = "perc95_ci"
#   ),
#   file = "processed_data/summary_tables/overall_nh_top20.csv",
#   row.names = FALSE
# )

### Boxplot ----
nh_top20_gg <- ps_filt |>
  microbiome::transform(transform = "compositional") |>
  ps2df() |>
  filter(Genus %in% overall_nh_top20$Genus) |>
  filter(Community2 == "Elderly") |>
  # filter(Abundance > 0) |>
  mutate(
    Genus = factor(Genus, levels = overall_nh_top20$Genus)
  ) |>
  mutate(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    .by = "Genus"
  ) |>
  ggplot(
    mapping = aes(
      x = Genus,
      y = Abundance * 100
    )
  ) +
  geom_point(alpha = 0.2, position = position_jitter(width = 0.25)) +
  geom_boxplot(alpha = 0.1, outliers = FALSE) +
  # stat_summary(fun = "mean", geom = "point", size = 3, color = "salmon", shape = 23, fun.args = list(na.rm = TRUE)) +
  geom_point(mapping = aes(y = mean_abundance * 100), size = 3, color = "salmon", shape = 23) +
  # geom_text(mapping = aes(y = mean_abundance * 100, label = round( mean_abundance * 100, 1)), size = 3, color = "black") +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.25, color = "salmon", position = position_nudge(x = -0.2)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, color = "lightgreen") +
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.25, color = "dodgerblue", position = position_nudge(x = 0.2)) +
  labs(
    x = "Genus",
    y = "Relative abundance",
    title = "Top 20 taxa from nursing homes",
    caption = "Boxplot = five-number summary (min, Q1, median, Q3, max);<br>salmon diamond = mean; salmon error bar = mean ± 1 SD;<br>green error bar = mean ± SE; blue error bar = mean ± 95% CI"
  ) +
  scale_y_continuous(
    transform = "sqrt",
    breaks = c(0, 4, 16, 36, 64, 100),
    limits = c(-10, 100),
    labels = scales::number_format(accuracy = 1, scale = 1, suffix = "%")
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = ggtext::element_markdown()
  )
nh_top20_gg

# ggsave(
#   plot = nh_top20_gg,
#   filename = "processed_data/figures/version4/nh_top20_boxplot.png",
#   units = "px", width = 1500, height = 1500, dpi = 320
# )

## Relative abundance of top 20 taxa from marginalized samples only ----
overall_ma_top20 <- ps_filt |>
  microbiome::transform(transform = "compositional") |>
  ps2df() |>
  filter(Community2 == "Marginalized") |>
  mutate(Abundance = Abundance * 100) |>
  # filter(Abundance > 0) |>
  summarize(
    avg_ra = mean(Abundance, na.rm = TRUE),
    med_ra = median(Abundance, na.rm = TRUE),
    min_ra = min(Abundance, na.rm = TRUE),
    max_ra = max(Abundance, na.rm = TRUE),
    sd_ra = sd(Abundance, na.rm = TRUE),
    se_ra = sd(Abundance, na.rm = TRUE) / sqrt(n()),
    perc95_ci = qt(0.975, df = n() - 1) * (se_ra),
    .by = c("Genus")
  ) |>
  arrange(desc(avg_ra)) |>
  # mutate_if(is.numeric, function(x){x*100}) |>
  head(20)
# overall_ma_top20 |> View()

# Save excel file ---
# write.csv(
#   overall_ma_top20 |>
#   dplyr::rename(
#     Mean = "avg_ra",
#     Median = "med_ra",
#     Min = "min_ra",
#     Max = "max_ra",
#     SD = "sd_ra",
#     SE = "se_ra",
#     `95% CI` = "perc95_ci"
#   ),
#   file = "processed_data/summary_tables/overall_ma_top20.csv",
#   row.names = FALSE
# )

### Boxplot ----
ma_top20_gg <-
  ps_filt |>
  microbiome::transform(transform = "compositional") |>
  ps2df() |>
  filter(Genus %in% overall_ma_top20$Genus) |>
  filter(Community2 == "Marginalized") |>
  # filter(Abundance > 0) |>
  mutate(
    Genus = factor(Genus, levels = overall_ma_top20$Genus)
  ) |>
  mutate(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    .by = "Genus"
  ) |>
  ggplot(
    mapping = aes(
      x = Genus,
      y = Abundance * 100
    )
  ) +
  geom_point(alpha = 0.2, position = position_jitter(width = 0.25)) +
  geom_boxplot(alpha = 0.1, outliers = FALSE) +
  # stat_summary(fun = "mean", geom = "point", size = 3, color = "salmon", shape = 23) +
  geom_point(mapping = aes(y = mean_abundance * 100), size = 3, color = "salmon", shape = 23) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.25, color = "salmon", position = position_nudge(x = -0.2)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, color = "lightgreen") +
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.25, color = "dodgerblue", position = position_nudge(x = 0.2)) +
  labs(
    x = "Genus",
    y = "Relative abundance",
    title = "Top 20 taxa from marginalized",
    caption = "Boxplot = five-number summary (min, Q1, median, Q3, max);<br>salmon diamond = mean; salmon error bar = mean ± 1 SD;<br>green error bar = mean ± SE; blue error bar = mean ± 95% CI"
  ) +
  scale_y_continuous(
    transform = "sqrt",
    breaks = c(0, 4, 16, 36, 64, 100),
    limits = c(-10, 100),
    labels = scales::number_format(accuracy = 1, scale = 1, suffix = "%")
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = ggtext::element_markdown()
  )
ma_top20_gg

# ggsave(
#   plot = ma_top20_gg,
#   filename = "processed_data/figures/version4/ma_top20_boxplot.png",
#   units = "px", width = 1500, height = 1500, dpi = 320
# )

## Relative abundance of top 20 taxa from minority samples only ----
overall_mi_top20 <- ps_filt |>
  microbiome::transform(transform = "compositional") |>
  ps2df() |>
  mutate(Abundance = Abundance * 100) |>
  filter(Community2 == "Minority") |>
  summarize(
    avg_ra = mean(Abundance, na.rm = TRUE),
    med_ra = median(Abundance, na.rm = TRUE),
    min_ra = min(Abundance, na.rm = TRUE),
    max_ra = max(Abundance, na.rm = TRUE),
    sd_ra = sd(Abundance, na.rm = TRUE),
    se_ra = sd(Abundance, na.rm = TRUE) / sqrt(n()),
    perc95_ci = qt(0.975, df = n() - 1) * (se_ra),
    .by = c("Genus")
  ) |>
  arrange(desc(avg_ra)) |>
  # mutate_if(is.numeric, function(x){x*100}) |>
  head(20)
# overall_mi_top20 |> View()

# Save excel file ---
# write.csv(
#   overall_mi_top20 |>
#   dplyr::rename(
#     Mean = "avg_ra",
#     Median = "med_ra",
#     Min = "min_ra",
#     Max = "max_ra",
#     SD = "sd_ra",
#     SE = "se_ra",
#     `95% CI` = "perc95_ci"
#   ),
#   file = "processed_data/summary_tables/overall_mi_top20.csv",
#   row.names = FALSE
# )

### Boxplot ----
mi_top20_gg <- ps_filt |>
  microbiome::transform(transform = "compositional") |>
  ps2df() |>
  filter(Genus %in% overall_mi_top20$Genus) |>
  filter(Community2 == "Minority") |>
  # filter(Abundance > 0) |>
  mutate(
    Genus = factor(Genus, levels = overall_mi_top20$Genus)
  ) |>
  mutate(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    .by = "Genus"
  ) |>
  ggplot(
    mapping = aes(
      x = Genus,
      y = Abundance * 100
    )
  ) +
  geom_point(alpha = 0.2, position = position_jitter(width = 0.25)) +
  geom_boxplot(alpha = 0.1, outliers = FALSE) +
  # stat_summary(fun = "mean", geom = "point", size = 3, color = "salmon", shape = 23) +
  geom_point(mapping = aes(y = mean_abundance * 100), size = 3, color = "salmon", shape = 23) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.25, color = "salmon", position = position_nudge(x = -0.2)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, color = "lightgreen") +
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.25, color = "dodgerblue", position = position_nudge(x = 0.2)) +
  labs(
    x = "Genus",
    y = "Relative abundance",
    title = "Top 20 taxa from minority",
    caption = "Boxplot = five-number summary (min, Q1, median, Q3, max);<br>salmon diamond = mean; salmon error bar = mean ± 1 SD;<br>green error bar = mean ± SE; blue error bar = mean ± 95% CI"
  ) +
  scale_y_continuous(
    transform = "sqrt",
    breaks = c(0, 4, 16, 36, 64, 100),
    limits = c(-10, 100),
    labels = scales::number_format(accuracy = 1, scale = 1, suffix = "%")
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = ggtext::element_markdown()
  )
mi_top20_gg

# ggsave(
#   plot = mi_top20_gg,
#   filename = "processed_data/figures/version4/mi_top20_boxplot.png",
#   units = "px", width = 1500, height = 1500, dpi = 320
# )

# Run LEfSe with percentage data instead of CPM ----
## Community comparisons across the covid status ----
com_res_perc <- run_lefse(
  ps_filt,
  group = "Community2",
  subgroup = "Covid",
  multigrp_strat = TRUE, # enable multi-class mode
  strict = "1",
  kw_cutoff = 0.05,
  wilcoxon_cutoff = 0.05,
  lda_cutoff = 4,
  norm = "TSS" # TSS = total sum scaling, aka percentage
)
com_res_perc
# No microbiome markers were found

## Covid-positive Community comparisons ----
cov_pos_res_perc <- run_lefse(
  ps_filt |> subset_samples(Covid == "Positive"),
  group = "Community2",
  # subgroup = "Covid",
  multigrp_strat = TRUE, # enable multi-class mode
  strict = "1",
  kw_cutoff = 0.05,
  wilcoxon_cutoff = 0.05,
  lda_cutoff = 4,
  norm = "TSS" # TSS = total sum scaling, aka percentage
)
cov_pos_res_perc
# No microbiome markers were found

## Covid-positive Community comparisons ----
cov_neg_res_perc <- run_lefse(
  ps_filt |> subset_samples(Covid == "Negative"),
  group = "Community2",
  # subgroup = "Covid",
  multigrp_strat = TRUE, # enable multi-class mode
  strict = "1",
  kw_cutoff = 0.05,
  wilcoxon_cutoff = 0.05,
  lda_cutoff = 4,
  norm = "TSS" # TSS = total sum scaling, aka percentage
)
cov_neg_res_perc
# No microbiome markers were found

# Run ANCOMBC2 with non-filtered data ----
library("ANCOMBC")
## Prepare data ----
tse <- mia::convertFromPhyloseq(ps |> subset_samples(Covid != "NA"))

## Run ANCOMBC2 ----
tse_ancom_res <- ancombc2(
  data = tse,
  assay_name = "counts",
  tax_level = "Genus",
  fix_formula = "Community2 + Covid",
  # rand_formula = "(1 | Site)",
  p_adj_method = "BH",
  pseudo_sens = TRUE,
  prv_cut = 0.001,
  lib_cut = 0,
  group = "Community2",
  struc_zero = TRUE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 10,
  verbose = TRUE,
  global = TRUE,
  pairwise = TRUE,
  dunnet = FALSE,
  trend = FALSE
)
tse_ancom_res$res |> View()
tse_ancom_res$res_pair |> View()
tse_ancom_res$res_global |> View()
