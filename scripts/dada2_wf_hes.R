# Script header ----
# Programmed by: Kairi Tanaka
# Programmed on: 2023-04-04
# Programmed to: Create a template to run dada2 pipeline
# Last modified by:
# Last modified on:
# Last modified to:

# Load packages ----
packages <- c(
  "tidyverse",
  "dada2",
  "xlsx"
)

lapply(
  packages,
  library,
  character.only = TRUE
)

# Switch working directory to the project folder ----
setwd("./sewage_health_equity/") # CHANGE this to the project directory

# Run DADA2 pipeline ----
## Specify path ----
path <- "./intermediate_data/16s_clean_fastq/" # CHANGE this to the directory containing the fastq files
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
# for (i in 1:length(qual_vec)) {
#   # Forward read quality pre QC
#   jpeg(
#     file = paste0("intermediate_data/dada2/F_qual_",i,".jpg"),
#     units = "in", width = 8, height = 8, res = 300
#   )
#   plotQualityProfile(fnFs[qual_vec[[i]]])
#   dev.off()
#
#   # Reverse read quality pre QC
#   jpeg(
#     file = paste0("intermediate_data/dada2/R_qual_",i,".jpg"),
#     units = "in", width = 8, height = 8, res = 300
#   )
#   plotQualityProfile(fnRs[qual_vec[[i]]])
#   dev.off()
# }

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
saveRDS(
  out,
  file = "intermediate_data/dada2/out.rds"
)

## Plot quality profiles of filtered reads ----
# Plot the quality scores of the sequencing read files via for loop
# for (i in 1:length(qual_vec)) {
#   # Forward read quality pre QC
#   jpeg(
#     file = paste0("intermediate_data/dada2/filtF_qual_",i,".jpg"),
#     units = "in", width = 8, height = 8, res = 300
#   )
#   plotQualityProfile(filtFs[qual_vec[[i]]])
#   dev.off()
#
#   # Reverse read quality pre QC
#   jpeg(
#     file = paste0("intermediate_data/dada2/filtR_qual_",i,".jpg"),
#     units = "in", width = 8, height = 8, res = 300
#   )
#   plotQualityProfile(filtRs[qual_vec[[i]]])
#   dev.off()
# }

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
jpeg(
  file = "intermediate_data/dada2/F_error_plot.jpg",
  units = "in", width = 9, height = 9, res = 300
)
plotErrors(
  errF,
  nominalQ = TRUE
)
dev.off()

jpeg(
  file = "intermediate_data/dada2/R_error_plot.jpg",
  units = "in", width = 9, height = 9, res = 300
)
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

# Inspect the merger dataframe from the first sample
# head(mergers[[1]])
# View(mergers[[1]])
# length(mergers)

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
sum(seqtab.nochim / sum(seqtab))
# ??% was chimeric sequences

# Save the data
saveRDS(
  seqtab.nochim,
  file = "intermediate_data/dada2/seqtab_nochim.rds"
)

# Save sequence table as xlsx file
write.xlsx(
  seqtab.nochim,
  file = "processed_data/summary_tables/dada2_sequence_table.xlsx"
)

## Track reads ----
getN <- function(x) sum(getUniques(x))

track <-
  cbind(
    out,
    sapply(dadaFs, getN),
    sapply(dadaRs, getN),
    sapply(mergers, getN),
    rowSums(seqtab.nochim)
  )
# If processing a single sample, remove the sapply calls:
# e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c(
  "input",
  "filtered",
  "denoisedF",
  "denoisedR",
  "merged",
  "nonchim"
)

rownames(track) <- sample.names

head(track)
# View(track)

# Save data
saveRDS(
  track,
  file = "intermediate_data/dada2/track.rds"
)

# Save n of reads tracking table as an xlsx file
write.xlsx(
  track,
  file = "processed_data/summary_tables/dada2_read_tracking.xlsx"
)

## Assign taxonomy ----
# Silva 138.2 database
silva_taxa <-
  assignTaxonomy(
    seqtab.nochim,
    "../databases/dada2/silva_nr99_v138.2_toGenus_trainset.fa.gz",
    tryRC = TRUE,
    multithread = TRUE,
    verbose = TRUE
  )

# Save as rds
saveRDS(
  silva_taxa,
  file = "intermediate_data/dada2/silva_taxa.rds"
)

# Save taxa table as xlsx
write.xlsx(
  silva_taxa,
  file = "processed_data/summary_tables/dada2_silva_taxa.xlsx"
)

# Alternatively, GreenGenes2 (comment out above and comment in below)
gg2_taxa <-
  assignTaxonomy(
    seqtab.nochim,
    "../databases/dada2/gg2_2024_09_toGenus_trainset.fa.gz",
    tryRC = TRUE,
    multithread = TRUE,
    verbose = TRUE
  )

# Save as rds
saveRDS(
  gg2_taxa,
  file = "intermediate_data/dada2/gg2_taxa.rds"
)

# Save taxa table as xlsx
write.xlsx(
  gg2_taxa,
  file = "processed_data/summary_tables/dada2_gg2_taxa.xlsx"
)

# Inspect the taxonomic assignments
# Copy the dataset for display purpose only
taxa.print <- silva_taxa
# taxa.print <- gg2_taxa

# Removing sequence rownames for DISPLAY ONLY
rownames(taxa.print) <- NULL

# Show the first few
# head(taxa.print)
# Or the whole taxa table
View(taxa.print)

# Merge taxa and sequence tables ----
taxseq_df <- silva_taxa |>
  data.frame() |>
  rownames_to_column(var = "ASV") |>
  merge(
    seqtab.nochim |>
      data.frame() |>
      rownames_to_column(var = "SampleID") |>
      pivot_longer(
        cols = -SampleID,
        names_to = "ASV",
        values_to = "n_reads"
      ),
    by = "ASV"
  )
taxseq_df |> View()

# Get taxonomic levels in dataset
tax_levels <- colnames(silva_taxa)

# Initiate a list to hold data
tax_summ_list <- list()

# Summarize data via for loop
for (i in 1:length(tax_levels)) {
  tax_summ_list[[i]] <- taxseq_df |>
    mutate(
      longKingdom = Kingdom,
      longPhylum = paste0(Kingdom, "; ", Phylum),
      longClass = paste0(Kingdom, "; ", Phylum, "; ", Class),
      longOrder = paste0(Kingdom, "; ", Phylum, "; ", Class, "; ", Order),
      longFamily = paste0(Kingdom, "; ", Phylum, "; ", Class, "; ", Order, "; ", Family),
      longGenus = paste0(Kingdom, "; ", Phylum, "; ", Class, "; ", Order, "; ", Family, "; ", Genus)
    ) |>
    select(!Kingdom:Genus) |>
    rename_at(vars(starts_with("long")), ~ str_remove(., "long")) |>
    summarize(
      sum_n_reads = sum(n_reads, na.rm = TRUE),
      .by = c((tax_levels[i]), SampleID)
    ) |>
    pivot_wider(
      names_from = "SampleID",
      values_from = "sum_n_reads"
    )
}

# Save as xlsx file
for (i in 1:length(tax_levels)) {
  if (i == 1) {
    write.xlsx(
      tax_summ_list[[i]],
      file = "processed_data/summary_tables/dada2_tax_summed_n_reads.xlsx",
      sheetName = tax_levels[i],
      append = FALSE
    )
  } else {
    write.xlsx(
      tax_summ_list[[i]],
      file = "processed_data/summary_tables/dada2_tax_summed_n_reads.xlsx",
      sheetName = tax_levels[i],
      append = TRUE
    )
  }
}
