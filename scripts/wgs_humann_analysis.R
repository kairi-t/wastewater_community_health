# Script header ----
# Programmed by: Kairi Tanaka
# Programmed on: 05-10-2024
# Porgrammed to: run MaAsLin2 on the genefamilies data
# Last modified by: Kairi Tanaka
# Last modified on: 06-17-2025
# Last modified to: visualize the assimilatory sulfate reduction pathways

# Install packages ----
# install.packages('tidyverse')
# if(!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install('Maaslin2')

# Load packages ----
packages <- c(
    "tidyverse",
    "Maaslin2",
    "cowplot",
    "ggpubr",
    "ggtext"
)

lapply(
    packages,
    library,
    character.only = TRUE
)

# Set working directory ----
setwd("~/Documents/LyuLab/sewage_health_equity")
getwd()

# Load data ----
## Uniref90 gene families ----
uniref90_cpm <- read_delim(
    file = "./intermediate_data/biobakery/humann/humann_combined_named_genefamilies_cpm.tsv",
    delim = "\t"
)
# uniref_cpm |> View()

## metacyc reaction ----
rxn_cpm <- read_delim(
    file = "./intermediate_data/biobakery/humann/humann_combined_rxn_cpm_named.tsv",
    delim = "\t"
)
# rxn_cpm |> View()

## Metacyc pathway ----
pathway_cpm <- read_delim(
    file = "./intermediate_data/biobakery/humann/humann_combined_pathabundance_cpm.tsv",
    delim = "\t"
)
# pathway_cpm |> View()

## EC ----
ec_cpm <- read_delim(
    file = "./intermediate_data/biobakery/humann/humann_combined_named_level4ec_cpm.tsv",
    delim = "\t"
)
# ec_cpm |> View()

## KO ----
ko_cpm <- read_delim(
    file = "./intermediate_data/biobakery/humann/humann_combined_named_ko_cpm.tsv",
    delim = "\t"
)
# ko_cpm |> View()

## eggNog  ----
eggnog_cpm <- read_delim(
    file = "./intermediate_data/biobakery/humann/humann_combined_eggnog_cpm.tsv",
    delim = "\t"
)
# eggnog_cpm |> View()

## pfam  ----
pfam_cpm <- read_delim(
    file = "./intermediate_data/biobakery/humann/humann_combined_named_pfam_cpm.tsv",
    delim = "\t"
)
# pfam_cpm |> View()

## go ----
go_cpm <- read_delim(
    file = "./intermediate_data/biobakery/humann/humann_combined_named_go_cpm.tsv",
    delim = "\t"
)
# go_cpm |> View()

# Clean data ----
## Create a function to clean humann dataframe ----
clean_humann_df <- function(df) {
    clean_df <-
        df |>
        rename(
            feature = contains("#")
        ) |>
        rename_with(
            .cols = where(is.numeric),
            .fn = function(x) {
                # Extract the sample names
                pre_name <- str_split_i(
                    string = x,
                    pattern = "_",
                    i = 1
                )

                # Refine samplennames
                chamber <- str_split_i(
                    string = pre_name,
                    pattern = "-",
                    i = 1
                )

                sample_id <- str_split_i(
                    string = pre_name,
                    pattern = "-",
                    i = 4
                )

                lane <- str_split_i(
                    string = pre_name,
                    pattern = "-",
                    i = 6
                )

                # Combine
                paste0(chamber, "_", sample_id, "_", lane)
            }
        ) |>
        # select(!(contains('D_'))) |>
        column_to_rownames(var = "feature") |>
        t() |>
        data.frame()
    return(clean_df)
}

## UniRef90 ----
cleaned_uniref90_cpm <- clean_humann_df(uniref90_cpm)
# cleaned_uniref90_cpm |> View()

## Metacyc reaction ----
cleaned_rxn_cpm <- clean_humann_df(rxn_cpm)
# cleaned_rxn_cpm |> View()

## Metacyc pathway ----
cleaned_path_cpm <- clean_humann_df(pathway_cpm)
# cleaned_ec_cpm |> View()

## EC ----
cleaned_ec_cpm <- clean_humann_df(ec_cpm)
# cleaned_ec_cpm |> View()

## KO ----
cleaned_ko_cpm <- clean_humann_df(ko_cpm)
# cleaned_ko_cpm |> View()

## eggNog ----
cleaned_eggnog_cpm <- clean_humann_df(eggnog_cpm)
# cleaned_eggnog_cpm |> View()

## pfam ----
cleaned_pfam_cpm <- clean_humann_df(pfam_cpm)
# cleaned_pfam_cpm |> View()

## GO ----
cleaned_go_cpm <- clean_humann_df(go_cpm)
# cleaned_go_cpm |> View()

## Prepare the input data for maaslin2 ----
### Create a function to select only the feature total and get rid of the stratified ----
select_feature_total <- function(df) {
    feature_total_df <-
        df |>
        select(!(contains(".g_") | contains(".unclassified")))
    return(feature_total_df)
}

### Uniref90 gene families ----
mal2_uniref_input <-
    cleaned_uniref_cpm |>
    select_feature_total()
# mal2_uniref_input

### Metacyc pathways ----
mal2_path_input <-
    cleaned_pathway_cpm |>
    select_feature_total()
# mal2_path_input

### Metacyc reactions ----
mal2_rxn_input <-
    cleaned_rxn_cpm |>
    select_feature_total()
# mal2_rxn_input

### EC ----
mal2_ec_input <-
    cleaned_ec_cpm |>
    select_feature_total()
# mal2_ec_input

### KO ----
mal2_ko_input <-
    cleaned_ko_cpm |>
    select_feature_total()
# mal2_ko_input

### pfam ----
mal2_pfam_input <-
    cleaned_pfam_cpm |>
    select_feature_total()
# mal2_pfam_input

### GO ----
mal2_go_input <-
    cleaned_go_cpm |>
    select_feature_total()
# mal2_go_input

# Create metadata ----
## Overall comparison ----
mal2_metadata <-
    data.frame(
        SampleID = rownames(cleaned_ec_cpm)
    ) |>
    separate(
        SampleID,
        into = c("Chamber", "Sample", "Lane"),
        sep = "_"
    ) |>
    mutate(
        SampleID = paste0(Chamber, "_", Sample, "_", Lane)
    ) |>
    column_to_rownames(var = "SampleID")
# filter(Chamber != "D")
mal2_metadata

## AB comparison ----
ab_mal2_metadata <-
    mal2_metadata |>
    filter(Chamber != "C")

## BC comparison ----
bc_mal2_metadata <-
    mal2_metadata |>
    filter(Chamber != "A")

# Run maaslin2 ----
# Create a list of metadata
meta_list <- list(
    "AllChamber" = mal2_metadata,
    "AB" = ab_mal2_metadata,
    "BC" = bc_mal2_metadata
)

# Make vector to loop through for the pairwise comparisons
cham_vec <- c("", "C_", "A_")

## Create a function to loop through the analysis ----
run_maaslin2 <- function(df, path) {
    # Make vector to loop through for the pairwise comparisons
    cham_vec <- c("", "C_", "A_")

    # Initiate a list to store results
    res_list <- list()

    # Loop through the maaslin2 runs for overall and pairwise comparisons
    for (i in 1:3) {
        if (i == 1) {
            print("i=", i)

            res_list[[i]] <- Maaslin2(
                input_data = df,
                input_metadata = meta_list[[i]],
                output = paste0(path, names(meta_list[i])),
                fixed_effects = c("Chamber"),
                random_effects = c("Sample", "Lane"),
                normalization = "NONE",
                transform = "NONE",
                cores = 10,
                reference = "Chamber,B",
                plot_scatter = FALSE
            )
        } else {
            print("i=", i)

            res_list[[i]] <- Maaslin2(
                input_data = df |> select(!(contains(cham_vec[i]))),
                input_metadata = meta_list[[i]],
                output = paste0(path, names(meta_list[i])),
                fixed_effects = c("Chamber"),
                random_effects = c("Sample", "Lane"),
                normalization = "NONE",
                transform = "NONE",
                cores = 10,
                reference = "Chamber,B",
                plot_scatter = FALSE
            )
        }
    }
    return(res_list)
}

## Uniref90 gene families ----
mal2_uniref_res_list <- run_maaslin2(
    mal2_uniref_input,
    path = "./intermediate_data/biobakery/humann/uniref90_mal2_"
)

## Pathways ----
mal2_path_res_list <- run_maaslin2(
    mal2_path_input,
    "./intermediate_data/biobakery/humann/path_mal2"
)

## Metacyc reaction ----
mal2_rxn_res_list <- run_maaslin2(
    mal2_rxn_input,
    "./intermediate_data/biobakery/humann/rxn_mal2_"
)

## EC ----
mal2_ec_res_list <- run_maaslin2(
    mal2_ec_input,
    "./intermediate_data/biobakery/humann/ec_mal2_"
)

## KO ----
mal2_ko_res_list <- run_maaslin2(
    mal2_ko_input,
    "./intermediate_data/biobakery/humann/ko_mal2_"
)

## pfam ----
mal2_pfam_res_list <- run_maaslin2(
    mal2_pfam_input,
    "./intermediate_data/biobakery/humann/pfam_mal2_"
)

## GO ----
mal2_go_res_list <- run_maaslin2(
    mal2_go_input,
    "./intermediate_data/biobakery/humann/go_mal2_"
)

# Prepare the output for the humann barplot ----
# The first row is sample ID
# the second set of N rows is metadata (N number of columns)
# The rest is the pathway data (where each row is the pathway and each column is the sample)
## Uniref90 gene families ----
uniref_humann_input <-
    mal2_metadata |>
    t() |>
    data.frame() |>
    rbind(
        named_uniref_cpm |>
            rename(
                gene_family = `# Gene Family`
            ) |>
            select(gene_family, contains("-RPKs")) |>
            rename_with(
                .cols = where(is.numeric),
                .fn = function(x) {
                    # Extract the sample names
                    pre_name <- str_split_i(
                        string = x,
                        pattern = "_",
                        i = 1
                    )

                    # Refine samplennames
                    chamber <- str_split_i(
                        string = pre_name,
                        pattern = "-",
                        i = 1
                    )

                    sample_id <- str_split_i(
                        string = pre_name,
                        pattern = "-",
                        i = 4
                    )

                    lane <- str_split_i(
                        string = pre_name,
                        pattern = "-",
                        i = 6
                    )

                    # Combine
                    paste0(chamber, "_", sample_id, "_", lane)
                }
            ) |>
            select(!(contains("D_"))) |>
            column_to_rownames(var = "gene_family")
    )
# uniref_humann_input |> View()

## Save table ----
write_delim(
    uniref_humann_input |> rownames_to_column(var = "Feature \\ SampleID"),
    file = "./intermediate_data/biobakery/humann/uniref90_mal2/uniref90_humann_input.tsv",
    delim = "\t"
)

## Metacyc reaction ----
rxn_humann_input <-
    mal2_metadata |>
    t() |>
    data.frame() |>
    rbind(
        named_rxn_cpm |>
            rename(
                gene_family = `# Gene Family`
            ) |>
            select(gene_family, contains("-RPKs")) |>
            rename_with(
                .cols = where(is.numeric),
                .fn = function(x) {
                    # Extract the sample names
                    pre_name <- str_split_i(
                        string = x,
                        pattern = "_",
                        i = 1
                    )

                    # Refine samplennames
                    chamber <- str_split_i(
                        string = pre_name,
                        pattern = "-",
                        i = 1
                    )

                    sample_id <- str_split_i(
                        string = pre_name,
                        pattern = "-",
                        i = 4
                    )

                    lane <- str_split_i(
                        string = pre_name,
                        pattern = "-",
                        i = 6
                    )

                    # Combine
                    paste0(chamber, "_", sample_id, "_", lane)
                }
            ) |>
            select(!(contains("D_"))) |>
            column_to_rownames(var = "gene_family")
    )
# rxn_humann_input |> View()

## Save table ----
write_delim(
    rxn_humann_input |> rownames_to_column(var = "Feature \\ SampleID"),
    file = "./intermediate_data/biobakery/humann/rxn_mal2/rxn_humann_input.tsv",
    delim = "\t"
)

## EC ----
ec_humann_input <-
    mal2_metadata |>
    t() |>
    data.frame() |>
    rbind(
        ec_cpm |>
            rename(
                EC = `# Gene Family`
            ) |>
            select(EC, contains("-RPKs")) |>
            rename_with(
                .cols = where(is.numeric),
                .fn = function(x) {
                    # Extract the sample names
                    pre_name <- str_split_i(
                        string = x,
                        pattern = "_",
                        i = 1
                    )

                    # Refine samplennames
                    chamber <- str_split_i(
                        string = pre_name,
                        pattern = "-",
                        i = 1
                    )

                    sample_id <- str_split_i(
                        string = pre_name,
                        pattern = "-",
                        i = 4
                    )

                    lane <- str_split_i(
                        string = pre_name,
                        pattern = "-",
                        i = 6
                    )

                    # Combine
                    paste0(chamber, "_", sample_id, "_", lane)
                }
            ) |>
            # select(!(contains('D_'))) |>
            column_to_rownames(var = "EC")
    )
# ec_humann_input |> View()

## Save table ----
write_delim(
    ec_humann_input |> rownames_to_column(var = "Feature \\ SampleID"),
    file = "./intermediate_data/biobakery/humann/ec_mal2/ec_humann_input.tsv",
    delim = "\t"
)

## Pathways ----
path_humann_input <-
    mal2_metadata |>
    t() |>
    data.frame() |>
    rbind(
        pathway_cpm |>
            rename(
                Pathway = `# Pathway`
            ) |>
            rename_with(
                .cols = where(is.numeric),
                .fn = function(x) {
                    # Extract the sample names
                    pre_name <- str_split_i(
                        string = x,
                        pattern = "_",
                        i = 1
                    )

                    # Refine samplennames
                    chamber <- str_split_i(
                        string = pre_name,
                        pattern = "-",
                        i = 1
                    )

                    sample_id <- str_split_i(
                        string = pre_name,
                        pattern = "-",
                        i = 4
                    )

                    lane <- str_split_i(
                        string = pre_name,
                        pattern = "-",
                        i = 6
                    )

                    # Combine
                    paste0(chamber, "_", sample_id, "_", lane)
                }
            ) |>
            select(!(contains("D_"))) |>
            column_to_rownames(var = "Pathway")
    )
# path_humann_input |> View()

## Save table ----
write_delim(
    path_humann_input |> rownames_to_column(var = "Feature \\ SampleID"),
    file = "./intermediate_data/biobakery/humann/pathway_mal2/path_humann_input.tsv",
    delim = "\t"
)

# Create figures for NCMS ----
## Establish theme for the barplot ----
font_size <- 25

bar_theme <-
    theme_set(
        theme_cowplot() +
            theme(
                panel.background = element_blank(),
                panel.grid.major.y = element_line(
                    linetype = "dashed",
                    color = "grey77"
                ),
                panel.grid.major.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_line(
                    color = "black"
                ),
                axis.line = element_line(
                    color = "black"
                ),
                axis.title.y = element_text(
                    size = font_size,
                    angle = 90
                ),
                axis.text.y = element_text(
                    color = "black",
                    size = font_size,
                ),
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                strip.background = element_blank(),
                strip.text = element_text(
                    color = "black",
                    margin = margin(2, 0, 2, 0, "pt"),
                    size = font_size
                ),
                legend.title = element_text(
                    size = font_size
                ),
                legend.text = ggtext::element_markdown(
                    size = font_size
                ),
            )
    )

## Create the metadata bar ----
meta_squares <-
    mal2_metadata |>
    rownames_to_column(var = "SampleID") |>
    mutate(
        Condition = factor(
            Chamber,
            labels = c("Hypoxic 5", "Hypoxic 1", "Anoxic", "Reductive"),
            levels = c("A", "D", "B", "C")
        ),
        conditions = case_when(
            Chamber == "A" ~ 1,
            Chamber == "D" ~ 2,
            Chamber == "B" ~ 3,
            Chamber == "C" ~ 4
        )
    ) |>
    ggplot(
        mapping = aes(
            x = reorder(SampleID, conditions),
            y = 1,
            fill = Condition
        )
    ) +
    geom_col(
        width = 1,
        color = "black"
    ) +
    geom_text(
        mapping = aes(
            label = Sample,
            x = SampleID,
            y = 0.5
        ),
        color = "white",
        inherit.aes = FALSE
    ) +
    # scale_fill_brewer(
    #   palette = 'RdBu',
    #   direction = -1
    # ) +
    scale_fill_manual(
        values = viridis::viridis(n = 4, begin = 0.1, end = 0.9, direction = -1)
    ) +
    labs(fill = "Redox") +
    theme_void() +
    theme(
        plot.margin = margin(0, 0, 0, 0, "pt"),
        legend.position = "bottom",
        legend.text.position = "right",
        legend.title = element_text(
            size = font_size, hjust = 0.5
        ),
        legend.title.position = "bottom",
        legend.text = element_text(
            size = font_size
        )
    )
meta_squares

# Create a function to plot bar figures ----
plot_ggbar <- function(df, ec, genus_filt, y_limits = c(0, 20)) {
    ggbar <-
        df |>
        select(!(contains("-CPM"))) |>
        rename_with(
            .cols = where(is.numeric),
            .fn = function(x) {
                # Extract the sample names
                pre_name <- str_split_i(
                    string = x,
                    pattern = "_",
                    i = 1
                )
                # Refine samplennames
                chamber <- str_split_i(
                    string = pre_name,
                    pattern = "-",
                    i = 1
                )
                sample_id <- str_split_i(
                    string = pre_name,
                    pattern = "-",
                    i = 4
                )
                lane <- str_split_i(
                    string = pre_name,
                    pattern = "-",
                    i = 6
                )
                # Combine
                paste0(chamber, "_", sample_id, "_", lane)
            }
        ) |>
        rename(ec = `# Gene Family`) |>
        pivot_longer(
            cols = where(is.numeric),
            names_to = "Sample",
            values_to = "Abundance_cpm"
        ) |>
        pivot_wider(
            names_from = "ec",
            values_from = "Abundance_cpm"
        ) |>
        column_to_rownames(var = "Sample") |>
        select(contains(ec)) |>
        select(contains("|")) |>
        rownames_to_column(var = "Sample") |>
        pivot_longer(
            cols = where(is.numeric),
            names_to = "ec",
            values_to = "Abundance_cpm"
        ) |>
        separate(
            col = "ec",
            into = c("EC", "Taxon"),
            sep = "\\|"
        ) |>
        mutate(
            Species = str_split_i(
                string = Taxon,
                pattern = ".s__",
                i = 2
            ),
            Genus = str_split_i(
                string = Species,
                pattern = "_",
                i = 1
            ),
            Species = str_replace(
                string = Species,
                pattern = "_",
                replacement = " "
            ),
            Species = paste0("*", Species, "*"),
            Chamber = str_split_i(
                strin = Sample,
                pattern = "_",
                i = 1
            ),
            conditions = case_when(
                Chamber == "A" ~ 1,
                Chamber == "D" ~ 2,
                Chamber == "B" ~ 3,
                Chamber == "C" ~ 4
            )
        ) |>
        filter(Genus %in% genus_filt) |>
        # filter(Chamber != 'D') |>
        ggplot(
            mapping = aes(
                x = reorder(Sample, conditions),
                y = Abundance_cpm,
            )
        ) +
        geom_col(
            mapping = aes(
                fill = Species
            ),
            width = 1,
            color = "black"
        ) +
        facet_wrap(~EC) +
        labs(
            y = "Counts (cpm)"
        ) +
        scale_y_continuous(
            limits = y_limits,
            expand = expansion(mult = c(0, 0.1))
        ) +
        # scale_fill_brewer(palette = 'Set2') +
        scale_fill_viridis_d(option = "turbo", begin = 0.1, end = 0.9) +
        theme(
            panel.grid.major.y = element_blank(),
            plot.margin = margin(0, 0, 0, 0, "pt"),
            legend.position = "inside",
            legend.position.inside = c(0.55, 0.8)
        )
    return(ggbar)
}

# Visualize Desulfovibrio pathways ----
## sulfate adenylyltransferase ----
sulf_adetrans_bar <-
    plot_ggbar(ec_cpm, ec = "2.7.7.4:", genus_filt = "Desulfovibrio")
sulf_adetrans_bar

ggarrange(
    sulf_adetrans_bar,
    meta_squares + theme(plot.margin = margin(-4.5, 0, 0, 0, "pt")),
    ncol = 1, align = "v", heights = c(1, 0.125)
) + bgcolor("white")

ggsave(
    filename = "figures/sewageProj/sulfate_adenylyltransferase.png",
    units = "in", width = 8, height = 8.5, dpi = 320
)
