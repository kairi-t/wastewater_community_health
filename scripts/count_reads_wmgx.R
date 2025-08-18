# Programmed by: Kairi Tanaka
# Programmed on: 2024-04-02
# Programmed to: Count the number of reads from the 16S sequencing data

# Load necessary libraries
library(ShortRead)
library(dplyr)

# Set working directory
setwd("/Volumes/zlyu3/raw_sequencing_data/WGS/sewage/10_08_2022_TSU/wgs_raw_data/")

# Set the folder path containing FASTQ files
folder_path <- "./unzip_files"
# folder_path <- "./gzip_files"

# Get a list of FASTQ files in the folder
fastq_files <- list.files(folder_path, pattern = "\\.fastq(.gz)?$", full.names = TRUE)

# Function to count reads in a single file
count_reads <- function(file_path) {
    reads <- readFastq(file_path)
    length(reads)
}

# Apply the function to each file and create a data frame
read_counts <- tibble(
    file_name = basename(fastq_files), # Extract file names (sample names)
    num_reads = sapply(fastq_files, count_reads)
)

print(read_counts)

# Save data
saveRDS(read_counts, file = "./processed_data/summary_tables/read_counts_wmgx_raw.rds")
xlsx::write.xlsx(
    read_counts,
    file = "./processed_data/summary_tables/read_counts_wmgx_raw.xlsx"
)

# Visualize
read_counts |>
    ggplot(
        mapping = aes(
            x = num_reads,
            y = file_name
        )
    ) +
    geom_col(
        fill = "white",
        color = "black"
    ) +
    geom_text(
        mapping = aes(
            label = num_reads,
            x = num_reads / 2
        )
    ) +
    scale_x_continuous(
        expand = expansion(mult = c(0, 0.1))
    ) +
    labs(
        y = "",
        x = "Number of reads"
    )

# Save plot
ggsave(
    filename = "./processed_data/figures/num_reads_wmgx.png",
    units = "in", width = 8, height = 24, dpi = 300
)
