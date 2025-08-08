library(readr)
library(dplyr)

# ---- single Kraken file ----
kraken_file <- "../../inputs-sequoia/kraken/1130_ESP_SOL_S235_L002.krk.txt"

# Read and process single file
kraken_single <- read_tsv(
  kraken_file,
  col_names = c(
    "percentage",
    "reads_clade",
    "reads_taxon",
    "rank_code",
    "taxon_ID",
    "taxon_name"
  ),
  show_col_types = FALSE
) %>%
  mutate(
    sample = sub("^(.*)_.*_.*$", "\\1", sub("\\..*", "", basename(kraken_file))),
    sample = gsub("_", "-", sample)  # match DRAGEN format
  ) %>%
  select(sample, everything())

# Filter for PMMoV (taxon_ID == 12239) and keep only relevant columns
kraken_pmmov <- kraken_single %>%
  filter(taxon_ID == 12239) %>%
  select(
    sample,
    percent.pmmov.kraken = percentage,
    pmmov.reads.kraken = reads_taxon
  )
