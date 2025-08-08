library(jsonlite)
library(dplyr)

# ---- inputs ----
file <- "../../inputs-sequoia/sequoia-vsp/0818-WOD-SOL.VSPv2.report.json"

# ---- helper ----
extract_microorganism_df <- function(microorganisms, sample_name) {
  df_na <- data.frame(
    sample = sample_name,
    name = NA_character_,
    alignedReadCount = NA_real_,
    ani = NA_real_,
    coverage = NA_real_,
    medianDepth = NA_real_,
    rpkm = NA_real_,
    stringsAsFactors = FALSE
  )
  
  if (length(microorganisms) > 0) {
    df_data <- do.call(rbind, lapply(microorganisms, function(micro) {
      data.frame(
        sample = sample_name,
        name = micro$name,
        alignedReadCount = micro$alignedReadCount,
        ani = micro$ani,
        coverage = micro$coverage,
        medianDepth = micro$medianDepth,
        rpkm = micro$rpkm,
        stringsAsFactors = FALSE
      )
    }))
    rbind(df_na, df_data)
  } else {
    df_na
  }
}

# ---- load one JSON ----
json_data <- fromJSON(file, simplifyDataFrame = FALSE)
microorganisms <- json_data$targetReport$microorganisms
sample_name <- sub("\\..*", "", basename(file))

df_single_sample <- extract_microorganism_df(microorganisms, sample_name)

df_single_sample$genus.name <- panel$Genus_Name[match(df_single_sample$name, panel$Reporting_Name)]
df_single_sample$strand.type <- panel$strand.type[match(df_single_sample$name, panel$Reporting_Name)]

df_single_sample <- df_single_sample %>% 
  filter(!is.na(name))

# ---- define spikes ----
bcov.spike <- "Human coronavirus OC43 (HCoV_OC43)"

tri.spike <- c("Human coronavirus OC43 (HCoV_OC43)", 
               "Human adenovirus F", 
               "Influenza A virus (H1N1)")

# ---- remove spike(s) ----
df_removed_spike <- df_single_sample %>% filter(!name %in% bcov.spike)

df_single_sample2 <- df_single_sample
df_single_sample2$relabund.rpkm <- as.numeric(df_single_sample$rpkm) / sum(as.numeric(df_single_sample$rpkm), na.rm = TRUE)

df_removed_spike$relabund.rpkm <- as.numeric(df_removed_spike$rpkm) / sum(as.numeric(df_removed_spike$rpkm), na.rm = TRUE)
