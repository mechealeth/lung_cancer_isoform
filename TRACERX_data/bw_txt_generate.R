suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(fst)
  library(readr)
})

# ---- inputs ----
mdata_fst <- "/QRISdata/Q7816/TRACERX_outs/20221014_transcriptomic_DATA/2022-10-18all_metadata.fst"

# Folder containing your per-sample CPM bigWigs
bw_dir <- "/QRISdata/Q7816/TRACERX_bw"   # <-- FIXED (no duplicated prefix)
bw_suffix <- "_Aligned.sortedByCoord.out.bam.CPM.bw"                

mdata <- read_fst(mdata_fst)

# ---- build group list from mdata$region ----
df <- mdata %>%
  transmute(
    sample = region,          # <-- sample ID used to name bigWig
    patient = patient,
    tumour_id = tumour_id,
    Histology = Histology,

    # define normal vs LN vs tumour using region string
    is_ln     = str_detect(sample, "_LN01$"),
    is_normal = str_detect(sample, "_N01$") & !is_ln,

    group = case_when(
      is_normal ~ "Normal",
      Histology %in% c("LUAD", "LUSC") ~ Histology,
      TRUE ~ "Other"
    ),

    bw = file.path(bw_dir, paste0(sample, bw_suffix))
  ) %>%
  filter(!is_ln) %>%              # exclude LN01 always
  filter(file.exists(bw))         # keep only bigWigs that exist

# ---- write lists (one bigWig path per line) ----
write_lines(df %>% filter(group == "LUAD")   %>% pull(bw), "LUAD.txt")
write_lines(df %>% filter(group == "LUSC")   %>% pull(bw), "LUSC.txt")
write_lines(df %>% filter(group == "Normal") %>% pull(bw), "Normal.txt")

# quick sanity check
df %>% count(group)
