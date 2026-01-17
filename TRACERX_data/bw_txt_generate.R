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


######
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(fst)
  library(readr)
})

# ---- inputs ----
mdata_fst <- "/QRISdata/Q7816/TRACERX_outs/20221014_transcriptomic_DATA/2022-10-18all_metadata.fst"
bw_dir    <- "/QRISdata/Q7816/TRACERX_bw"   # FIX: remove duplicated path
bw_suffix <- ".cpm.bw"                      # adjust if needed

mdata <- read_fst(mdata_fst)

# ---- 1) list all bigWigs that exist ----
bw_files <- list.files(bw_dir, pattern = "\\.bw$", full.names = TRUE)

bw_df <- tibble(
  bw = bw_files,
  bw_base = basename(bw_files)
) %>%
  # sample id = remove suffix (".cpm.bw" etc)
  mutate(sample = str_replace(bw_base, paste0(str_replace_all(bw_suffix, "\\.", "\\\\."), "$"), "")) %>%
  # if your filenames are like "..._Aligned.sortedByCoord.out.bam.cpm.bw",
  # and you want region like "CRUK0768_SU_N01", strip the STAR tail:
  mutate(region = str_replace(sample, "_Aligned.*$", "")) %>%
  mutate(
    patient = str_extract(region, "^CRUK\\d+"),
    is_ln     = str_detect(region, "_LN\\d+$"),          # LN01
    is_normal = str_detect(region, "_N\\d+$") & !is_ln   # N01
  ) %>%
  filter(!is_ln)  # always exclude LN

# ---- 2) join tumour histology when possible ----
# mdata region looks like "CRUK0005_SU_T1-R1" etc.
# We'll join by region string directly.
bw_df2 <- bw_df %>%
  left_join(
    mdata %>% select(region, Histology, tumour_id, patient),
    by = "region",
    suffix = c("", "_mdata")
  ) %>%
  mutate(
    group = case_when(
      is_normal ~ "Normal",
      Histology %in% c("LUAD", "LUSC") ~ Histology,
      TRUE ~ "Other"
    )
  )

# ---- 3) write lists ----
write_lines(bw_df2 %>% filter(group == "Normal") %>% pull(bw), "Normal.txt")
#write_lines(bw_df2 %>% filter(group == "LUAD")   %>% pull(bw), "LUAD.txt")
#write_lines(bw_df2 %>% filter(group == "LUSC")   %>% pull(bw), "LUSC.txt")
#write_lines(bw_df2 %>% filter(group == "Other")  %>% pull(bw), "Other.txt")
