# Goal: To remove centromere and big blacklisted regions

# INIT --------------------------------------------------------------------
source("0_annotate_samples.R")

# Load CNV ----------------------------------------------------------------

# if ABSOLUTE
# cnv.pattern=".AllelicCapSeg_PP_CCF_fit_v3.ABSOLUTE.segtab.txt"
# if non collapsed tau
cnv.pattern=".acs.ccf.tsv"

ccf.cnv.path='../data/alleliccapseg_pp_ccf_fit_v3_collapsed_out/'

ccf.cnv <- list.files(ccf.cnv.path, pattern = cnv.pattern)
ccf.dlist <- lapply(paste0(ccf.cnv.path, ccf.cnv), read_tsv)
names(ccf.dlist) <- str_remove(ccf.cnv, cnv.pattern)
#if non collapsed tau
ccf.dlist <- lapply(ccf.dlist, function(x) { x$gender <- as.character(x$gender); x})

# Load exclusion list -----------------------------------------------------

# default is that 3MB are removed around centromere but gatk has a larger file. 
# Will only remove large segments here to avoid slicing bed file into thousands of pieces

exclude <- read_tsv("../0_utils/CNV_and_centromere_blacklist.hg19.bed", col_names = c("Chr", "Start", "End"))
exclude$Size <- exclude$End-exclude$Start +1
exclude.large.only <- exclude %>% filter(Size >= 2.9E6)
write_tsv(exclude.large.only, "../0_utils//CNV_and_centromere_blacklist_3MB.bed", col_names = FALSE)

# Assemble CNV ------------------------------------------------------------

# purity and ploidy are extracted from the clinical sample sheet (absolute) and then written again in acs ccf output. avoid repeat
ccf.df <- rbindlist(ccf.dlist, idcol = "Pair") %>% 
  relocate(Pair) %>%
  # select(-ploidy, -purity) %>% 
  inner_join(select(clinicaldata, c(-ploidy, -purity)), by=c("Pair"="entity:sample_id"))

ccf.df.clean <- ccf.df %>%
  mutate(Chromosome = factor(Chromosome, levels=1:23, labels=c( as.character(1:22), "X")),
         Start=Start.bp, 
         End=End.bp, 
         Avg_Ploidy = purity * ploidy + 2 * (1 - purity),
         Segment_Mean=1/purity*(round(Avg_Ploidy)/2*tau-2*(1-purity)))
# formula adapted from Carter 2012 absolute paper.
# Can also use columns CN_minor / CN_major instead.

# Subtract CNV / Blacklisted centromeres ----------------------------------

ccf.bed <- ccf.df.clean %>%
  select(Chromosome, Start, End, Pair, Segment_Mean)
write_delim(ccf.bed, "data/tmp.cnv.2.bed", delim = "\t", col_names=FALSE)

# TODO how to return line?
system('bedtools subtract -a data/tmp.cnv.2.bed -b ../0_utils/CNV_and_centromere_blacklist_3MB.bed > data/tmp.cnv.2.regions.excluded.bed')

CNV.bed.clean <- read_tsv("data/tmp.cnv.2.regions.excluded.bed",
                          col_names = c("Chromosome", "Start", "End", "entity:sample_id", "Segment_Mean"))

# CNV.data.clean <- CNV.bed.clean %>% left_join(clinicaldata, by=c("entity:sample_id"))

# save as a seg file
CNV.bed.clean %>% 
  relocate(`entity:sample_id`) %>% 
  relocate(Segment_Mean, .after = last_col()) %>%
  write_tsv("data/Sept9_reprocessed.cnv.revisions3.seg")