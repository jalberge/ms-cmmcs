source("../2_wgs/0_annotate_samples.R")
dir.create("figures", showWarnings = FALSE)

# STEP 1
# Rescue correct protein change column from total MAF file
last.maf <- read_maf("../data/AllRevisions.txt")

# STEP 2
# Aggregate single phylogic maf files

# STEP 2.1 matched BM PB
ldfr.matched <- list.files("../data/phylogic-bm-pb-output/", "mut_ccfs.txt", full.names = TRUE) %>%
  map(fread, sep = "\t", stringsAsFactors = FALSE, verbose = FALSE, 
      data.table = TRUE, showProgress = TRUE, header = TRUE, fill = TRUE, 
      skip = "Patient_ID", quote = "") %>% 
  bind_rows()
# STEP 2.2 londitufinal
ldfr.longitudinal <- list.files("../data/phylogic-longitudinal-output/", "mut_ccfs.txt", full.names = TRUE) %>%
  map(fread, sep = "\t", stringsAsFactors = FALSE, verbose = FALSE, 
      data.table = TRUE, showProgress = TRUE, header = TRUE, fill = TRUE, 
      skip = "Patient_ID", quote = "") %>% 
  bind_rows()
# STEP 2.3 matched and serial
ldfr.matched.longitudinal <- list.files("../data/phylogic-serial-and-bm-pb-output/", "mut_ccfs.txt", full.names = TRUE) %>%
  map(fread, sep = "\t", stringsAsFactors = FALSE, verbose = FALSE, 
      data.table = TRUE, showProgress = TRUE, header = TRUE, fill = TRUE, 
      skip = "Patient_ID", quote = "") %>% 
  bind_rows()

# STEP 3
# fix protein change
# Connor has a tool, too, for uploading into Terra directly - See One Note script
last.maf.protein.change <- last.maf %>%
  filter(Protein_Change!="") %>%
  select(Chromosome, Start_position, Reference_Allele, Tumor_Seq_Allele2, Protein_Change, COSMIC_overlapping_mutations, COSMIC_total_alterations_in_gene) %>%
  distinct()

# 3.1
ldfr.fixed.protein.change.matched <- ldfr.matched %>% 
  mutate(Chromosome = as.character(Chromosome)) %>%
  select(-Protein_change) %>% 
  left_join(last.maf.protein.change, by=c("Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele"="Tumor_Seq_Allele2")) %>%
  relocate(Protein_Change, COSMIC_overlapping_mutations, COSMIC_total_alterations_in_gene, .after=Tumor_Seq_Allele)

# 3.2
ldfr.fixed.protein.change.longitudinal <- ldfr.longitudinal %>% 
  mutate(Chromosome = as.character(Chromosome)) %>%
  select(-Protein_change) %>% 
  left_join(last.maf.protein.change, by=c("Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele"="Tumor_Seq_Allele2")) %>%
  relocate(Protein_Change, COSMIC_overlapping_mutations, COSMIC_total_alterations_in_gene, .after=Tumor_Seq_Allele)

# 3.3
ldfr.fixed.protein.change.matched.longitudinal <- ldfr.matched.longitudinal %>% 
  mutate(Chromosome = as.character(Chromosome)) %>%
  select(-Protein_change) %>% 
  left_join(last.maf.protein.change, by=c("Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele"="Tumor_Seq_Allele2")) %>%
  relocate(Protein_Change, COSMIC_overlapping_mutations, COSMIC_total_alterations_in_gene, .after=Tumor_Seq_Allele)

#save all
write_tsv(ldfr.fixed.protein.change.matched, "../data/phylogic_aggregation_snv_Sept9.tsv")
write_tsv(ldfr.fixed.protein.change.longitudinal, "../data/phylogic_aggregation_snv_Mar14_longitudinal.tsv")
write_tsv(ldfr.fixed.protein.change.matched.longitudinal, "../data/phylogic_aggregation_snv_Mar14_matched_longitudinal.tsv")
