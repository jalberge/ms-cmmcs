---
title: "Extract post-clustering CCF from phylogicNDT"
author: "Ankit Dutta and JB Alberge"
date: February 16, 2022
output: 
  html_document:
    toc: true
    toc_depth: 2
    df_print: kable
    keep_md: true
---

Analyzing initial cmmcs pairs.


```r
source("../2_wgs/0_annotate_samples.R")
dir.create("figures")

paired.participants <- paste0("CTF", 
                              c("001", "002", "003", 
                                "012", "013", "015", "016",
                                "017", "018", "019", 
                                "021", "023", "025", 
                                "034", "036", "046", "047"))
```

# Pre-processing

Annotate variants and pivot data with 1 set of column per tissue. Here the difference between samples is tissue in CTCs, BMPCs


```r
# ldfr.fixed.protein.change <- read_tsv("../data/phylogic_aggregation_snv_Feb14.tsv")
ldfr.fixed.protein.change <- read_tsv("../data/phylogic_aggregation_snv_Mar14.tsv")
```

```
## Rows: 195216 Columns: 125
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr   (9): Patient_ID, Sample_ID, Hugo_Symbol, Reference_Allele, Tumor_Seq_A...
## dbl (113): Chromosome, Start_position, COSMIC_total_alterations_in_gene, t_r...
## lgl   (3): Sample_Alias, Allelic_CN_minor, Allelic_CN_major
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
ldfr.fixed.protein.change.clinical.annot <- ldfr.fixed.protein.change %>%
  inner_join(clinicaldata, by=c("Sample_ID"="entity:sample_id")) %>%
  mutate(Reclassification=case_when(
    Variant_Classification %in% non.synonymous ~ "Non-Silent",
    Variant_Classification == "Intron" ~ "Intron",
    TRUE ~ "Silent"
  ))

flat.matched.maf <- ldfr.fixed.protein.change.clinical.annot %>%
  filter(isRef==TRUE) %>%
  pivot_wider(id_cols = c(Patient_ID, Stage, Hugo_Symbol, Chromosome, Start_position, Reference_Allele, Tumor_Seq_Allele, Protein_Change, Variant_Classification, Variant_Type, COSMIC_overlapping_mutations, COSMIC_total_alterations_in_gene), names_from = tissue, values_from = c(starts_with("clust_ccf"), preDP_ccf_mean))
```

# Pre-clustering SNV density plot

Annotate post clustering


```r
plot.ccf.density(p.id = "CTF013",
                 maf = flat.matched.maf,
                 maf.labels = TRUE,
                 arrow = FALSE,
                 clustering = c("pre"),
                 n.bins = 30,
                 hex.color = "darkblue",
                 genes.of.interest=genes.of.interest,
                 min.COSMIC.total.alterations.in.gene=10000,
                 non.silent.only=TRUE,
                 save = TRUE,
                 plot = TRUE)
```

```
## Found 2 events to annotate out of 4825 total variants.
```

![](1_post_clustering_SNVs_matched_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

<!-- # Post-clustering CCF and confidence interval -->

<!-- Using old maf since no forcecalling there! -->

<!-- ```{r} -->
<!-- last.maf <- read_maf("../data/Paper2022_v2_filtered_Mar9.maf") -->
<!-- annot.maf <- last.maf %>% inner_join(clinicaldata, by=c("Tumor_Sample_Barcode")) %>% filter(isRef==TRUE) -->

<!-- annot.maf$matchedID <- paste0(annot.maf$Hugo_Symbol, " ", annot.maf$Protein_Change, " ", annot.maf$participant) -->

<!-- # # flat.matched.maf -->
<!-- # # prep.ccf.matched <- short.CTF.pair.maf %>%  -->
<!-- prep.ccf.matched <- flat.matched.maf %>% -->
<!--   rowwise() %>% -->
<!--   mutate(ID=paste0(Hugo_Symbol, " ", Protein_Change, " ", Patient_ID), -->
<!--          twolinesID=paste0(Hugo_Symbol, " ", Protein_Change, "\n", Patient_ID)) -->

<!-- # TODO make this code more efficient -->
<!-- tiny.prep.ccf.matched <-   prep.ccf.matched %>% -->
<!--   filter(Hugo_Symbol %in% genes.of.interest) %>% -->
<!--   mutate(is_in_CMMCs = ID %in% subset(annot.maf, tissue=="CMMCs" & participant==Patient_ID)$matchedID, -->
<!--          is_in_BMPCs = ID %in% subset(annot.maf, tissue=="BMPCs" & participant==Patient_ID)$matchedID, -->
<!--          CMMC_V = ifelse((!is_in_CMMCs & is_in_BMPCs), "*", ""), -->
<!--          BMPC_V = ifelse((!is_in_BMPCs & is_in_CMMCs), "*", "")) -->

<!-- tiny.prep.ccf.matched %>% -->
<!--   filter((clust_ccf_mean_BMPCs<.9|clust_ccf_mean_CMMCs<.9)) %>% # non clonal mutations -->
<!--   ungroup() %>% -->
<!--   summarise(Med_BMPCs=median(clust_ccf_mean_BMPCs),  -->
<!--             Med_CMMCs=median(clust_ccf_mean_CMMCs), -->
<!--             Q25_BMPCs=quantile(clust_ccf_mean_BMPCs, .25), -->
<!--             Q25_CMMCs=quantile(clust_ccf_mean_CMMCs, .25), -->
<!--             Q75_BMPCs=quantile(clust_ccf_mean_BMPCs, .75), -->
<!--             Q75_CMMCs=quantile(clust_ccf_mean_CMMCs, .75)) -->
<!-- # View(prep.ccf.matched %>% select(1:10, is_in_CMMCs, is_in_BMPCs)) -->

<!-- prep.ccf.matched %>% -->
<!--   filter((clust_ccf_mean_BMPCs<.9|clust_ccf_mean_CMMCs<.9)) %>% # non clonal mutations -->
<!--   filter(Hugo_Symbol %in% genes.of.interest & Variant_Classification %in% c("Intron", "Silent")) %>%# non clonal mutations   -->
<!--   ungroup() %>% -->
<!--   summarise(Med_BMPCs=median(clust_ccf_mean_BMPCs),  -->
<!--             Med_CMMCs=median(clust_ccf_mean_CMMCs), -->
<!--             Q25_BMPCs=quantile(clust_ccf_mean_BMPCs, .25), -->
<!--             Q25_CMMCs=quantile(clust_ccf_mean_CMMCs, .25), -->
<!--             Q75_BMPCs=quantile(clust_ccf_mean_BMPCs, .75), -->
<!--             Q75_CMMCs=quantile(clust_ccf_mean_CMMCs, .75)) -->
<!-- # View(prep.ccf.matched %>% select(1:10, is_in_CMMCs, is_in_BMPCs)) -->
<!-- ``` -->


<!-- ## Comparison of CCF post-clustering -->

<!-- ```{r} -->
<!-- paired.mut <- tiny.prep.ccf.matched %>% -->
<!--   ggplot(aes(x=fct_reorder(twolinesID, clust_ccf_mean_CMMCs))) + -->
<!--   # ggplot(aes(x=fct_reorder(paste0(Hugo_Symbol, "\n", Patient_ID), clust_ccf_mean_CMMCs))) +  -->
<!--   geom_bar(aes(y=clust_ccf_mean_BMPCs), fill=paired.pals["BMPCs"], stat = "identity", position=position_nudge(x=0.2), width=0.4) + -->
<!--   geom_point(aes(y=clust_ccf_mean_BMPCs), fill=paired.pals["BMPCs"], position=position_nudge(x=0.2), size=0.3, color="#444444") + -->
<!--   # geom_point(aes(y=preDP_ccf_mean_BMPCs), fill=paired.pals["BMPCs"], position=position_nudge(x=0.2), size=1, shape=2, color="#3B3B3B") + -->
<!--   geom_linerange(aes(ymin=clust_ccf_CI_low_BMPCs, ymax=clust_ccf_CI_high_BMPCs), size=0.5, position=position_nudge(x=0.2), color="#444444") + -->
<!--   geom_text(aes(label=CMMC_V), y=0.02, position=position_nudge(x=0.2), size=3) + -->
<!--   geom_bar(aes(y=clust_ccf_mean_CMMCs), fill=paired.pals["CMMCs"], stat = "identity", position=position_nudge(x=-0.2), width=0.4) + -->
<!--   geom_point(aes(y=clust_ccf_mean_CMMCs), fill=paired.pals["CMMCs"], position=position_nudge(x=-0.2), size=0.3) + -->
<!--   # geom_point(aes(y=preDP_ccf_mean_CMMCs), fill=paired.pals["CMMCs"], position=position_nudge(x=-0.2), size=1, shape=2, color="#3B3B3B") + -->
<!--   geom_linerange(aes(ymin=clust_ccf_CI_low_CMMCs, ymax=clust_ccf_CI_high_CMMCs), size=0.5, position=position_nudge(x=-0.2), color="#444444") + -->
<!--   geom_text(aes(label=BMPC_V), y=0.02, position=position_nudge(x=-0.2), size=3) + -->
<!--   scale_y_continuous(labels=scales::percent, limits = c(0, 1), expand=expansion(add=c(0, 0.05))) + -->
<!--   # scale_color_manual(values=paired.pals) + -->
<!--   # scale_fill_manual(values=paired.pals) + -->
<!--   labs(y="Cancer cell fraction", x="", fill="") + -->
<!--   theme_bw() + -->
<!--   # coord_flip() + -->
<!--   theme(panel.grid = element_blank(),  -->
<!--         text = element_text(size=7), -->
<!--         # aspect.ratio = .5,  -->
<!--         axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5), -->
<!--         legend.position = "none") -->

<!-- paired.summarized.mut <- prep.ccf.matched %>% -->
<!--   filter((clust_ccf_mean_BMPCs<.9 | clust_ccf_mean_CMMCs<.9)) %>% # non clonal mutations -->
<!--   filter(Hugo_Symbol %in% genes.of.interest & Variant_Classification %in% c("Intron", "Silent")) %>%# non clonal mutations -->
<!--   # filter(Variant_Classification %nin% non.synonymous) %>%  -->
<!--   pivot_longer(cols=c(clust_ccf_mean_BMPCs, clust_ccf_mean_CMMCs), names_to = "tissue", values_to = "Mean_CCF") %>% -->
<!--   mutate(tissue=str_extract(tissue, "(CMMCs|BMPCs)")) %>% -->
<!--   ggplot(aes(x=tissue, y=Mean_CCF, fill=tissue)) + -->
<!--   geom_boxplot(size=0.3, color="#3B3B3B") + -->
<!--   scale_y_continuous(labels=scales::percent, limits = c(0, 1), expand=expansion(add=c(0, 0.05))) + -->
<!--   stat_compare_means(paired = TRUE,  label = "p.signif", tip.length = 0, symnum.args = symnum.args, label.x = 1.5, size=3) + -->
<!--   scale_fill_manual(values=paired.pals) + -->
<!--   labs(y="Cancer cell fraction", x="", fill="") + -->
<!--   theme_bw() + -->
<!--   theme(panel.grid = element_blank(),  -->
<!--         text = element_text(size=7), -->
<!--         # aspect.ratio = 2,  -->
<!--         axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5), -->
<!--         legend.position = "none") -->

<!-- tinypaired.summarized.mut <- prep.ccf.matched %>% -->
<!--   # filter((clust_ccf_mean_BMPCs<.9|clust_ccf_mean_CMMCs<.9)) %>% # non clonal mutations -->
<!--   filter(Hugo_Symbol %in% genes.of.interest & Variant_Classification %in% non.synonymous) %>% -->
<!--   pivot_longer(cols=c(clust_ccf_mean_BMPCs, clust_ccf_mean_CMMCs), names_to = "tissue", values_to = "Mean_CCF") %>% -->
<!--   mutate(tissue=str_extract(tissue, "(CMMCs|BMPCs)")) %>% -->
<!--   ggplot(aes(x=tissue, y=Mean_CCF, fill=tissue)) + -->
<!--   geom_boxplot(size=0.3, color="#3B3B3B") + -->
<!--   # geom_line(aes(group=ID)) + -->
<!--   scale_y_continuous(labels=scales::percent, limits = c(0, 1), expand=expansion(add=c(0, 0.05))) + -->
<!--   stat_compare_means(paired = TRUE, label = "p.signif", tip.length = 0, symnum.args = symnum.args, label.x = 1.5, size=3) + -->
<!--   scale_fill_manual(values=paired.pals) + -->
<!--   labs(y="Cancer cell fraction", x="", fill="") + -->
<!--   facet_wrap(~Stage) + -->
<!--   theme_bw() + -->
<!--   theme(panel.grid = element_blank(),  -->
<!--         text = element_text(size=7), -->
<!--         # aspect.ratio = 2,  -->
<!--         axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5), -->
<!--         legend.position = "none") -->

<!-- driver.plot <- plot_grid(paired.mut, tinypaired.summarized.mut, paired.summarized.mut, nrow = 1, align = "h", axis = "blrt", rel_widths = c(4, 1, 1)) -->

<!-- driver.plot -->

<!-- ggsave("figures/paired-mut-drivers-row.pdf", driver.plot, width = 7, height = 3) -->
<!-- ggsave("figures/paired-mut-drivers-row.png", driver.plot, width = 7, height = 3) -->


<!-- ``` -->

<!-- ## Comparison of CCF PRE-clustering -->

<!-- ```{r} -->
<!-- paired.mut <- tiny.prep.ccf.matched %>% -->
<!--   ggplot(aes(x=fct_reorder(twolinesID, preDP_ccf_mean_CMMCs))) + -->
<!--   # ggplot(aes(x=fct_reorder(paste0(Hugo_Symbol, "\n", Patient_ID), preDP_ccf_mean_CMMCs))) +  -->
<!--   geom_bar(aes(y=preDP_ccf_mean_BMPCs), fill=paired.pals["BMPCs"], stat = "identity", position=position_nudge(x=0.2), width=0.4) + -->
<!--   geom_point(aes(y=preDP_ccf_mean_BMPCs), fill=paired.pals["BMPCs"], position=position_nudge(x=0.2), size=0.3, color="#444444") + -->
<!--   # geom_point(aes(y=preDP_ccf_mean_BMPCs), fill=paired.pals["BMPCs"], position=position_nudge(x=0.2), size=1, shape=2, color="#3B3B3B") + -->
<!--   # geom_linerange(aes(ymin=clust_ccf_CI_low_BMPCs, ymax=clust_ccf_CI_high_BMPCs), size=0.5, position=position_nudge(x=0.2), color="#444444") + -->
<!--   geom_text(aes(label=CMMC_V), y=0.02, position=position_nudge(x=0.2), size=3) + -->
<!--   geom_bar(aes(y=preDP_ccf_mean_CMMCs), fill=paired.pals["CMMCs"], stat = "identity", position=position_nudge(x=-0.2), width=0.4) + -->
<!--   geom_point(aes(y=preDP_ccf_mean_CMMCs), fill=paired.pals["CMMCs"], position=position_nudge(x=-0.2), size=0.3) + -->
<!--   # geom_point(aes(y=preDP_ccf_mean_CMMCs), fill=paired.pals["CMMCs"], position=position_nudge(x=-0.2), size=1, shape=2, color="#3B3B3B") + -->
<!--   # geom_linerange(aes(ymin=clust_ccf_CI_low_CMMCs, ymax=clust_ccf_CI_high_CMMCs), size=0.5, position=position_nudge(x=-0.2), color="#444444") + -->
<!--   geom_text(aes(label=BMPC_V), y=0.02, position=position_nudge(x=-0.2), size=3) + -->
<!--   scale_y_continuous(labels=scales::percent, limits = c(0, 1), expand=expansion(add=c(0, 0.05))) + -->
<!--   # scale_color_manual(values=paired.pals) + -->
<!--   # scale_fill_manual(values=paired.pals) + -->
<!--   labs(y="Cancer cell fraction", x="", fill="") + -->
<!--   theme_bw() + -->
<!--   # coord_flip() + -->
<!--   theme(panel.grid = element_blank(),  -->
<!--         text = element_text(size=7), -->
<!--         # aspect.ratio = .5,  -->
<!--         axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5), -->
<!--         legend.position = "none") -->

<!-- paired.summarized.mut <- prep.ccf.matched %>% -->
<!--   filter((preDP_ccf_mean_BMPCs<.9|preDP_ccf_mean_CMMCs<.9)) %>% # non clonal mutations -->
<!--   filter(Hugo_Symbol %in% genes.of.interest & Variant_Classification %in% c("Intron", "Silent")) %>%# non clonal mutations -->
<!--   pivot_longer(cols=c(preDP_ccf_mean_BMPCs, preDP_ccf_mean_CMMCs), names_to = "tissue", values_to = "Mean_CCF") %>% -->
<!--   mutate(tissue=str_extract(tissue, "(CMMCs|BMPCs)")) %>% -->
<!--   ggplot(aes(x=tissue, y=Mean_CCF, fill=tissue)) + -->
<!--   geom_boxplot(size=0.3, color="#3B3B3B") + -->
<!--   scale_y_continuous(labels=scales::percent, limits = c(0, 1), expand=expansion(add=c(0, 0.05))) + -->
<!--   stat_compare_means(paired = TRUE,  label = "p.signif", tip.length = 0, symnum.args = symnum.args, label.x = 1.5, size=3) + -->
<!--   scale_fill_manual(values=paired.pals) + -->
<!--   labs(y="Cancer cell fraction", x="", fill="") + -->
<!--   theme_bw() + -->
<!--   theme(panel.grid = element_blank(),  -->
<!--         text = element_text(size=7), -->
<!--         # aspect.ratio = 2,  -->
<!--         axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5), -->
<!--         legend.position = "none") -->

<!-- tinypaired.summarized.mut <- tiny.prep.ccf.matched %>% -->
<!--   filter((preDP_ccf_mean_BMPCs<.9|preDP_ccf_mean_CMMCs<.9)) %>% # non clonal mutations -->
<!--   pivot_longer(cols=c(preDP_ccf_mean_BMPCs, preDP_ccf_mean_CMMCs), names_to = "tissue", values_to = "Mean_CCF") %>% -->
<!--   mutate(tissue=str_extract(tissue, "(CMMCs|BMPCs)")) %>% -->
<!--   ggplot(aes(x=tissue, y=Mean_CCF, fill=tissue)) + -->
<!--   geom_boxplot(size=0.3, color="#3B3B3B") + -->
<!--   geom_line(aes(group=ID)) + -->
<!--   scale_y_continuous(labels=scales::percent, limits = c(0, 1), expand=expansion(add=c(0, 0.05))) + -->
<!--   stat_compare_means(paired = TRUE, label = "p.signif", tip.length = 0, symnum.args = symnum.args, label.x = 1.5, size=3) + -->
<!--   scale_fill_manual(values=paired.pals) + -->
<!--   labs(y="Cancer cell fraction", x="", fill="") + -->
<!--   theme_bw() + -->
<!--   theme(panel.grid = element_blank(),  -->
<!--         text = element_text(size=7), -->
<!--         # aspect.ratio = 2,  -->
<!--         axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5), -->
<!--         legend.position = "none") -->

<!-- driver.plot <- plot_grid(paired.mut, tinypaired.summarized.mut, paired.summarized.mut, nrow = 1, align = "h", axis = "blrt", rel_widths = c(4, 1, 1)) -->

<!-- driver.plot -->

<!-- ggsave("figures/paired-mut-drivers-row.pdf", driver.plot, width = 7, height = 3) -->
<!-- ggsave("figures/paired-mut-drivers-row.png", driver.plot, width = 7, height = 3) -->


<!-- ``` -->

