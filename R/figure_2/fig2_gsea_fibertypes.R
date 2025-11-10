source(here::here("R/Library.R"))

#open gsea_bp files of fibertype differences
f_fibertypes_bp <- vroom::vroom(here::here("results/gsea_bp_female_fibertypes.csv")) %>%
    mutate(sex = "female") %>%
    mutate(regulated = case_when(
        enrichmentScore > 0 ~ "mhc2",
        enrichmentScore < 0 ~ "mhc1"))

m_fibertypes_bp <- vroom::vroom(here::here("results/gsea_bp_male_fibertypes.csv")) %>%
    mutate(sex = "male") %>%
    mutate(regulated = case_when(
        enrichmentScore > 0 ~ "mhc2",
        enrichmentScore < 0 ~ "mhc1"))



#merge BP data frames
fibertypes_bp <- rbind(f_fibertypes_bp, m_fibertypes_bp)
fibertypes_bp$regulated <- factor(fibertypes_bp$regulated, levels = c("mhc1", "mhc2"))

#write.csv(fibertypes_bp,
#here::here("results/gsea_bp_fibertypes.csv"))


#Enrichment plot of BP for fibertype differences
gsea <- fibertypes_bp %>%
    dplyr::filter(Description %in% c(
        "striated muscle cell development",
        "proton transmembrane transport",
        "NADH regeneration",
        "mitochondrial gene expression",
        "cytoplasmic translation",
        "ATP metabolic process",
        "lipid catabolic process",
        "fatty acid metabolic process"
    )) |>
    ggplot(aes(x = sex, y = Description, color = qvalue, size = abs(NES))) +
    geom_point() +
    facet_wrap(~regulated, labeller = labeller(regulated = c("mhc1" = "Type I", "mhc2" = "Type II"))) +
    scale_x_discrete(breaks = c("female", "male"),
                     labels = c("Female", "Male" )) +
    theme_bw() +
    theme(text = ggplot2::element_text(size = 7),
          strip.text = ggplot2::element_text(size = 8),
          legend.key.size = ggplot2::unit(4, units = "mm")) +
    xlab("") +
    ylab("") +
    labs(color = "FDR", size = "NES")

#ggsave(plot = gsea, here::here('figures/figure_2/gsea_fibertypes.pdf'), height = 60, width = 120, units = "mm")


# LINEAR MIXED MODEL OF ENRICHED TERMS ------------------------------------

#load pre data
df_pre <- vroom::vroom(here::here("data/data_long_keywords.csv")) %>%
    filter(trial == "pre")

#Add keywords and GO terms
annotations <- vroom::vroom(here::here("data/keywords.csv")) %>%
    dplyr::rename_with(snakecase::to_snake_case) %>%
    dplyr::select(c("gene_names", "keywords", "gene_ontology_biological_process", "gene_ontology_cellular_component", "gene_ontology_molecular_function")) %>%
    dplyr::rename(gobp = gene_ontology_biological_process,
                  gocc = gene_ontology_cellular_component,
                  gomf = gene_ontology_molecular_function,
                  protein = gene_names) %>%
    dplyr::mutate(protein = gsub("\\ .*","", protein)) %>%
    dplyr::mutate(protein = make.names(protein, unique=TRUE), protein)

#create dataframe of mean expression of each protein pre
df_pre_mean <- df_pre %>%
    group_by(protein, group, fibertype, sex) %>%
    summarise(mean_expression = mean(expression, na.rm = TRUE)) %>%
    ungroup() %>%
    merge(annotations, by="protein", all.x = TRUE)

#create dataframes for enriched terms
df_translation <- df_pre_mean %>%
    filter(grepl("cytoplasmic translation",
                 gobp, ignore.case = TRUE) |
               protein %in% c("YBX1", "PKM", "RPL19", "RPL12", "RPS19", "FAU", "RPS8", "RPL13", "RPL18A", "RPL24", "RPL28", "RPS18", "RPS28", "RPL26L1", "RPS5", "RPL22", "RPS11", "YBX3",
                          "RPS25", "DENR", "CSDE1", "PAIP1", "RPL38", "RPLP2", "EIF3H", "DRG1", "RPL14", "RPL23", "RPS24", "RPS23", "RPL27A", "EIF2S2", "RPS29", "EIF4A2", "RPS6",
                          "RPL11", "EIF4G1", "RPL23A", "RPS26", "CNBP", "EIF2B3", "RPL4", "SH3BGRL", "EIF5", "RPS7", "DHX9", "RPL32", "RPL10", "EIF3L", "RPL5", "EIF4A1", "RPS3A",
                          "NCK1", "EIF3F", "SYNCRIP", "MIF4GD", "RPS4X"))

df_mito <- df_pre_mean %>%
    filter(grepl("mitochondrial gene expression",
                 gobp, ignore.case = TRUE) |
               protein %in% c(
        "C1QBP", "RCC1L", "SARS2", "MRPL3", "MRPL40", "MRPS25", "MRPL39", "MRPS34", "MRPL38", "MRPS10", "QRSL1", "YARS2", "MRPS9", "MRPL4", "MRPS18B", "DAP3", "TBRG4", "GFM1", "MRPL37",
        "FASTKD2", "MRPS28", "TRMT10C", "LRPPRC", "MRPS27", "MRPL15", "PRKAA1", "MRPS7", "IARS2", "MRPS23", "MRPL45", "MRPL46", "PUS1", "MTIF2", "NDUFA7", "TACO1", "TSFM", "HSD17B10",
        "MRPS22", "NGRN", "MRPS17", "MRPL44", "PTCD3", "COA3", "TUFM", "MRPL48", "TFAM", "DARS2", "MRPS30", "MRPS11", "MTRF1L", "TFB2M", "MRPS35", "LARS2", "UQCC2"
    ))

df_fa <- df_pre_mean %>%
    filter(grepl("fatty acid metabolic process",
                 gobp, ignore.case = TRUE) |
               protein %in% c(
        "PRKAA1", "THNSL2", "GPX1", "FABP3", "NDUFAB1", "ACSL1", "PCCB", "ETFDH", "SIRT2", "MGLL",
        "PCCA", "ACSF2", "ACADM", "CPT2", "ACADSB", "PTGR3", "HSD17B10", "CYGB", "MMUT", "RGN",
        "ACAA1", "ECHS1", "ACADVL", "ACAD9", "ECI1", "MLYCD", "PDK2", "HADH", "MGST3", "OXSM",
        "ABHD5", "HADHA", "ACOT2", "HADHB", "HSD17B8", "ACAT1", "GGT5", "AUH", "ACADS", "ECH1",
        "CRAT", "DECR1", "CPT1B", "AKR1C2", "ETFA", "PLIN5", "ECI2", "ETFB", "ACAA2", "ACSS1",
        "AKR1C3", "HTD2", "ACOT1", "CES2"
    ))



#linear mixed model of difference in cytoplasmic translation
lmm_translation <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = df_translation, REML = FALSE)
effects_translation <- anova(lmm_translation)

#linear mixed model of difference in mitochondrial gene expression
lmm_mito <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = df_mito, REML = FALSE)
effects_mito <- anova(lmm_mito)

#linear mixed model of difference in mitochondrial gene expression
lmm_fa <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = df_fa, REML = FALSE)
effects_fa <- anova(lmm_fa)

