source(here::here("R/Library.R"))

#load pre data
df_pre <- readRDS(here::here("data/data_long_keywords.rds")) %>%
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

# CYTOSOLIC RIBOSOME ------------------------------------------

# Make fibertype a factor with type II as reference
df_pre_mean$fibertype <- factor(df_pre_mean$fibertype, levels = c("mhc2", "mhc1"))

#filter for cytosolic ribosome
ribo_proteins <- df_pre_mean %>%
    filter(
        grepl("cytosolic ribosome", gocc, ignore.case = TRUE))

#linear mixed model of difference in intermediate filament protein abundance
lmm_ribo <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = ribo_proteins, REML = FALSE)
effects_ribo <- anova(lmm_ribo)

#calculate means for each fibertype for each sex
means_ribo <- emmeans(lmm_ribo, ~ sex | fibertype)
print(means_ribo)

#run linear model of sex for each fibertype
lm_ribo_sex <- contrast(means_ribo, method = "pairwise", by = "fibertype")
print(lm_ribo_sex)

#run linear model of fibertype for each sex
lm_ribo_fibertype <- contrast(means_ribo, method = "pairwise", by = "sex")
print(lm_ribo_fibertype)

# MITOCHONDRIA ------------------------------------------------------------

#filter for mitochondria
mito_proteins <- df_pre_mean %>%
    filter(
        grepl("mitochondrion", gocc, ignore.case = TRUE))

#linear mixed model of difference in intermediate filament protein abundance
lmm_mito <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = mito_proteins, REML = FALSE)
effects_mito <- anova(lmm_mito)

#calculate means for each fibertype for each sex
means_mito <- emmeans(lmm_mito, ~ sex | fibertype)
print(means_mito)

#run linear model of sex for each fibertype
lm_mito_sex <- contrast(means_mito, method = "pairwise", by = "fibertype")
print(lm_mito_sex)

#run linear model of fibertype for each sex
lm_mito_fibertype <- contrast(means_mito, method = "pairwise", by = "sex")
print(lm_mito_fibertype)


#Compute emmeans from linear mixed model
emm_ribo <- emmeans(lm_ribo_fibertype, ~ sex | sex) %>%
    as.data.frame() %>%
    mutate(cc = "Ribosome") %>%
    mutate(
        fibertype = case_when(
            estimate < 0 ~ "type1",
            estimate > 0 ~ "type2"))


emm_mito <- emmeans(lm_mito_fibertype, ~ sex | sex) %>%
    as.data.frame() %>%
    mutate(cc = "Mitochondria") %>%
    mutate(
        fibertype = case_when(
            estimate < 0 ~ "type1",
            estimate > 0 ~ "type2"))

emm_cc <- rbind(emm_ribo, emm_mito)


#define text for figure
p_cc <- tibble(
    sex = c("female", "female",
            "male", "male"),
    x = c(1, 2,
          1, 2),
    y = c(0.03, 0.20,
          0.03, 0.24),
    label = c("p<0.001", "p<0.001",
              "p<0.001", "p<0.001")
)

##MITOCHONDRIA AND RIBOSOME FIGURE##
cc_fig <- ggplot(emm_cc, aes(x = cc, y = estimate, fill = fibertype)) +
    geom_col(position = position_dodge(width = 0.9),
             width = 0.9, color = NA, alpha = 0.7) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, position = position_dodge(width = 0.9),
                  linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    scale_fill_manual(values=c("#440154FF", "#67CC5CFF"))+
    scale_x_discrete(labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    facet_wrap(~sex, labeller = labeller(sex = c("female" = "Females", "male" = "Males"))) +
    #scale_y_continuous(limits = c(-0.5, 2)) +
    #geom_segment(data = brackets_cc, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_cc, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab(NULL) +
    ylab("Log2fold difference (type II - type I)") +
    ggtitle("Mitochondria and ribosomes")

ggsave(plot = cc_fig, here::here('figures/figure_2/mito_ribo_figure.pdf'), height = 60, width = 90, units = "mm")


# GLYCOLYSIS --------------------------------------------------------------
#filter for glycolysis
glyco_proteins <- df_pre_mean %>%
    filter(
        grepl("glycolysis", gobp, ignore.case = TRUE))

#linear mixed model of difference in intermediate filament protein abundance
lmm_glyco <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = glyco_proteins, REML = FALSE)
effects_glyco <- anova(lmm_glyco)

#calculate means for each fibertype for each sex
means_glyco <- emmeans(lmm_glyco, ~ sex | fibertype)
print(means_glyco)

#run linear model of sex for each fibertype
lm_glyco_sex <- contrast(means_glyco, method = "pairwise", by = "fibertype")
print(lm_glyco_sex)

#run linear model of fibertype for each sex
lm_glyco_fibertype <- contrast(means_glyco, method = "pairwise", by = "sex")
print(lm_glyco_fibertype)


# FATTY ACID METABOLISM ---------------------------------------------------

#filter for fatty acid metabolism
fa_proteins <- df_pre_mean %>%
    filter(
        grepl("fatty acid metabolic process", gobp, ignore.case = TRUE))

#linear mixed model of difference in intermediate filament protein abundance
lmm_fa <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = fa_proteins, REML = FALSE)
effects_fa <- anova(lmm_fa)

#calculate means for each fibertype for each sex
means_fa <- emmeans(lmm_fa, ~ sex | fibertype)
print(means_fa)

#run linear model of sex for each fibertype
lm_fa_sex <- contrast(means_fa, method = "pairwise", by = "fibertype")
print(lm_fa_sex)

#run linear model of fibertype for each sex
lm_fa_fibertype <- contrast(means_fa, method = "pairwise", by = "sex")
print(lm_fa_fibertype)


#Compute emmeans from linear mixed model
emm_glyco <- emmeans(lm_glyco_fibertype, ~ sex | sex) %>%
    as.data.frame() %>%
    mutate(bp = "Glycolysis") %>%
    mutate(
        fibertype = case_when(
            estimate < 0 ~ "type1",
            estimate > 0 ~ "type2"))


emm_fa <- emmeans(lm_fa_fibertype, ~ sex | sex) %>%
    as.data.frame() %>%
    mutate(bp = "Fatty acid metabolism") %>%
    mutate(
        fibertype = case_when(
            estimate < 0 ~ "type1",
            estimate > 0 ~ "type2"))

emm_bp <- rbind(emm_glyco, emm_fa)


#define text for figure
p_bp <- tibble(
    sex = c("female", "female",
            "male", "male"),
    x = c(1, 2,
          1, 2),
    y = c(0.1, 0.9,
          0.1, 1.1),
    label = c("p<0.001", "p<0.001",
              "p<0.001", "p<0.001")
)

##FIGURE OF GLYCOLYSIS AND FATTY ACID METABOLISM##
bp_fig <- ggplot(emm_bp, aes(x = bp, y = estimate, fill = fibertype)) +
    geom_col(position = position_dodge(width = 0.9),
             width = 0.9, color = NA, alpha = 0.7) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, position = position_dodge(width = 0.9),
                  linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    scale_fill_manual(values=c("#440154FF", "#67CC5CFF"))+
    scale_x_discrete(labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    facet_wrap(~sex, labeller = labeller(sex = c("female" = "Females", "male" = "Males"))) +
    #scale_y_continuous(limits = c(-0.5, 2)) +
    #geom_segment(data = brackets_bp, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_bp, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab(NULL) +
    ylab("Log2fold difference (type II - type I)") +
    ggtitle("Glycolysis and fatty acid metabolism")

ggsave(plot = bp_fig, here::here('figures/figure_2/glyco_fa_figure.pdf'), height = 60, width = 90, units = "mm")

