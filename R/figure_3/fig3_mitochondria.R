source(here::here("R/Library.R"))

#load data long form with keywords
df_long <- readRDS(here::here("data/data_long_keywords.rds"))

#create dataframe of mean expression of each protein pre
df_pre_mean <- df_long %>%
    filter(trial == "pre") %>%
    group_by(protein, group, fibertype, sex) %>%
    summarise(mean_expression = mean(expression, na.rm = TRUE)) %>%
    ungroup()

#Load log2fold difference between sexes
log2fd_sex <- vroom::vroom(here::here("data/results_sex_keywords.csv"))

#Load mitocarta
mitocarta <- read_excel(here::here('data-raw/mitocarta.xls'))%>%
    dplyr::select('symbol', 'pathways') %>%
    dplyr::rename(protein=symbol)

# Make sex a factor with male as reference
df_pre_mean$sex <- factor(df_pre_mean$sex, levels = c("male", "female"))

# MITOCHONDRIAL COMPLEXES -------------------------------------------------

#Add mitocarta to df
df_pre_mean <- df_pre_mean %>%
    merge(mitocarta, by="protein", all.x = T) %>%
    dplyr::rename(mito = pathways)

ci <- df_pre_mean %>%
    dplyr::filter(grepl('CI subunits', mito)) %>%
    dplyr::mutate(complex = "CI")

cii <- df_pre_mean %>%
    dplyr::filter(grepl('CII subunits', mito)) %>%
    dplyr::mutate(complex = "CII")

ciii <- df_pre_mean %>%
    dplyr::filter(grepl('CIII subunits', mito)) %>%
    dplyr::mutate(complex = "CIII")

civ <- df_pre_mean %>%
    dplyr::filter(grepl('CIV subunits', mito)) %>%
    dplyr::mutate(complex = "CIV")

cv <- df_pre_mean %>%
    dplyr::filter(grepl('CV subunits', mito)) %>%
    dplyr::mutate(complex = "CV")

mito_df <- rbind(ci, cii, ciii, civ, cv)

##CI##
#linear mixed model
lmm_ci <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = ci, REML = FALSE)
fixed_effects_ci <- anova(lmm_ci)

#calculate means for each fibertype for each sex
means_ci <- emmeans(lmm_ci, ~ sex | fibertype)
print(means_ci)

#run linear model of sex for each fibertype
lm_ci_sex <- contrast(means_ci, method = "pairwise", by = "fibertype")
print(lm_ci_sex)

#run linear model of fibertype for each sex
lm_ci_fibertype <- contrast(means_ci, method = "pairwise", by = "sex")
print(lm_ci_fibertype)

##CII###
#linear mixed model
lmm_cii <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = cii, REML = FALSE)
fixed_effects_cii <- anova(lmm_cii)

#calculate means for each fibertype for each sex
means_cii <- emmeans(lmm_cii, ~ sex | fibertype)
print(means_cii)

#run linear model of sex for each fibertype
lm_cii_sex <- contrast(means_cii, method = "pairwise", by = "fibertype")
print(lm_cii_sex)

#run linear model of fibertype for each sex
lm_cii_fibertype <- contrast(means_cii, method = "pairwise", by = "sex")
print(lm_cii_fibertype)

##CIII##
#linear mixed model
lmm_ciii <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = ciii, REML = FALSE)
fixed_effects_ciii <- anova(lmm_ciii)

#calculate means for each fibertype for each sex
means_ciii <- emmeans(lmm_ciii, ~ sex | fibertype)
print(means_ciii)

#run linear model of sex for each fibertype
lm_ciii_sex <- contrast(means_ciii, method = "pairwise", by = "fibertype")
print(lm_ciii_sex)

#run linear model of fibertype for each sex
lm_ciii_fibertype <- contrast(means_ciii, method = "pairwise", by = "sex")
print(lm_ciii_fibertype)

##CIV##
#linear mixed model
lmm_civ <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = civ, REML = FALSE)
fixed_effects_civ <- anova(lmm_civ)


#calculate means for each fibertype for each sex
means_civ <- emmeans(lmm_civ, ~ sex | fibertype)
print(means_civ)

#run linear model of sex for each fibertype
lm_civ_sex <- contrast(means_civ, method = "pairwise", by = "fibertype")
print(lm_civ_sex)

#run linear model of fibertype for each sex
lm_civ_fibertype <- contrast(means_civ, method = "pairwise", by = "sex")
print(lm_civ_fibertype)

##CV##
#linear mixed model
lmm_cv <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = cv, REML = FALSE)
fixed_effects_cv <- anova(lmm_cv)

#calculate means for each fibertype for each sex
means_cv <- emmeans(lmm_cv, ~ sex | fibertype)
print(means_cv)

#run linear model of sex for each fibertype
lm_cv_sex <- contrast(means_cv, method = "pairwise", by = "fibertype")
print(lm_cv_sex)

#run linear model of fibertype for each sex
lm_cv_fibertype <- contrast(means_cv, method = "pairwise", by = "sex")
print(lm_cv_fibertype)

##FIGURE OF SEX DIFFERENCE IN MITOCHONDRIAL COMPLEXES##

#Add mitocarta to log2fold difference data
log2fd_sex <- log2fd_sex %>%
    merge(mitocarta, by="protein", all.x = T) %>%
    dplyr::rename(mito = pathways)

ci_log2fd <- log2fd_sex %>%
    dplyr::filter(grepl('CI subunits', mito)) %>%
    dplyr::mutate(complex = "CI")

cii_log2fd <- log2fd_sex %>%
    dplyr::filter(grepl('CII subunits', mito)) %>%
    dplyr::mutate(complex = "CII")

ciii_log2fd <- log2fd_sex %>%
    dplyr::filter(grepl('CIII subunits', mito)) %>%
    dplyr::mutate(complex = "CIII")

civ_log2fd <- log2fd_sex %>%
    dplyr::filter(grepl('CIV subunits', mito)) %>%
    dplyr::mutate(complex = "CIV")

cv_log2fd <- log2fd_sex %>%
    dplyr::filter(grepl('CV subunits', mito)) %>%
    dplyr::mutate(complex = "CV")

mito_log2fd <- rbind(ci_log2fd, cii_log2fd, ciii_log2fd, civ_log2fd, cv_log2fd)

#define text for figure
p_mito_complex <- tibble(
    fibertype = c("mhc1", "mhc2",
                  "mhc1", "mhc2",
                  "mhc1", "mhc2",
                  "mhc1", "mhc2",
                  "mhc1", "mhc2"),
    x = c(1, 1,
          2, 2,
          3, 3,
          4, 4,
          5, 5),
    y = c(0.7, 0.95,
          0.9, 1.05,
          0.5, 0.8,
          0.85, 0.8,
          0.55, 0.65),
    label = c("p=0.060", "p<0.001",
              "p<0.001", "p<0.001",
              "p<0.001", "p<0.001",
              "p<0.001", "p<0.001",
              "p<0.001", "p<0.001")
)

#mito complexes box plot
mito_complexes <- mito_log2fd %>%
    ggplot(aes(x = complex, y = logFC, fill = fibertype)) +
    geom_violin(trim = TRUE, width = 1, linewidth = 0.5, alpha = 0.5) +
    geom_boxplot(width = 0.25, color = "black", fill = "white", alpha = 0.5, outlier.size = 1.5, outlier.stroke = 0) +
    geom_jitter(size=1.5, width=0, height=0, alpha = 0.5,  stroke = 0)+
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values=c("#440154FF", "#67CC5CFF"))+
    scale_x_discrete(labels = c("CI" = "Complex \n I", "CII" = "Complex \n II", "CIII" = "Complex \n III", "CIV" = "Complex \n IV", "CV" = "Complex \n V")) +
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
    facet_wrap(~fibertype, labeller = labeller(fibertype = c("mhc1" = "Type I", "mhc2" = "Type II"))) +
    scale_y_continuous(limits = c(-1, 1)) +
    geom_text(data = p_mito_complex, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("") +
    ylab("Log2fold difference (male - female)") +
    ggtitle("Mitochondrial complexes")

#mito complexes bar plot
#Compute emmeans from linear mixed model
emm_ci <- emmeans(lm_ci_sex, ~ fibertype | fibertype) %>%
    as.data.frame() %>%
    dplyr::mutate(complex = "CI")

emm_cii <- emmeans(lm_cii_sex, ~ fibertype | fibertype) %>%
    as.data.frame() %>%
    dplyr::mutate(complex = "CII")

emm_ciii <- emmeans(lm_ciii_sex, ~ fibertype | fibertype) %>%
    as.data.frame() %>%
    dplyr::mutate(complex = "CIII")

emm_civ <- emmeans(lm_civ_sex, ~ fibertype | fibertype) %>%
    as.data.frame() %>%
    dplyr::mutate(complex = "CIV")

emm_cv <- emmeans(lm_cv_sex, ~ fibertype | fibertype) %>%
    as.data.frame() %>%
    dplyr::mutate(complex = "CV")

emm_mito <- rbind(emm_ci, emm_cii, emm_ciii, emm_civ, emm_cv)

#define text for figure
p_mito_complex <- tibble(
    fibertype = c("mhc1", "mhc2",
                  "mhc1", "mhc2",
                  "mhc1", "mhc2",
                  "mhc1", "mhc2",
                  "mhc1", "mhc2"),
    x = c(1, 1,
          2, 2,
          3, 3,
          4, 4,
          5, 5),
    y = c(0.25, 0.30,
          0.65, 0.70,
          0.45, 0.50,
          0.50, 0.50,
          0.35, 0.40),
    label = c("p=0.060", "p<0.001",
              "p<0.001", "p<0.001",
              "p<0.001", "p<0.001",
              "p<0.001", "p<0.001",
              "p<0.001", "p<0.001")
)

##FIGURE##
mito_complexes <- emm_mito %>%
    ggplot(aes(x = complex, y = estimate, fill = fibertype)) +
    geom_col(position = position_dodge(width = 0.9),
             width = 0.9, color = NA, alpha = 0.7) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, position = position_dodge(width = 0.9),
                  linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    scale_fill_manual(values=c("#440154FF", "#67CC5CFF"))+
    scale_x_discrete(labels = c("CI" = "Complex \n I", "CII" = "Complex \n II", "CIII" = "Complex \n III", "CIV" = "Complex \n IV", "CV" = "Complex \n V")) +
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
    facet_wrap(~fibertype, labeller = labeller(fibertype = c("mhc1" = "Type I", "mhc2" = "Type II"))) +
    scale_y_continuous(limits = c(-0.2, 0.8)) +
    geom_text(data = p_mito_complex, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("") +
    ylab("Log2fold difference (male - female)") +
    ggtitle("Mitochondrial complexes")

ggsave(plot = mito_complexes, here::here('figures/figure_3/mito_complexes.pdf'), height = 45, width = 110, units = "mm")

##MITORIBOSOME SUBUNITS##

#filter for mitoribosome subunits
mito_ribo <- df_pre_mean %>%
    dplyr::filter(grepl("MRPS|MRPL", protein)) %>%
    dplyr::mutate(subunit = case_when(str_detect(protein, "^MRPS") ~ "MRPS",
                                      TRUE ~ "MRPL"))

mrpl <- mito_ribo %>%
    filter(subunit == "MRPL")

mrps <- mito_ribo %>%
    filter(subunit == "MRPS")

##LARGE SUBUNIT##
#linear mixed model
lmm_mrpl <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = mrpl, REML = FALSE)
fixed_effects_mrpl <- anova(lmm_mrpl)

#calculate means for each fibertype for each sex
means_mrpl <- emmeans(lmm_mrpl, ~ sex | fibertype)
print(means_mrpl)

#run linear model of sex for each fibertype
lm_mrpl_sex <- contrast(means_mrpl, method = "pairwise", by = "fibertype")
print(lm_mrpl_sex)

#run linear model of fibertype for each sex
lm_mrpl_fibertype <- contrast(means_mrpl, method = "pairwise", by = "sex")
print(lm_mrpl_fibertype)

##SMALL SUBUNIT##
#linear mixed model
lmm_mrps <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = mrps, REML = FALSE)
fixed_effects_mrps <- anova(lmm_mrps)

#calculate means for each fibertype for each sex
means_mrps <- emmeans(lmm_mrps, ~ sex | fibertype)
print(means_mrps)

#run linear model of sex for each fibertype
lm_mrps_sex <- contrast(means_mrps, method = "pairwise", by = "fibertype")
print(lm_mrps_sex)

#run linear model of fibertype for each sex
lm_mrps_fibertype <- contrast(means_mrps, method = "pairwise", by = "sex")
print(lm_mrps_fibertype)

#filter log2fd for mitoribosome subunits
mito_ribo_log2fd <- log2fd_sex %>%
    dplyr::filter(grepl("MRPS|MRPL", protein)) %>%
    dplyr::mutate(subunit = case_when(str_detect(protein, "^MRPS") ~ "MRPS",
                                      TRUE ~ "MRPL"))

#define text for figure
p_mito_ribo <- tibble(
    fibertype = c("mhc1", "mhc2",
                  "mhc1", "mhc2"),
    x = c(1, 1,
          2, 2),
    y = c(0.7, 0.9,
          0.7, 0.8),
    label = c("p<0.001", "p<0.001",
              "p=0.001", "p=0.003")
)

#mito ribosome box plot
mitoribo_fig <- mito_ribo_log2fd %>%
    ggplot(aes(x = subunit, y = logFC, fill = fibertype)) +
    geom_violin(trim = TRUE, width = 1, linewidth = 0.5, alpha = 0.5) +
    geom_boxplot(width = 0.25, color = "black", fill = "white", alpha = 0.5, outlier.size = 1.5, outlier.stroke = 0) +
    geom_jitter(size=1.5, width=0, height=0, alpha = 0.5,  stroke = 0)+
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values=c("#440154FF", "#67CC5CFF"))+
    scale_x_discrete(labels = c("MRPL" = "39S \n subunit", "MRPS" = "28S \n subunit")) +
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
    facet_wrap(~fibertype, labeller = labeller(fibertype = c("mhc1" = "Type I", "mhc2" = "Type II"))) +
    scale_y_continuous(limits = c(-1, 1)) +
    geom_text(data = p_mito_ribo, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("") +
    ylab("Log2fold difference (male - female)") +
    ggtitle("Mitoribosome subunits")

#mito ribosome bar plot
#Compute means from linear mixed model
emm_mrpl <- emmeans(lm_mrpl_sex, ~ fibertype | fibertype) %>%
    as.data.frame() %>%
    dplyr::mutate(subunit = "MRPL")

emm_mrps <- emmeans(lm_mrps_sex, ~ fibertype | fibertype) %>%
    as.data.frame() %>%
    dplyr::mutate(subunit = "MRPS")

emm_mitoribo <- rbind(emm_mrpl, emm_mrps)

#define text for figure
p_mito_ribo <- tibble(
    fibertype = c("mhc1", "mhc2",
                  "mhc1", "mhc2"),
    x = c(1, 1,
          2, 2),
    y = c(0.4, 0.45,
          0.4, 0.4),
    label = c("p<0.001", "p<0.001",
              "p=0.001", "p=0.003")
)

##FIGURE##
mitoribo_fig <- emm_mitoribo %>%
    ggplot(aes(x = subunit, y = estimate, fill = fibertype)) +
    geom_col(position = position_dodge(width = 0.9),
             width = 0.9, color = NA, alpha = 0.7) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, position = position_dodge(width = 0.9),
                  linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    scale_fill_manual(values=c("#440154FF", "#67CC5CFF"))+
    scale_x_discrete(labels = c("MRPL" = "39S \n subunit", "MRPS" = "28S \n subunit")) +
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
    facet_wrap(~fibertype, labeller = labeller(fibertype = c("mhc1" = "Type I", "mhc2" = "Type II"))) +
    scale_y_continuous(limits = c(-0.25, 0.5)) +
    geom_text(data = p_mito_ribo, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("") +
    ylab("Log2fold difference (male - female)") +
    ggtitle("Mitoribosome subunits")

ggsave(plot = mitoribo_fig, here::here('figures/figure_3/mito_ribosome.pdf'), height = 45, width = 60, units = "mm")

# CS ACTIVITY -------------------------------------------------------------

#load data
cs_df <- vroom::vroom(here::here("data-raw/additional_data.csv"))

#filter for pre
cs_df_pre <- cs_df %>%
    filter(trial == "pre")

#Linear mixed model
lmm_cs <- lmer(cs_activity ~ trial * sex + (1 | id), data = cs_df, REML = FALSE)
fixed_effects_cs <- anova(lmm_cs)

#calculate means for each trial for each sex
means_cs <- emmeans(lmm_cs, ~ sex | trial)
print(means_cs)

#run linear model of sex for each trial
lm_cs_sex <- contrast(means_cs, method = "pairwise", by = "trial", adjust = "none")
print(lm_cs_sex)

#define brackets for figure
brackets_cs <- tibble(
    x = c(1, 1, 2),
    xend = c(2, 1, 2),
    y = c(55, 55, 55),
    yend = c(55, 45, 50)
)

#define p-value for figure
p_cs <- tibble(
    x = 1.5,
    y = 58,
    label = "p=0.043"
)

#Box plot of CS activity pre
cs_fig <- cs_df_pre %>%
    ggplot(aes(x = sex, y = cs_activity, fill = sex)) +
    geom_violin(trim = TRUE, width=1, linewidth = 0.5, alpha=0.5)+
    geom_boxplot(width=0.25, color="black", fill="white", alpha=0.5, outlier.size = 2.5, outlier.stroke = 0)+
    geom_jitter(size=2.5, width=0, aes(color = sex), alpha = 0.5, stroke=0)+
    scale_fill_manual(values=c("#000000", "#FF7518"))+
    scale_color_manual(values = c("#000000", "#FF7518")) +
    scale_x_discrete(labels=c("Females", "Males"))+
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill=NA, linewidth = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        text = element_text(size = 6),
        axis.text.x= element_text(color="black", size = 6),
        axis.text.y= element_text(color="black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    )+
    scale_y_continuous(limits = c(10, 60)) +
    geom_segment(data = brackets_cs, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_cs, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("none") +
    ylab("Citrate synthase activity \n (µmol/g/min)") +
    ggtitle("Citrate synthase activity")


#Bar plot of CS activity pre
cs_fig <- cs_df_pre %>%
    ggplot(aes(x = sex, y = cs_activity, fill = sex)) +
    stat_summary(fun = mean, geom = "bar",
                 position = position_dodge(width = 0.9),
                 width = 0.9, color = NA, alpha = 0.8) +
    geom_jitter(size=2, width=0, aes(color = sex), alpha = 0.5, stroke=0)+
    scale_fill_manual(values=c("#000000", "#FF7518"))+
    scale_color_manual(values = c("#000000", "#FF7518")) +
    scale_x_discrete(labels=c("Females", "Males"))+
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill=NA, linewidth = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        text = element_text(size = 6),
        axis.text.x= element_text(color="black", size = 6),
        axis.text.y= element_text(color="black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    )+
    coord_cartesian(ylim = c(10, 60)) +
    geom_segment(data = brackets_cs, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_cs, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("none") +
    ylab("Citrate synthase activity \n (µmol/g/min)") +
    ggtitle("Citrate synthase activity")

ggsave(plot = cs_fig, here::here('figures/figure_3/cs_activity_pre.pdf'), height = 45, width = 45, units = "mm")
