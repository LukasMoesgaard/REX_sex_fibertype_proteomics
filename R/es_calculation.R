
#load data
df <- vroom::vroom(here::here("data/data_log2_normalized.csv"))

df_pre <- df %>%
    select(contains("_pre"))

df_post <- df %>%
    select(contains("_post"))

# Calculate effect size for each protein
effect_sizes <- df_post %>%
    mutate(
        effect_size = (rowMeans(df_post, na.rm = TRUE) - rowMeans(df_pre, na.rm = TRUE)) /
            apply(df_pre, 1, sd, na.rm = TRUE)
    ) %>%
    select(effect_size)

# Combine with protein identifiers
es_df <- bind_cols(df %>% select(1) %>% rename(protein = 1),
                   effect_sizes)

#load results
results_training <- vroom::vroom(here::here("results/trial_main_effect.csv"))

#Filter significant regulated proteins
sig_proteins <- results_training %>%
    filter(q < 0.05) %>%
    pull(protein)

# Filter the effect size results to only include these proteins
es_filtered <- es_df %>%
    filter(protein %in% sig_proteins)

