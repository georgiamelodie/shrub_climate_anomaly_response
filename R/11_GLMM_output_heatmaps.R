##############################################################################
# Alnus and Salix shrub anatomical response to climate
# ISR xylem anomalies - heatmap of anomaly GLMM outputs

library(here)

source(here("R/utils_packages.R"))
load_project_packages()

source(here("R/00_paths_setup.R")) 

data_dir    <- here("data")
figures_dir <- here("figures")
output_dir  <- here("output")

if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

###
# Read *_coeffs.csv files from output directory only

files <- list.files(
  path = output_dir,
  pattern = "_coeffs\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No *_coeffs.csv files found in output directory.")
}

effects_all <- purrr::map_dfr(files, function(f) {
  
  dat <- readr::read_csv(f, show_col_types = FALSE)
  fname <- basename(f)
  
  taxa <- case_when(
    str_detect(fname, regex("^Salix", ignore_case = TRUE)) ~ "Salix",
    str_detect(fname, regex("^Alnus", ignore_case = TRUE)) ~ "Alnus",
    TRUE ~ NA_character_
  )
  
  model <- case_when(
    str_detect(fname, regex("LWBR", ignore_case = TRUE)) ~ "LWBR",
    str_detect(fname, regex("LWFR", ignore_case = TRUE)) ~ "LWFR",
    str_detect(fname, regex("EWFR", ignore_case = TRUE)) ~ "EWFR",
    TRUE ~ NA_character_
  )
  
  year_type <- dplyr::case_when(
    str_detect(fname, regex("_prev(_|\\b)|prevyear|legacy|preceding", ignore_case = TRUE)) ~ "prev",
    str_detect(fname, regex("_current(_|\\b)|current", ignore_case = TRUE)) ~ "current",
    TRUE ~ "current"  # fallback if old files lack tokens
  )
  
  dat %>%
    mutate(
      file = fname,
      taxa = taxa,
      model = model,
      year_type = year_type
    )
})





# Select which model set (current vs prevyear) to use for each taxon × anomaly. determined from previous AICc comparison
best_model_key <- tibble::tribble(
  ~taxa,   ~model,  ~year_type_keep,
  "Salix", "EWFR",  "prev",
  "Salix", "LWFR",  "prev",
  "Salix", "LWBR",  "prev",
  "Alnus", "EWFR",  "prev",
  "Alnus", "LWFR",  "current",
  "Alnus", "LWBR",  "current"
)

# Keep only rows from the chosen model set 
effects_all <- effects_all %>%
  inner_join(best_model_key, by = c("taxa", "model")) %>%
  filter(year_type == year_type_keep) %>%
  dplyr::select(-year_type_keep)


#keep only significant predictors and remove intercept
effects_all <- effects_all %>%
  filter(
    !str_detect(term, regex("Intercept", ignore_case = TRUE)),
    p.value < 0.05
  )


#term names into month + predictor type and prev/current 
month_order <- c("May", "Jun", "Jul", "Aug", "Sep")


effects_all <- effects_all %>%
  mutate(
    term_prev = str_detect(term, "_prev"),
    
    month_num = str_extract(term, "month[5-9]") %>%
      str_remove("month") %>%
      as.integer(),
    
    Month = case_when(
      month_num == 5 ~ "May",
      month_num == 6 ~ "Jun",
      month_num == 7 ~ "Jul",
      month_num == 8 ~ "Aug",
      month_num == 9 ~ "Sep",
      TRUE ~ NA_character_
    ),
    
    Predictor_type = case_when(
      str_detect(term, regex("precip", ignore_case = TRUE)) ~ "Precip",
      str_detect(term, regex("CDD", ignore_case = TRUE)) ~ "CDD",
      str_detect(term, regex("GDD", ignore_case = TRUE)) ~ "GDD",
      str_detect(term, regex("^Age_s$", ignore_case = TRUE)) ~ "ring age",
      str_detect(term, regex("lat", ignore_case = TRUE)) ~ "latitude",
      TRUE ~ NA_character_
    ),
    
    #define prev/current from the term itself, not the file type
    Year_type = if_else(term_prev, "prev", "current"),
    
    Full_label = case_when(
      !is.na(Predictor_type) & !is.na(Month) & Year_type == "prev"    ~ paste("prev", Month, Predictor_type),
      !is.na(Predictor_type) & !is.na(Month) & Year_type == "current" ~ paste(Month, Predictor_type),
      !is.na(Predictor_type) &  is.na(Month)                          ~ Predictor_type,
      TRUE                                                            ~ term
    ),
    
    log_OR = estimate
  ) %>%
  # drop ring age and latitude before factoring Predictor_type
  filter(!Predictor_type %in% c("ring age", "latitude")) %>%
  mutate(
    Month = factor(Month, levels = month_order),
    Predictor_type = factor(Predictor_type, levels = c("Precip", "CDD", "GDD"))
  )






#effect ordering
label_levels <- effects_all %>%
  distinct(taxa, Year_type, Month, Predictor_type, Full_label) %>%
  mutate(
    is_month = !is.na(Month) & !is.na(Predictor_type),
    group_rank = case_when(
      is_month & Year_type == "prev"    ~ 1L,
      is_month & Year_type == "current" ~ 2L,
      TRUE                              ~ 3L
    ),
    Full_label_sort = as.character(Full_label)
  ) %>%
  arrange(group_rank, Month, Predictor_type, Full_label_sort) %>%
  pull(Full_label) %>%
  unique()

effects_all <- effects_all %>%
  mutate(
    Full_label = factor(Full_label, levels = label_levels, ordered = TRUE),
    model = factor(model, levels = c("EWFR", "LWBR", "LWFR")),
    taxa  = factor(taxa, levels = c("Salix", "Alnus"))
  ) %>%
  mutate(
    Full_label = droplevels(Full_label)
  )

effects_all <- effects_all %>%
  mutate(model = factor(model, levels = c("EWFR", "LWFR", "LWBR")))



#Plot heatmap (facet by taxa)
ggplot(effects_all, aes(x = Full_label, y = model, fill = log_OR)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(~ taxa, ncol = 2, scales = "free_x") +
  scale_fill_gradient2(
    low = "#5E3C99",
    mid = "white",
    high = "#1B7837",
    midpoint = 0,
    name = "log(Odds Ratio)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(colour = "black", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position = "bottom",
    legend.key.width = grid::unit(1.5, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )


ggsave(
  file.path(figures_dir, "GLMM_anomaly_heatmap.png"),
  width = 10,
  height = 5,
  dpi = 300
)








#split panels so well spaced
effects_salix <- effects_all %>% filter(taxa == "Salix")
effects_alnus <- effects_all %>% filter(taxa == "Alnus")

effects_salix <- effects_salix %>%
  mutate(model = factor(model, levels = c("EWFR", "LWFR", "LWBR")))

effects_alnus <- effects_alnus %>%
  mutate(model = factor(model, levels = c("EWFR", "LWFR", "LWBR")))


plot_heatmap <- function(dat, title) {
  ggplot(dat, aes(x = Full_label, y = model, fill = log_OR)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "#5E3C99",
      mid = "white",
      high = "#1B7837",
      midpoint = 0,
      name = "log(Odds Ratio)"
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1,colour = "black"),
      axis.text.y = element_text(colour = "black"),
      legend.position = "bottom"
    )
}



p_salix <- plot_heatmap(effects_salix, "Salix")
p_alnus <- plot_heatmap(effects_alnus, "Alnus")

p_salix 
ggsave(
  file.path(figures_dir, "Salix_GLMM_anomaly_heatmap.png"),
  plot = p_salix,
  width = 10, height = 5, dpi = 300
)

p_alnus
ggsave(
  file.path(figures_dir, "Alnus_GLMM_anomaly_heatmap.png"),
  plot = p_alnus,
  width = 10, height = 5, dpi = 300
)


