##############################################################################
# Alnus and Salix shrub anatomical response to climate
#########Plot of BR and FR occurrence in Alnus and Salix
library(here)

source(here("R/utils_packages.R"))
load_project_packages()

source(here("R/00_paths_setup.R")) 

data_dir    <- here("data")
figures_dir <- here("figures")
output_dir  <- here("output")

if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
if (!dir.exists(output_dir))  dir.create(output_dir, recursive = TRUE)



#####read in BR and FR datasets

dfsalix <- read.csv(
  file.path(data_dir, "salixBRFRdata.csv"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "Na", "N/A")
)

dfalnus <- read.csv(
  file.path(data_dir, "alnusBRFRdata.csv"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "Na", "N/A")
)

# Force: year numeric, everything else character
dfsalix <- dfsalix %>% mutate(year = as.integer(year), across(-year, as.character))
dfalnus <- dfalnus %>% mutate(year = as.integer(year), across(-year, as.character))





#function to process alnus and salix BRs and FRs
process_species <- function(df, species_name) {
  
  long <- df %>%
    pivot_longer(
      cols = -year,
      names_to = "sample",
      values_to = "anomaly_raw"
    ) %>%
    mutate(anomaly_raw = str_trim(anomaly_raw)) %>%
    separate_rows(anomaly_raw, sep = ",") %>%
    mutate(
      anomaly_raw = str_trim(anomaly_raw),
      
      # 1. classify blue/frost (full or partial)
      anomaly = case_when(
        str_detect(anomaly_raw, regex("BLW", ignore_case = TRUE)) ~ "blue",
        str_detect(anomaly_raw, regex("FLW", ignore_case = TRUE)) ~ "frost",
        TRUE ~ NA_character_
      ),
      
      # 2. classify full vs partial
      anomaly_strength = case_when(
        str_detect(anomaly_raw, regex("^pBLW$|^pFLW$", ignore_case = TRUE)) ~ "partial",
        str_detect(anomaly_raw, regex("BLW|FLW", ignore_case = TRUE)) ~ "full",
        TRUE ~ NA_character_
      ),
      
      species = species_name
    )
  
  spans <- long %>%
    filter(!is.na(anomaly_raw) & anomaly_raw != "") %>%
    group_by(sample, species) %>%
    summarise(
      first_year = min(year, na.rm = TRUE),
      last_year  = max(year, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(last_year)) %>%
    mutate(sample_order = factor(sample, levels = sample))
  
  list(long = long, spans = spans)
}


#process
salix_data <- process_species(dfsalix, "Salix")
alnus_data <- process_species(dfalnus, "Alnus")

long_all  <- bind_rows(salix_data$long,  alnus_data$long)
spans_all <- bind_rows(salix_data$spans, alnus_data$spans)



#create combined factor (blue_full, blue_partial, frost_full, frost_partial)
long_all <- long_all %>%
  mutate(anom_class = paste(anomaly, anomaly_strength, sep = "_"))

#colours
lifespan_colours <- c(
  "Salix" = "#2196F3",   # blue
  "Alnus" = "#D55E00"    # orange
)

# shapes
shape_vals <- c(
  blue_full      = 22,  # filled square
  blue_partial   = 0,   # empty square
  frost_full     = 23,  # filled diamond
  frost_partial  = 5    # empty diamond
)

# sizes
size_vals <- c(
  blue_full = 1.5,
  blue_partial = 1,
  frost_full = 2,
  frost_partial = 1.5
)

# fills (NA = hollow symbols)
fill_vals <- c(
  blue_full      = "#264bff",
  blue_partial   = NA,
  frost_full     = 'lightblue',
  frost_partial  = NA
)

# x-axis breaks at multiples of 5
x_breaks <- seq(
  from = floor(min(long_all$year, na.rm = TRUE) / 5) * 5,
  to   = ceiling(max(long_all$year, na.rm = TRUE) / 5) * 5,
  by   = 5
)


#plot

brfr_plot <- ggplot() + 
  # Lifespan bars
  geom_segment(
    data = spans_all,
    aes(
      x = first_year, xend = last_year,
      y = sample_order, yend = sample_order,
      colour = species
    ),
    linewidth = 1.1
  ) +
  scale_colour_manual(values = lifespan_colours, guide = "none") +
  
  #full anomalies
  geom_point(
    data = long_all %>% filter(anomaly_strength == "full" & !is.na(anom_class)),
    aes(
      x = year, y = sample,
      shape = anom_class, fill = anom_class, size = anom_class
    ),
    colour = "black",
    stroke = 0.4
  ) +
  
  #partial anomalies
  geom_point(
    data = long_all %>% filter(anomaly_strength == "partial" & !is.na(anom_class)),
    aes(
      x = year, y = sample,
      shape = anom_class, fill = anom_class, size = anom_class
    ),
    colour = "black",
    stroke = 0.4
  ) +
  
  scale_shape_manual(
    values = shape_vals,
    name = "Anomaly type",
    labels = c(
      blue_full = "LWBR",
      blue_partial = "pLWBR",
      frost_full = "LWFR",
      frost_partial = "pLWFR"
    )
  ) +
  scale_fill_manual(
    values = fill_vals,
    name = "Anomaly type",
    labels = c(
      blue_full = "LWBR",
      blue_partial = "pLWBR",
      frost_full = "LWFR",
      frost_partial = "pLWFR"
    )
  ) +
  scale_size_manual(values = size_vals, guide = "none") +
  
  facet_grid(
    rows = vars(species),
    scales = "free_y",
    space = "free_y"
  ) +
  
  scale_x_continuous(breaks = x_breaks) +
  
  labs(x = "Year", y = "Sample ID") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 13, face = "bold"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    legend.position = "right",
    legend.justification = c(5, 0.5),
    legend.margin = margin(r = 0, l = -5)
  )

# Print
brfr_plot

#save plot
ggsave(file.path(figures_dir, "LWBRLWFR_plot.jpg"), brfr_plot, width = 9, height = 10, dpi = 300)




#####Plot of EWBR, EWFR (supp. materials plot)
# function to process alnus and salix BRs and FRs (earlywood: BEW/FEW)
process_species <- function(df, species_name) {
  
  long <- df %>%
    pivot_longer(
      cols = -year,
      names_to = "sample",
      values_to = "anomaly_raw"
    ) %>%
    mutate(anomaly_raw = str_trim(anomaly_raw)) %>%
    separate_rows(anomaly_raw, sep = ",") %>%
    mutate(
      anomaly_raw = str_trim(anomaly_raw),
      
      # 1) classify blue vs frost (earlywood codes)
      anomaly = case_when(
        str_detect(anomaly_raw, regex("BEW", ignore_case = TRUE)) ~ "blue",
        str_detect(anomaly_raw, regex("FEW", ignore_case = TRUE)) ~ "frost",
        TRUE ~ NA_character_
      ),
      
      # 2) classify full vs partial (earlywood codes)
      anomaly_strength = case_when(
        str_detect(anomaly_raw, regex("^pBEW$|^pFEW$", ignore_case = TRUE)) ~ "partial",
        str_detect(anomaly_raw, regex("BEW|FEW", ignore_case = TRUE)) ~ "full",
        TRUE ~ NA_character_
      ),
      
      species = species_name
    )
  
  spans <- long %>%
    filter(!is.na(anomaly_raw) & anomaly_raw != "") %>%
    group_by(sample, species) %>%
    summarise(
      first_year = min(year, na.rm = TRUE),
      last_year  = max(year, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(last_year)) %>%
    mutate(sample_order = factor(sample, levels = sample))
  
  list(long = long, spans = spans)
}

# process
salix_data <- process_species(dfsalix, "Salix")
alnus_data <- process_species(dfalnus, "Alnus")

long_all  <- bind_rows(salix_data$long,  alnus_data$long)
spans_all <- bind_rows(salix_data$spans, alnus_data$spans)

# combined factor (blue_full, blue_partial, frost_full, frost_partial)
long_all <- long_all %>%
  mutate(anom_class = paste(anomaly, anomaly_strength, sep = "_"))

# colours
lifespan_colours <- c(
  "Salix" = "#2196F3",
  "Alnus" = "#D55E00"
)

# shapes
shape_vals <- c(
  blue_full      = 22,  # filled square
  blue_partial   = 0,   # empty square
  frost_full     = 23,  # filled diamond
  frost_partial  = 5    # empty diamond
)

# sizes
size_vals <- c(
  blue_full = 1.5,
  blue_partial = 1,
  frost_full = 2,
  frost_partial = 1.5
)

# fills
fill_vals <- c(
  blue_full      = "#264bff",
  blue_partial   = NA,
  frost_full     = "lightblue",
  frost_partial  = NA
)

# x-axis breaks at multiples of 5
x_breaks <- seq(
  from = floor(min(long_all$year, na.rm = TRUE) / 5) * 5,
  to   = ceiling(max(long_all$year, na.rm = TRUE) / 5) * 5,
  by   = 5
)

# plot (earlywood)
brfr_plot <- ggplot() +
  geom_segment(
    data = spans_all,
    aes(
      x = first_year, xend = last_year,
      y = sample_order, yend = sample_order,
      colour = species
    ),
    linewidth = 1.1
  ) +
  scale_colour_manual(values = lifespan_colours, guide = "none") +
  
  geom_point(
    data = long_all %>% filter(anomaly_strength == "full" & !is.na(anom_class)),
    aes(x = year, y = sample, shape = anom_class, fill = anom_class, size = anom_class),
    colour = "black",
    stroke = 0.4
  ) +
  geom_point(
    data = long_all %>% filter(anomaly_strength == "partial" & !is.na(anom_class)),
    aes(x = year, y = sample, shape = anom_class, fill = anom_class, size = anom_class),
    colour = "black",
    stroke = 0.4
  ) +
  
  scale_shape_manual(
    values = shape_vals,
    name = "Anomaly type",
    labels = c(
      blue_full = "EWBR",
      blue_partial = "pEWBR",
      frost_full = "EWFR",
      frost_partial = "pEWFR"
    )
  ) +
  scale_fill_manual(
    values = fill_vals,
    name = "Anomaly type",
    labels = c(
      blue_full = "EWBR",
      blue_partial = "pEWBR",
      frost_full = "EWFR",
      frost_partial = "pEWFR"
    )
  ) +
  scale_size_manual(values = size_vals, guide = "none") +
  
  facet_grid(
    rows = vars(species),
    scales = "free_y",
    space = "free_y"
  ) +
  
  scale_x_continuous(breaks = x_breaks) +
  labs(x = "Year", y = "Sample ID") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 13, face = "bold"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    legend.position = "right",
    legend.justification = c(5, 0.5),
    legend.margin = margin(r = 0, l = -5)
  )

brfr_plot

##save plot
ggsave(file.path(figures_dir, "EWBREFWFR_plot.jpg"),
       brfr_plot, width = 9, height = 10, dpi = 300)





# Frequency subplots for LWBR and LWFR
# Filter display by minimum sample depth 

#Function: build frequency table with frequency filter based on sample depth
make_anom_freq <- function(df, start_year, samp_dep = 5) {
  
  # Long including 0s (for denominator + measured span)
  long_all <- df %>%
    pivot_longer(cols = -year, names_to = "sample", values_to = "anomaly_raw") %>%
    mutate(anomaly_raw = str_trim(anomaly_raw))
  
  # Denominator: active samples (measured) per year
  # blank = missing; "0" counts as active
  active_samples <- long_all %>%
    mutate(anomaly_raw = na_if(anomaly_raw, "")) %>%
    filter(!is.na(anomaly_raw), year >= start_year) %>%
    group_by(sample) %>%
    summarise(first_year = min(year), last_year = max(year), .groups = "drop") %>%
    rowwise() %>%
    mutate(year = list(seq(first_year, last_year))) %>%
    unnest(year) %>%
    filter(year >= start_year) %>%
    group_by(year) %>%
    summarise(n_active = n_distinct(sample), .groups = "drop")
  
  # Measurement span (for completing years + x-axis limits)
  min_year_measured <- min(active_samples$year, na.rm = TRUE)
  max_year_measured <- max(active_samples$year, na.rm = TRUE)
  all_years <- min_year_measured:max_year_measured
  
  # Anomaly-only long (drop 0/blank, split codes)
  long <- long_all %>%
    mutate(
      anomaly_raw = na_if(anomaly_raw, "0"),
      anomaly_raw = na_if(anomaly_raw, "")
    ) %>%
    filter(!is.na(anomaly_raw), year >= start_year) %>%
    separate_rows(anomaly_raw, sep = ",") %>%
    mutate(
      anomaly_raw = str_trim(anomaly_raw),
      anomaly = case_when(
        str_detect(anomaly_raw, regex("^p?BLW$", ignore_case = TRUE)) ~ "LWBR",
        str_detect(anomaly_raw, regex("^p?FLW$", ignore_case = TRUE)) ~ "LWFR",
        TRUE ~ NA_character_
      ),
      type = case_when(
        str_detect(anomaly_raw, regex("^pBLW$|^pFLW$", ignore_case = TRUE)) ~ "partial",
        str_detect(anomaly_raw, regex("BLW$|FLW$", ignore_case = TRUE)) ~ "full",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(anomaly)) %>%
    distinct(year, sample, anomaly, type)
  
  #distinct anomalous samples per anomaly + type + year
  anom_counts <- long %>%
    group_by(anomaly, type, year) %>%
    summarise(n_anom = n_distinct(sample), .groups = "drop")
  
  # Frequency table (complete across full measured span)
  anom_freq <- anom_counts %>%
    left_join(active_samples, by = "year") %>%
    mutate(freq = n_anom / n_active) %>%
    group_by(anomaly, type) %>%
    complete(year = all_years, fill = list(n_anom = 0, n_active = NA, freq = 0)) %>%
    ungroup()
  
  # Filter display: only show when sample depth >= samp_dep
  anom_freq <- anom_freq %>%
    mutate(freq_filt = if_else(n_active >= samp_dep, freq, NA_real_))
  
  list(anom_freq = anom_freq, active_samples = active_samples)
}



# Build frequency tables (Alnus + Salix)

samp_dep <- 5

# Alnus: start year for analysis
al_start <- 1915
al_out <- make_anom_freq(dfalnus, start_year = al_start, samp_dep = samp_dep)
anom_freq_al <- al_out$anom_freq

# Salix: start year for analysis
sa_start <- 1915
sa_out <- make_anom_freq(dfsalix, start_year = sa_start, samp_dep = samp_dep)
anom_freq_sa <- sa_out$anom_freq


# X-axis settings
# Alnus x-axis
x_min_al <- min(al_out$active_samples$year, na.rm = TRUE)
x_max_al <- max(al_out$active_samples$year, na.rm = TRUE)
x_breaks_al <- seq(floor(x_min_al / 5) * 5, floor(x_max_al / 5) * 5, by = 5)

# Salix x-axis (force to 1915 for neatness)
x_min_sa <- 1915
x_max_sa <- max(sa_out$active_samples$year, na.rm = TRUE)
x_breaks_sa <- seq(floor(x_min_sa / 5) * 5, floor(x_max_sa / 5) * 5, by = 5)



#plot - Alnus 

# Fig 1 LWBR and pLWBR (Alnus)
al_lwbr <- ggplot(
  filter(anom_freq_al, anomaly == "LWBR"),
  aes(x = year, y = freq_filt, fill = type, colour = type)
) +
  geom_col(width = 0.9) +
  scale_fill_manual(
    values = c("full" = "blue", "partial" = "white"),
    labels = c("full" = "LWBR", "partial" = "pLWBR"),
    guide = "legend"
  ) +
  scale_colour_manual(
    values = c("full" = "blue", "partial" = "blue"),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.25),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_continuous(
    limits = c(x_min_al, x_max_al),
    breaks = x_breaks_al
  ) +
  labs(
    x = NULL,
    y = "Proportion",
    title = "LWBR"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title         = element_text(size = 12, face = "bold"),
    axis.title.y       = element_text(size = 12),
    panel.border       = element_blank(),
    axis.line          = element_line(colour = "black", size = 0.5)
  )

# Fig 2 LWFR and pLWFR (Alnus)
al_lwfr <- ggplot(
  filter(anom_freq_al, anomaly == "LWFR"),
  aes(x = year, y = freq_filt, fill = type, colour = type)
) +
  geom_col(width = 0.9) +
  scale_fill_manual(
    values = c("full" = "lightblue", "partial" = "white"),
    labels = c("full" = "LWFR", "partial" = "pLWFR"),
    guide = "legend"
  ) +
  scale_colour_manual(
    values = c("full" = "lightblue", "partial" = "lightblue"),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.25),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_continuous(
    limits = c(x_min_al, x_max_al),
    breaks = x_breaks_al
  ) +
  labs(
    x = NULL,
    y = "Proportion",
    title = "LWFR"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title         = element_text(size = 12, face = "bold"),
    axis.title.y       = element_text(size = 12),
    panel.border       = element_blank(),
    axis.line          = element_line(colour = "black", size = 0.5)
  )



#plot — Salix 

# Fig 1 LWBR and pLWBR (Salix)
sa_lwbr <- ggplot(
  filter(anom_freq_sa, anomaly == "LWBR"),
  aes(x = year, y = freq_filt, fill = type, colour = type)
) +
  geom_col(width = 0.9) +
  scale_fill_manual(
    values = c("full" = "blue", "partial" = "white"),
    labels = c("full" = "LWBR", "partial" = "pLWBR"),
    guide = "legend"
  ) +
  scale_colour_manual(
    values = c("full" = "blue", "partial" = "blue"),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.25),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_continuous(
    limits = c(x_min_sa, x_max_sa),
    breaks = x_breaks_sa
  ) +
  labs(
    x = NULL,
    y = "Proportion",
    title = "LWBR"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title         = element_text(size = 12, face = "bold"),
    axis.title.y       = element_text(size = 12),
    panel.border       = element_blank(),
    axis.line          = element_line(colour = "black", size = 0.5)
  )

# Fig 2 LWFR and pLWFR (Salix)
sa_lwfr <- ggplot(
  filter(anom_freq_sa, anomaly == "LWFR"),
  aes(x = year, y = freq_filt, fill = type, colour = type)
) +
  geom_col(width = 0.9) +
  scale_fill_manual(
    values = c("full" = "lightblue", "partial" = "white"),
    labels = c("full" = "LWFR", "partial" = "pLWFR"),
    guide = "legend"
  ) +
  scale_colour_manual(
    values = c("full" = "lightblue", "partial" = "lightblue"),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.25),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_continuous(
    limits = c(x_min_sa, x_max_sa),
    breaks = x_breaks_sa
  ) +
  labs(
    x = NULL,
    y = "Proportion",
    title = "LWFR"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title         = element_text(size = 12, face = "bold"),
    axis.title.y       = element_text(size = 12),
    panel.border       = element_blank(),
    axis.line          = element_line(colour = "black", size = 0.5)
  )


#print plots
al_lwbr
al_lwfr
sa_lwbr
sa_lwfr


# Save plots
ggsave(file.path(figures_dir, paste0("Alnus_LWBR_freq_minsampdepth", samp_dep, ".png")),
       al_lwbr, width = 8, height = 1.5, dpi = 300, bg = "white")

ggsave(file.path(figures_dir, paste0("Alnus_LWFR_freq_minsampdepth", samp_dep, ".png")),
       al_lwfr, width = 8, height = 1.5, dpi = 300, bg = "white")

ggsave(file.path(figures_dir, paste0("Salix_LWBR_freq_minsampdepth", samp_dep, ".png")),
       sa_lwbr, width = 8, height = 1.5, dpi = 300, bg = "white")

ggsave(file.path(figures_dir, paste0("Salix_LWFR_freq_minsampdepth", samp_dep, ".png")),
       sa_lwfr, width = 8, height = 1.5, dpi = 300, bg = "white")



