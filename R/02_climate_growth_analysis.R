##############################################################################
# Alnus and Salix shrub anatomical response to climate
# Growth–climate relationships 
library(here)

source(here("R/utils_packages.R"))
load_project_packages()

source(here("R/00_paths_setup.R")) 

# Project dirs
data_dir    <- here("data")
output_dir  <- here("output")
figures_dir <- here("figures")

if (!dir.exists(output_dir))  dir.create(output_dir)
if (!dir.exists(figures_dir)) dir.create(figures_dir)


##########CRU data
#read in CRU tmp data from dat file

##########CRU data
#read in CRU tmp data from dat file
cru_tmp_file <- file.path(data_dir, "CRU4_tmp.dat")
readLines(cru_tmp_file, n = 29)
CRUT <- read.table(cru_tmp_file, header = FALSE, skip = 26)

head(CRUT)

#Assign column headers
colnames(CRUT) <- c("Year", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
head(CRUT)


# Filter the data for the years 1969-2023
CRUT <- CRUT %>%
  filter(Year >= 1969 & Year <= 2023)



# Make first Year column into rownames for detrending
CRUTd <- data.frame(CRUT, row.names = 1)
CRUTd


#detrending climate data   #difference TRUE for temp, FALSE for precip
dclim <- detrend(CRUTd, method = "Spline", nyrs = 10 * length(CRUTd), difference = TRUE)
dclim
colnames(dclim) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
head(dclim)

# Convert rownames to a column
dclim <- dclim %>%
  rownames_to_column(var = "Year")

# Ensure all values are numeric
dclim <- dclim %>%
  mutate_if(is.character, as.numeric)

dclim

write.csv(dclim, file.path(output_dir, "dclim.csv"), row.names = FALSE)





#######precipitation

##########ERA precip data

#read in climate from dat file
era_prcp_file <- file.path(data_dir, "ERA5_prcp_daily.dat")
readLines(era_prcp_file, n = 25)
ERA5P <- read.table(era_prcp_file, header = FALSE, skip = 24)

head(ERA5P)

# Assign column names
colnames(ERA5P) <- c("date", "precip_mm")

# Convert the 'date' column from YYYYMMDD format to a standard Date format
ERA5P$date <- as.Date(as.character(ERA5P$date), format = "%Y%m%d")

# Check the structure of the data
head(ERA5P)



# Extract year, month, and day from the 'time' column, padding single-digit month and day with leading zero
ERA5P <- ERA5P %>%
  mutate(Year = year(date),
         month = sprintf("%02d", month(date)),
         day = sprintf("%02d", day(date))) %>%
  relocate(Year, month, day, .after = date)


# Filter data for the range 1948-2023
ERA5P_filtered <- ERA5P %>%
  filter(Year >= 1948 & Year <= 2023)


# Plot 
ggplot(ERA5P_filtered, aes(x = date, y = precip_mm)) +
  geom_line(color = "blue") +
  labs(title = "precip 1973-2023)",
       x = "Date",
       y = "precip (mm)") +
  theme_minimal()




# Aggregate to monthly totals and reshape wide with months as 1–12
precip <- ERA5P_filtered %>%
  group_by(Year, month) %>%
  summarise(precip_mm = sum(precip_mm), .groups = "drop") %>%
  pivot_wider(
    names_from  = month,
    values_from = precip_mm
  ) %>%
  arrange(Year)

# Check result
str(precip)
head(precip)

#Assign column headers
colnames(precip) <- c("Year", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
head(precip)


# convert all month columns to numeric 
precip <- precip %>%
  mutate(across(-Year, ~as.numeric(.)))

# Check structure
str(precip)

# Rename columns to numeric 1-12
colnames(precip) <- c("Year", 1:12)

str(precip)






ERA5P_monthly <- aggregate(precip_mm~Year+month, ERA5P_filtered, sum)

ERA5P_yearly <- aggregate(precip_mm ~ Year, ERA5P_filtered, sum)



# pivot_wider: convert "01" -> integer 1, then use month.name
ERA5P_wide <- ERA5P_monthly %>%
  pivot_wider(
    names_from = month,
    values_from = precip_mm,
    names_glue = "{month.abb[as.integer(month)]}_precip",
    values_fill = 0
  )

head(ERA5P_wide)
# Ensure all values are numeric
precip <- ERA5P_wide %>%
  mutate_if(is.character, as.numeric)

precip

# Make first Year column into rownames for detrending
ERA5P_wided <- data.frame(ERA5P_wide, row.names = 1)
ERA5P_wided


#detrending climate data   #difference TRUE for temp, FALSE for precip
dprecip <- detrend(ERA5P_wided, method = "Spline", nyrs = 10 * length(ERA5P_wided), difference = FALSE)
dprecip
colnames(dprecip) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
head(dprecip)

# Convert rownames to a column
dprecip <- dprecip %>%
  rownames_to_column(var = "Year")

# Ensure all values are numeric
dprecip <- dprecip %>%
  mutate_if(is.character, as.numeric)

dprecip


##########load or compute detrended alnus and salix chronologies

read_crn_any_csv <- function(file, prefer = "vsc") {
  df <- read.csv(file, check.names = FALSE)
  
  names_lower <- tolower(names(df))
  
  # case 1: two-col year/vsc export
  if (ncol(df) == 2 && all(c("year", "vsc") %in% names_lower)) {
    df <- setNames(df, names_lower)
    out <- data.frame(vsc = as.numeric(df$vsc))
    rownames(out) <- as.integer(df$year)
    return(out)
  }
  
  # case 2: first col is year, remaining columns are series
  rownames(df) <- df[[1]]
  df <- df[, -1, drop = FALSE]
  df[] <- lapply(df, as.numeric)
  
  # if a preferred column exists (e.g., vsc), return only that
  if (!is.null(prefer) && prefer %in% colnames(df)) {
    df <- df[, prefer, drop = FALSE]
  }
  
  df
}


saspl10SV <- read_crn_any_csv(file.path(output_dir, "salix_variance_stabilised_chron.csv"))
alspl10SV <- read_crn_any_csv(file.path(output_dir, "alnus_variance_stabilised_chron.csv"))


############Temperature growth analysis 

#correlation analysis  
####for alnus
ac <- dcc(chrono = alspl10SV, climate = dclim, selection = -1:9,   #selection indicates previous year months (e.g-6, and current years months (e.g. 7))
          method = "correlation", dynamic = "static", win_size = 20, win_offset = 1, start_last = FALSE,
          timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE) #this is the main function in treeclim
coef <- coef(ac) #model coefficients 
plot(ac) #plot the model coefficients

summary(ac)


####for salix
sc <- dcc(chrono = saspl10SV, climate = dclim, selection = -1:9,   #selection indicates previous year months (e.g-6, and current years months (e.g. 7))
          method = "correlation", dynamic = "static", win_size = 20, win_offset = 1, start_last = FALSE,
          timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE) #this is the main function in treeclim
coef <- coef(sc) #model coefficients 
plot(sc) #plot the model coefficients

summary(sc)


#############save output as table

# Extract summary table from sc
ac_summary <- summary(ac)

# Add row names as a column (Month info)
ac_summary <- ac_summary %>%
  tibble::rownames_to_column(var = "Variable")

# Save to CSV
write.csv(ac_summary, file.path(output_dir, "Alnus_static_tempdcc_summary.csv"), row.names = FALSE)


# Extract summary table from ac
sc_summary <- summary(sc)

# Add row names as a column (Month info)
sc_summary <- sc_summary %>%
  tibble::rownames_to_column(var = "Variable")

# Save to CSV
write.csv(sc_summary, file.path(output_dir, "Salix_static_dcc_summary.csv"), row.names = FALSE)


#################plot alnus (ac)  and salix (sc)  dcc output

library(ggplot2)
library(dplyr)
library(stringr)

#function to process static dcc object
process_dcc_summary <- function(dcc_obj, species_name) {
  coef_summary <- summary(dcc_obj)  # Extract summary table
  varnames <- rownames(dcc_obj$coef)  # Original row names
  
  coef_summary <- coef_summary %>%
    mutate(
      VarString = varnames,
      TimePeriod = case_when(
        str_detect(VarString, "prev") ~ "prev.",
        str_detect(VarString, "curr") ~ "cur.",
        TRUE ~ "Other"
      ),
      MonthClean = str_to_title(month),
      MonthLabel = paste(TimePeriod, MonthClean),
      Significant = significant,
      SigLabel = if_else(Significant, "TRUE", "FALSE"),
      Species = species_name
    )
  
  month_levels <- c(
    paste("prev.", c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")),
    paste("cur.", c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  )
  
  coef_summary$MonthLabel <- factor(coef_summary$MonthLabel, levels = month_levels)
  coef_summary$SigLabel <- factor(coef_summary$SigLabel, levels = c("TRUE", "FALSE"))
  
  return(coef_summary)
}


# Plot function 
plot_combined_species <- function(combined_df,
                                  save_path = NULL,
                                  width = 8,
                                  height = 6,
                                  dpi = 300,
                                  unify_y = FALSE,
                                  y_step = 0.4,
                                  symmetric = FALSE  # set TRUE to force +/- same limit around 0
) {
  # compute shared y-axis limits across both species
  if (unify_y) {
    # Use CI bounds if available; fall back to coef range
    y_min_raw <- if ("ci_lower" %in% names(combined_df)) min(combined_df$ci_lower, na.rm = TRUE) else min(combined_df$coef, na.rm = TRUE)
    y_max_raw <- if ("ci_upper" %in% names(combined_df)) max(combined_df$ci_upper, na.rm = TRUE) else max(combined_df$coef, na.rm = TRUE)
    
    if (symmetric) {
      y_ext <- max(abs(c(y_min_raw, y_max_raw)))
      y_min_raw <- -y_ext
      y_max_raw <-  y_ext
    }
    
    # Round to nearest step outward
    y_min <- floor(y_min_raw / y_step) * y_step
    y_max <- ceiling(y_max_raw / y_step) * y_step
    y_breaks <- seq(y_min, y_max, by = y_step)
  }
  
  p <- ggplot(combined_df, aes(x = MonthLabel, y = coef, color = coef)) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, linetype = SigLabel), width = 0.2, alpha = 0.8) +
    geom_point(aes(shape = SigLabel), size = 4, stroke = 1) +
    scale_color_gradient2(low = "#ca0020", mid = "white", high = "#0571b0", midpoint = 0) +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), name = "Significance") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), name = "Significance") +
    labs(
      x = "Month",
      y = "Correlation Coefficient",
      color = "Coefficient"
    ) +
    facet_wrap(~Species, ncol = 1, scales = if (unify_y) "fixed" else "free_y") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10),
      plot.title = element_blank()
    )
  
  if (unify_y) {
    p <- p + scale_y_continuous(limits = c(y_min, y_max), breaks = y_breaks)
  }
  
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = p, width = width, height = height, dpi = dpi)
    message("Plot saved to: ", save_path)
  }
  
  p
}


# Process
alnus_summary <- process_dcc_summary(ac, "Alnus")
salix_summary <- process_dcc_summary(sc, "Salix")

# Combine
combined_df <- bind_rows(alnus_summary, salix_summary) %>%
  mutate(Species = factor(Species, levels = c("Alnus", "Salix")))


# Plot 
plot_static <- plot_combined_species(combined_df, unify_y = TRUE, y_step = 0.4)
plot_static

# Save plot
ggsave(file.path(figures_dir, "Alnus_Salix_staticDCC.png"),
       plot_static, width = 10, height = 6, dpi = 300)




#################temperature growth moving windows plots

#correlation analyses for each species
###salix##
smc <- dcc(chrono = saspl10SV, climate = dclim, selection = 6:7,   #selection indicates previous year months (e.g-6, and current years months (e.g. 7))
           method = "correlation", dynamic = "moving", win_size = 20, win_offset = 1, start_last = FALSE,
           timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE) #this is the main function in treeclim
coef <- coef(smc) #model coefficients 
plot(smc) #plot the model coefficients


summary(smc)

###alnus###
amc <- dcc(chrono = alspl10SV, climate = dclim, selection = 6:7,   #selection indicates previous year months (e.g-6, and current years months (e.g. 7))
           method = "correlation", dynamic = "moving", win_size = 20, win_offset = 1, start_last = FALSE,
           timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE) #this is the main function in treeclim
coef <- coef(amc) #model coefficients 
plot(amc) #plot the model coefficients

summary(amc)


############## plots of temp-growth moving window correlations
plot_dcc_points_ci <- function(smc,
                               title = NULL) {
  # Extract data as before
  coef_df <- smc$coef$coef
  ci_lower_df <- smc$coef$ci_lower
  ci_upper_df <- smc$coef$ci_upper
  signif_df <- smc$coef$significant
  
  pivot_coef_data <- function(df, value_name) {
    df %>%
      rownames_to_column(var = "Month") %>%
      pivot_longer(-Month, names_to = "Window", values_to = value_name)
  }
  
  coef_long <- pivot_coef_data(coef_df, "Correlation")
  ci_lower_long <- pivot_coef_data(ci_lower_df, "CI_Lower")
  ci_upper_long <- pivot_coef_data(ci_upper_df, "CI_Upper")
  signif_long <- pivot_coef_data(signif_df, "Significant")
  
  combined <- coef_long %>%
    left_join(ci_lower_long, by = c("Month", "Window")) %>%
    left_join(ci_upper_long, by = c("Month", "Window")) %>%
    left_join(signif_long, by = c("Month", "Window"))
  
  combined <- combined %>%
    mutate(
      MonthType = case_when(
        str_detect(Month, "\\.curr\\.") ~ "Current",
        str_detect(Month, "\\.prev\\.") ~ "Previous",
        TRUE ~ "Other"
      ),
      Month = str_replace(Month, "X\\d+\\.(curr|prev)\\.", ""),
      Month = str_to_title(Month),
      Month = paste(MonthType, Month),
      # Factor for linetype and shape based on significance
      SigLineType = if_else(Significant, "solid", "dashed"),
      SigShape = if_else(Significant, 16, 1) # filled circle if sig, hollow if not
    )
  
  ggplot(combined, aes(x = Window, y = Correlation, color = Correlation)) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper, linetype = SigLineType), width = 0.3, alpha = 0.8) +
    geom_point(aes(shape = SigShape), size = 3, stroke = 1) +
    facet_wrap(~ Month, ncol = 1, scales = "free_y") +
    scale_color_gradient2(low = "#ca0020", mid = "white", high = "#0571b0", midpoint = 0) +
    scale_linetype_identity() +    # Use linetypes as is
    scale_shape_identity() +       # Use shapes as is
    labs(title = title, x = "Window", y = "Correlation") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5)
    )
}

plot_dcc_points_ci(smc)




 ########## temp-growth moving-window dcc object

process_moving_dcc <- function(mc, species_name = "Species", label_style = "short") {
  coef_df     <- mc$coef$coef
  ci_lower_df <- mc$coef$ci_lower
  ci_upper_df <- mc$coef$ci_upper
  signif_df   <- mc$coef$significant
  
  pivot_coef_data <- function(df, value_name) {
    df %>%
      rownames_to_column(var = "Month") %>%
      pivot_longer(-Month, names_to = "Window", values_to = value_name)
  }
  
  coef_long     <- pivot_coef_data(coef_df,     "Correlation")
  ci_lower_long <- pivot_coef_data(ci_lower_df, "CI_Lower")
  ci_upper_long <- pivot_coef_data(ci_upper_df, "CI_Upper")
  signif_long   <- pivot_coef_data(signif_df,   "Significant")
  
  df <- coef_long %>%
    left_join(ci_lower_long, by = c("Month","Window")) %>%
    left_join(ci_upper_long, by = c("Month","Window")) %>%
    left_join(signif_long,   by = c("Month","Window"))
  
  #### clean Window labels
  df <- df %>%
    mutate(
      Window = gsub("^X", "", Window),
      Window = gsub("\\.", "-", Window)
    )
  
  ##### month label parts
  df <- df %>%
    mutate(
      MonthType = case_when(
        str_detect(Month, "\\.curr\\.") ~ ifelse(label_style == "short", "cur.",  "Current"),
        str_detect(Month, "\\.prev\\.") ~ ifelse(label_style == "short", "prev.", "Previous"),
        TRUE ~ "Other"
      ),
      MonthName = str_replace(Month, "X\\d+\\.(curr|prev)\\.", ""),
      MonthName = str_to_title(MonthName),
      MonthLabel = paste(MonthType, MonthName),
      # YOUR significance labels preserved
      SigLabel = factor(if_else(Significant, "TRUE", "FALSE"), levels = c("TRUE","FALSE")),
      Species  = species_name
    )
  
  df
}



#########plots 
create_moving_species_plot <- function(data,
                                       title = NULL,
                                       every_n = 2,          # show every nth window label
                                       legend_title = "Significance",
                                       legend_text_size = 8,
                                       legend_title_size = 10,
                                       drop_unused = FALSE   # keep blank space for missing windows
) {
  ######  window order
  window_info <- data %>%
    distinct(Window) %>%
    mutate(
      Window = trimws(Window),
      start_yr = as.integer(sub("-.*", "", Window)),
      end_yr   = as.integer(sub(".*-", "", Window))
    ) %>%
    arrange(start_yr, end_yr)
  
  window_levels <- window_info$Window
  
  ######### ordered factor
  data <- data %>%
    mutate(Window = factor(Window, levels = window_levels))
  
  # breaks every_n
  breaks_vec <- window_levels[seq(1, length(window_levels), by = every_n)]
  
  #### calculate y-axis limits
  y_min <- floor(min(data$Correlation, na.rm = TRUE) * 4) / 4  # round down to nearest 0.25
  y_max <- ceiling(max(data$Correlation, na.rm = TRUE) * 4) / 4 # round up to nearest 0.25
  
  ####### build plot 
  p <- ggplot(data, aes(x = Window, y = Correlation, color = Correlation)) +
    geom_errorbar(
      aes(ymin = CI_Lower, ymax = CI_Upper, linetype = SigLabel),
      width = 0.3, alpha = 0.8
    ) +
    geom_point(
      aes(shape = SigLabel),
      size = 3, stroke = 1
    ) +
    facet_grid(MonthLabel ~ Species, scales = "fixed") +  # scales fixed for consistent y-axis
    scale_color_gradient2(
      low = "#ca0020", mid = "white", high = "#0571b0", midpoint = 0
    ) +
    scale_linetype_manual(
      values = c("TRUE" = "solid", "FALSE" = "dashed"),
      name = legend_title
    ) +
    scale_shape_manual(
      values = c("TRUE" = 16, "FALSE" = 1),
      name = legend_title
    ) +
    scale_x_discrete(breaks = breaks_vec, drop = drop_unused) +
    scale_y_continuous(
      breaks = seq(y_min, y_max, by = 0.25),
      limits = c(y_min, y_max+0.25)
    ) +
    labs(
      title = title,
      x = "Window",
      y = "Correlation",
      color = "Coefficient"
    ) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      strip.text   = element_text(face = "bold"),
      legend.text  = element_text(size = legend_text_size),
      legend.title = element_text(size = legend_title_size),
      plot.title   = element_text(hjust = 0.5)
    )
  
  p
}



##########use alnus and salix dcc outputs
#process each
alnus_moving <- process_moving_dcc(amc, "Alnus", label_style = "short")
salix_moving <- process_moving_dcc(smc, "Salix", label_style = "short")


#### order

moving_df <- bind_rows(alnus_moving, salix_moving) %>%
  mutate(
    Species   = factor(Species, levels = c("Alnus", "Salix")),
    MonthLabel = factor(MonthLabel, levels = c("cur. Jun", "cur. Jul"))
  )


#### plot, show every other label
moving_plot <- create_moving_species_plot(moving_df, title = NULL, every_n = 2)

print(moving_plot)







##Precipitation

###############precipitation correlation analysis 

####for alnus
acp <- dcc(chrono = alspl10SV, climate = dprecip, selection = -6:8,   #selection indicates previous year months (e.g-6, and current years months (e.g. 7))
           method = "correlation", dynamic = "static", win_size = 20, win_offset = 1, start_last = FALSE,
           timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE) #this is the main function in treeclim
coef <- coef(acp) #model coefficients 
plot(acp) #plot the model coefficients
summary(acp)


####for salix
scp <- dcc(chrono = saspl10SV, climate = dprecip, selection = -6:8,   #selection indicates previous year months (e.g-6, and current years months (e.g. 7))
           method = "correlation", dynamic = "static", win_size = 20, win_offset = 1, start_last = FALSE,
           timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE) #this is the main function in treeclim
coef <- coef(scp) #model coefficients 
plot(scp) #plot the model coefficients
summary(scp)



#############save output as table

# Extract summary table from acp
acp_summary <- summary(acp)

# Add row names as a column (Month)
acp_summary <- acp_summary %>%
  tibble::rownames_to_column(var = "Variable")

# Save to CSV
write.csv(acp_summary, file.path(output_dir, "Alnus_static_precipdcc_summary.csv"), row.names = FALSE)



# Extract summary table from scp
scp_summary <- summary(scp)

# Add row names as a column (Month)
scp_summary <- scp_summary %>%
  tibble::rownames_to_column(var = "Variable")

# Save to CSV
write.csv(sc_summary, file.path(output_dir, "Salix_static_tempdcc_summary.csv"), row.names = FALSE)


#################plot alnus (acp)  and salix (scp)  dcc output
####process function
process_dcc_summary <- function(dcc_obj, species_name) {
  coef_summary <- summary(dcc_obj)  # Extract summary table
  varnames <- rownames(dcc_obj$coef)  # Original row names
  
  coef_summary <- coef_summary %>%
    mutate(
      VarString = varnames,
      TimePeriod = case_when(
        str_detect(VarString, "prev") ~ "prev.",
        str_detect(VarString, "curr") ~ "cur.",
        TRUE ~ "Other"
      ),
      MonthClean = str_to_title(month),
      MonthLabel = paste(TimePeriod, MonthClean),
      Significant = significant,
      SigLabel = if_else(Significant, "TRUE", "FALSE"),
      Species = species_name
    )
  
  month_levels <- c(
    paste("prev.", c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")),
    paste("cur.", c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  )
  
  coef_summary$MonthLabel <- factor(coef_summary$MonthLabel, levels = month_levels)
  coef_summary$SigLabel <- factor(coef_summary$SigLabel, levels = c("TRUE", "FALSE"))
  
  return(coef_summary)
}


#### plot function
plot_combined_species <- function(combined_df,
                                  save_path = NULL,
                                  width = 8,
                                  height = 6,
                                  dpi = 300,
                                  unify_y = FALSE,
                                  y_step = 0.4,
                                  symmetric = FALSE  # set TRUE to force +/- same limit around 0
) {
  # compute shared y-axis limits across both species
  if (unify_y) {
    # Use CI bounds if available; fall back to coef range
    y_min_raw <- if ("ci_lower" %in% names(combined_df)) min(combined_df$ci_lower, na.rm = TRUE) else min(combined_df$coef, na.rm = TRUE)
    y_max_raw <- if ("ci_upper" %in% names(combined_df)) max(combined_df$ci_upper, na.rm = TRUE) else max(combined_df$coef, na.rm = TRUE)
    
    if (symmetric) {
      y_ext <- max(abs(c(y_min_raw, y_max_raw)))
      y_min_raw <- -y_ext
      y_max_raw <-  y_ext
    }
    
    # Round to nearest step outward
    y_min <- floor(y_min_raw / y_step) * y_step
    y_max <- ceiling(y_max_raw / y_step) * y_step
    y_breaks <- seq(y_min, y_max, by = y_step)
  }
  
  p <- ggplot(combined_df, aes(x = MonthLabel, y = coef, color = coef)) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, linetype = SigLabel), width = 0.2, alpha = 0.8) +
    geom_point(aes(shape = SigLabel), size = 4, stroke = 1) +
    scale_color_gradient2(low = "#ca0020", mid = "white", high = "#0571b0", midpoint = 0) +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), name = "Significance") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), name = "Significance") +
    labs(
      x = "Month",
      y = "Correlation Coefficient",
      color = "Coefficient"
    ) +
    facet_wrap(~Species, ncol = 1, scales = if (unify_y) "fixed" else "free_y") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10),
      plot.title = element_blank()
    )
  
  if (unify_y) {
    p <- p + scale_y_continuous(limits = c(y_min, y_max), breaks = y_breaks)
  }
  
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = p, width = width, height = height, dpi = dpi)
    message("Plot saved to: ", save_path)
  }
  
  p
}




####process
alnusp_summary <- process_dcc_summary(acp, "Alnus")
salixp_summary <- process_dcc_summary(scp, "Salix")


#combine
combined_df <- bind_rows(alnusp_summary, salixp_summary) %>%
  mutate(Species = factor(Species, levels = c("Alnus", "Salix")))



#### plot
plot_static <- plot_combined_species(combined_df, unify_y = TRUE, y_step = 0.4)
plot_static

ggsave(file.path(figures_dir, "Salix_Alnus_precipstaticDCC.png"),
       plot_static, width = 10, height = 6, dpi = 300)






################# precipitation-growth moving windows plots
#correlation analyses for each species
###salix
smcp <- dcc(chrono = saspl10SV, climate = dprecip, selection = -10:8,   #selection indicates previous year months (e.g-6, and current years months (e.g. 7))
            method = "correlation", dynamic = "moving", win_size = 20, win_offset = 1, start_last = FALSE,
            timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE) #this is the main function in treeclim
coef <- coef(smcp) #model coefficients 
plot(smcp) #plot the model coefficients
summary(smcp)


###alnus
amcp <- dcc(chrono = alspl10SV, climate = dprecip, selection = -10:8,   #selection indicates previous year months (e.g-6, and current years months (e.g. 7))
            method = "correlation", dynamic = "moving", win_size = 20, win_offset = 1, start_last = FALSE,
            timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE) #this is the main function in treeclim
coef <- coef(amcp) #model coefficients 
plot(amcp) #plot the model coefficients
summary(amcp)



#############save moving corr output as tables
# Extract summary table from acp
amcp_summary <- summary(amcp)

# Save to CSV
write.csv(amcp_summary, file.path(output_dir, "Alnus_moving_precipdcc_summary.csv"), row.names = FALSE)

# Extract summary table from ac
smcp_summary <- summary(smcp)


# Save to CSV
write.csv(smcp_summary, file.path(output_dir, "Salix_moving_precipdcc_summary.csv"), row.names = FALSE)

# Build combined moving-correlation dataframe for precip (mcc_df)
salix_precip_moving <- process_moving_dcc(smcp, "Salix", label_style = "short")
alnus_precip_moving <- process_moving_dcc(amcp, "Alnus", label_style = "short")

mcc_df <- bind_rows(alnus_precip_moving, salix_precip_moving) %>%
  mutate(
    Species = factor(Species, levels = c("Alnus", "Salix"))
  )





######## plot Salix moving correlations for significant months
# Plot Salix moving correlations for significant month: May only
plot_smcp_may <- function(smcp, title = "Salix") {
  
  # Extract coefficient tables
  coef_df     <- smcp$coef$coef
  ci_lower_df <- smcp$coef$ci_lower
  ci_upper_df <- smcp$coef$ci_upper
  signif_df   <- smcp$coef$significant
  
  # Convert to long format
  pivot_coef_data <- function(df, value_name) {
    df %>%
      tibble::rownames_to_column(var = "Month") %>%
      tidyr::pivot_longer(-Month, names_to = "Window", values_to = value_name)
  }
  
  coef_long     <- pivot_coef_data(coef_df, "Correlation")
  ci_lower_long <- pivot_coef_data(ci_lower_df, "CI_Lower")
  ci_upper_long <- pivot_coef_data(ci_upper_df, "CI_Upper")
  signif_long   <- pivot_coef_data(signif_df, "Significant")
  
  # Merge tables
  combined <- coef_long %>%
    dplyr::left_join(ci_lower_long, by = c("Month", "Window")) %>%
    dplyr::left_join(ci_upper_long, by = c("Month", "Window")) %>%
    dplyr::left_join(signif_long,   by = c("Month", "Window"))
  
  # Filter for May only
  combined <- combined %>%
    dplyr::filter(Month == "X3.curr.may") %>%
    dplyr::mutate(
      MonthLabel  = "cur. May",
      SigLineType = if_else(Significant, "solid", "dashed"),
      SigShape    = if_else(Significant, 16, 1)
    )
  
  # Plot
  ggplot(combined, aes(x = Window, y = Correlation, color = Correlation)) +
    geom_errorbar(
      aes(ymin = CI_Lower, ymax = CI_Upper, linetype = SigLineType),
      width = 0.3, alpha = 0.8
    ) +
    geom_point(aes(shape = SigShape), size = 3, stroke = 1) +
    facet_wrap(~ MonthLabel, ncol = 1, scales = "free_y") +
    scale_color_gradient2(low = "#ca0020", mid = "white", high = "#0571b0", midpoint = 0) +
    scale_linetype_identity() +
    scale_shape_identity() +
    labs(title = title, x = "Window", y = "Correlation") +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      strip.text   = element_text(face = "bold"),
      plot.title   = element_text(hjust = 0.5)
    )
}

# Run and save
p_smcp <- plot_smcp_may(smcp, title = "Salix")
print(p_smcp)

ggsave(
  filename = file.path(figures_dir, "Salix_precipMCC_may.png"),
  plot     = p_smcp,
  width    = 10,
  height   = 3,
  dpi      = 300
)



######### plot Alnus moving correlations for significant months

###### filter and plot only previous November and current August
plot_amcp_nov_aug <- function(amcp, title = "Alnus") {
  # Extract coefficient tables
  coef_df <- amcp$coef$coef
  ci_lower_df <- amcp$coef$ci_lower
  ci_upper_df <- amcp$coef$ci_upper
  signif_df <- amcp$coef$significant
  
  ###### convert to long format
  pivot_coef_data <- function(df, value_name) {
    df %>%
      tibble::rownames_to_column(var = "Month") %>%
      tidyr::pivot_longer(-Month, names_to = "Window", values_to = value_name)
  }
  
  coef_long <- pivot_coef_data(coef_df, "Correlation")
  ci_lower_long <- pivot_coef_data(ci_lower_df, "CI_Lower")
  ci_upper_long <- pivot_coef_data(ci_upper_df, "CI_Upper")
  signif_long <- pivot_coef_data(signif_df, "Significant")
  
  ####  merge tables
  combined <- coef_long %>%
    dplyr::left_join(ci_lower_long, by = c("Month", "Window")) %>%
    dplyr::left_join(ci_upper_long, by = c("Month", "Window")) %>%
    dplyr::left_join(signif_long, by = c("Month", "Window"))
  
  ##### filter for  target months
  combined <- combined %>%
    dplyr::filter(Month %in% c("X3.prev.nov", "X3.curr.aug")) %>%
    dplyr::mutate(
      MonthLabel = dplyr::recode(Month,
                                 "X3.prev.nov" = "prev. November",
                                 "X3.curr.aug" = "cur. August"),
      SigLineType = if_else(Significant, "solid", "dashed"),
      SigShape = if_else(Significant, 16, 1)
    )
  
  #### plot
  ggplot(combined, aes(x = Window, y = Correlation, color = Correlation)) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper, linetype = SigLineType),
                           width = 0.3, alpha = 0.8) +
    geom_point(aes(shape = SigShape), size = 3, stroke = 1) +
    facet_wrap(~ MonthLabel, ncol = 1, scales = "free_y") +
    scale_color_gradient2(low = "#ca0020", mid = "white", high = "#0571b0", midpoint = 0) +
    scale_linetype_identity() +
    scale_shape_identity() +
    labs(title = title, x = "Window", y = "Correlation") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5)
    )
}

#### run
p_amcp <- plot_amcp_nov_aug(amcp, title = "Alnus")

print(p_amcp)
p_amcp
ggsave(
  filename = file.path(figures_dir, "Alnus_precipMCC.png"),
  plot     = p_amcp,
  width    = 10,
  height   = 6,
  dpi      = 300
)







###### plot Alnus and Salix precip moving corr together


# so we filter using labels
salix_may <- dplyr::filter(mcc_df, Species == "Salix", MonthLabel == "cur. May")
alnus_aug <- dplyr::filter(mcc_df, Species == "Alnus", MonthLabel == "cur. Aug")
alnus_nov <- dplyr::filter(mcc_df, Species == "Alnus", MonthLabel == "prev. Nov")

# Shared colour limits across all panels (symmetric around 0)
col_lims <- range(mcc_df$Correlation, na.rm = TRUE)
col_lims <- c(-max(abs(col_lims)), max(abs(col_lims)))

# Shared x-axis breaks (every 2nd window label, same as temp code)
window_info <- mcc_df %>%
  dplyr::distinct(Window) %>%
  dplyr::mutate(
    Window   = trimws(Window),
    start_yr = as.integer(sub("-.*", "", Window)),
    end_yr   = as.integer(sub(".*-", "", Window))
  ) %>%
  dplyr::arrange(start_yr, end_yr)

window_levels <- window_info$Window
breaks_vec <- window_levels[seq(1, length(window_levels), by = 2)]

# Single-panel plotting function using temperature-style SigLabel
plot_mcc_single <- function(df, title, col_lims, breaks_vec) {
  # ensure Window is ordered consistently across panels
  df <- df %>% dplyr::mutate(Window = factor(Window, levels = window_levels))
  
  ggplot(df, aes(x = Window, y = Correlation, color = Correlation)) +
    geom_errorbar(
      aes(ymin = CI_Lower, ymax = CI_Upper, linetype = SigLabel),
      width = 0.3, alpha = 0.8
    ) +
    geom_point(aes(shape = SigLabel), size = 3, stroke = 1) +
    scale_color_gradient2(
      low = "#ca0020", mid = "white", high = "#0571b0", midpoint = 0,
      limits = col_lims
    ) +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = "none") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), guide = "none") +
    scale_x_discrete(breaks = breaks_vec) +
    labs(title = title, x = NULL, y = "Correlation") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(hjust = 0.5, size = 8, face = "bold")
    )
}

p_salix_may <- plot_mcc_single(salix_may, "cur. May", col_lims, breaks_vec)
p_alnus_aug <- plot_mcc_single(alnus_aug, "cur. Aug", col_lims, breaks_vec)
p_alnus_nov <- plot_mcc_single(alnus_nov, "prev. Nov", col_lims, breaks_vec)

# Layout:
# Row 1: Alnus Nov | blank
# Row 2: Alnus Aug | Salix May
p_combo <- (p_alnus_nov | patchwork::plot_spacer()) /
  (p_alnus_aug | p_salix_may) +
  plot_layout(
    heights = c(1, 1),
    widths  = c(1.2, 0.8),
    guides  = "collect"
  ) &
  theme(legend.position = "right")

p_combo

ggsave(
  filename = file.path(figures_dir, "Salix_Alnus_precipMCC_panels.png"),
  plot     = p_combo,
  width    = 12, height = 6, dpi = 300
)






##########age cohort split alnus and salix climate-growth relationships

#Read data
al_spline10_trunc <- csv2rwl(fname = file.path(output_dir, "al_spline10_trunc.csv"))


ring_al_age_matrix <- function(rwl_mat) {
  out <- rwl_mat * NA_real_
  for (j in seq_len(ncol(rwl_mat))) {
    idx <- which(!is.na(rwl_mat[, j]))
    if (length(idx) > 0) out[idx, j] <- seq_along(idx)
  }
  out
}

al_age_mat <- ring_al_age_matrix(al_spline10_trunc)
al_age_mat


###create old and young age cohorts
al_young <- al_spline10_trunc
al_young[al_age_mat > 20] <- NA
al_young <- al_young[rowSums(!is.na(al_young)) >= 5, ]

al_old <- al_spline10_trunc
al_old[al_age_mat < 21] <- NA
al_old <- al_old[rowSums(!is.na(al_old)) >= 5, ]


al_spline10_trunc_young <- chron.stabilized(al_young, winLength = 20)
al_spline10_trunc_old   <- chron.stabilized(al_old,   winLength = 20)


alnus_young_crn    <- al_spline10_trunc_young[, "vsc", drop = FALSE]
alnus_old_crn      <- al_spline10_trunc_old[,   "vsc", drop = FALSE]

stopifnot(is.numeric(as.numeric(rownames(alnus_young_crn))))


dclim <- read.csv(file.path(output_dir, "dclim.csv"), check.names = FALSE)

dclim
clim_mat <- dclim %>%
  dplyr::select(Year, `1`:`12`) %>%
  dplyr::mutate(dplyr::across(`1`:`12`, as.numeric)) %>%
  dplyr::filter(stats::complete.cases(dplyr::across(`1`:`12`))) %>%
  tibble::column_to_rownames("Year")


al_clim_young <- clim_mat[rownames(alnus_young_crn), , drop = FALSE]
al_clim_old   <- clim_mat[rownames(alnus_old_crn),   , drop = FALSE]


al_clim_old <- al_clim_old %>%
  tibble::rownames_to_column(var = "Year")

al_clim_young <- al_clim_young %>%
  tibble::rownames_to_column(var = "Year")


# Ensure all values are numeric
al_clim_young <- al_clim_young %>%
  dplyr::mutate_if(is.character, as.numeric)
al_clim_young

al_clim_old <- al_clim_old %>%
  dplyr::mutate_if(is.character, as.numeric)
al_clim_old


al_clim_young <- al_clim_young %>%
  dplyr::mutate(Year = as.numeric(Year)) %>%
  dplyr::filter(is.finite(Year))

al_clim_old <- al_clim_old %>%
  dplyr::mutate(Year = as.numeric(Year)) %>%
  dplyr::filter(is.finite(Year))



###static correlations
al_dcc_young <- dcc(
  chrono = alnus_young_crn,
  climate = al_clim_young,
  selection = -1:9,
  method = "correlation",
  dynamic = "static",
  win_size = 20,
  win_offset = 1, start_last = FALSE,
  timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE) 

summary(al_dcc_young)


al_dcc_old <- dcc(
  chrono = alnus_old_crn,
  climate = al_clim_old,
  selection = -1:9,
  method = "correlation",
  dynamic = "static",
  win_size = 20,
  win_offset = 1, start_last = FALSE,
  timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE) 

summary(al_dcc_old)



#moving correlations
al_mc_young <- dcc(
  chrono = alnus_young_crn,
  climate = al_clim_young,
  selection = 6:7,
  method = "correlation",
  dynamic = "moving",
  win_size = 20, win_offset = 1, start_last = FALSE,
  ci = 0.05, boot = "std", sb = FALSE
)

al_mc_old <- dcc(
  chrono = alnus_old_crn,
  climate = al_clim_old,
  selection = 6:7,
  method = "correlation",
  dynamic = "moving",
  win_size = 20, win_offset = 1, start_last = FALSE,
  ci = 0.05, boot = "std", sb = FALSE
)

# climate complete
stopifnot(all(complete.cases(al_clim_young[, as.character(1:12)])))
stopifnot(all(complete.cases(al_clim_old[,   as.character(1:12)])))

# year overlap with chronologies
range(al_clim_young$Year); range(as.numeric(rownames(alnus_young_crn)))
range(al_clim_old$Year);   range(as.numeric(rownames(alnus_old_crn)))

# number of years used in each model
summary(al_dcc_young)$data$timespan
# and after you run dcc_old:
summary(al_dcc_old)$data$timespan




###plot common coverage for both age cohorts
keep_yrs <- 1988:2023

#restrict young chronology to keep_yrs
alnus_young_crn_88 <- alnus_young_crn[rownames(alnus_young_crn) %in% keep_yrs, , drop = FALSE]

#restrict young climate to keep_yrs AND keep the correct format for treeclim
al_clim_young_88 <- al_clim_young %>%
  dplyr::filter(Year %in% keep_yrs) %>%
  dplyr::select(Year, `1`:`12`) %>%
  dplyr::mutate(dplyr::across(`1`:`12`, as.numeric)) %>%
  dplyr::filter(stats::complete.cases(dplyr::select(., `1`:`12`)))

#enforce intersection (important if any years dropped by complete.cases)
common_yrs <- intersect(as.numeric(rownames(alnus_young_crn_88)), al_clim_young_88$Year)

alnus_young_crn_88 <- alnus_young_crn_88[as.character(common_yrs), , drop = FALSE]
al_clim_young_88   <- al_clim_young_88 %>% dplyr::filter(Year %in% common_yrs)

#checks
stopifnot(nrow(alnus_young_crn_88) == nrow(al_clim_young_88))
stopifnot(all(al_clim_young_88$Year == as.numeric(rownames(alnus_young_crn_88))))

#moving correlation (Jun/Jul)
al_mc_young_88 <- dcc(
  chrono = alnus_young_crn_88,
  climate = al_clim_young_88,
  selection = 6:7,
  method = "correlation", dynamic = "moving",
  win_size = 20, win_offset = 1, start_last = FALSE,
  ci = 0.05, boot = "std", sb = FALSE
)

plot(al_mc_young_88)
summary(al_mc_young_88)


###create long dataframe
dcc_moving_to_long <- function(dcc_obj, age_class) {
  
  pivot_coef_data <- function(df, value_name) {
    df %>%
      rownames_to_column(var = "Month") %>%
      pivot_longer(-Month, names_to = "Window", values_to = value_name)
  }
  
  coef_long     <- pivot_coef_data(dcc_obj$coef$coef,        "Correlation")
  ci_lower_long <- pivot_coef_data(dcc_obj$coef$ci_lower,    "CI_Lower")
  ci_upper_long <- pivot_coef_data(dcc_obj$coef$ci_upper,    "CI_Upper")
  signif_long   <- pivot_coef_data(dcc_obj$coef$significant, "Significant")
  
  coef_long %>%
    left_join(ci_lower_long, by = c("Month", "Window")) %>%
    left_join(ci_upper_long, by = c("Month", "Window")) %>%
    left_join(signif_long,   by = c("Month", "Window")) %>%
    mutate(
      AgeClass    = age_class,
      Window      = factor(Window, levels = unique(Window)),
      SigLineType = if_else(Significant, "solid", "dashed"),
      SigShape    = if_else(Significant, 16, 1)
    )
}


###create age cohort dataframes
al_young_df <- dcc_moving_to_long(al_mc_young_88, "Young (≤20), 1988+")
al_old_df   <- dcc_moving_to_long(al_mc_old,      "Old (≥21)")

al_moving_age_df <- bind_rows(al_young_df, al_old_df) %>%
  mutate(AgeClass = factor(AgeClass, levels = c("Young (≤20), 1988+", "Old (≥21)")))



###plot
plot_moving_jun_jul_by_age <- function(df, title = NULL, add_trend = TRUE) {
  
  plot_df <- df %>%
    filter(Month %in% c("X3.curr.jun", "X3.curr.jul")) %>%
    mutate(
      MonthLabel = recode(Month,
                          "X3.curr.jun" = "June",
                          "X3.curr.jul" = "July"),
      MonthLabel = factor(MonthLabel, levels = c("June", "July"))  # controls row order
    )
  
  p <- ggplot(plot_df, aes(x = Window, y = Correlation, color = Correlation)) +
    geom_errorbar(
      aes(ymin = CI_Lower, ymax = CI_Upper, linetype = SigLineType),
      width = 0.3, alpha = 0.8
    ) +
    geom_point(
      aes(shape = SigShape),
      size = 3, stroke = 1
    ) +
    facet_grid(MonthLabel ~ AgeClass) +
    scale_color_gradient2(low = "#ca0020", mid = "white", high = "#0571b0", midpoint = 0) +
    scale_linetype_identity() +
    scale_shape_identity() +
    labs(title = title, x = "Window", y = "Correlation") +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      strip.text   = element_text(face = "bold"),
      plot.title   = element_text(hjust = 0.5)
    )
  
  if (add_trend) {
    p <- p +
      geom_smooth(
        data = plot_df %>%
          mutate(WindowIndex = as.numeric(Window)),
        aes(x = WindowIndex, y = Correlation, group = interaction(MonthLabel, AgeClass)),
        inherit.aes = FALSE,
        method = "lm",
        se = FALSE,
        linetype = "dashed",
        linewidth = 0.5,
        color = "black"
      )
  }
  
  p
}

al_moving_age_plot <- plot_moving_jun_jul_by_age(
  al_moving_age_df,
  title = "Alnus",
  add_trend = TRUE
)

print(al_moving_age_plot)
ggsave(file.path(figures_dir, "Alnus_agecohortsMCC.png"),
       al_moving_age_plot, width = 10, height = 6, dpi = 300)




##################Salix

#Read data
sa_spline10_trunc <- csv2rwl(fname = file.path(output_dir, "sa_spline10_trunc.csv"))

ring_sa_age_matrix <- function(rwl_mat) {
  out <- rwl_mat * NA_real_
  for (j in seq_len(ncol(rwl_mat))) {
    idx <- which(!is.na(rwl_mat[, j]))
    if (length(idx) > 0) out[idx, j] <- seq_along(idx)
  }
  out
}

sa_age_mat <- ring_sa_age_matrix(sa_spline10_trunc)
sa_age_mat


###create old and young age cohorts
sa_young <- sa_spline10_trunc
sa_young[sa_age_mat > 20] <- NA
sa_young <- sa_young[rowSums(!is.na(sa_young)) >= 5, ]

sa_old <- sa_spline10_trunc
sa_old[sa_age_mat < 21] <- NA
sa_old <- sa_old[rowSums(!is.na(sa_old)) >= 5, ]


sa_saspl10SV_young <- chron.stabilized(sa_young, winLength = 20)
sa_saspl10SV_old   <- chron.stabilized(sa_old,   winLength = 20)

sa_salix_young_crn <- sa_saspl10SV_young[, "vsc", drop = FALSE]
sa_salix_old_crn   <- sa_saspl10SV_old[,   "vsc", drop = FALSE]


stopifnot(is.numeric(as.numeric(rownames(sa_salix_young_crn))))


sa_dclim <- read.csv(file.path(output_dir, "dclim.csv"), check.names = FALSE)
sa_dclim
sa_clim_mat <- sa_dclim %>%
  dplyr::select(Year, `1`:`12`) %>%
  dplyr::mutate(dplyr::across(`1`:`12`, as.numeric)) %>%
  dplyr::filter(stats::complete.cases(dplyr::across(`1`:`12`))) %>%
  tibble::column_to_rownames("Year")


sa_clim_young <- sa_clim_mat[rownames(sa_salix_young_crn), , drop = FALSE]
sa_clim_old   <- sa_clim_mat[rownames(sa_salix_old_crn),   , drop = FALSE]

# Convert rownames to a column
sa_clim_young <- sa_clim_young %>%
  rownames_to_column(var = "Year")

sa_clim_old <- sa_clim_old %>%
  rownames_to_column(var = "Year")


# Ensure all values are numeric
sa_clim_young <- sa_clim_young %>%
  mutate_if(is.character, as.numeric)
sa_clim_young

sa_clim_old <- sa_clim_old %>%
  mutate_if(is.character, as.numeric)
sa_clim_old


al_clim_young <- al_clim_young %>%
  dplyr::mutate(Year = as.numeric(Year)) %>%
  dplyr::filter(is.finite(Year))

al_clim_old <- al_clim_old %>%
  dplyr::mutate(Year = as.numeric(Year)) %>%
  dplyr::filter(is.finite(Year))


###static correlations
sa_dcc_young <- dcc(
  chrono = sa_salix_young_crn,
  climate = sa_clim_young,
  selection = -1:9,
  method = "correlation",
  dynamic = "static",
  win_size = 20,
  win_offset = 1, start_last = FALSE,
  timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE) 

summary(sa_dcc_young)

sa_dcc_old <- dcc(
  chrono = sa_salix_old_crn,
  climate = sa_clim_old,
  selection = -1:9,
  method = "correlation",
  dynamic = "static",
  win_size = 20,
  win_offset = 1, start_last = FALSE,
  timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE) 

summary(sa_dcc_old)





#moving correlations
sa_mc_young <- dcc(
  chrono = sa_salix_young_crn,
  climate = sa_clim_young,
  selection = 6:7,
  method = "correlation",
  dynamic = "moving",
  win_size = 20, win_offset = 1, start_last = FALSE,
  ci = 0.05, boot = "std", sb = FALSE
)

sa_mc_old <- dcc(
  chrono = sa_salix_old_crn,
  climate = sa_clim_old,
  selection = 6:7,
  method = "correlation",
  dynamic = "moving",
  win_size = 20, win_offset = 1, start_last = FALSE,
  ci = 0.05, boot = "std", sb = FALSE
)

# climate complete
stopifnot(all(complete.cases(sa_clim_young[, as.character(1:12)])))
stopifnot(all(complete.cases(sa_clim_old[,   as.character(1:12)])))

# year overlap with chronologies
range(sa_clim_young$Year); range(as.numeric(rownames(sa_salix_young_crn)))
range(sa_clim_old$Year);   range(as.numeric(rownames(sa_salix_old_crn)))

# number of years used in each model
summary(sa_dcc_young)$data$timespan
# and after you run dcc_old:
summary(sa_dcc_old)$data$timespan





###plot common coverage for both age cohorts
keep_yrs <- 1993:2023

#restrict young chronology to keep_yrs
sa_salix_young_crn_93 <- sa_salix_young_crn[rownames(sa_salix_young_crn) %in% keep_yrs, , drop = FALSE]

#restrict young climate to keep_yrs AND keep the correct format for treeclim
sa_clim_young_93 <- sa_clim_young %>%
  dplyr::filter(Year %in% keep_yrs) %>%
  dplyr::select(Year, `1`:`12`) %>%
  dplyr::mutate(dplyr::across(`1`:`12`, as.numeric)) %>%
  dplyr::filter(stats::complete.cases(dplyr::select(., `1`:`12`)))

#enforce intersection (important if any years dropped by complete.cases)
common_yrs <- intersect(as.numeric(rownames(sa_salix_young_crn_93)), sa_clim_young_93$Year)

sa_salix_young_crn_93 <- sa_salix_young_crn_93[as.character(common_yrs), , drop = FALSE]
sa_clim_young_93      <- sa_clim_young_93 %>% dplyr::filter(Year %in% common_yrs)

# checks
stopifnot(nrow(sa_salix_young_crn_93) == nrow(sa_clim_young_93))
stopifnot(all(sa_clim_young_93$Year == as.numeric(rownames(sa_salix_young_crn_93))))

#moving correlation (Jun/Jul)
sa_mc_young_93 <- dcc(
  chrono = sa_salix_young_crn_93,
  climate = sa_clim_young_93,
  selection = 6:7,
  method = "correlation", dynamic = "moving",
  win_size = 20, win_offset = 1, start_last = FALSE,
  ci = 0.05, boot = "std", sb = FALSE
)

plot(sa_mc_young_93)
summary(sa_mc_young_93)


dcc_moving_to_long <- function(dcc_obj, age_class) {
  
  pivot_coef_data <- function(df, value_name) {
    df %>%
      rownames_to_column(var = "Month") %>%
      pivot_longer(-Month, names_to = "Window", values_to = value_name)
  }
  
  coef_long     <- pivot_coef_data(dcc_obj$coef$coef,        "Correlation")
  ci_lower_long <- pivot_coef_data(dcc_obj$coef$ci_lower,    "CI_Lower")
  ci_upper_long <- pivot_coef_data(dcc_obj$coef$ci_upper,    "CI_Upper")
  signif_long   <- pivot_coef_data(dcc_obj$coef$significant, "Significant")
  
  coef_long %>%
    left_join(ci_lower_long, by = c("Month", "Window")) %>%
    left_join(ci_upper_long, by = c("Month", "Window")) %>%
    left_join(signif_long,   by = c("Month", "Window")) %>%
    mutate(
      AgeClass    = age_class,
      Window      = factor(Window, levels = unique(Window)),
      SigLineType = if_else(Significant, "solid", "dashed"),
      SigShape    = if_else(Significant, 16, 1)
    )
}

sa_young_df <- dcc_moving_to_long(sa_mc_young_93, "Young (≤20), 1993+")
sa_old_df   <- dcc_moving_to_long(sa_mc_old,      "Old (≥21)")

sa_moving_age_df <- bind_rows(sa_young_df, sa_old_df) %>%
  mutate(AgeClass = factor(AgeClass, levels = c("Young (≤20), 1993+", "Old (≥21)")))



plot_moving_jun_jul_by_age <- function(df, title = NULL, add_trend = TRUE) {
  
  plot_df <- df %>%
    filter(Month %in% c("X3.curr.jun", "X3.curr.jul")) %>%
    mutate(
      MonthLabel = recode(Month,
                          "X3.curr.jun" = "June",
                          "X3.curr.jul" = "July"),
      MonthLabel = factor(MonthLabel, levels = c("June", "July"))  # controls row order
    )
  
  p <- ggplot(plot_df, aes(x = Window, y = Correlation)) +
    geom_errorbar(
      aes(ymin = CI_Lower, ymax = CI_Upper, linetype = SigLineType),
      width = 0.3, alpha = 0.5
    ) +
    geom_point(
      aes(shape = SigShape, color = Correlation),
      size = 3, stroke = 1
    ) +
    facet_grid(MonthLabel ~ AgeClass) +
    scale_color_gradient2(low = "#ca0020", mid = "white", high = "#0571b0", midpoint = 0) +
    scale_linetype_identity() +
    scale_shape_identity() +
    labs(title = title, x = "Window", y = "Correlation") +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      strip.text   = element_text(face = "bold"),
      plot.title   = element_text(hjust = 0.5)
    )
  
  if (add_trend) {
    p <- p +
      geom_smooth(
        data = plot_df %>%
          mutate(WindowIndex = as.numeric(Window)),
        aes(x = WindowIndex, y = Correlation, group = interaction(MonthLabel, AgeClass)),
        inherit.aes = FALSE,
        method = "lm",
        se = FALSE,
        linetype = "dashed",
        linewidth = 0.5,
        color = "black"
      )
  }
  
  p
}



sa_moving_age_plot <- plot_moving_jun_jul_by_age(
  sa_moving_age_df,
  title = "Salix",
  add_trend = TRUE
)

print(sa_moving_age_plot)
ggsave(file.path(figures_dir, "Salix_agecohortsMCC.png"),
       sa_moving_age_plot, width = 10, height = 6, dpi = 300)



