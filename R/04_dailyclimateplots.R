##############################################################################
# Alnus and Salix shrub anatomical response to climate
#########Climate plots - cold degree days, cumulative growing degree days, cumulative growing degree days, precipitation
library(here)

source(here("R/utils_packages.R"))
load_project_packages()

source(here("R/00_paths_setup.R")) 

#project data folders
data_dir        <- here("data")
anomalies_dir   <- here("data", "anomalies")
rwi_dir         <- here("data", "rwi")
meta_dir        <- here("data", "sample_metadata")
crosswalk_dir   <- here("data", "id_crosswalk")
climate_dir     <- here("data", "climate_raw")
derived_dir     <- here("data", "derived")

#Create figures folder if needed
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir)
}

#global plot settings
theme_set(theme_minimal(base_size = 13))


###plots of cold degree days CDD and growing degree days GDD for highlight years

###climate data
tmin_file <- file.path(climate_dir, "ERA5_tmin_daily.dat")
readLines(tmin_file, n = 25)
ERA5TM <- read.table(tmin_file, header = FALSE, skip = 24)

#check correct dataframe is made
head(ERA5TM)
colnames(ERA5TM) <- c("Date", "Temp")
head(ERA5TM)

# Convert the date column to Date type
ERA5TM$Date <- ymd(ERA5TM$Date)
head(ERA5TM)

# Extract year, month, and day from the 'time' column, padding single-digit month and day with leading zero
ERA5TM <- ERA5TM %>%
  mutate(year = year(Date),
         month = sprintf("%02d", month(Date)),
         day = sprintf("%02d", day(Date))) %>%
  relocate(year, month, day, .after = Date)
head(ERA5TM)

#Filter data for the range 1970-2022
ETM <- ERA5TM %>%
  filter(year >= 1950 & year <= 2022)
head(ETM)


# Convert Month column to numeric
ETM$month <- as.numeric(ETM$month)
ETM$day <- as.numeric(ETM$day)
ETM$Date <- as.Date(ETM$Date)

# Plot the data with 'time' as the x-axis and 'temp' as the y-axis
ggplot(ETM, aes(x = Date, y = Temp)) +
  geom_line(color = "darkorchid4") +
  labs(title = "Daily minimum temperature (1950-2022)",
       x = "Date",
       y = "Temp") +
  theme_bw(base_size = 15)



###calculate julian dates and restrict to growing season window
ETM$julian <- yday(ETM$Date)  
head(ETM)
tail(ETM)

ETMm <- aggregate(Temp ~ julian, ETM, FUN = mean)
head(ETMm)
tail(ETMm)
colnames(ETMm)[2]<-"mTemp"

ETM<-merge(ETM,ETMm,by="julian")
ETM<-ETM[which(ETM$julian>120 & ETM$julian<280), ]

# Calculate the mean daily temperature for each julian day across all years
mean_temp_all <- ETM %>%
  dplyr::group_by(julian) %>%
  dplyr::summarize(mean_temp = mean(Temp, na.rm = TRUE), .groups = "drop")

head(mean_temp_all)




####CDD

###make the CDD shaded area 
make_cdd_ribbon <- function(data, year_value, base_temp = 5,
                            shade_start = 152, shade_end = 245) {
  
  hd <- data %>%
    dplyr::filter(year == year_value) %>%
    dplyr::arrange(julian) %>%
    dplyr::mutate(below = Temp < base_temp)
  
  if (nrow(hd) < 2) stop("Not enough data for this year.")
  
  ribbon_list <- vector("list", nrow(hd))
  k <- 0
  
  for (i in 1:(nrow(hd) - 1)) {
    
    t1 <- hd$Temp[i];   t2 <- hd$Temp[i + 1]
    x1 <- hd$julian[i]; x2 <- hd$julian[i + 1]
    b1 <- hd$below[i];  b2 <- hd$below[i + 1]
    
    # keep point if below threshold AND inside shading window
    if (b1 && x1 > shade_start && x1 <= shade_end) {
      k <- k + 1
      ribbon_list[[k]] <- data.frame(julian = x1, Temp = t1)
    }
    
    # add interpolated crossing point if threshold crossed
    if (b1 != b2) {
      if (!isTRUE(all.equal(t2, t1))) {
        f <- (base_temp - t1) / (t2 - t1)
        x_cross <- x1 + f * (x2 - x1)
        
        if (x_cross > shade_start && x_cross <= shade_end) {
          k <- k + 1
          ribbon_list[[k]] <- data.frame(julian = x_cross, Temp = base_temp)
        }
      }
    }
  }
  
  # last point
  if (hd$below[nrow(hd)] &&
      hd$julian[nrow(hd)] > shade_start &&
      hd$julian[nrow(hd)] <= shade_end) {
    
    k <- k + 1
    ribbon_list[[k]] <- data.frame(
      julian = hd$julian[nrow(hd)],
      Temp = hd$Temp[nrow(hd)]
    )
  }
  
  ribbon_df <- dplyr::bind_rows(ribbon_list[seq_len(k)]) %>%
    dplyr::arrange(julian) %>%
    dplyr::mutate(
      ymin = Temp,
      ymax = base_temp,
      run_id = cumsum(c(0, diff(julian) > 1.1))
    )
  
  list(hd = hd, ribbon_df = ribbon_df)
}



###plot function for highlight year
plot_cdd_highlight_year <- function(ETM, mean_temp_all,
                                    highlight_year,
                                    base_temp = 5,
                                    outfile = NULL,
                                    width = 8, height = 5, dpi = 300) {
  
  obj <- make_cdd_ribbon(ETM, highlight_year, base_temp)
  hd <- obj$hd
  ribbon_df <- obj$ribbon_df
  
  p <- ggplot() +
    geom_line(data = ETM,
              aes(x = julian, y = Temp, group = year),
              color = "grey80") +
    
    geom_ribbon(
      data = ribbon_df,
      aes(x = julian, ymin = ymin, ymax = ymax, group = run_id),
      fill = "darkorchid4", alpha = 0.3, colour = NA
    ) +
    
    geom_line(data = hd,
              aes(x = julian, y = Temp),
              color = "darkorchid4") +
    
    geom_line(data = mean_temp_all,
              aes(x = julian, y = mean_temp),
              colour = "red", linetype = "dashed") +
    
    geom_hline(yintercept = base_temp,
               linetype = "dashed", color = "green") +
    geom_hline(yintercept = 0,
               linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(153, 183, 214, 245),
               linetype = "dashed", color = "grey") +
    
    annotate("text", x = 170, y = -24, label = "Jun") +
    annotate("text", x = 200, y = -24, label = "Jul") +
    annotate("text", x = 230, y = -24, label = "Aug") +
    
    labs(
      title = paste("Daily ERA5 Tmin -", highlight_year),
      x = "day of year",
      y = "Tmin (°C)"
    ) +
    
    scale_y_continuous(breaks = seq(-30, 40, by = 5)) +
    scale_x_continuous(breaks = seq(0, 366, by = 50)) +
    
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 8),
      legend.position = c(0.9, 0.9),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 10)
    ) +
    guides(linetype = guide_legend(title = NULL))
  
  print(p)
  
  if (!is.null(outfile)) {
    ggsave(outfile, p, width = width, height = height, dpi = dpi)
  }
  
  invisible(p)
}


###CDD plots for years used in publication figures

###saving to /figures using relative paths)

#build output path in repository
save_plot_path <- function(prefix, year, ext = "png") {
  file.path(figures_dir, paste0(prefix, "_highlightyear", year, "_plot.", ext))
}

# Choose output format
out_ext <- "png"   

#1974
plot_cdd_highlight_year(
  ETM, mean_temp_all, 1974,
  outfile = save_plot_path("CDD", 1974, out_ext),
  width = 8, height = 5, dpi = 300
)


#1994
plot_cdd_highlight_year(
  ETM, mean_temp_all, 1994,
  outfile = save_plot_path("CDD", 1994, out_ext),
  width = 8, height = 5, dpi = 300
)

#1996
plot_cdd_highlight_year(
  ETM, mean_temp_all, 1996,
  outfile = save_plot_path("CDD", 1996, out_ext),
  width = 8, height = 5, dpi = 300
)


#2000
plot_cdd_highlight_year(
  ETM, mean_temp_all, 2000,
  outfile = save_plot_path("CDD", 2000, out_ext),
  width = 8, height = 5, dpi = 300
)

#2002
plot_cdd_highlight_year(
  ETM, mean_temp_all, 2002,
  outfile = save_plot_path("CDD", 2002, out_ext),
  width = 8, height = 5, dpi = 300
)

#2012
plot_cdd_highlight_year(
  ETM, mean_temp_all, 2012,
  outfile = save_plot_path("CDD", 2012, out_ext),
  width = 8, height = 5, dpi = 300
)

#2015
plot_cdd_highlight_year(
  ETM, mean_temp_all, 2015,
  outfile = save_plot_path("CDD", 2015, out_ext),
  width = 8, height = 5, dpi = 300
)

#2017
plot_cdd_highlight_year(
  ETM, mean_temp_all, 2017,
  outfile = save_plot_path("CDD", 2017, out_ext),
  width = 8, height = 5, dpi = 300
)

#2018
plot_cdd_highlight_year(
  ETM, mean_temp_all, 2018,
  outfile = save_plot_path("CDD", 2018, out_ext),
  width = 8, height = 5, dpi = 300
)







####
########
###########
###############GDD - Growing Degree Days (days above 5C) plot for highlight year

#Load ERA5 t2m temperature data
t2m_file <- file.path(climate_dir, "ERA5_t2m_daily.dat")
readLines(t2m_file, n = 25)

ERA5T2M <- read.table(t2m_file, header = FALSE, skip = 22)
colnames(ERA5T2M) <- c("Date", "Temp")  # daily mean 2m temperature
ERA5T2M$Date <- ymd(ERA5T2M$Date)

# Extract year, month, day, julian day
ERA5T2M <- ERA5T2M %>%
  mutate(
    year = year(Date),
    month = month(Date),
    day = day(Date),
    julian = yday(Date)
  ) %>%
  relocate(year, month, day, julian, .after = Date)

# Filter for chosen years window (example: 1950-2022)
ET2M <- ERA5T2M %>% filter(year >= 1950 & year <= 2022)


ET2M$julian <- yday(ET2M$Date)  
head(ET2M)
tail(ET2M)

ET2Mm <- aggregate(Temp ~ julian, ET2M, FUN = mean)
head(ET2Mm)
tail(ET2Mm)
colnames(ET2Mm)[2]<-"mTemp"

ET2M<-merge(ET2M,ET2Mm,by="julian")
ET2M<-ET2M[which(ET2M$julian>120 & ET2M$julian<280), ]


#####plot code for years used in publication figures 1984, 2001, 2015
#shading for GDD
make_threshold_ribbon <- function(data, year_value, base_temp = 5) {
  hd <- data %>%
    dplyr::filter(year == year_value) %>%
    dplyr::arrange(julian) %>%
    dplyr::mutate(above = Temp > base_temp)
  
  if (nrow(hd) < 2) stop("Not enough data for this year.")
  
  ribbon_list <- vector("list", nrow(hd))  # preallocate-ish
  
  k <- 0
  for (i in 1:(nrow(hd) - 1)) {
    t1 <- hd$Temp[i];   t2 <- hd$Temp[i + 1]
    x1 <- hd$julian[i]; x2 <- hd$julian[i + 1]
    a1 <- hd$above[i];  a2 <- hd$above[i + 1]
    
    # keep point i if above
    if (a1) {
      k <- k + 1
      ribbon_list[[k]] <- data.frame(julian = x1, Temp = t1)
    }
    
    # add interpolated crossing point if it crosses the threshold
    if (a1 != a2) {
      # guard against divide-by-zero (flat line)
      if (!isTRUE(all.equal(t2, t1))) {
        f <- (base_temp - t1) / (t2 - t1)
        x_cross <- x1 + f * (x2 - x1)
        k <- k + 1
        ribbon_list[[k]] <- data.frame(julian = x_cross, Temp = base_temp)
      }
    }
  }
  
  # add last point if last day is above
  if (hd$above[nrow(hd)]) {
    k <- k + 1
    ribbon_list[[k]] <- data.frame(julian = hd$julian[nrow(hd)], Temp = hd$Temp[nrow(hd)])
  }
  
  ribbon_df <- dplyr::bind_rows(ribbon_list[seq_len(k)]) %>%
    dplyr::arrange(julian) %>%
    dplyr::mutate(
      ymin = base_temp,
      ymax = Temp,
      run_id = cumsum(c(0, diff(julian) > 1.1))
    )
  
  list(hd = hd, ribbon_df = ribbon_df)
}

###plot code for highlight year
plot_highlight_year <- function(ET2M, ET2Mm, highlight_year,
                                base_temp = 5,
                                outfile = NULL,
                                width = 8, height = 5, dpi = 300) {
  
  obj <- make_threshold_ribbon(ET2M, year_value = highlight_year, base_temp = base_temp)
  hd <- obj$hd
  ribbon_df <- obj$ribbon_df
  
  p <- ggplot() +
    geom_line(data = ET2M, aes(x = julian, y = Temp, group = year), color = "grey80") +
    
    geom_ribbon(
      data = ribbon_df,
      aes(x = julian, ymin = ymin, ymax = ymax, group = run_id),
      fill = "darkgreen", alpha = 0.3, colour = NA
    ) +
    
    geom_line(data = hd, aes(x = julian, y = Temp), color = "darkgreen") +
    
    geom_line(data = ET2Mm, aes(x = julian, y = mTemp),
                       colour = "red", linetype = "dashed") +
    
    geom_hline(yintercept = base_temp, linetype = "dashed", color = "green") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(153, 183, 214, 245), linetype = "dashed", color = "grey") +
    
    annotate("text", x = 170, y = min(ET2M$Temp, na.rm = TRUE) - 1, label = "Jun") +
    annotate("text", x = 200, y = min(ET2M$Temp, na.rm = TRUE) - 1, label = "Jul") +
    annotate("text", x = 230, y = min(ET2M$Temp, na.rm = TRUE) - 1, label = "Aug") +
    
    annotate("text", x = 130, y = max(ET2M$Temp, na.rm = TRUE) - 1,
                      label = paste(highlight_year), color = "black", size = 4, hjust = 0) +
    
    labs(x = "Day of year", y = "T2m (°C)") +
    scale_y_continuous(breaks = seq(-30, 40, by = 5)) +
    scale_x_continuous(limits = c(120, 280), breaks = seq(120, 280, by = 25)) +
    
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(size = 8),
      legend.position = c(0.9, 0.9),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      plot.title = element_blank()
    ) +
    guides(linetype = guide_legend(title = NULL))
  
  print(p)
  
  if (!is.null(outfile)) {
    ggsave(outfile, p, width = width, height = height, dpi = dpi)
  }
  
  invisible(p)
}



###plot and save for years used in publication figure

#1984
plot_highlight_year(
  ET2M, ET2Mm, 1984,
  base_temp = 5,
  outfile = file.path(figures_dir, "GDD_highlightyear1984_plot.png"),
  width = 10,
  height = 5,
  dpi = 300
)

#2001
plot_highlight_year(
  ET2M, ET2Mm, 2001,
  base_temp = 5,
  outfile = file.path(figures_dir, "GDD_highlightyear2001_plot.jpg"),
  width = 8,
  height = 5,
  dpi = 300
)

#2015
plot_highlight_year(
  ET2M, ET2Mm, 2015,
  base_temp = 5,
  outfile = file.path(figures_dir, "GDD_highlightyear2015_plot.jpg"),
  width = 8,
  height = 5,
  dpi = 300
)


##########
#############
#################

####### cGDD - cumulative growing degree days plot for highlight year


ERA5T2M <- read.table(t2m_file, header = FALSE, skip = 22)
colnames(ERA5T2M) <- c("Date", "Temp")  # daily mean 2m temperature
ERA5T2M$Date <- ymd(ERA5T2M$Date)

#Extract year, month, day, julian day
ERA5T2M <- ERA5T2M %>%
  mutate(
    year = year(Date),
    month = month(Date),
    day = day(Date),
    julian = yday(Date)
  ) %>%
  relocate(year, month, day, julian, .after = Date)

#Filter for chosen years window (example: 1950-2022)
ET2M <- ERA5T2M %>% filter(year >= 1950 & year <= 2022)

#Set base temperature for GDD calculation
base_temp <- 5   
ET2M <- ET2M %>%
  mutate(GDD = pmax(Temp - base_temp, 0))   # only accumulate if Temp > base_temp


#Calculate cumulative GDD for each year
ET2M <- ET2M %>%
  group_by(year) %>%
  arrange(julian) %>%
  mutate(cumGDD = cumsum(GDD)) %>%
  ungroup()

#Calculate mean daily GDD
mean_gdd <- ET2M %>%
  group_by(julian) %>%
  summarise(mean_cumGDD = mean(cumGDD, na.rm = TRUE), .groups = "drop")



###plot for highlight year

plot_cGDD_highlight_year <- function(ET2M, mean_gdd,
                                     highlight_year,
                                     outfile = NULL,
                                     width = 8,
                                     height = 5,
                                     dpi = 300) {
  
  highlight_data <- ET2M %>%
    dplyr::filter(year == highlight_year)
  
  p <- ggplot() +
    
    # Grey lines for all years
    geom_line(data = ET2M,
              aes(x = julian, y = cumGDD, group = year),
              color = "grey80") +
    
    # Highlighted year
    geom_line(data = highlight_data,
              aes(x = julian, y = cumGDD),
              color = "darkgreen", linewidth = 1) +
    
    # Mean cumulative GDD
    geom_line(data = mean_gdd,
              aes(x = julian, y = mean_cumGDD),
              color = "red", linetype = "dashed") +
    
    # Vertical lines for months
    geom_vline(xintercept = c(153, 183, 214, 245),
               linetype = "dashed", color = "grey") +
    
    annotate("text", x = 170,
             y = max(ET2M$cumGDD, na.rm = TRUE) * 0.9,
             label = "Jun") +
    annotate("text", x = 200,
             y = max(ET2M$cumGDD, na.rm = TRUE) * 0.9,
             label = "Jul") +
    annotate("text", x = 230,
             y = max(ET2M$cumGDD, na.rm = TRUE) * 0.9,
             label = "Aug") +
    
    annotate("text",
             x = 50,
             y = max(ET2M$cumGDD, na.rm = TRUE) * 0.9,
             label = highlight_year,
             size = 6,
             color = "black",
             fontface = "plain",
             hjust = 1) +
    
    labs(x = NULL, y = "Cumulative GDD") +
    scale_x_continuous(breaks = seq(0, 366, by = 50)) +
    theme_minimal(base_size = 13)
  
  print(p)
  
  if (!is.null(outfile)) {
    ggsave(outfile, p,
           width = width,
           height = height,
           dpi = dpi)
  }
  
  invisible(p)
}



#Highlight years for publication
highlight_years <- c(1974, 1996, 2002)

for (yr in highlight_years) {
  
  plot_cGDD_highlight_year(
    ET2M,
    mean_gdd,
    highlight_year = yr,
    outfile = file.path(figures_dir,
                        paste0("cGDD_highlightyear", yr, "_plot.png")),
    width = 8,
    height = 5,
    dpi = 300
  )
}









##########
#####Precipitation for highlight anomaly years

#precipitation - using ERA as CRU lacks data for the region
#read in and convert ERA precip data

precip_file <- file.path(climate_dir, "ERA5_prcp_daily.dat")
readLines(precip_file, n = 25)
ERA5P <- read.table(precip_file, header = FALSE, skip = 22)

head(ERA5P)
colnames(ERA5P) <- c("Date", "precip")
head(ERA5P)


# Convert the date column to Date type
ERA5P$Date <- ymd(ERA5P$Date)
head(ERA5P)


# Extract year, month, and day
ERA5P <- ERA5P %>%
  mutate(
    Year = year(Date),
    Month = sprintf("%02d", month(Date)),
    Day = sprintf("%02d", day(Date))
  ) %>%
  relocate(Year, Month, Day, .after = Date)

str(ERA5P$Date)


# Filter data for the range 1948-2023
ERA5P_filtered <- ERA5P %>%
  filter(Year >= 1948 & Year <= 2023)



# Aggregate to monthly totals and reshape wide with months as 1–12 
precip <- ERA5P_filtered %>%
  group_by(Year, Month) %>%
  summarise(precip = sum(precip), .groups = "drop") %>%
  pivot_wider( names_from = Month, values_from = precip ) %>%
  arrange(Year)


# Check result
str(precip)
head(precip)

#Assign column headers
colnames(precip) <- c("Year", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
head(precip)


# convert all month columns to numeric to prevent list-columns
precip <- precip %>%
  mutate(across(-Year, ~as.numeric(.)))

# Check structure
str(precip)

# Rename columns to numeric 1-12
colnames(precip) <- c("Year", 1:12)

str(precip)



ERA5P_monthly <- aggregate(precip~Year+Month, ERA5P_filtered, sum)

ERA5P_yearly <- aggregate(precip ~ Year, ERA5P_filtered, sum)



# pivot_wider
ERA5P_wide <- ERA5P_monthly %>%
  pivot_wider(
    names_from = Month,
    values_from = precip,
    names_glue = "{month.abb[as.integer(Month)]}_precip",
    values_fill = 0
  )

head(ERA5P_wide)
# Ensure all values are numeric
precip <- ERA5P_wide %>%
  mutate_if(is.character, as.numeric)

precip



#Long format
precip_long <- precip %>%
  pivot_longer(
    cols = ends_with("precip"),
    names_to = "Month_raw",
    values_to = "precip"
  ) %>%
  mutate(
    # Extract 3-letter month abbreviation
    Month = substr(Month_raw, 1, 3),
    # Convert to ordered factor
    Month = factor(Month, levels = month.abb, ordered = TRUE)
  )

#monthly precip (1950–2023)
climatology <- precip_long %>%
  filter(Year >= 1950, Year <= 2023) %>%
  group_by(Month) %>%
  summarise(mean_precip = mean(precip, na.rm = TRUE), .groups = "drop")



##Precipitation plots

# Make sure figures folder exists (safe on fresh clone)
if (!dir.exists(figures_dir)) dir.create(figures_dir)

# Plot + save function (one year per plot)
plot_precip_highlight_year <- function(precip_long, climatology,
                                       highlight_year,
                                       outfile = NULL,
                                       width = 6,
                                       height = 4,
                                       dpi = 300) {
  
  yearly_data <- precip_long %>% filter(Year == highlight_year)
  if (nrow(yearly_data) == 0) stop("No precipitation data for that year: ", highlight_year)
  
  p <- ggplot() +
    geom_col(data = yearly_data,
             aes(x = Month, y = precip),
             fill = "skyblue", alpha = 0.7) +
    geom_line(data = climatology,
              aes(x = Month, y = mean_precip, group = 1),
              color = "red", linewidth = 1) +
    geom_point(data = climatology,
               aes(x = Month, y = mean_precip),
               color = "red", size = 2) +
    labs(
      title = highlight_year,
      subtitle = "red line: 1950–2023 average",
      y = "precipitation (mm)",
      x = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.subtitle = element_text(size = 12)
    )
  
  print(p)
  
  if (!is.null(outfile)) {
    ggsave(outfile, p, width = width, height = height, dpi = dpi)
  }
  
  invisible(p)
}

# Choose highlight years 
highlight_years <- c(1994, 2015, 2017)  # <- edit this list as needed

for (yr in highlight_years) {
  plot_precip_highlight_year(
    precip_long = precip_long,
    climatology = climatology,
    highlight_year = yr,
    outfile = file.path(figures_dir, paste0("precip_highlightyear", yr, "_plot.png")),
    width = 6,
    height = 4,
    dpi = 300
  )
}


