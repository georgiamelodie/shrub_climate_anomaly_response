##############################################################################
# Alnus and Salix shrub anatomical response to climate
# ISR xylem anomalies - LWBR thermal threshold comparison between Alnus and Salix:
# Paired year bootstrap of difference in AIC-determined thermal threshold T*
library(here)

source(here("R/utils_packages.R"))
load_project_packages()

source(here("R/00_paths_setup.R"))

# project data folders
data_dir        <- here("data")
anomalies_dir   <- here("data", "anomalies")
rwi_dir         <- here("data", "rwi")
meta_dir        <- here("data", "sample_metadata")
crosswalk_dir   <- here("data", "id_crosswalk")
climate_dir     <- here("data", "climate_raw")
derived_dir     <- here("data", "derived")

if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
if (!dir.exists(output_dir))  dir.create(output_dir, recursive = TRUE)



#load wood anatomy anomaly data

raw_alnus <- read.csv(
  file.path(anomalies_dir, "alnusBRFRdata.csv"),
  check.names = FALSE,
  na.strings = c("", "NA")
)

raw_salix <- read.csv(
  file.path(anomalies_dir, "salixBRFRdata.csv"),
  check.names = FALSE,
  na.strings = c("", "NA")
)

rings_alnus <- raw_alnus
rings_salix <- raw_salix



# W326 has degraded inner rings and is excluded from Salix LWBR modelling
if ("W326" %in% colnames(rings_salix)) {
  rings_salix <- rings_salix %>%
    dplyr::select(-W326)
}

stopifnot(!"W326" %in% colnames(rings_salix))


#Yearly LWBR counts only

make_annual_LWBR <- function(rings) {
  
  rings %>%
    mutate(
      across(-year, as.character)
    ) %>%
    tidyr::pivot_longer(
      cols = -year,
      names_to = "SampleID",
      values_to = "raw_value"
    ) %>%
    filter(
      !is.na(raw_value),
      trimws(raw_value) != ""
    ) %>%
    mutate(
      raw_value = trimws(raw_value),
      tBLW = as.integer(
        grepl("\\b(BLW|pBLW)\\b", raw_value)
      )
    ) %>%
    filter(year >= 1950) %>%
    group_by(Year = year) %>%
    summarise(
      k = sum(tBLW),
      n = n(),
      .groups = "drop"
    )
}


yr_alnus <- make_annual_LWBR(rings_alnus)
yr_salix <- make_annual_LWBR(rings_salix)

summary(yr_alnus$n)
summary(yr_salix$n)

range(yr_alnus$Year)
range(yr_salix$Year)





#read daily Tmin
tmin_path <- file.path(climate_dir, "ERA5_tmin_daily.dat")
if (!file.exists(tmin_path)) stop("Missing file: ", tmin_path)

tmin_daily <- read.table(tmin_path, header = FALSE, skip = 24)
colnames(tmin_daily) <- c("date", "tmin")

tmin_daily <- tmin_daily %>%
  mutate(
    date  = as.Date(as.character(date), format = "%Y%m%d"),
    Year  = year(date),
    Month = month(date)
  )

#function: August cold exposure below threshold T
# E_aug(T) = Σ max(T − Tmin, 0) over August days (degree days)
make_E_aug <- function(T, dat = tmin_daily) {
  dat %>%
    filter(Month == 8) %>%
    group_by(Year) %>%
    summarise(
      E_aug = sum(pmax(T - tmin, 0), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(T = T)
}



###find T range for august
aug_tmin <- tmin_daily %>% filter(Month == 8)

Tmin_obs <- min(aug_tmin$tmin, na.rm = TRUE)
Tmax_obs <- max(aug_tmin$tmin, na.rm = TRUE)

c(Tmin_obs, Tmax_obs)


# Extend observed range by 1 degree C for candidate T* thresholds
T_low  <- floor(Tmin_obs) - 1
T_high <- ceiling(Tmax_obs) + 1

T_low
T_high


Ts <- seq(T_low, T_high, by = 0.5)  
E_aug_grid <- bind_rows(lapply(Ts, make_E_aug, dat = tmin_daily))

# E_aug_grid has: Year, E_aug, T
dplyr::glimpse(E_aug_grid)


# Merge yearly LWBR counts with threshold temps

alnus_thresh <- yr_alnus %>%
  left_join(
    E_aug_grid,
    by = "Year"
  )

salix_thresh <- yr_salix %>%
  left_join(
    E_aug_grid,
    by = "Year"
  )


##models to find T*

fit_one_T <- function(Tval, dat) {
  
  d <- dat %>%
    filter(T == Tval)
  
  m <- glm(
    cbind(k, n - k) ~ E_aug,
    family = binomial,
    data = d
  )
  
  tibble(
    T = Tval,
    AIC = AIC(m),
    beta = coef(m)[["E_aug"]]
  )
}


# Create T theshold (T*) for each taxon

AIC_alnus <- bind_rows(
  lapply(
    sort(unique(alnus_thresh$T)),
    fit_one_T,
    dat = alnus_thresh
  )
) %>%
  mutate(
    deltaAIC = AIC - min(AIC, na.rm = TRUE)
  )

AIC_salix <- bind_rows(
  lapply(
    sort(unique(salix_thresh$T)),
    fit_one_T,
    dat = salix_thresh
  )
) %>%
  mutate(
    deltaAIC = AIC - min(AIC, na.rm = TRUE)
  )


Tstar_alnus <- AIC_alnus$T[
  which.min(AIC_alnus$AIC)
]

Tstar_salix <- AIC_salix$T[
  which.min(AIC_salix$AIC)
]


Tstar_alnus
Tstar_salix


support_alnus <- AIC_alnus %>%
  filter(deltaAIC <= 2)

support_salix <- AIC_salix %>%
  filter(deltaAIC <= 2)

support_alnus
support_salix




# ensure common comparison of years
# common grid for taxon comparison.
# grid covers combined delta-AIC <= 2 support ranges, +/-  1 degree C.

T_compare_low <- min(
  support_alnus$T,
  support_salix$T
) - 1

T_compare_high <- max(
  support_alnus$T,
  support_salix$T
) + 1

Ts_compare <- seq(
  T_compare_low,
  T_compare_high,
  by = 0.5
)

Ts_compare <- intersect(
  Ts_compare,
  Ts
)

Ts_compare


#bootstrapping of the taxon difference in T* with paired years.
# calendar years are sampled ONCE per bootstrap replicate, with same year then used for both taxa due to shared regional climate history.

common_years <- intersect(
  sort(unique(yr_alnus$Year)),
  sort(unique(yr_salix$Year))
)

length(common_years)
range(common_years)


# restrict both taxa to the same calendar-year period
alnus_compare <- alnus_thresh %>%
  filter(
    Year %in% common_years,
    T %in% Ts_compare
  )

salix_compare <- salix_thresh %>%
  filter(
    Year %in% common_years,
    T %in% Ts_compare
  )


# T* for one dataset
get_Tstar <- function(dat, Ts_use) {
  
  aic_tbl <- bind_rows(
    lapply(
      Ts_use,
      function(Tval) {
        
        d <- dat %>%
          filter(T == Tval)
        
        m <- tryCatch(
          glm(
            cbind(k, n - k) ~ E_aug,
            family = binomial,
            data = d
          ),
          error = function(e) NULL
        )
        
        if (is.null(m)) {
          return(
            tibble(
              T = Tval,
              AIC = NA_real_
            )
          )
        }
        
        tibble(
          T = Tval,
          AIC = AIC(m)
        )
      }
    )
  ) %>%
    filter(is.finite(AIC))
  
  if (nrow(aic_tbl) == 0) {
    return(NA_real_)
  }
  
  aic_tbl$T[
    which.min(aic_tbl$AIC)
  ]
}



# Observed T* using the same common years and threshold grid
# used in the paired bootstrap
Tstar_alnus_compare <- get_Tstar(
  alnus_compare,
  Ts_compare
)

Tstar_salix_compare <- get_Tstar(
  salix_compare,
  Ts_compare
)

Delta_Tstar_observed <-
  Tstar_salix_compare - Tstar_alnus_compare

Tstar_alnus_compare
Tstar_salix_compare
Delta_Tstar_observed



# resample a sequence of years 
resample_years <- function(dat, boot_yrs) {
  
  bind_rows(
    lapply(
      boot_yrs,
      function(y) {
        dat %>%
          filter(Year == y)
      }
    )
  )
}


set.seed(123)

B <- 1000

boot_results <- data.frame(
  replicate = seq_len(B),
  Tstar_Alnus = NA_real_,
  Tstar_Salix = NA_real_,
  Delta_Tstar = NA_real_
)


for (b in seq_len(B)) {
  
  # sample calendar years once
  boot_yrs <- sample(
    common_years,
    length(common_years),
    replace = TRUE
  )
  
  # same resampled years passed to both taxa
  boot_alnus <- resample_years(
    alnus_compare,
    boot_yrs
  )
  
  boot_salix <- resample_years(
    salix_compare,
    boot_yrs
  )
  
  T_alnus_b <- get_Tstar(
    boot_alnus,
    Ts_compare
  )
  
  T_salix_b <- get_Tstar(
    boot_salix,
    Ts_compare
  )
  
  boot_results$Tstar_Alnus[b] <- T_alnus_b
  boot_results$Tstar_Salix[b] <- T_salix_b
  
  boot_results$Delta_Tstar[b] <-
    T_salix_b - T_alnus_b
}





#summary

boot_results_clean <- boot_results %>%
  filter(
    !is.na(Tstar_Alnus),
    !is.na(Tstar_Salix),
    !is.na(Delta_Tstar)
  )


nrow(boot_results_clean)

delta_quantiles <- quantile(
  boot_results_clean$Delta_Tstar,
  probs = c(
    0.025,
    0.10,
    0.50,
    0.90,
    0.975
  ),
  type = 1,
  na.rm = TRUE
)

delta_quantiles


prop_salix_higher <- mean(
  boot_results_clean$Delta_Tstar > 0,
  na.rm = TRUE
)

prop_equal <- mean(
  boot_results_clean$Delta_Tstar == 0,
  na.rm = TRUE
)

prop_salix_lower <- mean(
  boot_results_clean$Delta_Tstar < 0,
  na.rm = TRUE
)


paired_boot_summary <- data.frame(
  observed_Tstar_Alnus =
    Tstar_alnus_compare,
  
  observed_Tstar_Salix =
    Tstar_salix_compare,
  
  observed_difference =
    Delta_Tstar_observed,
  
  boot_median_difference =
    unname(delta_quantiles["50%"]),
  
  boot_95_low =
    unname(delta_quantiles["2.5%"]),
  
  boot_95_high =
    unname(delta_quantiles["97.5%"]),
  
  boot_80_low =
    unname(delta_quantiles["10%"]),
  
  boot_80_high =
    unname(delta_quantiles["90%"]),
  
  proportion_Salix_higher =
    prop_salix_higher,
  
  proportion_equal =
    prop_equal,
  
  proportion_Salix_lower =
    prop_salix_lower
)

paired_boot_summary


write.csv(
  boot_results_clean,
  file.path(
    output_dir,
    "LWBR_taxon_paired_bootstrap_Tstar.csv"
  ),
  row.names = FALSE
)

write.csv(
  paired_boot_summary,
  file.path(
    output_dir,
    "LWBR_taxon_paired_bootstrap_summary.csv"
  ),
  row.names = FALSE
)


# plot spread of between-taxa T threshold deltaT*.

delta_plot <- ggplot(
  boot_results_clean,
  aes(x = Delta_Tstar)
) +
  geom_histogram(
    binwidth = 0.5,
    boundary = -0.25,
    closed = "left",
    fill = "grey70",
    colour = "black"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = Delta_Tstar_observed,
    linewidth = 1
  ) +
  labs(
    x = expression(
      Delta*T^"*" == T[Salix]^"*" - T[Alnus]^"*"~degree*C
    ),
    y = "Bootstrap frequency"
  ) +
  theme_classic()

delta_plot


ggsave(
  file.path(
    figures_dir,
    "LWBR_taxon_paired_bootstrap_deltaTstar.png"
  ),
  delta_plot,
  width = 7,
  height = 5,
  dpi = 300
)

paired_boot_summary

