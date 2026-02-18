##############################################################################
# Alnus and Salix shrub anatomical response to climate
# ISR Salix xylem anomalies - LWBR Temp threshold
library(here)

source(here("R/utils_packages.R"))
load_project_packages()

source(here("R/00_paths_setup.R")) 

data_dir    <- here("data")
figures_dir <- here("figures")
output_dir  <- here("output")

if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
if (!dir.exists(output_dir))  dir.create(output_dir, recursive = TRUE)


#load wood anatomy anomaly data
raw_salix <- read.csv(
  file.path(data_dir, "salixBRFRdata.csv"),
  check.names = FALSE,
  na.strings = c("", "NA")
)

names(raw_salix)
str(raw_salix)

# Use raw as "rings" (years as rows, samples as columns)
rings_salix <- raw_salix


#load sample latitudes
lat_df_salix <- read.csv(
  file.path(data_dir, "ISR_salix_subset_samples.csv")
)

# #create long dataframe
df_long_salix <- data.frame(
  Year = integer(),
  SampleID = character(),
  raw_value = character(),
  present = logical(),
  FLW = integer(), pFLW = integer(),
  FEW = integer(), pFEW = integer(),
  BLW = integer(), pBLW = integer(),
  BEW = integer(), pBEW = integer(),
  stringsAsFactors = FALSE
)

for (sample in colnames(rings_salix)[-1]) {
  for (i in 1:nrow(rings_salix)) {
    
    year <- rings_salix$year[i]
    value <- rings_salix[i, sample]
    
    # presence if entry exists
    pres <- !(is.na(value) || trimws(as.character(value)) == "")
    
    value_chr <- if (pres) trimws(as.character(value)) else NA_character_
    
    FLW  <- as.integer(pres && grepl("\\bFLW\\b",  value_chr))
    pFLW <- as.integer(pres && grepl("\\bpFLW\\b", value_chr))
    FEW  <- as.integer(pres && grepl("\\bFEW\\b",  value_chr))
    pFEW <- as.integer(pres && grepl("\\bpFEW\\b", value_chr))
    BLW  <- as.integer(pres && grepl("\\bBLW\\b",  value_chr))
    pBLW <- as.integer(pres && grepl("\\bpBLW\\b", value_chr))
    BEW  <- as.integer(pres && grepl("\\bBEW\\b",  value_chr))
    pBEW <- as.integer(pres && grepl("\\bpBEW\\b", value_chr))
    
    df_long_salix <- rbind(df_long_salix, data.frame(
      Year = year,
      SampleID = sample,
      raw_value = value_chr,
      present = pres,
      FLW = FLW, pFLW = pFLW,
      FEW = FEW, pFEW = pFEW,
      BLW = BLW, pBLW = pBLW,
      BEW = BEW, pBEW = pBEW
    ))
  }
}



library(dplyr)

df_long_salix <- df_long_salix %>%
  arrange(SampleID, Year) %>%
  group_by(SampleID) %>%
  mutate(
    Age = ifelse(present, cumsum(present), NA_integer_)
  ) %>%
  ungroup()


df_long_salix <- df_long_salix %>% filter(present)


# merge latitude info
df_long_salix <- df_long_salix %>%
  left_join(lat_df_salix, by = c("SampleID" = "ID"))

#anomaly counts
df_long_salix <- df_long_salix %>%
  mutate(
    tBLW = as.integer((BLW + pBLW) > 0),
    tBEW = as.integer((BEW + pBEW) > 0),
    tFLW = as.integer((FLW + pFLW) > 0),
    tFEW = as.integer((FEW + pFEW) > 0)
  )

# filter years 
df_long_salix <- df_long_salix %>%
  filter(Year >= 1950)



# Number of shrubs contributing per year
n_by_year <- df_long_salix %>%
  filter(Year >= 1950) %>%
  count(Year, name = "n")

summary(n_by_year$n)
plot(n_by_year$Year, n_by_year$n, type = "b",
     ylab = "Number of shrubs with rings",
     xlab = "Year")



yr_dat <- df_long_salix %>%
  filter(Year >= 1950) %>%
  group_by(Year) %>%
  summarise(
    k = sum(tBLW, na.rm=TRUE),
    n = n(),
    .groups="drop"
  )



#read daily Tmin
tmin_path <- file.path(data_dir, "ERA5_tmin_daily.dat")
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


#Tmin and max plus and minus 1 for brackets for T* thresholds
T_low  <- floor(min(aug_tmin$tmin, na.rm=TRUE)) - 1
T_high <- ceiling(max(aug_tmin$tmin, na.rm=TRUE)) + 1

T_low
T_high


Ts <- seq(T_low, T_high, by = 0.5)  
E_aug_grid <- bind_rows(lapply(Ts, make_E_aug, dat = tmin_daily))

# E_aug_grid has: Year, E_aug, T
dplyr::glimpse(E_aug_grid)


#combine with E_aug for all candidate thresholds
yr_thresh <- yr_dat %>%
  left_join(E_aug_grid, by = c("Year" = "Year"))

head(yr_thresh)



fit_one_T <- function(Tval, dat){
  d <- dat %>% filter(T == Tval)
  
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

AIC_tbl <- bind_rows(
  lapply(unique(yr_thresh$T), fit_one_T, dat = yr_thresh)
)

# Inspect
head(AIC_tbl)



###temperature threshold without lowest AIC - Tstar
T_star <- AIC_tbl$T[which.min(AIC_tbl$AIC)]
T_star


AIC_tbl <- AIC_tbl %>%
  mutate(deltaAIC = AIC - min(AIC, na.rm = TRUE))

# Candidate support range
subset(AIC_tbl, deltaAIC <= 2)



#bootstrapping
set.seed(123)

boot_Tstar <- function(dat){
  yrs <- sort(unique(dat$Year))
  boot_yrs <- sample(yrs, length(yrs), replace = TRUE)
  
  boot_dat <- dplyr::bind_rows(lapply(boot_yrs, function(y){
    dat %>% dplyr::filter(Year == y)
  }))
  
  AIC_boot <- dplyr::bind_rows(
    lapply(unique(boot_dat$T), fit_one_T, dat = boot_dat)
  )
  
  AIC_boot$T[which.min(AIC_boot$AIC)]
}



B <- 1000
T_boot <- replicate(B, boot_Tstar(yr_thresh))

quantile(T_boot, probs = c(0.025, 0.5, 0.975))



hist(T_boot, breaks = 30, col = "grey",
     main = "Bootstrap distribution of T*",
     xlab = "Threshold temperature (°C)")
abline(v = T_star, col = "red", lwd = 2)



quantile(T_boot, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
quantile(T_boot, probs = c(0.1, 0.9), na.rm = TRUE)




########constrained bootstrap to dAIC

boot_support_interval <- function(dat, Ts){
  yrs <- unique(dat$Year)
  boot_yrs <- sample(yrs, length(yrs), replace = TRUE)
  
  boot_dat <- dplyr::bind_rows(lapply(boot_yrs, function(y){
    dat %>% dplyr::filter(Year == y)
  }))
  
  aic_tbl <- dplyr::bind_rows(lapply(Ts, function(Tval){
    d <- boot_dat %>% dplyr::filter(T == Tval)
    if (nrow(d) == 0) return(tibble::tibble(T = Tval, AIC = NA_real_))
    m <- glm(cbind(k, n - k) ~ E_aug, family = binomial, data = d)
    tibble::tibble(T = Tval, AIC = AIC(m))
  })) %>%
    dplyr::filter(!is.na(AIC))
  
  if (nrow(aic_tbl) == 0) return(c(T_best = NA, T_lo = NA, T_hi = NA))
  
  aic_tbl <- aic_tbl %>%
    dplyr::mutate(deltaAIC = AIC - min(AIC, na.rm = TRUE))
  
  supported <- aic_tbl %>% dplyr::filter(deltaAIC <= 2)
  if (nrow(supported) == 0) supported <- aic_tbl
  
  c(
    T_best = aic_tbl$T[which.min(aic_tbl$AIC)],
    T_lo   = min(supported$T, na.rm = TRUE),
    T_hi   = max(supported$T, na.rm = TRUE)
  )
}



set.seed(123)
B <- 1000

#calculate constrained Ts derived from dAIC
AIC_tbl <- AIC_tbl %>%
  mutate(deltaAIC = AIC - min(AIC, na.rm = TRUE))

supported_T <- AIC_tbl %>% filter(deltaAIC <= 2)

T_sup_low  <- min(supported_T$T, na.rm = TRUE)
T_sup_high <- max(supported_T$T, na.rm = TRUE)

bracket <- 1     # degrees C 
step    <- 0.5   # to match Ts step

Ts_used <- seq(T_sup_low - bracket, T_sup_high + bracket, by = step)
Ts_used <- intersect(Ts_used, sort(unique(yr_thresh$T)))

boot_mat <- replicate(B, boot_support_interval(yr_thresh, Ts_used))
boot_df <- as.data.frame(t(boot_mat))

# Point estimate distribution
quantile(boot_df$T_best, c(0.025, 0.5, 0.975), na.rm=TRUE)

# Support-interval uncertainty 
quantile(boot_df$T_lo, c(0.1, 0.5, 0.9), na.rm=TRUE)
quantile(boot_df$T_hi, c(0.1, 0.5, 0.9), na.rm=TRUE)



#constrained bootstrap to dAIC CI central range (10-90%)
T_rangecentral <- quantile(boot_df$T_best, c(0.1,0.5, 0.9), na.rm = TRUE)
T_range_df1 <- data.frame(
  quantile = names(T_rangecentral),
  T_best   = as.numeric(T_rangecentral),
  row.names = NULL
)

write.csv(
  T_range_df1,
  file.path(output_dir, "Salix_Tstar_bootstrap_range_10_90.csv"),
  row.names = FALSE
)
T_rangecentral
###constrained bootstrap to dAIC CI range (2.5-97.5%)
T_range_full <- quantile(boot_df$T_best, c(0.025, 0.5, 0.975), na.rm = TRUE)
T_range_df2 <- data.frame(
  quantile = names(T_range_full),
  T_best   = as.numeric(T_range_full),
  row.names = NULL
)

write.csv(
  T_range_df2,
  file.path(output_dir, "Salix_Tstar_bootstrap_range_2.5_97.5.csv"),
  row.names = FALSE
)
T_range_full



library(dplyr)
library(ggplot2)

#fit model at Tstar and generate prediction plot and E50 (p=0.5)  
#E50: cumulative cumulative cold exposure below Tstar associated with a predicted LWBR probability of 0.5.
#T_star <- 4

#subset data to the chosen threshold
d_star <- yr_thresh %>% filter(T == T_star)

# fit the binomial model at T*
m_star <- glm(cbind(k, n - k) ~ E_aug, family = binomial, data = d_star)

#prediction grid over observed E_aug range
newdat <- data.frame(
  E_aug = seq(min(d_star$E_aug, na.rm = TRUE),
              max(d_star$E_aug, na.rm = TRUE),
              length.out = 200)
)

pr <- predict(m_star, newdata = newdat, type = "link", se.fit = TRUE)

newdat <- newdat %>%
  mutate(
    p_fit  = plogis(pr$fit),
    p_low  = plogis(pr$fit - 1.96 * pr$se.fit),
    p_high = plogis(pr$fit + 1.96 * pr$se.fit)
  )

#compute E_aug where predicted p = 0.5
b0 <- coef(m_star)[1]
b1 <- coef(m_star)[2]
E50 <- -b0 / b1

#####E50: cumulative cumulative cold exposure below Tstar associated with a predicted LWBR probability of 0.5.
E50
