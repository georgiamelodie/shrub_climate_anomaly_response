##############################################################################
# Alnus and Salix shrub anatomical response to climate
# ISR Alnus xylem anomalies - EWFR GLMMs

library(here)

source(here("R/utils_packages.R"))
load_project_packages()

source(here("R/0_paths_setup.R")) 

data_dir    <- here("data")
figures_dir <- here("figures")
output_dir  <- here("output")

if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
if (!dir.exists(output_dir))  dir.create(output_dir, recursive = TRUE)


library(utils)
library(dplR)
library(tidyr)
library(readr)
library(ggplot2)
library(corrplot)
library(lme4)
library(MuMIn)
library(MASS)
library(dplyr)
library(lubridate)
library(stringr)
library(broom.mixed)
library(tibble)

#load wood anatomy anomaly data
raw <- read.csv(
  file.path(data_dir, "alnusBRFRdata.csv"),
  check.names = FALSE,
  na.strings = c("", "NA")
)

names(raw)

# Use raw as "rings" (years as rows, samples as columns)
rings <- raw


#load sample latitudes
lat_df <- read.csv(
  file.path(data_dir, "ISR_alnus_subset_samples.csv")
)


# #create long dataframe
df_long <- data.frame(
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

for (sample in colnames(rings)[-1]) {
  for (i in 1:nrow(rings)) {
    
    year <- rings$year[i]
    value <- rings[i, sample]
    
    # presence: TRUE only if the original entry is genuinely present
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
    
    df_long <- rbind(df_long, data.frame(
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

df_long <- df_long %>%
  arrange(SampleID, Year) %>%
  group_by(SampleID) %>%
  mutate(
    Age = ifelse(present, cumsum(present), NA_integer_)
  ) %>%
  ungroup()


df_long <- df_long %>% filter(present)


# Merge latitude info
df_long <- df_long %>%
  left_join(lat_df, by = c("SampleID" = "ID"))

#anomaly counts
df_long <- df_long %>%
  mutate(
    tBLW = as.integer((BLW + pBLW) > 0),
    tBEW = as.integer((BEW + pBEW) > 0),
    tFLW = as.integer((FLW + pFLW) > 0),
    tFEW = as.integer((FEW + pFEW) > 0)
  )

# Filter years 
df_long <- df_long %>%
  filter(Year >= 1950)







####climate 
#function for climate data parametetrs 
prepare_monthly_climate <- function(filepath,
                                    varname,
                                    value_name = "value",
                                    calc_GDD = FALSE,
                                    calc_CDD = FALSE,
                                    precip = FALSE) {
  #read file
  dat <- read.table(filepath, header = FALSE, skip = 24)
  colnames(dat) <- c("date", value_name)
  dat$date <- as.Date(as.character(dat$date), format = "%Y%m%d")
  
  # add year, month
  dat <- dat %>%
    mutate(year = year(date),
           month = month(date))
  
  #compute summaries
  if (precip) {
    dat_month <- dat %>%
      group_by(year, month) %>%
      summarise(!!value_name := sum(.data[[value_name]], na.rm = TRUE),
                .groups = "drop")
  } else {
    dat_month <- dat %>%
      group_by(year, month) %>%
      summarise(
        mean_value = mean(.data[[value_name]], na.rm = TRUE),
        GDD = if (calc_GDD) sum(pmax(.data[[value_name]] - 5, 0), na.rm = TRUE) else NA,
        CDD = if (calc_CDD) sum(pmax(5 - .data[[value_name]], 0), na.rm = TRUE) else NA,
        .groups = "drop"
      )
  }
  
  #pivot to wide format
  if (precip) {
    dat_wide <- dat_month %>%
      pivot_wider(
        id_cols = year,
        names_from = month,
        values_from = !!sym(value_name),
        names_glue = paste0(varname, "_month{month}")
      )
  } else {
    dat_wide <- dat_month %>%
      pivot_wider(
        id_cols = year,
        names_from = month,
        values_from = c(mean_value, GDD, CDD),
        names_glue = paste0(varname, "_{.value}_month{month}")
      )
  }
  
  #add lagged version (previous year)
  dat_lag <- dat_wide %>%
    mutate(year = year + 1) %>%
    rename_with(~ paste0(., "_prev"), -year)
  
  #combine current and lagged
  dat_combined <- full_join(dat_wide, dat_lag, by = "year")
  return(dat_combined)
}

#load climate data 
clim_t2m <- prepare_monthly_climate(
  filepath = file.path(data_dir,"ERA5_t2m_daily.dat"),
  varname  = "t2m",
  value_name = "t2m",
  calc_GDD = TRUE,
  calc_CDD = TRUE
)

clim_tmin <- prepare_monthly_climate(
  filepath = file.path(data_dir,"ERA5_tmin_daily.dat"),
  varname  = "tmin",
  value_name = "tmin",
  calc_GDD = TRUE,
  calc_CDD = TRUE
)

clim_precip <- prepare_monthly_climate(
  filepath = file.path(data_dir,"ERA5_prcp_daily.dat"),
  varname  = "precip",
  value_name = "precip_mm",
  precip = TRUE
)

# Merge all climate data
clim_list <- list(clim_t2m, clim_tmin, clim_precip)
climate_all <- Reduce(function(x, y) merge(x, y, by = "year", all = TRUE), clim_list)


# merge all data
df_full <- df_long %>%
  left_join(climate_all, by = c("Year" = "year"))

df_full <- df_full %>%
  mutate(lat = as.numeric(lat),
         long = as.numeric(long))



#scale variables incl ring age
df_full_scaled <- df_full %>%
  mutate(
    across(
      .cols = c(lat, Age, matches("^(t2m_|tmin_|.*GDD.*|.*CDD.*|precip_month)")),
      .fns = ~ {
        s <- sd(.x, na.rm = TRUE)
        
        # if all NA or constant, return value
        if (is.na(s) || s == 0) return(.x)
        
        sc <- scale(.x)
        out <- as.numeric(sc)
        
        attr(out, "scaled:center") <- attr(sc, "scaled:center")
        attr(out, "scaled:scale")  <- attr(sc, "scaled:scale")
        out
      },
      .names = "{.col}_s"
    )
  )


df_mod <- df_full_scaled %>%
  mutate(SampleID = factor(SampleID))


nrow(df_full_scaled)
nrow(df_mod)
attr(df_mod$tmin_CDD_month8_s, "scaled:center")









###GLMMs
#hypothesis-driven model testing 

###test CDD versus GDD
thermal_vars <- c(
  "tFEW", "SampleID", "lat_s", "Age_s",
  # CDD (current year)
  "tmin_CDD_month5_s", "tmin_CDD_month6_s",
  "tmin_CDD_month7_s", "tmin_CDD_month8_s",
  # GDD (current year)
  "t2m_GDD_month5_s", "t2m_GDD_month6_s",
  "t2m_GDD_month7_s", "t2m_GDD_month8_s"
)

df_thermal <- df_mod %>%
  dplyr::select(dplyr::any_of(thermal_vars)) %>%
  tidyr::drop_na()


###models

#baseline
m0_null <- glmer(
  tFEW ~ lat_s + Age_s + (1 | SampleID),
  data = df_thermal, family = binomial,
  control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e5))
)

#mCDD
m_CDD <- glmer(
  tFEW ~ lat_s + Age_s +
    tmin_CDD_month5_s + tmin_CDD_month6_s +
    tmin_CDD_month7_s + tmin_CDD_month8_s +
    (1 | SampleID),
  data = df_thermal, family = binomial,
  control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e5))
)

#mGDD
m_GDD <- glmer(
  tFEW ~ lat_s + Age_s +
    t2m_GDD_month5_s + t2m_GDD_month6_s +
    t2m_GDD_month7_s + t2m_GDD_month8_s +
    (1 | SampleID),
  data = df_thermal, family = binomial,
  control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e5))
)


sel_thermal <- model.sel(
  list(null = m0_null, CDD = m_CDD, GDD = m_GDD),
  rank = "AICc"
)
print(sel_thermal)

supported_thermal <- get.models(sel_thermal, subset = delta <= 2)
names(supported_thermal)
#"CDD"  
###cold stress / CDD dominant thermal mechanism




##test CDD  and precip

##EWFR: early growing season climate inform test models

vars_current <- c(
  "tFEW","SampleID","lat_s","Age_s",
  # CDD
  "tmin_CDD_month5_s","tmin_CDD_month6_s","tmin_CDD_month7_s",
  # precip
  "precip_month5_s","precip_month6_s"
)

df_current <- df_mod %>%
  dplyr::select(dplyr::any_of(vars_current)) %>%
  tidyr::drop_na()


# glmm control
glmm_ctrl <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 5e5)
)

# baseline model
m0_null <- glmer(
  tFEW ~ lat_s + Age_s + (1 | SampleID),
  data = df_current,
  family = binomial,
  control = glmm_ctrl
)

# M1: Early growing season cold stress (May-June)
m1_earlyCDD <- glmer(
  tFEW ~ lat_s + Age_s +
    tmin_CDD_month5_s + tmin_CDD_month6_s +
    (1 | SampleID),
  data = df_current,
  family = binomial,
  control = glmm_ctrl
)

# M2: Mid growing season cold events (July)
m2_julyCDD <- glmer(
  tFEW ~ lat_s + Age_s +
    tmin_CDD_month7_s +
    (1 | SampleID),
  data = df_current,
  family = binomial,
  control = glmm_ctrl
)

# M3: Early growing season precipitation (May-June)
m3_earlyPrecip <- glmer(
  tFEW ~ lat_s + Age_s +
    precip_month5_s + precip_month6_s +
    (1 | SampleID),
  data = df_current,
  family = binomial,
  control = glmm_ctrl
)

# M4: Combined early growing season CDD and precip (May-June)
m4_earlyCDDprecip <- glmer(
  tFEW ~ lat_s + Age_s +
    tmin_CDD_month5_s + tmin_CDD_month6_s +
    precip_month5_s + precip_month6_s +
    (1 | SampleID),
  data = df_current,
  family = binomial,
  control = glmm_ctrl
)


library(MuMIn)
options(na.action = "na.fail")

cand_set_current <- list(
  null            = m0_null,
  early_CDD       = m1_earlyCDD,
  july_CDD        = m2_julyCDD,
  early_precip    = m3_earlyPrecip,
  earlyCDD_precip = m4_earlyCDDprecip
)

sel_current <- model.sel(cand_set_current, rank = "AICc")
print(sel_current)

supported_current <- get.models(sel_current, subset = delta <= 2)
names(supported_current)



#"early_CDD"       "earlyCDD_precip"

best_current <- get.models(sel_current, 1)[[1]]
best_current


# R2
r.squaredGLMM(m1_earlyCDD)

###check overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = ratio, rdf = rdf, p = pval)
}

overdisp_fun(best_current)

###to check the sign moderate overdispersion, 
library(DHARMa)
res <- simulateResiduals(best_current, n = 1000)
plot(res)
testDispersion(res)
testZeroInflation(res)
#no significant overdispersion from simulation residual diagnostics

#tidy format effect summaries
tidy_best <- broom.mixed::tidy(best_current, effects = "fixed", conf.int = TRUE)
print(tidy_best)





####coeff table with ORs
best_coef <- coef(summary(best_current)) %>%
  as.data.frame() %>%
  rownames_to_column("term")

# Optional if using averaged models
if (!"Estimate" %in% names(best_coef)) {
  est_col <- setdiff(names(best_coef), "term")[1]
  best_coef <- best_coef %>% dplyr::rename(Estimate = all_of(est_col))
}

# Confidence intervals for fixed effects 
best_ci <- confint(best_current, parm = "beta_", method = "Wald") %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  dplyr::rename(conf.low = `2.5 %`, conf.high = `97.5 %`)

# Combine 
tidy_best <- best_coef %>%
  dplyr::left_join(best_ci, by = "term") %>%
  dplyr::rename(
    estimate  = Estimate,
    std.error = `Std. Error`,
    z.value   = `z value`,
    p.value   = `Pr(>|z|)`
  ) %>%
  dplyr::mutate(
    odds_ratio   = exp(estimate),
    conf.low.OR  = exp(conf.low),
    conf.high.OR = exp(conf.high)
  ) %>%
  dplyr::select(
    term, estimate, std.error, z.value, p.value,
    conf.low, conf.high,
    odds_ratio, conf.low.OR, conf.high.OR
  )


tidy_best
write.csv(
  tidy_best,
  file.path(output_dir, "Alnus_EWFR_model_current_coeffs.csv"),
  row.names = FALSE
)




######Testing preceding year impact

###### EWFR: testing preceding-year (legacy) climate effects
##current-year EWFR best models predictors: May and June CDD, or May and  June CDD and May and June precip

library(dplyr)
library(tidyr)
library(lme4)
library(MuMIn)

# variables
vars_legacy <- c(
  "tFEW","SampleID","lat_s","Age_s",
  # current year from first model set
  "tmin_CDD_month5_s","tmin_CDD_month6_s",
  "precip_month5_s","precip_month6_s",
  # preceding year 
  "tmin_CDD_month8_prev_s","tmin_CDD_month9_prev_s",
  "precip_month8_prev_s","precip_month9_prev_s"
)

df_legacy <- df_mod %>%
  dplyr::select(dplyr::any_of(vars_legacy)) %>%
  tidyr::drop_na()

# glmm control 
glmm_ctrl <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 5e5)
)

#baseline 
m0_null_legacy <- glmer(
  tFEW ~ lat_s + Age_s + (1 | SampleID),
  data = df_legacy,
  family = binomial,
  control = glmm_ctrl
)

#current year CDD
m_current_earlyCDD <- glmer(
  tFEW ~ lat_s + Age_s +
    tmin_CDD_month5_s + tmin_CDD_month6_s +
    (1 | SampleID),
  data = df_legacy,
  family = binomial,
  control = glmm_ctrl
)



#preceding year

#P1: preceding year Aug cDD
m_prevAugCDD <- glmer(
  tFEW ~ lat_s + Age_s +
    tmin_CDD_month5_s + tmin_CDD_month6_s +
    tmin_CDD_month8_prev_s +
    (1 | SampleID),
  data = df_legacy,
  family = binomial,
  control = glmm_ctrl
)

#P2: preceding year Sep precip
m_prevSepPrecip <- glmer(
  tFEW ~ lat_s + Age_s +
    tmin_CDD_month5_s + tmin_CDD_month6_s +
    precip_month9_prev_s +
    (1 | SampleID),
  data = df_legacy,
  family = binomial,
  control = glmm_ctrl
)

#P3: add both preceding
m_prevAugCDD_prevSepPrecip <- glmer(
  tFEW ~ lat_s + Age_s +
    tmin_CDD_month5_s + tmin_CDD_month6_s +
    tmin_CDD_month8_prev_s + precip_month9_prev_s +
    (1 | SampleID),
  data = df_legacy,
  family = binomial,
  control = glmm_ctrl
)



#model selection
options(na.action = "na.fail")

legacy_set <- list(
  null                       = m0_null_legacy,
  current_earlyCDD           = m_current_earlyCDD,
  prevAugCDD                 = m_prevAugCDD,
  prevSepPrecip              = m_prevSepPrecip,
  prevAugCDD_prevSepPrecip   = m_prevAugCDD_prevSepPrecip
)

sel_legacy <- model.sel(legacy_set, rank = "AICc")
print(sel_legacy)

supported_legacy <- get.models(sel_legacy, subset = delta <= 2)
names(supported_legacy)
# "prevAugCDD"      "prevAugCDD_prevSepPrecip"
 

summary(m_prevAugCDD)
summary(m_prevAugCDD_prevSepPrecip)




nrow(df_current)
nrow(df_legacy)

# event rate
sum(df_current$tFEW == 1); mean(df_current$tFEW == 1)
sum(df_legacy$tFEW == 1);  mean(df_legacy$tFEW == 1)

# how many rows overlap?
nrow(dplyr::semi_join(df_current, df_legacy, by = c("SampleID","tFEW","lat_s","Age_s")))




best_legacy <- get.models(sel_legacy, 1)[[1]]
best_legacy
#prevAugCDD


# R2
r.squaredGLMM(m_prevAugCDD)

###check overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = ratio, rdf = rdf, p = pval)
}

overdisp_fun(best_legacy)

###to check the sign moderate overdispersion, 
library(DHARMa)
res <- simulateResiduals(best_legacy, n = 1000)
plot(res)
testDispersion(res)
testZeroInflation(res)
#no significant overdispersion from simulation residual diagnostics

#tidy format effect summaries
tidy_best <- broom.mixed::tidy(best_legacy, effects = "fixed", conf.int = TRUE)
print(tidy_best)






####coeff table with ORs
best_coef <- coef(summary(best_legacy)) %>%
  as.data.frame() %>%
  rownames_to_column("term")

# Optional if using averaged models
if (!"Estimate" %in% names(best_coef)) {
  est_col <- setdiff(names(best_coef), "term")[1]
  best_coef <- best_coef %>% dplyr::rename(Estimate = all_of(est_col))
}

# Confidence intervals for fixed effects 
best_ci <- confint(best_legacy, parm = "beta_", method = "Wald") %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  dplyr::rename(conf.low = `2.5 %`, conf.high = `97.5 %`)

# Combine 
tidy_legacy <- best_coef %>%
  dplyr::left_join(best_ci, by = "term") %>%
  dplyr::rename(
    estimate  = Estimate,
    std.error = `Std. Error`,
    z.value   = `z value`,
    p.value   = `Pr(>|z|)`
  ) %>%
  dplyr::mutate(
    odds_ratio   = exp(estimate),
    conf.low.OR  = exp(conf.low),
    conf.high.OR = exp(conf.high)
  ) %>%
  dplyr::select(
    term, estimate, std.error, z.value, p.value,
    conf.low, conf.high,
    odds_ratio, conf.low.OR, conf.high.OR
  )


tidy_legacy
write.csv(
  tidy_legacy,
  file.path(output_dir, "Alnus_EWFR_prevyear_model_coeffs.csv"),
  row.names = FALSE
)





