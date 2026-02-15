##############################################################################
# Alnus and Salix shrub anatomical response to climate
# ISR Salix xylem anomalies - LWBR GLMMs

library(here)

source(here("R/utils_packages.R"))
load_project_packages()

source(here("R/00_paths_setup.R")) 
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
  file.path(data_dir, "salixBRFRdata.csv"),
  check.names = FALSE,
  na.strings = c("", "NA")
)

names(raw)
str(raw)

# Use raw as "rings" (years as rows, samples as columns)
rings <- raw

# remove curtailed sample (W326) before processing
if ("W326" %in% colnames(rings)) {
  rings <- rings %>% dplyr::select(-W326)
}
stopifnot(!"W326" %in% colnames(rings))

#load sample latitudes
lat_df <- read.csv(
  file.path(data_dir, "ISR_salix_subset_samples.csv")
)

lat_df <- lat_df %>% filter(ID != "W326")

#create long dataframe
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
    
    #presence: TRUE if ring present and measured for year
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

# Filter years for when climate data available
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
  "tBLW", "SampleID", "lat_s", "Age_s",
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
  tBLW ~ lat_s + Age_s + (1 | SampleID),
  data = df_thermal, family = binomial,
  control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e5))
)

#mCDD
m_CDD <- glmer(
  tBLW ~ lat_s + Age_s +
    tmin_CDD_month5_s + tmin_CDD_month6_s +
    tmin_CDD_month7_s + tmin_CDD_month8_s +
    (1 | SampleID),
  data = df_thermal, family = binomial,
  control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e5))
)

#mGDD
m_GDD <- glmer(
  tBLW ~ lat_s + Age_s +
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

###cold stress / CDD dominant thermal mechanism




##test CDD  and precip

vars_current <- c(
  "tBLW","SampleID","lat_s","Age_s",
  #CDD
  "tmin_CDD_month5_s","tmin_CDD_month6_s","tmin_CDD_month7_s",
  "tmin_CDD_month8_s","tmin_CDD_month9_s",
  #precip
  "precip_month6_s","precip_month7_s","precip_month8_s"
)

df_current <- df_mod %>%
  dplyr::select(dplyr::any_of(vars_current)) %>%
  tidyr::drop_na()


# glmm control
glmm_ctrl <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 5e5)
)

#baseline model
m0_null <- glmer(
  tBLW ~ lat_s + Age_s + (1 | SampleID),
  data = df_current,
  family = binomial,
  control = glmm_ctrl
)

#M1: Early growing season cold stress 
# May-June minimum-temperature cold exposure
m1_earlyCDD <- glmer(
  tBLW ~ lat_s + Age_s +
    tmin_CDD_month5_s + tmin_CDD_month6_s +
    (1 | SampleID),
  data = df_current,
  family = binomial,
  control = glmm_ctrl
)

#M2: Mid growing season cold events (Jul)
m2_midCDD <- glmer(
  tBLW ~ lat_s + Age_s +
    tmin_CDD_month7_s +
    (1 | SampleID),
  data = df_current,
  family = binomial,
  control = glmm_ctrl
)

#M3: Late growing season cold stress (Aug, Sep) 
# Tmin
m3_lateCDD <- glmer(
  tBLW ~ lat_s + Age_s +
    tmin_CDD_month8_s + tmin_CDD_month9_s +
    (1 | SampleID),
  data = df_current,
  family = binomial,
  control = glmm_ctrl
)

#M4: growing season exposure (June - Aug)
m4_totalCDD <- glmer(
  tBLW ~ lat_s + Age_s +
    tmin_CDD_month6_s + tmin_CDD_month7_s +
    tmin_CDD_month8_s + 
    (1 | SampleID),
  data = df_current,
  family = binomial,
  control = glmm_ctrl
)

#M5: precipitation constraint (Jun-Aug)
m5_precip <- glmer(
  tBLW ~ lat_s + Age_s +
    precip_month6_s + precip_month7_s + precip_month8_s +
    (1 | SampleID),
  data = df_current,
  family = binomial,
  control = glmm_ctrl
)

#M6: Combined CDD and precip
m6_CDDprecip <- glmer(
  tBLW ~ lat_s + Age_s +
    tmin_CDD_month7_s + tmin_CDD_month8_s +
    precip_month7_s + precip_month8_s +
    (1 | SampleID),
  data = df_current,
  family = binomial,
  control = glmm_ctrl
)

#M7: early and late growing season cold
m7_earlylateCDD <- glmer(
  tBLW ~ lat_s + Age_s +
    tmin_CDD_month5_s +
    tmin_CDD_month8_s + 
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
  mid_CDD         = m2_midCDD,
  late_CDD        = m3_lateCDD,
  season_CDD      = m4_totalCDD,
  precip          = m5_precip,
  lateCDD_precip  = m6_CDDprecip,
  earlyLate_CDD   = m7_earlylateCDD
)

sel_current <- model.sel(cand_set_current, rank = "AICc")
print(sel_current)

supported_current <- get.models(sel_current, subset = delta <= 2)
names(supported_current)


best_current <- get.models(sel_current, 1)[[1]]
best_current

# R2 of model
r.squaredGLMM(m7_earlylateCDD)


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
  file.path(output_dir, "Salix_LWBR_current_model_coeffs.csv"),
  row.names = FALSE
)


effect_table <- tidy_best %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(
    Predictor = dplyr::recode(term,
                              "lat_s" = "Latitude",
                              "Age_s" = "Ring age",
                              "tmin_CDD_month5_s" = "May CDD",
                              "tmin_CDD_month8_s" = "August CDD"
    ),
    `Odds ratio (95% CI)` = sprintf("%.2f (%.2f-%.2f)", odds_ratio, conf.low.OR, conf.high.OR)
  ) %>%
  dplyr::select(Predictor, `Odds ratio (95% CI)`)

effect_table

effect_table <- effect_table %>%
  mutate(`Odds ratio (95% CI)` = gsub("\\.00", "", `Odds ratio (95% CI)`))
effect_table


write.csv(
  effect_table,
  file.path(output_dir, "Salix_LWBR_current_model_OR.csv"),
  row.names = FALSE
)






######Testing preceding year impact

vars_legacy <- c(
  "tBLW","SampleID", "lat_s","Age_s", 
  "tmin_CDD_month5_s", "tmin_CDD_month8_s", 
  "tmin_CDD_month7_prev_s", "tmin_CDD_month8_prev_s", 
  "tmin_CDD_month9_prev_s" ) 

df_legacy <- df_mod %>% 
  dplyr::select(dplyr::any_of(vars_legacy)) %>% 
  tidyr::drop_na()



# Baseline (no climate)
m0_null_legacy <- glmer(
  tBLW ~ lat_s + Age_s + (1 | SampleID),
  data = df_legacy,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 5e5))
)

##current year May and August CDD only
m_late_CDD <- glmer(
  tBLW ~ lat_s + Age_s +
    tmin_CDD_month5_s + tmin_CDD_month8_s +
    (1 | SampleID),
  data = df_legacy,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 5e5))
)

#previous year CDD
m_prevCDD <- glmer(
  tBLW ~ lat_s + Age_s +
    tmin_CDD_month7_prev_s +
    tmin_CDD_month8_prev_s +
    tmin_CDD_month9_prev_s +
    (1 | SampleID),
  data = df_legacy,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 5e5))
)

#carry over
m_carryoverCDD <- glmer(
  tBLW ~ lat_s + Age_s +
    tmin_CDD_month5_s + tmin_CDD_month8_s +
    tmin_CDD_month8_prev_s +
    tmin_CDD_month9_prev_s +
    (1 | SampleID),
  data = df_legacy,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 5e5))
)




library(MuMIn)
options(na.action = "na.fail")

legacy_set <- list(
  late_CDD      = m_late_CDD,
  prevCDD       = m_prevCDD,
  carryoverCDD  = m_carryoverCDD,
  null          = m0_null_legacy
)

sel_legacy <- model.sel(legacy_set, rank = "AICc")
print(sel_legacy)

# Supported models
supported_legacy <- get.models(sel_legacy, subset = delta <= 2)
names(supported_legacy)

#carryoverCDD optimum model
best_legacy <- m_carryoverCDD

tidy_bestlegacy <- broom.mixed::tidy(best_legacy, effects = "fixed", conf.int = TRUE)
print(tidy_bestlegacy)


# R2 of legacy model
r.squaredGLMM(m_carryoverCDD)


###check overdispersion
overdisp_fun(best_legacy)






#####
####legacy coeff table with ORs
best_legacy <- m_carryoverCDD

#### coefficients table with ORs (legacy model)
legacy_coef <- coef(summary(best_legacy)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term")

# Safety fallback (usually not needed for a single glmer model)
if (!"Estimate" %in% names(legacy_coef)) {
  est_col <- setdiff(names(legacy_coef), "term")[1]
  legacy_coef <- legacy_coef %>% dplyr::rename(Estimate = all_of(est_col))
}

# Confidence intervals for FIXED effects only
legacy_ci <- confint(best_legacy, parm = "beta_", method = "Wald") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  dplyr::rename(conf.low = `2.5 %`, conf.high = `97.5 %`)

# Combine into tidy table + odds ratios
tidy_legacy <- legacy_coef %>%
  dplyr::left_join(legacy_ci, by = "term") %>%
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
  file.path(output_dir, "Salix_LWBR_prevyear_model_coeffs.csv"),
  row.names = FALSE
)






legacy_effect_table <- tidy_bestlegacy %>%
  # keep fixed effects only
  dplyr::filter(effect == "fixed", term != "(Intercept)") %>%
  
  # compute odds ratios
  dplyr::mutate(
    odds_ratio   = exp(estimate),
    conf.low.OR  = exp(conf.low),
    conf.high.OR = exp(conf.high)
  ) %>%
  
  # relabel predictors
  dplyr::mutate(
    Predictor = dplyr::recode(term,
                              "lat_s" = "Latitude",
                              "Age_s" = "Ring age",
                              "tmin_CDD_month5_s" = "May CDD (current year)",
                              "tmin_CDD_month8_s" = "August CDD (current year)",
                              "tmin_CDD_month8_prev_s" = "August CDD (previous year)",
                              "tmin_CDD_month9_prev_s" = "September CDD (previous year)"
    ),
    `Odds ratio (95% CI)` =
      sprintf("%.2f (%.2f-%.2f)", odds_ratio, conf.low.OR, conf.high.OR)
  ) %>%
  
  # keep only presentation columns
  dplyr::select(Predictor, `Odds ratio (95% CI)`)

legacy_effect_table



write.csv(
  legacy_effect_table,
  file.path(output_dir, "Salix_LWBR_prevyear_model_OR.csv"),
  row.names = FALSE
)






##########
library(dplyr)
library(ggplot2)

# best-supported current-year model
best_current <- get.models(sel_current, 1)[[1]]

#  scale parameters from the unscaled columns 
scale_params_current <- list(
  m5_mean = mean(df_mod$tmin_CDD_month5, na.rm = TRUE),
  m5_sd   = sd(df_mod$tmin_CDD_month5, na.rm = TRUE),
  m8_mean = mean(df_mod$tmin_CDD_month8, na.rm = TRUE),
  m8_sd   = sd(df_mod$tmin_CDD_month8, na.rm = TRUE)
)

make_prediction_data_current <- function(var, orig_seq, sp) {
  # Hold other predictors at 0 (mean on scaled scale)
  nd <- data.frame(
    lat_s = 0,
    Age_s = 0,
    tmin_CDD_month5_s = 0,
    tmin_CDD_month8_s = 0,
    tmin_CDD_month9_s = 0   # <- required because best_current includes this term
  )
  
  scaled_seq <- switch(
    var,
    "tmin_CDD_month5_s" = (orig_seq - sp$m5_mean) / sp$m5_sd,
    "tmin_CDD_month8_s" = (orig_seq - sp$m8_mean) / sp$m8_sd
  )
  
  nd <- nd[rep(1, length(orig_seq)), ]
  nd[[var]] <- scaled_seq
  nd
}

seqs_current <- list(
  tmin_CDD_month5_s = seq(min(df_mod$tmin_CDD_month5, na.rm = TRUE),
                          max(df_mod$tmin_CDD_month5, na.rm = TRUE),
                          length.out = 100),
  tmin_CDD_month8_s = seq(min(df_mod$tmin_CDD_month8, na.rm = TRUE),
                          max(df_mod$tmin_CDD_month8, na.rm = TRUE),
                          length.out = 100)
)

pred_list_current <- lapply(names(seqs_current), function(v) {
  nd <- make_prediction_data_current(v, seqs_current[[v]], scale_params_current)
  
  pr <- predict(best_current, newdata = nd, type = "link", se.fit = TRUE, re.form = NA)
  
  data.frame(
    orig_x   = seqs_current[[v]],
    fit      = plogis(pr$fit),
    lower    = plogis(pr$fit - 1.96 * pr$se.fit),
    upper    = plogis(pr$fit + 1.96 * pr$se.fit),
    variable = v
  )
})

all_preds_current <- bind_rows(pred_list_current)

#LWBR curent year climate prediction plot
ggplot(all_preds_current, aes(x = orig_x, y = fit)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = "LWBR"),
    alpha = 0.25
  ) +
  geom_line(
    aes(color = "LWBR"),
    linewidth = 1
  ) +
  facet_wrap(
    ~ variable,
    scales = "free_x",
    strip.position = "bottom",
    labeller = as_labeller(c(
      tmin_CDD_month5_s = "May CDD",
      tmin_CDD_month8_s = "Aug CDD"
    ))
  ) +
  labs(
    x = NULL,
    y = "Predicted probability",
    fill  = "Ring Type",
    color = "Ring Type",
    title = "Salix LWBR - 95% CI"
  ) +
  scale_color_manual(values = c("LWBR" = "#386CB0")) +
  scale_fill_manual(values  = c("LWBR" = "#386CB0")) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      margin = margin(t = 6, b = 1)
    ),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(
      margin = margin(t = 2, b = 2)
    ),
    legend.position = "none"
  )


ggsave(
  file.path(figures_dir, "Salix_LWBR_current_climate_prediction.png"),
  width = 6,
  height = 4,
  dpi = 300
)





##for graphical abstract
ggplot(all_preds_current, aes(x = orig_x, y = fit)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = "LWBR"),
    alpha = 0.25
  ) +
  geom_line(
    aes(color = "LWBR"),
    linewidth = 1
  ) +
  facet_wrap(
    ~ variable,
    scales = "free_x",
    strip.position = "bottom",
    labeller = as_labeller(c(
      tmin_CDD_month5_s = "May CDD",
      tmin_CDD_month8_s = "Aug CDD"
    ))
  ) +
  labs(
    x = NULL,
    y = "Predicted probability",
    fill  = "Ring Type",
    color = "Ring Type",
    title = "Salix LWBR - 95% CI"
  ) +
  scale_color_manual(values = c("LWBR" = "#386CB0")) +
  scale_fill_manual(values  = c("LWBR" = "#386CB0")) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),      # remove all gridlines
    
    axis.line = element_line(           # add x and y axis lines
      linewidth = 0.6,
      colour = "black"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      margin = margin(t = 6, b = 15)
    ),
    
    axis.text.y = element_text(
      margin = margin(r = 6)
    ),
    
    strip.placement = "outside",
    legend.position = "none"
  )


ggsave(
  file.path(figures_dir, "graphicalabstract_Salix_LWBR_current_climate_prediction.png"),
  width = 6,
  height = 4,
  dpi = 300
)



