##############################################################################
# Alnus and Salix shrub anatomical response to climate
# ISR chronologies (Salix and Alnus)

library(here)

source(here("R/utils_packages.R"))
load_project_packages()

source(here("R/00_paths_setup.R")) 

# Project folders (repo root detected by .Rproj / here)
data_dir    <- here("data")
figures_dir <- here("figures")
output_dir  <- here("output")   # optional, for CSV outputs

# Create folders if missing (safe on fresh clone)
if (!dir.exists(figures_dir)) dir.create(figures_dir)
if (!dir.exists(output_dir))  dir.create(output_dir)

#Read data
sa_file <- file.path(data_dir, "ISR_salix_RWI.csv")
sa <- csv2rwl(fname = sa_file)


# Crossdating 
sa_rho <- interseries.cor(sa, prewhiten = TRUE, method = "spearman")

ro_tbl_sa <- data.frame(
  series = rownames(sa_rho),
  ro = sa_rho[,1]
)

#check only series with ro>=0.2
ro_tbl_sa$ro_round1 <- round(ro_tbl_sa$ro, 1)
keep_series <- ro_tbl_sa$series[ro_tbl_sa$ro_round1 >= 0.2]

sa_qc <- sa[, keep_series]


write.csv(
  ro_tbl_sa %>% mutate(keep = series %in% keep_series),
  file.path(output_dir, "salix_series_ro_screen.csv"),
  row.names = FALSE
)


#Detrend: 10-yr spline
sa_spline10 <- detrend(sa_qc, method = "Spline", n = 10)

#Truncate to years with >=5 series
sa_spline10_trunc <- sa_spline10[rowSums(!is.na(sa_spline10)) >= 5, ]
write.csv(sa_spline10_trunc, file.path(output_dir, "sa_spline10_trunc.csv"))

#Residual chronology 
salix_crn_res <- chron(sa_spline10_trunc, prefix = "", biweight = TRUE, prewhiten = TRUE)
write.csv(salix_crn_res, file.path(output_dir, "salix_residual_chron_spline10_trunc.csv"))

#Variance-stabilised chronology 
saspl10SV <- chron.stabilized(sa_spline10_trunc, winLength = 20)
# saspl10SV is not crn -  store year and vsc
saspl10SV_out <- data.frame(
  year = as.numeric(rownames(saspl10SV)),
  vsc  = saspl10SV$vsc
)
write.csv(saspl10SV_out, file.path(output_dir, "salix_variance_stabilised_chron.csv"), row.names = FALSE)

#Chronology strength metrics (EPS/rbar etc.)
salix_stats <- rwi.stats(sa_spline10_trunc)
write.csv(salix_stats, file.path(output_dir, "salix_rwi_stats_spline10_trunc.csv"))
salix_stats

#SSS
salix_sss <- sss(sa_spline10_trunc)
salix_mean_sss <- mean(salix_sss, na.rm = TRUE)

sss_out_sa <- data.frame(year = as.numeric(names(salix_sss)), sss = as.numeric(salix_sss))
write.csv(sss_out_sa, file.path(output_dir, "salix_sss_by_year.csv"), row.names = FALSE)


# view residual chronology
plot(salix_crn_res, add.spline = TRUE, nyrs = 10)
abline(h = 1, lty = 2)


#Chronology + SSS plot with 0.8 SSS threshold boundary 
ss <- salix_sss
yr <- as.numeric(names(ss))
vsc <- saspl10SV[as.character(yr), "vsc"]
sa.crn <- saspl10SV


png(file.path(figures_dir, "salix_chronology_SSS_spline10.png"),
    width = 10, height = 5, units = "in", res = 300)



def.par <- par(no.readonly = TRUE)
par(mar = c(2, 2, 2, 2), mgp = c(1.1, 0.1, 0), tcl = 0.25, xaxs='i')

# Base plot (RWI)
plot(yr, vsc, type = "n", xlab = "Year", ylab = "RWI", axes = FALSE)

# Axis limits from current device
lims  <- par("usr")
ymin  <- lims[3]
ymax  <- lims[4]

# Identify stretches where SSS <= 0.8
below_cutoff <- ss <= 0.8
rle_below    <- rle(below_cutoff)
ends         <- cumsum(rle_below$lengths)
starts       <- c(1, head(ends, -1) + 1)

#data frame of below-threshold runs
bad_runs <- data.frame(
  start_i = starts[rle_below$values],
  end_i   = ends[rle_below$values]
)
bad_runs$end_year <- yr[bad_runs$end_i]

# Shade below-threshold intervals and draw vertical line
for (i in seq_len(nrow(bad_runs))) {
  x_start <- yr[bad_runs$start_i[i]]
  x_end   <- bad_runs$end_year[i]
  
  # grey shade
  polygon(
    x = c(x_start, x_end, x_end, x_start),
    y = c(ymin, ymin, ymax, ymax),
    col = "grey80", border = NA
  )
  
  # vertical line 
  abline(v = x_end, col = "blue", lty = 2, lwd = 1)
}

# Add RWI data and axes
abline(h = 1, lwd = 1.5)
lines(yr, vsc, col = "grey50")
lines(yr, caps(vsc, nyrs = 10), col = "red", lwd = 1.5)
axis(1); axis(2); axis(3)

#Add end boundary label
label_y <- mean(range(vsc, na.rm = TRUE))  # midpoint of RWI values
for (i in seq_len(nrow(bad_runs))) {
  text(bad_runs$end_year[i] + 1.3, label_y,  # shift slightly right of the line
       labels = bad_runs$end_year[i],
       col = "blue", cex = 0.8, srt = 90, pos = 3)  # pos=4 means right of x coord
}


#Overlay SSS
par(new = TRUE)
plot(yr, ss, type = "l", xlab = "", ylab = "", axes = FALSE,
     col = "blue", ylim = c(0, 1))
abline(h = 0.80, col = "blue", lty = "dashed")
axis(4, at = pretty(ss))
mtext("SSS", side = 4, line = 1.1, lwd = 1.5)
box()

par(def.par)


dev.off()











##############################################################################
# ISR alnus chronology 

#Read data
al_file <- file.path(data_dir, "ISR_alnus_RWI.csv")
al <- csv2rwl(fname = al_file)


# Crossdating 
al_rho <- interseries.cor(al, prewhiten = TRUE, method = "spearman")

ro_tbl_al <- data.frame(
  series = rownames(al_rho),
  ro = al_rho[,1]
)

#check only series with ro>=0.2
ro_tbl_al$ro_round1 <- round(ro_tbl_al$ro, 1)
keep_series <- ro_tbl_al$series[ro_tbl_al$ro_round1 >= 0.2]

al_qc <- al[, keep_series]


write.csv(
  ro_tbl_al %>% mutate(keep = series %in% keep_series),
  file.path(output_dir, "alnus_series_ro_screen.csv"),
  row.names = FALSE
)

#Detrend: 10-yr spline
al_spline10 <- detrend(al_qc, method = "Spline", n = 10)

#Truncate to years with >=5 series
al_spline10_trunc <- al_spline10[rowSums(!is.na(al_spline10)) >= 5, ]
write.csv(
  al_spline10_trunc,
  file.path(output_dir, "al_spline10_trunc.csv")
)


#Residual chronology 
alnus_crn_res <- chron(al_spline10_trunc, prefix = "", biweight = TRUE, prewhiten = TRUE)
write.csv(
  alnus_crn_res,
  file.path(output_dir, "alnus_residual_chron_spline10_trunc.csv")
)


#Variance-stabilised chronology 
alspl10SV <- chron.stabilized(al_spline10_trunc, winLength = 20)
# alspl10SV is not crn -  store year and vsc
alspl10SV_out <- data.frame(
  year = as.numeric(rownames(alspl10SV)),
  vsc  = alspl10SV$vsc
)
alspl10SV_out
write.csv(
  alspl10SV_out,
  file.path(output_dir, "alnus_variance_stabilised_chron.csv"),
  row.names = FALSE
)

#Chronology strength metrics (EPS/rbar etc.)
alnus_stats <- rwi.stats(al_spline10_trunc)
write.csv(
  alnus_stats,
  file.path(output_dir, "alnus_rwi_stats_spline10_trunc.csv")
)

alnus_stats

#SSS
alnus_sss <- sss(al_spline10_trunc)
alnus_mean_sss <- mean(alnus_sss, na.rm = TRUE)

sss_out_al <- data.frame(year = as.numeric(names(alnus_sss)), sss = as.numeric(alnus_sss))
write.csv(
  sss_out_al,
  file.path(output_dir, "alnus_sss_by_year.csv"),
  row.names = FALSE
)


#view residual chronology
plot(alnus_crn_res, add.spline = TRUE, nyrs = 10)
abline(h = 1, lty = 2)




#Chronology + SSS plot with 0.8 SSS threshold boundary 
ss <- alnus_sss
yr <- as.numeric(names(ss))
vsc <- alspl10SV[as.character(yr), "vsc"]
al.crn <- alspl10SV


png(file.path(figures_dir, "alnus_chronology_SSS_spline10.png"),
    width = 10, height = 5, units = "in", res = 300)


def.par <- par(no.readonly = TRUE)
par(mar = c(2, 2, 2, 2), mgp = c(1.1, 0.1, 0), tcl = 0.25, xaxs='i')

# Base plot (RWI)
plot(yr, vsc, type = "n", xlab = "Year", ylab = "RWI", axes = FALSE)

# Axis limits from current device
lims  <- par("usr")
ymin  <- lims[3]
ymax  <- lims[4]

# Identify stretches where SSS <= 0.8
below_cutoff <- ss <= 0.8
rle_below    <- rle(below_cutoff)
ends         <- cumsum(rle_below$lengths)
starts       <- c(1, head(ends, -1) + 1)

#data frame of below-threshold runs
bad_runs <- data.frame(
  start_i = starts[rle_below$values],
  end_i   = ends[rle_below$values]
)
bad_runs$end_year <- yr[bad_runs$end_i]

# Shade below-threshold intervals and draw vertical line
for (i in seq_len(nrow(bad_runs))) {
  x_start <- yr[bad_runs$start_i[i]]
  x_end   <- bad_runs$end_year[i]
  
  # grey shade
  polygon(
    x = c(x_start, x_end, x_end, x_start),
    y = c(ymin, ymin, ymax, ymax),
    col = "grey80", border = NA
  )
  
  # vertical line 
  abline(v = x_end, col = "blue", lty = 2, lwd = 1)
}

# Add RWI data and axes
abline(h = 1, lwd = 1.5)
lines(yr, vsc, col = "grey50")
lines(yr, caps(vsc, nyrs = 10), col = "red", lwd = 1.5)
axis(1); axis(2); axis(3)

#Add end boundary label
label_y <- mean(range(vsc, na.rm = TRUE))  # midpoint of RWI values
for (i in seq_len(nrow(bad_runs))) {
  text(bad_runs$end_year[i] + 1.3, label_y,  # shift slightly right of the line
       labels = bad_runs$end_year[i],
       col = "blue", cex = 0.8, srt = 90, pos = 3)  # pos=4 means right of x coord
}


#Overlay SSS
par(new = TRUE)
plot(yr, ss, type = "l", xlab = "", ylab = "", axes = FALSE,
     col = "blue", ylim = c(0, 1))
abline(h = 0.80, col = "blue", lty = "dashed")
axis(4, at = pretty(ss))
mtext("SSS", side = 4, line = 1.1, lwd = 1.5)
box()

par(def.par)

dev.off()
