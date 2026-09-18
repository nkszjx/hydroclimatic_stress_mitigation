
# Robustness checks: Alternative functional forms

rm(list=ls())
library(parallel)
library(lmtest)
library(DescTools)
library(foreign)
library(Matrix)
library(lfe)
library(magrittr)
library(naniar)
library(dplyr)
library(plotly)
library(zoo)
library(readxl)
library(mice)
library(rio)
library(car)
library(AER)
library(lme4)
library(ggplot2)
library(brms)
library(mgcv)
library(future)
library(plm)
library(data.table)
library(fixest)
library(clubSandwich)
library(ggpubr)
library(webshot)
library(htmlwidgets)
library(scales)
library(scatterplot3d)
library(extrafont)
library(showtext)
library(future.apply)
library(stargazer)
library(corrplot)
library(stringr)
library(patchwork)
library(cowplot)
library(splines)   

###########################################################################
#                   Two-way interactive quadratic regression
###########################################################################

setwd("XXX") # please set your folder path
merged_dt <- fread("Dataset_matched_buffer10km.csv")

# Absolute latitude to unify hemispheres
merged_dt$Latitude <- abs(merged_dt$Latitude)


###########################################################################
#                     Part 1: FVC ~ Water deficit
###########################################################################
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude',
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population',
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)),
                                                 .SDcols = cols_to_center]
merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))

knots_wd <- quantile(merged_dt_WaterDeficit$WaterDeficit, probs = c(0.10,0.50,0.90), na.rm = TRUE)

model_urban_rcs <- felm(
  FVC ~ 1 +
    ns(WaterDeficit, knots = knots_wd) +
    urban +
    ns(WaterDeficit, knots = knots_wd):urban +
    Seasonalityindex +
    AverageTemperature +
    TemperatureRange +
    WindSpeed +
    Elevation +
    Latitude +
    Precipitation +
    SoilMoisture +
    SoilPH +
    GDP_per_capita_PPP + HDI + Population +
    ImperviousSurface + HumanSettlement + CityArea
  | CityID + year | 0 | CityID,
  data = merged_dt_WaterDeficit
)
summary(model_urban_rcs)

output_file <- "FVC_urban_in_out_RCS.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_rcs))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# Compute effects, marginal effects and confidence intervals
P_min <- quantile(merged_dt_WaterDeficit$WaterDeficit, 0.001, na.rm = TRUE)
P_max <- quantile(merged_dt_WaterDeficit$WaterDeficit, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)
rcs_mat_pred <- ns(P_values, knots = knots_wd)
n_spline_basis <- ncol(rcs_mat_pred)
coef_all <- coef(model_urban_rcs)
rcs_names <- paste0("ns(WaterDeficit, knots = knots_wd)", 1:n_spline_basis)
rcs_inter_names <- paste0(rcs_names, ":urbanin")
beta_urban_fac <- coef_all["urbanin"]                   
beta_nonurban   <- coef_all[rcs_names]                  
beta_urban_incr <- coef_all[rcs_inter_names]            
beta_urban      <- beta_nonurban + beta_urban_incr 
effect_non_urban <- as.vector(rcs_mat_pred %*% beta_nonurban)
effect_urban     <- as.vector(beta_urban_fac + rcs_mat_pred %*% beta_urban)
vcov_mat <- vcov(model_urban_rcs)
vcov_nonurb <- vcov_mat[rcs_names, rcs_names]
var_nonurb <- rowSums((rcs_mat_pred %*% vcov_nonurb) * rcs_mat_pred)
se_nonurb <- sqrt(var_nonurb)
ci_nonurb_lower <- effect_non_urban - 1.96*se_nonurb
ci_nonurb_upper <- effect_non_urban + 1.96*se_nonurb
vcov_urb_sub <- vcov_mat[c("urbanin", rcs_names, rcs_inter_names), c("urbanin", rcs_names, rcs_inter_names)]
n_pred <- length(P_values)
var_urb <- numeric(n_pred)
for(i in 1:n_pred){
  g <- c(1, rcs_mat_pred[i, ], rcs_mat_pred[i, ])
  var_urb[i] <- as.vector(t(g) %*% vcov_urb_sub %*% g)
}
se_urb <- sqrt(var_urb)
ci_urb_lower <- effect_urban - 1.96*se_urb
ci_urb_upper <- effect_urban + 1.96*se_urb

effect_df <- rbind(
  data.frame(
    WaterDeficit = P_values,
    Region = "Non‑urban",
    Effect = effect_non_urban,
    SE = se_nonurb,
    CI_Lower = ci_nonurb_lower,
    CI_Upper = ci_nonurb_upper
  ),
  data.frame(
    WaterDeficit = P_values,
    Region = "Urban",
    Effect = effect_urban,
    SE = se_urb,
    CI_Lower = ci_urb_lower,
    CI_Upper = ci_urb_upper
  )
)

urban_colors <- c(
  "Urban" = "#FF7A3C",
  "Non‑urban" = "#00A850"
)
urban_colors_hist <- c(
  "in" = "#FF7A3C",
  "out" = "#00A850"
)

# Plot effect curve
effect_plot <- ggplot(effect_df, aes(x = WaterDeficit, y = Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.2, color = NA
  ) +
  scale_linetype_manual(values = c("Non‑urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non‑urban" = 0.4, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "Water deficit",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title.y = element_text(size = 6),
    legend.position = c(0.3, 0.3),
    legend.direction = "vertical",
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) +
  guides(
    color = guide_legend(
      keywidth = 0.9,
      keyheight = 0.5,
      override.aes = list(
        linetype = c("Non‑urban"="dashed", "Urban"="solid"),
        linewidth = c("Non‑urban" = 0.4, "Urban" = 0.4)
      )
    )
  )

# Plot distribution histogram
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill = urban), binwidth =  5,  alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600)) +
  scale_y_continuous(limits = c(0.0, 0.007),
                     breaks = seq(0.0, 0.007, by = 0.002),
                     labels = scales::number_format(accuracy = 0.001)) +
  labs(x = "Water deficit (mm)",
       y = "") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title.x = element_text(size = 6),
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist)

effect_waterdeficit_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_waterdeficit_FVC)



###########################################################################
#                     Part 2: NPP ~ Water deficit
###########################################################################
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude',
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population',
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)),
                                                 .SDcols = cols_to_center]
merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))

knots_wd <- quantile(merged_dt_WaterDeficit$WaterDeficit, probs = c(0.10,0.50,0.90), na.rm = TRUE)

model_urban_rcs <- felm(
  NPP ~ 1 +
    ns(WaterDeficit, knots = knots_wd) +
    urban +
    ns(WaterDeficit, knots = knots_wd):urban +
    Seasonalityindex +
    AverageTemperature +
    TemperatureRange +
    WindSpeed +
    Elevation +
    Latitude +
    Precipitation +
    SoilMoisture +
    SoilPH +
    GDP_per_capita_PPP + HDI + Population +
    ImperviousSurface + HumanSettlement + CityArea
  | CityID + year | 0 | CityID,
  data = merged_dt_WaterDeficit
)
summary(model_urban_rcs)

output_file <- "NPP_urban_in_out_RCS.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_rcs))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# Compute effects, marginal effects and confidence intervals
P_min <- quantile(merged_dt_WaterDeficit$WaterDeficit, 0.001, na.rm = TRUE)
P_max <- quantile(merged_dt_WaterDeficit$WaterDeficit, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)
rcs_mat_pred <- ns(P_values, knots = knots_wd)
n_spline_basis <- ncol(rcs_mat_pred)
coef_all <- coef(model_urban_rcs)
rcs_names <- paste0("ns(WaterDeficit, knots = knots_wd)", 1:n_spline_basis)
rcs_inter_names <- paste0(rcs_names, ":urbanin")
beta_urban_fac <- coef_all["urbanin"]                   
beta_nonurban   <- coef_all[rcs_names]                  
beta_urban_incr <- coef_all[rcs_inter_names]            
beta_urban      <- beta_nonurban + beta_urban_incr      
effect_non_urban <- as.vector(rcs_mat_pred %*% beta_nonurban)
effect_urban     <- as.vector(beta_urban_fac + rcs_mat_pred %*% beta_urban)
vcov_mat <- vcov(model_urban_rcs)
vcov_nonurb <- vcov_mat[rcs_names, rcs_names]
var_nonurb <- rowSums((rcs_mat_pred %*% vcov_nonurb) * rcs_mat_pred)
se_nonurb <- sqrt(var_nonurb)
ci_nonurb_lower <- effect_non_urban - 1.96*se_nonurb
ci_nonurb_upper <- effect_non_urban + 1.96*se_nonurb
vcov_urb_sub <- vcov_mat[c("urbanin", rcs_names, rcs_inter_names), c("urbanin", rcs_names, rcs_inter_names)]
n_pred <- length(P_values)
var_urb <- numeric(n_pred)
for(i in 1:n_pred){
  g <- c(1, rcs_mat_pred[i, ], rcs_mat_pred[i, ])
  var_urb[i] <- as.vector(t(g) %*% vcov_urb_sub %*% g)
}
se_urb <- sqrt(var_urb)
ci_urb_lower <- effect_urban - 1.96*se_urb
ci_urb_upper <- effect_urban + 1.96*se_urb


effect_df <- rbind(
  data.frame(
    WaterDeficit = P_values,
    Region = "Non‑urban",
    Effect = effect_non_urban,
    SE = se_nonurb,
    CI_Lower = ci_nonurb_lower,
    CI_Upper = ci_nonurb_upper
  ),
  data.frame(
    WaterDeficit = P_values,
    Region = "Urban",
    Effect = effect_urban,
    SE = se_urb,
    CI_Lower = ci_urb_lower,
    CI_Upper = ci_urb_upper
  )
)

urban_colors <- c(
  "Urban" = "#FF7A3C",
  "Non‑urban" = "#00A850"
)
urban_colors_hist <- c(
  "in" = "#FF7A3C",
  "out" = "#00A850"
)

# Plot effect curve
effect_plot <- ggplot(effect_df, aes(x = WaterDeficit, y = Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.2, color = NA
  ) +
  scale_linetype_manual(values = c("Non‑urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non‑urban" = 0.4, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "Water deficit",
    y = expression(Effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title.y = element_text(size = 6),
    legend.position = c(0.3, 0.3),
    legend.direction = "vertical",
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) +
  guides(
    color = guide_legend(
      keywidth = 0.9,
      keyheight = 0.5,
      override.aes = list(
        linetype = c("Non‑urban"="dashed", "Urban"="solid"),
        linewidth = c("Non‑urban" = 0.4, "Urban" = 0.4)
      )
    )
  )

# Plot distribution histogram
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill = urban), binwidth =  5,  alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600)) +
  scale_y_continuous(limits = c(0.0, 0.007),
                     breaks = seq(0.0, 0.007, by = 0.002),
                     labels = scales::number_format(accuracy = 0.001)) +
  labs(x = "Water deficit (mm)",
       y = "") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title.x = element_text(size = 6),
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist)

effect_waterdeficit_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_waterdeficit_NPP)




###########################################################################
#                     Part 3: FVC ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)),
                                                 .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))

knots_wd <- quantile(merged_dt_Seasonalityindex$Seasonalityindex, probs = c(0.10,0.50,0.90), na.rm = TRUE)

model_urban_rcs <- felm(
  FVC ~ 1 +
    ns(Seasonalityindex, knots = knots_wd) +
    urban +
    ns(Seasonalityindex, knots = knots_wd):urban +
    WaterDeficit +
    AverageTemperature +
    TemperatureRange +
    WindSpeed +
    Elevation +
    Latitude +
    Precipitation +
    SoilMoisture +
    SoilPH +
    GDP_per_capita_PPP + HDI + Population +
    ImperviousSurface + HumanSettlement + CityArea
  | CityID + year | 0 | CityID,
  data = merged_dt_Seasonalityindex
)
summary(model_urban_rcs)

output_file <- "FVC_urban_in_out_RCS.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_rcs))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# Compute effects, marginal effects and confidence intervals
P_min <- quantile(merged_dt_Seasonalityindex$Seasonalityindex, 0.001, na.rm = TRUE)
P_max <- quantile(merged_dt_Seasonalityindex$Seasonalityindex, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)
rcs_mat_pred <- ns(P_values, knots = knots_wd)
n_spline_basis <- ncol(rcs_mat_pred)
coef_all <- coef(model_urban_rcs)
rcs_names <- paste0("ns(Seasonalityindex, knots = knots_wd)", 1:n_spline_basis)
rcs_inter_names <- paste0(rcs_names, ":urbanin")
beta_urban_fac <- coef_all["urbanin"]                   
beta_nonurban   <- coef_all[rcs_names]                  
beta_urban_incr <- coef_all[rcs_inter_names]            
beta_urban      <- beta_nonurban + beta_urban_incr      
effect_non_urban <- as.vector(rcs_mat_pred %*% beta_nonurban)
effect_urban     <- as.vector(beta_urban_fac + rcs_mat_pred %*% beta_urban)
vcov_mat <- vcov(model_urban_rcs)
vcov_nonurb <- vcov_mat[rcs_names, rcs_names]
var_nonurb <- rowSums((rcs_mat_pred %*% vcov_nonurb) * rcs_mat_pred)
se_nonurb <- sqrt(var_nonurb)
ci_nonurb_lower <- effect_non_urban - 1.96*se_nonurb
ci_nonurb_upper <- effect_non_urban + 1.96*se_nonurb
vcov_urb_sub <- vcov_mat[c("urbanin", rcs_names, rcs_inter_names), c("urbanin", rcs_names, rcs_inter_names)]
n_pred <- length(P_values)
var_urb <- numeric(n_pred)
for(i in 1:n_pred){
  g <- c(1, rcs_mat_pred[i, ], rcs_mat_pred[i, ])
  var_urb[i] <- as.vector(t(g) %*% vcov_urb_sub %*% g)
}
se_urb <- sqrt(var_urb)
ci_urb_lower <- effect_urban - 1.96*se_urb
ci_urb_upper <- effect_urban + 1.96*se_urb

effect_df <- rbind(
  data.frame(
    Seasonalityindex = P_values,
    Region = "Non‑urban",
    Effect = effect_non_urban,
    SE = se_nonurb,
    CI_Lower = ci_nonurb_lower,
    CI_Upper = ci_nonurb_upper
  ),
  data.frame(
    Seasonalityindex = P_values,
    Region = "Urban",
    Effect = effect_urban,
    SE = se_urb,
    CI_Lower = ci_urb_lower,
    CI_Upper = ci_urb_upper
  )
)

urban_colors <- c(
  "Urban" = "#FF7A3C",
  "Non‑urban" = "#00A850"
)
urban_colors_hist <- c(
  "in" = "#FF7A3C",
  "out" = "#00A850"
)

# Plot effect curve
effect_plot <- ggplot(effect_df, aes(x = Seasonalityindex, y = Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.2, color = NA
  ) +
  scale_linetype_manual(values = c("Non‑urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non‑urban" = 0.4, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "Rainfall seasonality",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title.y = element_text(size = 6),
    legend.position = c(0.3, 0.3),
    legend.direction = "vertical",
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) +
  guides(
    color = guide_legend(
      keywidth = 0.9,
      keyheight = 0.5,
      override.aes = list(
        linetype = c("Non‑urban"="dashed", "Urban"="solid"),
        linewidth = c("Non‑urban" = 0.4, "Urban" = 0.4)
      )
    )
  )

# Plot distribution histogram
hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill = urban), binwidth =  0.0001,  alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02)) +
  labs(x = "Rainfall seasonality",
       y = "") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title.x = element_text(size = 6),
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist)

effect_Seasonalityindex_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_Seasonalityindex_FVC)



###########################################################################
#                     Part 4: NPP ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude',
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population',
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)),
                                                 .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))

knots_wd <- quantile(merged_dt_Seasonalityindex$Seasonalityindex, probs = c(0.10,0.50,0.90), na.rm = TRUE)

model_urban_rcs <- felm(
  NPP ~ 1 +
    ns(Seasonalityindex, knots = knots_wd) +
    urban +
    ns(Seasonalityindex, knots = knots_wd):urban +
    WaterDeficit +
    AverageTemperature +
    TemperatureRange +
    WindSpeed +
    Elevation +
    Latitude +
    Precipitation +
    SoilMoisture +
    SoilPH +
    GDP_per_capita_PPP + HDI + Population +
    ImperviousSurface + HumanSettlement + CityArea
  | CityID + year | 0 | CityID,
  data = merged_dt_Seasonalityindex
)
summary(model_urban_rcs)

output_file <- "NPP_urban_in_out_RCS.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_rcs))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# Compute effects, marginal effects and confidence intervals
P_min <- quantile(merged_dt_Seasonalityindex$Seasonalityindex, 0.001, na.rm = TRUE)
P_max <- quantile(merged_dt_Seasonalityindex$Seasonalityindex, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)
rcs_mat_pred <- ns(P_values, knots = knots_wd)
n_spline_basis <- ncol(rcs_mat_pred)
coef_all <- coef(model_urban_rcs)
rcs_names <- paste0("ns(Seasonalityindex, knots = knots_wd)", 1:n_spline_basis)
rcs_inter_names <- paste0(rcs_names, ":urbanin")
beta_urban_fac <- coef_all["urbanin"]                   
beta_nonurban   <- coef_all[rcs_names]                  
beta_urban_incr <- coef_all[rcs_inter_names]            
beta_urban      <- beta_nonurban + beta_urban_incr      
effect_non_urban <- as.vector(rcs_mat_pred %*% beta_nonurban)
effect_urban     <- as.vector(beta_urban_fac + rcs_mat_pred %*% beta_urban)
vcov_mat <- vcov(model_urban_rcs)
vcov_nonurb <- vcov_mat[rcs_names, rcs_names]
var_nonurb <- rowSums((rcs_mat_pred %*% vcov_nonurb) * rcs_mat_pred)
se_nonurb <- sqrt(var_nonurb)
ci_nonurb_lower <- effect_non_urban - 1.96*se_nonurb
ci_nonurb_upper <- effect_non_urban + 1.96*se_nonurb
vcov_urb_sub <- vcov_mat[c("urbanin", rcs_names, rcs_inter_names), c("urbanin", rcs_names, rcs_inter_names)]
n_pred <- length(P_values)
var_urb <- numeric(n_pred)
for(i in 1:n_pred){
  g <- c(1, rcs_mat_pred[i, ], rcs_mat_pred[i, ])
  var_urb[i] <- as.vector(t(g) %*% vcov_urb_sub %*% g)
}
se_urb <- sqrt(var_urb)
ci_urb_lower <- effect_urban - 1.96*se_urb
ci_urb_upper <- effect_urban + 1.96*se_urb

effect_df <- rbind(
  data.frame(
    Seasonalityindex = P_values,
    Region = "Non‑urban",
    Effect = effect_non_urban,
    SE = se_nonurb,
    CI_Lower = ci_nonurb_lower,
    CI_Upper = ci_nonurb_upper
  ),
  data.frame(
    Seasonalityindex = P_values,
    Region = "Urban",
    Effect = effect_urban,
    SE = se_urb,
    CI_Lower = ci_urb_lower,
    CI_Upper = ci_urb_upper
  )
)

urban_colors <- c(
  "Urban" = "#FF7A3C",
  "Non‑urban" = "#00A850"
)
urban_colors_hist <- c(
  "in" = "#FF7A3C",
  "out" = "#00A850"
)

# Plot effect curve
effect_plot <- ggplot(effect_df, aes(x = Seasonalityindex, y = Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.2, color = NA
  ) +
  scale_linetype_manual(values = c("Non‑urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non‑urban" = 0.4, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "Rainfall seasonality",
    y = expression(Effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title.y = element_text(size = 6),
    legend.position = c(0.3, 0.6),
    legend.direction = "vertical",
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) +
  guides(
    color = guide_legend(
      keywidth = 0.9,
      keyheight = 0.5,
      override.aes = list(
        linetype = c("Non‑urban"="dashed", "Urban"="solid"),
        linewidth = c("Non‑urban" = 0.4, "Urban" = 0.4)
      )
    )
  )

# Plot distribution histogram
hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill = urban), binwidth =  0.0001,  alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02)) +
  labs(x = "Rainfall seasonality",
       y = "") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title.x = element_text(size = 6),
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist)

effect_Seasonalityindex_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_Seasonalityindex_NPP)



################################################ Composite figure ###############################################################
p_final <- cowplot::plot_grid(
  effect_waterdeficit_FVC,
  effect_waterdeficit_NPP,
  effect_Seasonalityindex_FVC,
  effect_Seasonalityindex_NPP,
  nrow = 2,
  ncol = 2,
  labels = c("a","b","c","d"),
  label_size = 10,
  scale = 1,
  label_x = 0.02,
  label_y = 1.05,
  align = "hv",
  rel_widths = c(1,1),
  rel_heights = c(1,1)
) +
  theme(plot.margin = margin(25, 5, 5, 5))

# 输出
print(p_final)
ggsave("Supplementary_Fig7.pdf",
       plot = p_final, width = 130, height = 130, units = "mm", dpi = 600)








