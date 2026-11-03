

# Marginal effect on FVC and NPP

rm(list=ls())  
library(parallel)
library(lmtest)
library(DescTools)
library(foreign)
library(Matrix)
library(lfe)  
library(magrittr)
library(margins)
library(naniar)
library(dplyr)
library(plotly)
library(zoo)
library(readxl)
library(mice)
library(rio)
library(orca)
library(car)
library(AER)
library(lme4)
library(ggplot2)
library(brms)
library(mgcv)
library(future)
library(plm)
library(dplyr)
library(data.table)
library(fixest)
library(clubSandwich) 
library(ggpubr)
library(marginaleffects)
library(webshot)
library(htmlwidgets)
library(scales)
library(scatterplot3d)
library(extrafont)
library(showtext)
library(future)
library(future.apply)
library(stargazer)
library(data.table)
library(corrplot)
library(stringr)
library(patchwork)
library(cowplot)
library(patchwork)
library(scales)
library(ggrastr)  
library(scattermore)



###########################################################################
#                   Two-way interactive quadratic regression
###########################################################################

setwd("XXX") # please set your folder path
merged_dt <- fread("Dataset_matched_buffer10km.csv")

merged_dt$Latitude <- abs(merged_dt$Latitude)
merged_dt$WaterDeficit_2 <- merged_dt$WaterDeficit * merged_dt$WaterDeficit
merged_dt$Seasonalityindex_2 <- merged_dt$Seasonalityindex * merged_dt$Seasonalityindex



###########################################################################
#                     FVC ~ Water deficit
###########################################################################
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]
merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))

model_urban <- felm(FVC ~ 1+WaterDeficit +
                      WaterDeficit_2 +
                      urban + 
                      Seasonalityindex +
                      AverageTemperature + 
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Precipitation + 
                      SoilMoisture +
                      SoilPH +
                      WaterDeficit:urban +
                      WaterDeficit_2:urban +
                      GDP_per_capita_PPP + HDI + Population + 
                      ImperviousSurface + HumanSettlement + CityArea
                    | CityID + year | 0 | CityID, data = merged_dt_WaterDeficit)  
summary(model_urban)

output_file <- "FVC_urban_in_out2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

coeffs <- coef(model_urban)
beta_urban <- coeffs["urbanin"] 
beta1 <- coeffs["WaterDeficit"]              
beta2 <- coeffs["WaterDeficit_2"]             
beta1_city <- coeffs["WaterDeficit:urbanin"]  
beta2_city <- coeffs["WaterDeficit_2:urbanin"]
beta1_total <- beta1 + beta1_city
beta2_total <- beta2 + beta2_city
P_min <- quantile(merged_dt_WaterDeficit$WaterDeficit, 0.001, na.rm = TRUE)
P_max <- quantile(merged_dt_WaterDeficit$WaterDeficit, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200) 
effect_non_urban <- beta1 * P_values + beta2 * (P_values^2)  
effect_urban <- beta_urban + beta1_total * P_values + beta2_total * (P_values^2)  
vcov_matrix <- vcov(model_urban)
cov_vars_non_urban <- c("WaterDeficit", "WaterDeficit_2")
vcov_non_urban <- vcov_matrix[cov_vars_non_urban, cov_vars_non_urban]
cov_vars_urban <- c("WaterDeficit", "WaterDeficit_2", "WaterDeficit:urbanin", "WaterDeficit_2:urbanin")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]


marginal_non_urban <- beta1 + 2 * beta2 * P_values  
marginal_urban <- beta1_total + 2 * beta2_total * P_values  
d_marg_beta1 <- 1  
d_marg_beta2 <- 2 * P_values  
var_marginal_non_urban <- (d_marg_beta1^2) * vcov_non_urban["WaterDeficit", "WaterDeficit"] +
  (d_marg_beta2^2) * vcov_non_urban["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_non_urban["WaterDeficit", "WaterDeficit_2"]
se_marginal_non_urban <- sqrt(var_marginal_non_urban)
ci_marginal_non_urban_lower <- marginal_non_urban - 1.96 * se_marginal_non_urban
ci_marginal_non_urban_upper <- marginal_non_urban + 1.96 * se_marginal_non_urban
d_marg_beta1_u <- 1  
d_marg_beta2_u <- 2 * P_values  
var_marginal_urban <- (d_marg_beta1_u^2) * vcov_urban["WaterDeficit", "WaterDeficit"] +
  (d_marg_beta2_u^2) * vcov_urban["WaterDeficit_2", "WaterDeficit_2"] +
  (d_marg_beta1_u^2) * vcov_urban["WaterDeficit:urbanin", "WaterDeficit:urbanin"] +
  (d_marg_beta2_u^2) * vcov_urban["WaterDeficit_2:urbanin", "WaterDeficit_2:urbanin"] +
  2 * d_marg_beta1_u^2 * vcov_urban["WaterDeficit", "WaterDeficit:urbanin"] +  
  2 * d_marg_beta2_u^2 * vcov_urban["WaterDeficit_2", "WaterDeficit_2:urbanin"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2:urbanin"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit:urbanin", "WaterDeficit_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit:urbanin", "WaterDeficit_2:urbanin"]  
se_marginal_urban <- sqrt(var_marginal_urban)
ci_marginal_urban_lower <- marginal_urban - 1.96 * se_marginal_urban
ci_marginal_urban_upper <- marginal_urban + 1.96 * se_marginal_urban

marginal_df <- rbind(
  data.frame(
    WaterDeficit = P_values,
    Region = "Non-urban",
    Marginal_Effect = marginal_non_urban,
    SE = se_marginal_non_urban,
    CI_Lower = ci_marginal_non_urban_lower,
    CI_Upper = ci_marginal_non_urban_upper
  ),
  data.frame(
    WaterDeficit = P_values,
    Region = "Urban",
    Marginal_Effect = marginal_urban,
    SE = se_marginal_urban,
    CI_Lower = ci_marginal_urban_lower,
    CI_Upper = ci_marginal_urban_upper
  )
)

urban_colors <- c(
  "Urban" = "#D62728", 
  "Non-urban" = "#1F77B4"
)
urban_colors_hist <- c(
  "in" = "#D62728", 
  "out" = "#1F77B4"
)



hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill = urban), binwidth =  5,  alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600)) +
  scale_y_continuous(limits = c(0.0, 0.007),
                     breaks = seq(0.0, 0.007, by = 0.002),
                     labels = scales::number_format(accuracy = 0.001)) +
  labs(x = "Water deficit (mm)", 
       y = "",
       subtitle = paste0("")) +
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.title.x = element_text(size = 15), 
    plot.margin = margin(t=-1.30,b=0.5,l=0.15,r=0.2, unit = "cm"), 
    axis.text.y = element_text(color = "black", size = 14),  
    axis.text.x = element_text(color = "black", size = 14),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist) 

print(hist_plot)

# Plot marginal effect curve
marginal_plot <- ggplot(marginal_df, aes(x = WaterDeficit, y = Marginal_Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha=0.3, color = NA
  ) +
  labs(
    x = "WaterDeficit",
    y = "Marginal effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  scale_linetype_manual(values = c("Non-urban"="solid", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 1, "Urban" = 1), guide="none") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600)) +
  scale_y_continuous(limits = c(-0.0004, 0.0002),
                     breaks = seq(-0.0004, 0.0002, by = 0.0002),
                     labels = scales::number_format(accuracy = 0.0001)) +
  
  theme(
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 14),  
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),   
    axis.title.y = element_text(size = 15),
    legend.position = c(0.7, 0.9), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 14, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) 

print(marginal_plot)

marginal_waterdeficit_FVC <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(marginal_waterdeficit_FVC)



###########################################################################
#                     NPP ~ Water deficit
###########################################################################
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]
merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))

model_urban <- felm(NPP ~ 1 + WaterDeficit +
                      WaterDeficit_2 +
                      urban + 
                      Seasonalityindex +
                      AverageTemperature + 
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Precipitation + 
                      SoilMoisture +
                      SoilPH +
                      WaterDeficit:urban +
                      WaterDeficit_2:urban +
                      GDP_per_capita_PPP + HDI + Population + 
                      ImperviousSurface + HumanSettlement + CityArea
                    | CityID + year | 0 | CityID, data = merged_dt_WaterDeficit)  
summary(model_urban)

output_file <- "NPP_urban_in_out2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

coeffs <- coef(model_urban)
beta_urban <- coeffs["urbanin"]
beta1 <- coeffs["WaterDeficit"]              
beta2 <- coeffs["WaterDeficit_2"]             
beta1_city <- coeffs["WaterDeficit:urbanin"]  
beta2_city <- coeffs["WaterDeficit_2:urbanin"]
beta1_total <- beta1 + beta1_city
beta2_total <- beta2 + beta2_city
P_min <- quantile(merged_dt_WaterDeficit$WaterDeficit, 0.001, na.rm = TRUE)
P_max <- quantile(merged_dt_WaterDeficit$WaterDeficit, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200) 
effect_non_urban <- beta1 * P_values + beta2 * (P_values^2)  
effect_urban <- beta_urban + beta1_total * P_values + beta2_total * (P_values^2)  
vcov_matrix <- vcov(model_urban)
cov_vars_non_urban <- c("WaterDeficit", "WaterDeficit_2")
vcov_non_urban <- vcov_matrix[cov_vars_non_urban, cov_vars_non_urban]
cov_vars_urban <- c("WaterDeficit", "WaterDeficit_2", "WaterDeficit:urbanin", "WaterDeficit_2:urbanin")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]

marginal_non_urban <- beta1 + 2 * beta2 * P_values  
marginal_urban <- beta1_total + 2 * beta2_total * P_values  
d_marg_beta1 <- 1  
d_marg_beta2 <- 2 * P_values  
var_marginal_non_urban <- (d_marg_beta1^2) * vcov_non_urban["WaterDeficit", "WaterDeficit"] +
  (d_marg_beta2^2) * vcov_non_urban["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_non_urban["WaterDeficit", "WaterDeficit_2"]
se_marginal_non_urban <- sqrt(var_marginal_non_urban)
ci_marginal_non_urban_lower <- marginal_non_urban - 1.96 * se_marginal_non_urban
ci_marginal_non_urban_upper <- marginal_non_urban + 1.96 * se_marginal_non_urban
d_marg_beta1_u <- 1  
d_marg_beta2_u <- 2 * P_values  
var_marginal_urban <- (d_marg_beta1_u^2) * vcov_urban["WaterDeficit", "WaterDeficit"] +
  (d_marg_beta2_u^2) * vcov_urban["WaterDeficit_2", "WaterDeficit_2"] +
  (d_marg_beta1_u^2) * vcov_urban["WaterDeficit:urbanin", "WaterDeficit:urbanin"] +
  (d_marg_beta2_u^2) * vcov_urban["WaterDeficit_2:urbanin", "WaterDeficit_2:urbanin"] +
  2 * d_marg_beta1_u^2 * vcov_urban["WaterDeficit", "WaterDeficit:urbanin"] +  
  2 * d_marg_beta2_u^2 * vcov_urban["WaterDeficit_2", "WaterDeficit_2:urbanin"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2:urbanin"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit:urbanin", "WaterDeficit_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit:urbanin", "WaterDeficit_2:urbanin"]  
se_marginal_urban <- sqrt(var_marginal_urban)
ci_marginal_urban_lower <- marginal_urban - 1.96 * se_marginal_urban
ci_marginal_urban_upper <- marginal_urban + 1.96 * se_marginal_urban

marginal_df <- rbind(
  data.frame(
    WaterDeficit = P_values,
    Region = "Non-urban",
    Marginal_Effect = marginal_non_urban,
    SE = se_marginal_non_urban,
    CI_Lower = ci_marginal_non_urban_lower,
    CI_Upper = ci_marginal_non_urban_upper
  ),
  data.frame(
    WaterDeficit = P_values,
    Region = "Urban",
    Marginal_Effect = marginal_urban,
    SE = se_marginal_urban,
    CI_Lower = ci_marginal_urban_lower,
    CI_Upper = ci_marginal_urban_upper
  )
)

urban_colors <- c(
  "Urban" = "#D62728", 
  "Non-urban" = "#1F77B4"
)
urban_colors_hist <- c(
  "in" = "#D62728", 
  "out" = "#1F77B4"
)


hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill = urban), binwidth =  5,  alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600)) +
  scale_y_continuous(limits = c(0.0, 0.007),
                     breaks = seq(0.0, 0.007, by = 0.002),
                     labels = scales::number_format(accuracy = 0.001)) +
  labs(x = "Water deficit (mm)", 
       y = "",
       subtitle = paste0("")) +
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.title.x = element_text(size = 15), 
    plot.margin = margin(t=-1.30,b=0.5,l=0.15,r=0.2, unit = "cm"), 
    axis.text.y = element_text(color = "black", size = 14),  
    axis.text.x = element_text(color = "black", size = 14),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist) 

print(hist_plot)

# Plot marginal effect curve
marginal_plot <- ggplot(marginal_df, aes(x = WaterDeficit, y = Marginal_Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha=0.3, color = NA
  ) +
  labs(
    x = "WaterDeficit",
    y = expression(Marginal~effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  scale_linetype_manual(values = c("Non-urban"="solid", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 1, "Urban" = 1), guide="none") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600)) +
  theme(
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 14),  
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),      
    axis.title.y = element_text(size = 15),
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 14, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) 

print(marginal_plot)

marginal_waterdeficit_NPP <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(marginal_waterdeficit_NPP)




###########################################################################
#                     FVC ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                         .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))

model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      Seasonalityindex_2 +
                      urban +
                      AverageTemperature +  
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Precipitation +  
                      SoilMoisture +
                      SoilPH +
                      Seasonalityindex:urban +
                      Seasonalityindex_2:urban +
                      GDP_per_capita_PPP + HDI + Population +
                      ImperviousSurface + HumanSettlement + CityArea
                    |  CityID + year | 0 | CityID, data = merged_dt_Seasonalityindex)  
summary(model_urban)

output_file <- "FVC_urban_in_out2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

coeffs <- coef(model_urban)
beta_urban <- coeffs["urbanin"]
beta1 <- coeffs["Seasonalityindex"]             
beta2 <- coeffs["Seasonalityindex_2"]             
beta1_city <- coeffs["Seasonalityindex:urbanin"]  
beta2_city <- coeffs["Seasonalityindex_2:urbanin"]
beta1_total <- beta1 + beta1_city
beta2_total <- beta2 + beta2_city
P_min <- quantile(merged_dt_Seasonalityindex$Seasonalityindex, 0.001, na.rm = TRUE)
P_max <- quantile(merged_dt_Seasonalityindex$Seasonalityindex, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)  
effect_non_urban <- beta1 * P_values + beta2 * (P_values^2)  
effect_urban <- beta_urban + beta1_total * P_values + beta2_total * (P_values^2)  
vcov_matrix <- vcov(model_urban)
cov_vars_non_urban <- c("Seasonalityindex", "Seasonalityindex_2")
vcov_non_urban <- vcov_matrix[cov_vars_non_urban, cov_vars_non_urban]
cov_vars_urban <- c("Seasonalityindex", "Seasonalityindex_2", "Seasonalityindex:urbanin", "Seasonalityindex_2:urbanin")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]

marginal_non_urban <- beta1 + 2 * beta2 * P_values  
marginal_urban <- beta1_total + 2 * beta2_total * P_values  
d_marg_beta1 <- 1 
d_marg_beta2 <- 2 * P_values  
var_marginal_non_urban <- (d_marg_beta1^2) * vcov_non_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_marg_beta2^2) * vcov_non_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_non_urban["Seasonalityindex", "Seasonalityindex_2"]
se_marginal_non_urban <- sqrt(var_marginal_non_urban)
ci_marginal_non_urban_lower <- marginal_non_urban - 1.96 * se_marginal_non_urban
ci_marginal_non_urban_upper <- marginal_non_urban + 1.96 * se_marginal_non_urban
d_marg_beta1_u <- 1  
d_marg_beta2_u <- 2 * P_values  
var_marginal_urban <- (d_marg_beta1_u^2) * vcov_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_marg_beta2_u^2) * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  (d_marg_beta1_u^2) * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex:urbanin"] +
  (d_marg_beta2_u^2) * vcov_urban["Seasonalityindex_2:urbanin", "Seasonalityindex_2:urbanin"] +
  2 * d_marg_beta1_u^2 * vcov_urban["Seasonalityindex", "Seasonalityindex:urbanin"] +  
  2 * d_marg_beta2_u^2 * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2:urbanin"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2:urbanin"] + 
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex_2:urbanin"]  
se_marginal_urban <- sqrt(var_marginal_urban)
ci_marginal_urban_lower <- marginal_urban - 1.96 * se_marginal_urban
ci_marginal_urban_upper <- marginal_urban + 1.96 * se_marginal_urban

marginal_df <- rbind(
  data.frame(
    Seasonalityindex = P_values,
    Region = "Non-urban",
    Marginal_Effect = marginal_non_urban,
    SE = se_marginal_non_urban,
    CI_Lower = ci_marginal_non_urban_lower,
    CI_Upper = ci_marginal_non_urban_upper
  ),
  data.frame(
    Seasonalityindex = P_values,
    Region = "Urban",
    Marginal_Effect = marginal_urban,
    SE = se_marginal_urban,
    CI_Lower = ci_marginal_urban_lower,
    CI_Upper = ci_marginal_urban_upper
  )
)

urban_colors <- c(
  "Urban" = "#D62728", 
  "Non-urban" = "#1F77B4"
)
urban_colors_hist <- c(
  "in" = "#D62728", 
  "out" = "#1F77B4"
)


hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill = urban), binwidth =  0.0001,  alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02)) +
  labs(x = "Rainfall seasonality",  
       y = "",
       subtitle = paste0("")) +
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.title.x = element_text(size = 15),  
    plot.margin = margin(t=-1.30,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 14),  
    axis.text.x = element_text(color = "black", size = 14),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist)

print(hist_plot)

# Plot marginal effect curve
marginal_plot <- ggplot(marginal_df, aes(x = Seasonalityindex, y = Marginal_Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha=0.3, color = NA
  ) +
  labs(
    x = "Seasonalityindex",
    y = "Marginal effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  scale_linetype_manual(values = c("Non-urban"="solid", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 1, "Urban" = 1), guide="none") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02)) +
  theme(
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 14),  
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),      
    axis.title.y = element_text(size = 15),
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 14, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors)

print(marginal_plot)

marginal_Seasonalityindex_FVC <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(marginal_Seasonalityindex_FVC)


###########################################################################
#                     NPP ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude',
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population',
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)),
                                                         .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))

model_urban <- felm(NPP ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      Seasonalityindex_2 +
                      urban +
                      AverageTemperature +  
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Precipitation + 
                      SoilMoisture +
                      SoilPH +
                      Seasonalityindex:urban +
                      Seasonalityindex_2:urban +
                      
                      GDP_per_capita_PPP + HDI + Population +
                      ImperviousSurface + HumanSettlement + CityArea
                    |  CityID + year | 0 | CityID, data = merged_dt_Seasonalityindex)  
summary(model_urban)

output_file <- "NPP_urban_in_out2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

coeffs <- coef(model_urban)
beta_urban <- coeffs["urbanin"]
beta1 <- coeffs["Seasonalityindex"]             
beta2 <- coeffs["Seasonalityindex_2"]             
beta1_city <- coeffs["Seasonalityindex:urbanin"]  
beta2_city <- coeffs["Seasonalityindex_2:urbanin"]
beta1_total <- beta1 + beta1_city
beta2_total <- beta2 + beta2_city
P_min <- quantile(merged_dt_Seasonalityindex$Seasonalityindex, 0.001, na.rm = TRUE)
P_max <- quantile(merged_dt_Seasonalityindex$Seasonalityindex, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200) 
effect_non_urban <- beta1 * P_values + beta2 * (P_values^2)  
effect_urban <- beta_urban + beta1_total * P_values + beta2_total * (P_values^2)  
vcov_matrix <- vcov(model_urban)
cov_vars_non_urban <- c("Seasonalityindex", "Seasonalityindex_2")
vcov_non_urban <- vcov_matrix[cov_vars_non_urban, cov_vars_non_urban]
cov_vars_urban <- c("Seasonalityindex", "Seasonalityindex_2", "Seasonalityindex:urbanin", "Seasonalityindex_2:urbanin")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]

marginal_non_urban <- beta1 + 2 * beta2 * P_values  
marginal_urban <- beta1_total + 2 * beta2_total * P_values  
d_marg_beta1 <- 1 
d_marg_beta2 <- 2 * P_values  
var_marginal_non_urban <- (d_marg_beta1^2) * vcov_non_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_marg_beta2^2) * vcov_non_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_non_urban["Seasonalityindex", "Seasonalityindex_2"]
se_marginal_non_urban <- sqrt(var_marginal_non_urban)
ci_marginal_non_urban_lower <- marginal_non_urban - 1.96 * se_marginal_non_urban
ci_marginal_non_urban_upper <- marginal_non_urban + 1.96 * se_marginal_non_urban
d_marg_beta1_u <- 1  
d_marg_beta2_u <- 2 * P_values  
var_marginal_urban <- (d_marg_beta1_u^2) * vcov_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_marg_beta2_u^2) * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  (d_marg_beta1_u^2) * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex:urbanin"] +
  (d_marg_beta2_u^2) * vcov_urban["Seasonalityindex_2:urbanin", "Seasonalityindex_2:urbanin"] +
  2 * d_marg_beta1_u^2 * vcov_urban["Seasonalityindex", "Seasonalityindex:urbanin"] +  
  2 * d_marg_beta2_u^2 * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2:urbanin"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2:urbanin"] + 
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex_2:urbanin"]  
se_marginal_urban <- sqrt(var_marginal_urban)
ci_marginal_urban_lower <- marginal_urban - 1.96 * se_marginal_urban
ci_marginal_urban_upper <- marginal_urban + 1.96 * se_marginal_urban

marginal_df <- rbind(
  data.frame(
    Seasonalityindex = P_values,
    Region = "Non-urban",
    Marginal_Effect = marginal_non_urban,
    SE = se_marginal_non_urban,
    CI_Lower = ci_marginal_non_urban_lower,
    CI_Upper = ci_marginal_non_urban_upper
  ),
  data.frame(
    Seasonalityindex = P_values,
    Region = "Urban",
    Marginal_Effect = marginal_urban,
    SE = se_marginal_urban,
    CI_Lower = ci_marginal_urban_lower,
    CI_Upper = ci_marginal_urban_upper
  )
)

hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill = urban), binwidth =  0.0001,  alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02)) +
  labs(x = "Rainfall seasonality",  
       y = "",
       subtitle = paste0("")) +
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.title.x = element_text(size = 15),  
    plot.margin = margin(t=-1.30,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 14),  
    axis.text.x = element_text(color = "black", size = 14),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist)

print(hist_plot)

# Plot marginal effect curve
marginal_plot <- ggplot(marginal_df, aes(x = Seasonalityindex, y = Marginal_Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha=0.3, color = NA
  ) +
  labs(
    x = "Seasonalityindex",
    y = expression(Marginal~effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  scale_linetype_manual(values = c("Non-urban"="solid", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 1, "Urban" = 1), guide="none") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02)) +
  theme(
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 14),  
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),    
    axis.title.y = element_text(size = 15),
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 14, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) 

print(marginal_plot)

marginal_Seasonalityindex_NPP <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(marginal_Seasonalityindex_NPP)



crop_bottom <- function(plot, crop_ratio = 0.15) {
  ggdraw() +
    draw_plot(plot, x = 0, y = crop_ratio, width = 1, height = 1 - crop_ratio) +
    theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"))
}

p1_crop  <- crop_bottom(marginal_waterdeficit_FVC, crop_ratio = -0.02)
p2_crop  <- crop_bottom(marginal_waterdeficit_NPP, crop_ratio = -0.02)
p3_crop  <- crop_bottom(marginal_Seasonalityindex_FVC, crop_ratio = -0.05)
p4_crop  <- crop_bottom(marginal_Seasonalityindex_NPP, crop_ratio = -0.05)

fig_2x2 <- plot_grid(
  p1_crop, p2_crop,
  p3_crop, p4_crop,
  nrow = 2,
  ncol = 2,
  align = "hv",      
  axis = "tblr",      
  rel_widths = c(1,1),
  rel_heights = c(1,1),
  labels = c("a", "b", "c", "d"),
  label_size = 18,
  label_fontfamily = "Times",
  label_fontface = "bold",
  label_x = rep(0.08,4),  
  label_y = rep(1.04,4),  
  hjust = 0,
  vjust = 1
)

fig_2x2 <- ggdraw(fig_2x2) +
  theme(plot.margin = margin(5, 5, 5, 5, unit = "mm"))

print(fig_2x2)

ggsave(
  "./FigS3_marginal_2x2.pdf",
  plot = fig_2x2,
  width = 300, height = 260, units = "mm", dpi=300
)

