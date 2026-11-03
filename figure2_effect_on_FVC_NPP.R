

# Effect on FVC and NPP


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




# ============================================================
#    FVC × Water deficit
# ============================================================

setwd("XXX") # please set your folder path
merged_dt <- fread("Dataset_matched_buffer10km.csv")

merged_dt[, urban_plot := fcase(
  urban == "in",  "Urban",
  urban == "out", "Non-urban"
)]
merged_dt[, urban_plot := factor(urban_plot, levels = c("Urban", "Non-urban"))]
merged_dt_season <- merged_dt[Precipitation < 1000]
agg_cols <- c("FVC", "NPP", "WaterDeficit", "Seasonalityindex")
city_dt <- merged_dt[,
                     lapply(.SD, mean, na.rm = TRUE),
                     by = .(CityID, urban_plot),
                     .SDcols = agg_cols
]
city_dt_season <- merged_dt_season[,
                                   lapply(.SD, mean, na.rm = TRUE),
                                   by = .(CityID, urban_plot),
                                   .SDcols = agg_cols
]

pal_color <- c("Urban" = "#D62728", "Non-urban" = "#1F77B4")
pal_linetype <- c("Urban" = "solid",   "Non-urban" = "dashed")

p_a <- ggplot(city_dt, aes(x = WaterDeficit, y = FVC)) +
  geom_scattermore(aes(color = urban_plot),
                   pointsize = 5,    
                   pixels = c(1500,1125), 
                   interpolate = TRUE, 
  ) +
  geom_smooth(aes(linetype = urban_plot),
              color = "black", fill = NA,
              method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.6) +
  scale_color_manual(values = pal_color) +
  scale_linetype_manual(values = pal_linetype) +
  labs(x = "Water deficit (mm)", y = "FVC",
       color = "", linetype = "") +
  scale_x_continuous(limits = c(0, 2500), 
                     breaks = seq(0, 2400, by =600),
                     expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(0, 1.05),
                     breaks = seq(0, 1, by = 0.25),
                     labels = scales::number_format(accuracy = 0.01),
                     expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.text = element_text(color = "black", size = 12),  
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.title = element_text(size = 12),  
    legend.position = c(0.6, 1.0), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 1.5, keyheight = 0.4))

print(p_a)




# ============================================================
#    NPP × Water deficit
# ============================================================
p_b <- ggplot(city_dt, aes(x = WaterDeficit, y = NPP)) +
  geom_scattermore(aes(color = urban_plot),
                   pointsize = 5,    
                   pixels = c(1500,1125), 
                   interpolate = TRUE, 
  ) +
  geom_smooth(aes(linetype = urban_plot),
              color = "black", fill = NA,
              method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.6) +
  scale_color_manual(values = pal_color) +
  scale_linetype_manual(values = pal_linetype) +
  labs(x = "Water deficit (mm)",
       y = expression("NPP (g C m"^-2*" yr"^-1*")"),
       color = "", linetype = "") +
  scale_x_continuous(limits = c(0, 2500), 
                     breaks = seq(0, 2400, by =600),
                     expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(0, 1700),
                     breaks = seq(0, 1600, by = 400),
                     expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text = element_text(color = "black", size = 12),  
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.title = element_text(size = 12),  
    legend.position = "none", 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 1.5, keyheight = 0.4))

print(p_b)



# ============================================================
#    FVC × Rainfall seasonality
# ============================================================
p_c <- ggplot(city_dt_season, aes(x = Seasonalityindex, y = FVC)) +
  geom_scattermore(aes(color = urban_plot),
                   pointsize = 5,    
                   pixels = c(1500,1125), 
                   interpolate = TRUE, 
  ) +
  geom_smooth(aes(linetype = urban_plot),
              color = "black", fill = NA,
              method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.6) +
  scale_color_manual(values = pal_color) +
  scale_linetype_manual(values = pal_linetype) +
  labs(x = "Rainfall seasonality", y = "FVC",
       color = "", linetype = "") +
  scale_x_continuous(limits = c(0, 0.105), 
                     breaks = seq(0, 0.1, by =0.02),
                     expand = c(0.002, 0.002)) +
  scale_y_continuous(limits = c(0, 1.05),
                     breaks = seq(0, 1, by = 0.25),
                     labels = scales::number_format(accuracy = 0.01),
                     expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text = element_text(color = "black", size = 12),  
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.title = element_text(size = 12),  
    legend.position = "none", 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 1.5, keyheight = 0.4))

print(p_c)


# ============================================================
#    NPP × Rainfall seasonality
# ============================================================
p_d <- ggplot(city_dt_season, aes(x = Seasonalityindex, y = NPP)) +
  geom_scattermore(aes(color = urban_plot),
                   pointsize = 5,    
                   pixels = c(1500,1125), 
                   interpolate = TRUE, 
  ) +
  geom_smooth(aes(linetype = urban_plot),
              color = "black", fill = NA,
              method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.6) +
  scale_color_manual(values = pal_color) +
  scale_linetype_manual(values = pal_linetype) +
  labs(x = "Rainfall seasonality",
       y = expression("NPP (g C m"^-2*" yr"^-1*")"),
       color = "", linetype = "") +
  scale_x_continuous(limits = c(0, 0.105), 
                     breaks = seq(0, 0.1, by =0.02),
                     expand = c(0.002, 0.002)) +
  scale_y_continuous(limits = c(0, 1700),
                     breaks = seq(0, 1600, by = 400),
                     expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text = element_text(color = "black", size = 12),  
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.title = element_text(size = 12),  
    legend.position = "none", 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 1.5, keyheight = 0.4))

print(p_d)




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
d_beta1 <- P_values  
d_beta2 <- P_values^2 
var_effect_non_urban <- (d_beta1^2) * vcov_non_urban["WaterDeficit", "WaterDeficit"] +
  (d_beta2^2) * vcov_non_urban["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_beta1 * d_beta2 * vcov_non_urban["WaterDeficit", "WaterDeficit_2"]
se_effect_non_urban <- sqrt(var_effect_non_urban)
ci_effect_non_urban_lower <- effect_non_urban - 1.96 * se_effect_non_urban
ci_effect_non_urban_upper <- effect_non_urban + 1.96 * se_effect_non_urban
cov_vars_urban <- c("WaterDeficit", "WaterDeficit_2", "WaterDeficit:urbanin", "WaterDeficit_2:urbanin")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]
d_beta1_u <- P_values  
d_beta2_u <- P_values^2  
var_effect_urban <- (d_beta1_u^2) * vcov_urban["WaterDeficit", "WaterDeficit"] +
  (d_beta2_u^2) * vcov_urban["WaterDeficit_2", "WaterDeficit_2"] +
  (d_beta1_u^2) * vcov_urban["WaterDeficit:urbanin", "WaterDeficit:urbanin"] +
  (d_beta2_u^2) * vcov_urban["WaterDeficit_2:urbanin", "WaterDeficit_2:urbanin"] +
  2 * d_beta1_u^2 * vcov_urban["WaterDeficit", "WaterDeficit:urbanin"] +  
  2 * d_beta2_u^2 * vcov_urban["WaterDeficit_2", "WaterDeficit_2:urbanin"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2:urbanin"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit:urbanin", "WaterDeficit_2"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit:urbanin", "WaterDeficit_2:urbanin"]  
se_effect_urban <- sqrt(var_effect_urban)
ci_effect_urban_lower <- effect_urban - 1.96 * se_effect_urban
ci_effect_urban_upper <- effect_urban + 1.96 * se_effect_urban

effect_df <- rbind(
  data.frame(
    WaterDeficit = P_values,
    Region = "Non-urban",
    Effect = effect_non_urban,
    SE = se_effect_non_urban,
    CI_Lower = ci_effect_non_urban_lower,
    CI_Upper = ci_effect_non_urban_upper
  ),
  data.frame(
    WaterDeficit = P_values,
    Region = "Urban",
    Effect = effect_urban,
    SE = se_effect_urban,
    CI_Lower = ci_effect_urban_lower,
    CI_Upper = ci_effect_urban_upper
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


effect_plot <- ggplot(effect_df, aes(x = WaterDeficit, y = Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha=0.3, color = NA
  ) +
  scale_linetype_manual(values = c("Non-urban"="solid", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 1, "Urban" = 1), guide="none") +
  theme_classic() + 
  labs(
    x = "Water deficit",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.title.y = element_text(size = 13),  
    legend.position = c(0.3, 0.3), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) 

print(effect_plot)

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
    axis.title.x = element_text(size = 13), 
    plot.margin = margin(t=-1.30,b=0.5,l=0.15,r=0.2, unit = "cm"), 
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_text(color = "black", size = 12),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist) 

print(hist_plot)

effect_waterdeficit_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_waterdeficit_FVC)



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
d_beta1 <- P_values  
d_beta2 <- P_values^2 
var_effect_non_urban <- (d_beta1^2) * vcov_non_urban["WaterDeficit", "WaterDeficit"] +
  (d_beta2^2) * vcov_non_urban["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_beta1 * d_beta2 * vcov_non_urban["WaterDeficit", "WaterDeficit_2"]
se_effect_non_urban <- sqrt(var_effect_non_urban)
ci_effect_non_urban_lower <- effect_non_urban - 1.96 * se_effect_non_urban
ci_effect_non_urban_upper <- effect_non_urban + 1.96 * se_effect_non_urban
cov_vars_urban <- c("WaterDeficit", "WaterDeficit_2", "WaterDeficit:urbanin", "WaterDeficit_2:urbanin")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]
d_beta1_u <- P_values  
d_beta2_u <- P_values^2  
var_effect_urban <- (d_beta1_u^2) * vcov_urban["WaterDeficit", "WaterDeficit"] +
  (d_beta2_u^2) * vcov_urban["WaterDeficit_2", "WaterDeficit_2"] +
  (d_beta1_u^2) * vcov_urban["WaterDeficit:urbanin", "WaterDeficit:urbanin"] +
  (d_beta2_u^2) * vcov_urban["WaterDeficit_2:urbanin", "WaterDeficit_2:urbanin"] +
  2 * d_beta1_u^2 * vcov_urban["WaterDeficit", "WaterDeficit:urbanin"] +  
  2 * d_beta2_u^2 * vcov_urban["WaterDeficit_2", "WaterDeficit_2:urbanin"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2:urbanin"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit:urbanin", "WaterDeficit_2"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit:urbanin", "WaterDeficit_2:urbanin"]  
se_effect_urban <- sqrt(var_effect_urban)
ci_effect_urban_lower <- effect_urban - 1.96 * se_effect_urban
ci_effect_urban_upper <- effect_urban + 1.96 * se_effect_urban

effect_df <- rbind(
  data.frame(
    WaterDeficit = P_values,
    Region = "Non-urban",
    Effect = effect_non_urban,
    SE = se_effect_non_urban,
    CI_Lower = ci_effect_non_urban_lower,
    CI_Upper = ci_effect_non_urban_upper
  ),
  data.frame(
    WaterDeficit = P_values,
    Region = "Urban",
    Effect = effect_urban,
    SE = se_effect_urban,
    CI_Lower = ci_effect_urban_lower,
    CI_Upper = ci_effect_urban_upper
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

# Plot effect curve
effect_plot <- ggplot(effect_df, aes(x = WaterDeficit, y = Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha=0.3, color = NA
  ) +
  scale_linetype_manual(values = c("Non-urban"="solid", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 1, "Urban" = 1), guide="none") +
  theme_classic() + 
  labs(
    x = "Water deficit",
    y = expression(Effect~on~NPP~(g~C~m^-2~yr^-1)),
    y = "",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(color = "black", size = 0.5),   
    axis.title.y = element_text(size = 13),  
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 12, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) 

print(effect_plot)

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
    axis.title.x = element_text(size = 13), 
    plot.margin = margin(t=-1.30,b=0.5,l=0.15,r=0.2, unit = "cm"), 
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_text(color = "black", size = 12),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist) 

print(hist_plot)

effect_waterdeficit_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_waterdeficit_NPP)




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
d_beta1 <- P_values  
d_beta2 <- P_values^2  
var_effect_non_urban <- (d_beta1^2) * vcov_non_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_beta2^2) * vcov_non_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_beta1 * d_beta2 * vcov_non_urban["Seasonalityindex", "Seasonalityindex_2"]
se_effect_non_urban <- sqrt(var_effect_non_urban)
ci_effect_non_urban_lower <- effect_non_urban - 1.96 * se_effect_non_urban
ci_effect_non_urban_upper <- effect_non_urban + 1.96 * se_effect_non_urban
cov_vars_urban <- c("Seasonalityindex", "Seasonalityindex_2", "Seasonalityindex:urbanin", "Seasonalityindex_2:urbanin")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]
d_beta1_u <- P_values  
d_beta2_u <- P_values^2  
var_effect_urban <- (d_beta1_u^2) * vcov_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_beta2_u^2) * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  (d_beta1_u^2) * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex:urbanin"] +
  (d_beta2_u^2) * vcov_urban["Seasonalityindex_2:urbanin", "Seasonalityindex_2:urbanin"] +
  2 * d_beta1_u^2 * vcov_urban["Seasonalityindex", "Seasonalityindex:urbanin"] +  
  2 * d_beta2_u^2 * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2:urbanin"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2"] + 
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2:urbanin"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex_2"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex_2:urbanin"]  
se_effect_urban <- sqrt(var_effect_urban)
ci_effect_urban_lower <- effect_urban - 1.96 * se_effect_urban
ci_effect_urban_upper <- effect_urban + 1.96 * se_effect_urban

effect_df <- rbind(
  data.frame(
    Seasonalityindex = P_values,
    Region = "Non-urban",
    Effect = effect_non_urban,
    SE = se_effect_non_urban,
    CI_Lower = ci_effect_non_urban_lower,
    CI_Upper = ci_effect_non_urban_upper
  ),
  data.frame(
    Seasonalityindex = P_values,
    Region = "Urban",
    Effect = effect_urban,
    SE = se_effect_urban,
    CI_Lower = ci_effect_urban_lower,
    CI_Upper = ci_effect_urban_upper
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


# Plot effect curve
effect_plot <- ggplot(effect_df, aes(x = Seasonalityindex, y = Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha=0.3, color = NA
  ) +
  scale_linetype_manual(values = c("Non-urban"="solid", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 1, "Urban" = 1), guide="none") +
  theme_classic() + 
  labs(
    x = "Seasonality index",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(color = "black", size = 0.5),   
    axis.title.y = element_text(size = 13),  
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 12, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) 

print(effect_plot)

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
    axis.title.x = element_text(size = 13),  
    plot.margin = margin(t=-1.30,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_text(color = "black", size = 12),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist)

print(hist_plot)

effect_Seasonalityindex_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_Seasonalityindex_FVC)


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
d_beta1 <- P_values  
d_beta2 <- P_values^2  
var_effect_non_urban <- (d_beta1^2) * vcov_non_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_beta2^2) * vcov_non_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_beta1 * d_beta2 * vcov_non_urban["Seasonalityindex", "Seasonalityindex_2"]
se_effect_non_urban <- sqrt(var_effect_non_urban)
ci_effect_non_urban_lower <- effect_non_urban - 1.96 * se_effect_non_urban
ci_effect_non_urban_upper <- effect_non_urban + 1.96 * se_effect_non_urban
cov_vars_urban <- c("Seasonalityindex", "Seasonalityindex_2", "Seasonalityindex:urbanin", "Seasonalityindex_2:urbanin")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]
d_beta1_u <- P_values  
d_beta2_u <- P_values^2  
var_effect_urban <- (d_beta1_u^2) * vcov_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_beta2_u^2) * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  (d_beta1_u^2) * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex:urbanin"] +
  (d_beta2_u^2) * vcov_urban["Seasonalityindex_2:urbanin", "Seasonalityindex_2:urbanin"] +
  2 * d_beta1_u^2 * vcov_urban["Seasonalityindex", "Seasonalityindex:urbanin"] +  
  2 * d_beta2_u^2 * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2:urbanin"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2"] + 
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2:urbanin"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex_2"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex_2:urbanin"]  
se_effect_urban <- sqrt(var_effect_urban)
ci_effect_urban_lower <- effect_urban - 1.96 * se_effect_urban
ci_effect_urban_upper <- effect_urban + 1.96 * se_effect_urban

effect_df <- rbind(
  data.frame(
    Seasonalityindex = P_values,
    Region = "Non-urban",
    Effect = effect_non_urban,
    SE = se_effect_non_urban,
    CI_Lower = ci_effect_non_urban_lower,
    CI_Upper = ci_effect_non_urban_upper
  ),
  data.frame(
    Seasonalityindex = P_values,
    Region = "Urban",
    Effect = effect_urban,
    SE = se_effect_urban,
    CI_Lower = ci_effect_urban_lower,
    CI_Upper = ci_effect_urban_upper
  )
)


# Plot effect curve
effect_plot <- ggplot(effect_df, aes(x = Seasonalityindex, y = Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha=0.3, color = NA
  ) +
  scale_linetype_manual(values = c("Non-urban"="solid", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 1, "Urban" = 1), guide="none") +
  theme_classic() + 
  labs(
    x = "Seasonality index",
    y = expression(Effect~on~NPP~(g~C~m^-2~yr^-1)),
    y = "",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(color = "black", size = 0.5),   
    axis.title.y = element_text(size = 13),  
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 12, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) 

print(effect_plot)

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
    axis.title.x = element_text(size = 13),  
    plot.margin = margin(t=-1.30,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_text(color = "black", size = 12),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist)

print(hist_plot)

effect_Seasonalityindex_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_Seasonalityindex_NPP)



strip_margin <- function(p) {
  p + theme(plot.margin       = margin(0, 0, 0, 0, unit = "pt"),
            axis.ticks.length = unit(2, "pt"))
}

p1 <- strip_margin(p_a)                          # a
p2 <- strip_margin(p_b)                          # b
p3 <- strip_margin(p_c)                          # c
p4 <- strip_margin(p_d)                          # d
p5 <- strip_margin(effect_waterdeficit_FVC)      # e
p6 <- strip_margin(effect_waterdeficit_NPP)      # f
p7 <- strip_margin(effect_Seasonalityindex_FVC)  # g
p8 <- strip_margin(effect_Seasonalityindex_NPP)  # h

row_h <- c(1.0, 1.0)                 
col_w <- c(1.0, 1.0, 1.0, 1.0)    
make_row <- function(row_list, rel_widths = col_w) {
  plot_grid(
    plotlist   = row_list,
    nrow       = 1, ncol = 4,       
    align      = "h",
    axis       = "tb",
    rel_widths = rel_widths
  )
}
row1 <- make_row(list(p1, p2, p3, p4))   
row2 <- make_row(list(p5, p6, p7, p8))   

left_pad   <- 0.025
right_pad  <- 0.015
top_pad    <- 0.04
bottom_pad <- 0.00
row_gap    <- 0.008         
total_h <- 1 - top_pad - bottom_pad
h_norm  <- row_h / sum(row_h)
row1_h <- total_h * h_norm[1]
row2_h <- total_h * h_norm[2]
y_row1 <- bottom_pad + row2_h + row_gap   
y_row2 <- bottom_pad                      
total_w <- 1 - left_pad - right_pad
w_norm  <- col_w / sum(col_w)
col_gap <- 0.01
x_col1 <- left_pad
x_col2 <- x_col1 + total_w * w_norm[1] + col_gap
x_col3 <- x_col2 + total_w * w_norm[2] + col_gap
x_col4 <- x_col3 + total_w * w_norm[3] + col_gap
fig_2x4 <- ggdraw() +
  draw_plot(row2, x = left_pad, y = y_row2, width = total_w, height = row2_h) +
  draw_plot(row1, x = left_pad, y = y_row1, width = total_w, height = row1_h)
tag_x <- c(x_col1+0.01, x_col2+0.01, x_col3-0.01, x_col4-0.01) + 0.008
tag_y_row1 <- y_row1 + row1_h + 0.02
tag_y_row2 <- y_row2 + row2_h + 0.01
tag_x_all <- c(tag_x, tag_x)                       
tag_y_all <- c(rep(tag_y_row1, 4), rep(tag_y_row2, 4))
fig_2x4 <- fig_2x4 +
  draw_plot_label(
    label    = c("a","b","c","d","e","f","g","h"),
    x        = tag_x_all,
    y        = tag_y_all,
    size     = 14,
    fontface = "bold",
    family   = "Times",
    hjust    = 0, vjust = 1
  )

panel_w_mm <- 70; panel_h_mm <- 50
fig_w_mm <- panel_w_mm * 4 * 1.30 + 27
fig_h_mm <- panel_h_mm * (1 + 1.18) * 1.15 + 20

ggsave("./Fig2_2x4.pdf",
       plot   = fig_2x4,
       width  = fig_w_mm,
       height = fig_h_mm+50,
       units  = "mm",
       dpi    = 300)


