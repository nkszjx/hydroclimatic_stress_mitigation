
# The urban environment reshapes the nonlinear effects of water deficit on FVC and NPP

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


###########################################################################
#                   Two-way interactive quadratic regression
###########################################################################

setwd("XXX") # please set your folder path
merged_dt <- fread("Dataset_matched_buffer10km.csv")

# Absolute latitude to unify hemispheres
merged_dt$Latitude <- abs(merged_dt$Latitude)

# Multicollinearity test with VIF
lm_model <- lm(FVC ~  Precipitation + WaterDeficit + Seasonalityindex + AverageTemperature + TemperatureRange +  
                 WindSpeed + SoilMoisture +  SoilPH + 
                 Elevation + Latitude + 
                 GDP_per_capita_PPP + HDI + Population + 
                 ImperviousSurface + HumanSettlement + CityArea, 
               data = merged_dt)
vif_values <- vif(lm_model)
print(vif_values)

# Create quadratic terms for nonlinear effect estimation
merged_dt$WaterDeficit_2 <- merged_dt$WaterDeficit * merged_dt$WaterDeficit
merged_dt$Seasonalityindex_2 <- merged_dt$Seasonalityindex * merged_dt$Seasonalityindex



###########################################################################
#                     Part 1: FVC Regression & Visualization
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

# Compute effects, marginal effects and confidence intervals
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
  "Urban" = "#FF7A3C", 
  "Non-urban" = "#00A850"
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
    alpha=0.1, color = NA
  ) +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() + 
  labs(
    x = "Water deficit",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600), 
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
      keywidth = 0.6,
      keyheight = 0.5,
      override.aes = list(
        linetype = c("Non-urban"="dashed", "Urban"="solid"),
        linewidth = c("Non-urban" = 0.25, "Urban" = 0.4)
      )
    )
  )

  

print(effect_plot)

# Plot distribution histogram
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
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3), 
    axis.title.x = element_text(size = 6), 
    plot.margin = margin(t=-1.30,b=0.5,l=0.15,r=0.2, unit = "cm"), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_text(color = "black", size = 6),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist) 

print(hist_plot)

effect_waterdeficit_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_waterdeficit_FVC)




# Plot marginal effect curve
marginal_plot <- ggplot(marginal_df, aes(x = WaterDeficit, y = Marginal_Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha=0.1, color = NA
  ) +
  labs(
    x = "WaterDeficit",
    y = "Marginal effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600)) +
  scale_y_continuous(limits = c(-0.0004, 0.0002),
                     breaks = seq(-0.0004, 0.0002, by = 0.0002),
                     labels = scales::number_format(accuracy = 0.0001)) +
  
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.7, 0.9), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) +
  guides(
    color = guide_legend(
      keywidth = 0.6,
      keyheight = 0.5,
      override.aes = list(
        linetype = c("Non-urban"="dashed", "Urban"="solid"),
        linewidth = c("Non-urban" = 0.25, "Urban" = 0.4)
      )
    )
  )


print(marginal_plot)

marginal_waterdeficit_FVC <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))

print(marginal_waterdeficit_FVC)



# Plot effect curve --- legend
effect_plot_legend <- ggplot(effect_df, aes(x = WaterDeficit, y = Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha=0.1, color = NA
  ) +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), 
                        labels = c("Non-urban" = "Non-urban effect curve with 95% CI", 
                                   "Urban" = "Urban effect curve with 95% CI"),
                        guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() + 
  labs(
    x = "Water deficit",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  # geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(name = "Effect on FVC", limits = c(-0.4, 0.6),
                     breaks = seq(-0.4, 0.58, by = 0.2),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.line.y.left = element_blank(),
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.45, 1.03), 
    # legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 5, unit = "pt")),   
    legend.key.spacing.x = unit(30, "pt"),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = c(
    "Urban" = "gray50", 
    "Non-urban" = "gray50"
  ),
  labels = c("Non-urban" = "Non-urban effect curve with 95% CI", 
             "Urban" = "Urban effect curve with 95% CI")) +
  scale_fill_manual(values = c(
    "Urban" = "gray50", 
    "Non-urban" = "gray50"
  ),
  labels = c("Non-urban" = "Non-urban effect curve with 95% CI", 
             "Urban" = "Urban effect curve with 95% CI")) +
  guides(
    color = guide_legend(
      nrow = 1, ncol = 2,
      keywidth = 1,
      keyheight = 0.5,
      override.aes = list(
        linetype = c("Non-urban"="dashed", "Urban"="solid"),
        linewidth = c("Non-urban" = 0.25, "Urban" = 0.4)
      )
    )
  )


print(effect_plot_legend)






# Plot marginal effect curve --- legend
marginal_plot_legend <- ggplot(effect_df, aes(x = WaterDeficit, y = Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha=0.1, color = NA
  ) +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), 
                        labels = c("Non-urban" = "Non-urban marginal effect curve with 95% CI", 
                                   "Urban" = "Urban marginal effect curve with 95% CI"),
                        guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() + 
  labs(
    x = "Water deficit",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  # geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(name = "Effect on FVC", limits = c(-0.4, 0.6),
                     breaks = seq(-0.4, 0.58, by = 0.2),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.line.y.left = element_blank(),
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.45, 1.03), 
    # legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 5, unit = "pt")),   
    legend.key.spacing.x = unit(30, "pt"),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = c(
    "Urban" = "gray50", 
    "Non-urban" = "gray50"
  ),
  labels = c("Non-urban" = "Non-urban marginal effect curve with 95% CI", 
             "Urban" = "Urban marginal effect curve with 95% CI")) +
  scale_fill_manual(values = c(
    "Urban" = "gray50", 
    "Non-urban" = "gray50"
  ),
  labels = c("Non-urban" = "Non-urban marginal effect curve with 95% CI", 
             "Urban" = "Urban marginal effect curve with 95% CI")) +
  guides(
    color = guide_legend(
      nrow = 1, ncol = 2,
      keywidth = 1,
      keyheight = 0.5,
      override.aes = list(
        linetype = c("Non-urban"="dashed", "Urban"="solid"),
        linewidth = c("Non-urban" = 0.25, "Urban" = 0.4)
      )
    )
  )


print(marginal_plot_legend)


###########################################################################
#                     Part 2: NPP Regression & Visualization
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

# Compute effects, marginal effects and confidence intervals
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
  "Urban" = "#FF7A3C", 
  "Non-urban" = "#00A850"
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
    alpha=0.1, color = NA
  ) +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() + 
  labs(
    x = "Water deficit",
    y = expression(Effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600), 
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
      keywidth = 0.6,
      keyheight = 0.5,
      override.aes = list(
        linetype = c("Non-urban"="dashed", "Urban"="solid"),
        linewidth = c("Non-urban" = 0.25, "Urban" = 0.4)
      )
    )
  )


print(effect_plot)

# Plot distribution histogram
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
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3), 
    axis.title.x = element_text(size = 6), 
    plot.margin = margin(t=-1.30,b=0.5,l=0.15,r=0.2, unit = "cm"), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_text(color = "black", size = 6),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist) 

print(hist_plot)


effect_waterdeficit_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_waterdeficit_NPP)



# Plot marginal effect curve
marginal_plot <- ggplot(marginal_df, aes(x = WaterDeficit, y = Marginal_Effect, color = Region, linetype = Region, linewidth = Region)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha=0.1, color = NA
  ) +
  labs(
    x = "WaterDeficit",
    y = expression(Marginal~effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.7, 0.25), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) +
  guides(
    color = guide_legend(
      keywidth = 0.6,
      keyheight = 0.5,
      override.aes = list(
        linetype = c("Non-urban"="dashed", "Urban"="solid"),
        linewidth = c("Non-urban" = 0.25, "Urban" = 0.4)
      )
    )
  )


print(marginal_plot)

marginal_waterdeficit_NPP <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))

print(marginal_waterdeficit_NPP)





###########################################################################
#                   Three-way interactive quadratic regression
###########################################################################


###########################################################################
#                     Part 1: FVC Regression & Visualization
###########################################################################

# Climate zone
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]


merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))
merged_dt_WaterDeficit$ClimateZone <- factor(merged_dt_WaterDeficit$ClimateZone, levels = c("Temperate","Tropical","Cold","Arid"))


model_urban <- felm(FVC ~ 1 + 
                      WaterDeficit +                   
                      WaterDeficit_2 +  
                      urban + 
                      WaterDeficit:ClimateZone + 
                      WaterDeficit_2:ClimateZone + 
                      urban:ClimateZone +  
                      WaterDeficit:urban + 
                      WaterDeficit_2:urban + 
                      WaterDeficit:urban:ClimateZone +  
                      WaterDeficit_2:urban:ClimateZone + 
                      Seasonalityindex + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    | CityID + year | 0 | CityID, 
                    data = merged_dt_WaterDeficit)

summary(model_urban)

output_file <- "FVC_urban_in_out_climate.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)



# Compute effects, marginal effects and confidence intervals
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "WaterDeficit",
  "WaterDeficit_2",
  "urbanin",
  "WaterDeficit:ClimateZoneTropical",
  "WaterDeficit:ClimateZoneCold",
  "WaterDeficit:ClimateZoneArid",
  "WaterDeficit_2:ClimateZoneTropical",
  "WaterDeficit_2:ClimateZoneCold",
  "WaterDeficit_2:ClimateZoneArid",
  "urbanin:ClimateZoneTropical",
  "urbanin:ClimateZoneCold",
  "urbanin:ClimateZoneArid",
  "WaterDeficit:urbanin",
  "WaterDeficit_2:urbanin",
  "WaterDeficit:urbanin:ClimateZoneTropical",
  "WaterDeficit:urbanin:ClimateZoneCold",
  "WaterDeficit:urbanin:ClimateZoneArid",
  "WaterDeficit_2:urbanin:ClimateZoneTropical",
  "WaterDeficit_2:urbanin:ClimateZoneCold",
  "WaterDeficit_2:urbanin:ClimateZoneArid"
)

b <- beta_all[v_names]

climate_colors <- c(
  "Temperate" = "#2ecc71",
  "Tropical"  = "#e74c3c",
  "Arid"      = "#95a5a6",
  "Cold"      = "#3498db"
)

wd_seq <- seq(0, 2400, length.out = 200)
climates <- c("Temperate","Tropical","Cold","Arid")
plot_df <- data.frame()
me_plot_df <- data.frame()

for(cz in climates){
  for(wd in wd_seq){
    if(cz=="Temperate"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,0,
        0,0,0,
        0,0,0,
        0,0,0,0,0,
        0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,0,
        0,0,0,
        0,0,0,
        0,0,0,0,0,
        0,0,0
      )
    }else if(cz=="Tropical"){
      g_non <- c(
        wd, wd^2, 0,
        wd,0,0,
        wd^2,0,0,
        0,0,0,
        0,0,0,0,0,
        0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        1,0,0,
        2*wd,0,0,
        0,0,0,
        0,0,0,0,0,
        0,0,0
      )
    }else if(cz=="Cold"){
      g_non <- c(
        wd, wd^2, 0,
        0,wd,0,
        0,wd^2,0,
        0,0,0,
        0,0,0,0,0,
        0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,1,0,
        0,2*wd,0,
        0,0,0,
        0,0,0,0,0,
        0,0,0
      )
    }else if(cz=="Arid"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,wd,
        0,0,wd^2,
        0,0,0,
        0,0,0,0,0,
        0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,1,
        0,0,2*wd,
        0,0,0,
        0,0,0,0,0,
        0,0,0
      )
    }
    
    pred_non <- as.numeric(crossprod(g_non, b))
    se_non   <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
    me_non <- as.numeric(crossprod(gme_non, b))
    se_me_non <- sqrt(as.numeric(t(gme_non) %*% vcov_mat[v_names, v_names] %*% gme_non))
    me_ci_low_non <- me_non - 1.96*se_me_non
    me_ci_high_non <- me_non + 1.96*se_me_non
    
    
    if(cz=="Temperate"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,
        0,0,0,
        0,0,0,
        wd, wd^2,
        0,0,0,
        0,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,0,
        0,0,0,
        0,0,0,
        1, 2*wd,
        0,0,0,
        0,0,0
      )
    }else if(cz=="Tropical"){
      g_urb <- c(
        wd, wd^2, 1,
        wd,0,0,
        wd^2,0,0,
        1,0,0,
        wd, wd^2,
        wd,0,0,
        wd^2,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        1,0,0,
        2*wd,0,0,
        0,0,0,
        1, 2*wd,
        1,0,0,
        2*wd,0,0
      )
    }else if(cz=="Cold"){
      g_urb <- c(
        wd, wd^2, 1,
        0,wd,0,
        0,wd^2,0,
        0,1,0,
        wd, wd^2,
        0,wd,0,
        0,wd^2,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,1,0,
        0,2*wd,0,
        0,0,0,
        1, 2*wd,
        0,1,0,
        0,2*wd,0
      )
    }else if(cz=="Arid"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,wd,
        0,0,wd^2,
        0,0,1,
        wd, wd^2,
        0,0,wd,
        0,0,wd^2
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,1,
        0,0,2*wd,
        0,0,0,
        1, 2*wd,
        0,0,1,
        0,0,2*wd
      )
    }
    
    pred_urb <- as.numeric(crossprod(g_urb, b))
    se_urb   <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
    me_urb <- as.numeric(crossprod(gme_urb, b))
    se_me_urb <- sqrt(as.numeric(t(gme_urb) %*% vcov_mat[v_names, v_names] %*% gme_urb))
    me_ci_low_urb <- me_urb - 1.96*se_me_urb
    me_ci_high_urb <- me_urb + 1.96*se_me_urb
    
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                ClimateZone=cz,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                ClimateZone=cz,
                                LandType="Urban"))
    
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_non,
                                   ci_low=me_ci_low_non,
                                   ci_high=me_ci_high_non,
                                   ClimateZone=cz,
                                   LandType="Non-urban"))
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_urb,
                                   ci_low=me_ci_low_urb,
                                   ci_high=me_ci_high_urb,
                                   ClimateZone=cz,
                                   LandType="Urban"))
  }
}


climate_limits <- tibble(
  ClimateZone = c("Temperate", "Tropical", "Arid",  "Cold"),
  min_WD = c(
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone == "Temperate", ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone == "Tropical", ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone == "Arid", ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone == "Cold", ]$WaterDeficit, na.rm = TRUE)
  ),
  max_WD = c(
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone == "Temperate", ]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone == "Tropical", ]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone == "Arid", ]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone == "Cold", ]$WaterDeficit, na.rm = TRUE)
  )
)
print(climate_limits)

plot_df_clip <- plot_df %>%
  left_join(climate_limits, by = "ClimateZone") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

me_plot_df_clip <- me_plot_df %>%
  left_join(climate_limits, by = "ClimateZone") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)


# Plot effect curve
plot_df_clip$ClimateZone <- factor(plot_df_clip$ClimateZone, levels = c( "Tropical", "Temperate",  "Cold","Arid"))

effect_plot <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC,
                                        color = ClimateZone, fill=ClimateZone,
                                        group = interaction(ClimateZone, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"),  guide = "none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.8, 0.1),
                     breaks = seq(-0.8, 0.1, by = 0.2),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.2, 0.35), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=climate_colors) +
  scale_fill_manual(values=climate_colors) 

print(effect_plot)

# Plot distribution histogram
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=ClimateZone), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.011), breaks = seq(0, 0.011, by =0.004)) +
  labs(x = "Water deficit (mm)", y = "") +
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3), 
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 6),
    legend.position = "none"
  ) + 
  scale_fill_manual(values = climate_colors)

print(hist_plot)

effect_climate_waterdeficit_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_climate_waterdeficit_FVC)


# Plot marginal effect curve
me_plot_df_clip$ClimateZone <- factor(me_plot_df_clip$ClimateZone, levels = c( "Tropical", "Temperate",  "Cold","Arid"))

marginal_plot <- ggplot(me_plot_df_clip, aes(x = WaterDeficit, y = MarginalEffect,
                                       color = ClimateZone, fill=ClimateZone,
                                       group = interaction(ClimateZone, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = "Marginal effect on FVC",
    color = "", fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(-0.0010, 0.0005),
                     breaks = seq(-0.0010, 0.0005, by = 0.0005),
                     labels = scales::number_format(accuracy = 0.0001)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.78, 0.35), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=climate_colors)+
  scale_fill_manual(values=climate_colors)

print(marginal_plot)


marginal_climate_waterdeficit_FVC <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_climate_waterdeficit_FVC)





# Continent
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]


merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))
merged_dt_WaterDeficit$continent <- factor(merged_dt_WaterDeficit$continent, levels = c( "Asia",
                                                                                             "NorthAmerica",
                                                                                             "Africa",
                                                                                             "Europe",
                                                                                             "Oceania",
                                                                                             "SouthAmerica"))


model_urban <- felm(FVC ~ 1 + 
                      WaterDeficit + WaterDeficit_2 +
                      urban +
                      WaterDeficit:continent + 
                      WaterDeficit_2:continent +
                      urban:continent +
                      WaterDeficit:urban +
                      WaterDeficit_2:urban +
                      WaterDeficit:urban:continent +
                      WaterDeficit_2:urban:continent +
                      Seasonalityindex + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    | CityID + year | 0 | CityID, 
                    data = merged_dt_WaterDeficit)
summary(model_urban)

output_file <- "FVC_urban_in_out_continent.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)


# Compute effects, marginal effects and confidence intervals
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "WaterDeficit",
  "WaterDeficit_2",
  "urbanin",
  "WaterDeficit:continentNorthAmerica",
  "WaterDeficit:continentAfrica",
  "WaterDeficit:continentEurope",
  "WaterDeficit:continentOceania",
  "WaterDeficit:continentSouthAmerica",
  "WaterDeficit_2:continentNorthAmerica",
  "WaterDeficit_2:continentAfrica",
  "WaterDeficit_2:continentEurope",
  "WaterDeficit_2:continentOceania",
  "WaterDeficit_2:continentSouthAmerica",
  "urbanin:continentNorthAmerica",
  "urbanin:continentAfrica",
  "urbanin:continentEurope",
  "urbanin:continentOceania",
  "urbanin:continentSouthAmerica",
  "WaterDeficit:urbanin",
  "WaterDeficit_2:urbanin",
  "WaterDeficit:urbanin:continentNorthAmerica",
  "WaterDeficit:urbanin:continentAfrica",
  "WaterDeficit:urbanin:continentEurope",
  "WaterDeficit:urbanin:continentOceania",
  "WaterDeficit:urbanin:continentSouthAmerica",
  "WaterDeficit_2:urbanin:continentNorthAmerica",
  "WaterDeficit_2:urbanin:continentAfrica",
  "WaterDeficit_2:urbanin:continentEurope",
  "WaterDeficit_2:urbanin:continentOceania",
  "WaterDeficit_2:urbanin:continentSouthAmerica"
)

b <- beta_all[v_names]

wd_seq <- seq(0, 2400, length.out = 200)
continents <- c("Asia","NorthAmerica","Africa","Europe","Oceania","SouthAmerica")

plot_df     <- data.frame()
me_plot_df  <- data.frame()

for(ct in continents){
  for(wd in wd_seq){
    if(ct=="Asia"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(ct=="NorthAmerica"){
      g_non <- c(
        wd, wd^2, 0,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        1,0,0,0,0,
        2*wd,0,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(ct=="Africa"){
      g_non <- c(
        wd, wd^2, 0,
        0,wd,0,0,0,
        0,wd^2,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,1,0,0,0,
        0,2*wd,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(ct=="Europe"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,wd,0,0,
        0,0,wd^2,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,1,0,0,
        0,0,2*wd,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(ct=="Oceania"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,0,wd,0,
        0,0,0,wd^2,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,0,1,0,
        0,0,0,2*wd,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(ct=="SouthAmerica"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,0,0,wd,
        0,0,0,0,wd^2,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,0,0,1,
        0,0,0,0,2*wd,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
    me_non      <- as.numeric(crossprod(gme_non, b))
    se_me_non   <- sqrt(as.numeric(t(gme_non) %*% vcov_mat[v_names, v_names] %*% gme_non))
    me_ci_low_non  <- me_non - 1.96*se_me_non
    me_ci_high_non <- me_non + 1.96*se_me_non
    
    if(ct=="Asia"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        wd, wd^2,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        1, 2*wd,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(ct=="NorthAmerica"){
      g_urb <- c(
        wd, wd^2, 1,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
        1,0,0,0,0,
        wd, wd^2,
        wd,0,0,0,0,
        wd^2,0,0,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        1,0,0,0,0,
        2*wd,0,0,0,0,
        0,0,0,0,0,
        1, 2*wd,
        1,0,0,0,0,
        2*wd,0,0,0,0
      )
    }else if(ct=="Africa"){
      g_urb <- c(
        wd, wd^2, 1,
        0,wd,0,0,0,
        0,wd^2,0,0,0,
        0,1,0,0,0,
        wd, wd^2,
        0,wd,0,0,0,
        0,wd^2,0,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,1,0,0,0,
        0,2*wd,0,0,0,
        0,0,0,0,0,
        1, 2*wd,
        0,1,0,0,0,
        0,2*wd,0,0,0
      )
    }else if(ct=="Europe"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,wd,0,0,
        0,0,wd^2,0,0,
        0,0,1,0,0,
        wd, wd^2,
        0,0,wd,0,0,
        0,0,wd^2,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,1,0,0,
        0,0,2*wd,0,0,
        0,0,0,0,0,
        1, 2*wd,
        0,0,1,0,0,
        0,0,2*wd,0,0
      )
    }else if(ct=="Oceania"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,wd,0,
        0,0,0,wd^2,0,
        0,0,0,1,0,
        wd, wd^2,
        0,0,0,wd,0,
        0,0,0,wd^2,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,0,1,0,
        0,0,0,2*wd,0,
        0,0,0,0,0,
        1, 2*wd,
        0,0,0,1,0,
        0,0,0,2*wd,0
      )
    }else if(ct=="SouthAmerica"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,0,wd,
        0,0,0,0,wd^2,
        0,0,0,0,1,
        wd, wd^2,
        0,0,0,0,wd,
        0,0,0,0,wd^2
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,0,0,1,
        0,0,0,0,2*wd,
        0,0,0,0,0,
        1, 2*wd,
        0,0,0,0,1,
        0,0,0,0,2*wd
      )
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
    me_urb      <- as.numeric(crossprod(gme_urb, b))
    se_me_urb   <- sqrt(as.numeric(t(gme_urb) %*% vcov_mat[v_names, v_names] %*% gme_urb))
    me_ci_low_urb  <- me_urb - 1.96*se_me_urb
    me_ci_high_urb <- me_urb + 1.96*se_me_urb
    

    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                Continent=ct,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                Continent=ct,
                                LandType="Urban"))
 
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_non,
                                   ci_low=me_ci_low_non,
                                   ci_high=me_ci_high_non,
                                   Continent=ct,
                                   LandType="Non-urban"))
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_urb,
                                   ci_low=me_ci_low_urb,
                                   ci_high=me_ci_high_urb,
                                   Continent=ct,
                                   LandType="Urban"))
  }
}

continent_limits <- tibble(
  Continent = c("Asia","NorthAmerica","Africa","Europe","Oceania","SouthAmerica"),
  min_WD = c(
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Asia",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="NorthAmerica",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Africa",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Europe",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Oceania",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="SouthAmerica",]$WaterDeficit,na.rm=T)
  ),
  max_WD = c(
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Asia",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="NorthAmerica",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Africa",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Europe",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Oceania",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="SouthAmerica",]$WaterDeficit,na.rm=T)
  )
)
print(continent_limits)

plot_df_clip <- plot_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

me_plot_df_clip <- me_plot_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

plot_df_clip$Continent <- factor(plot_df_clip$Continent, 
                                 levels = c("Asia","NorthAmerica","Africa","Europe","Oceania","SouthAmerica"),
                                 labels = c("ASIA","NAM","AFR","EUR","OCE","SAM")) 
me_plot_df_clip$Continent <- factor(me_plot_df_clip$Continent,
                                    levels = c("Asia","NorthAmerica","Africa","Europe","Oceania","SouthAmerica"),
                                    labels = c("ASIA","NAM","AFR","EUR","OCE","SAM"))

continent_colors <- c(
  "ASIA" = "#3498db",
  "AFR" = "#95a5a6",
  "NAM" = "#2ecc71",
  "SAM" = "#4A1054",
  "EUR" = "#e74c3c",
  "OCE" = "#f39c12"
)

# Plot effect curve
plot_df_clip$Continent <- factor(plot_df_clip$Continent,
                                 levels = c("ASIA", "AFR", "NAM", "SAM", "EUR", "OCE"))

effect_plot <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC,
                                        color = Continent, fill=Continent,
                                        group = interaction(Continent, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"),  guide = "none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600),
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.32, 0.23), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6,
                               margin = margin(l=1, unit="pt")),
    legend.key.spacing.x = unit(1, "pt"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, ncol = 3,keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=continent_colors) +
  scale_fill_manual(values=continent_colors)

print(effect_plot)


# Plot distribution histogram
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=continent), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.023), breaks = seq(0, 0.023, by = 0.01)) +
  labs(x = "Water deficit (mm)", y = "") +
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3), 
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 6),
    legend.position = "none"
  ) + 
  scale_fill_manual(values = c("Asia" = "#3498db",
                               "Africa" = "#95a5a6",
                               "NorthAmerica" = "#2ecc71",
                               "SouthAmerica" = "#4A1054",
                               "Europe" = "#e74c3c",
                               "Oceania" = "#f39c12"
  ))

print(hist_plot)

effect_continent_waterdeficit_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_continent_waterdeficit_FVC)

# Plot marginal effect curve
me_plot_df_clip$Continent <- factor(me_plot_df_clip$Continent,
                                    levels = c("ASIA", "AFR", "NAM", "SAM", "EUR", "OCE"))

marginal_plot <- ggplot(me_plot_df_clip, aes(x = WaterDeficit, y = MarginalEffect,
                                             color = Continent, fill=Continent,
                                             group = interaction(Continent, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = "Marginal effect on FVC",
    color = "", fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(-0.0012, 0.001),
                     breaks = seq(-0.0012, 0.001, by = 0.0006),
                     labels = scales::number_format(accuracy = 0.0001)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.45, 0.25), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l=1, unit="pt")),  
    legend.key.spacing.x = unit(1, "pt"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=continent_colors)+
  scale_fill_manual(values=continent_colors)

print(marginal_plot)

marginal_continent_waterdeficit_FVC <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_continent_waterdeficit_FVC)




# City size
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]

merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))
merged_dt_WaterDeficit$CitySize <- factor(merged_dt_WaterDeficit$CitySize, levels = c("Small",
                                                                                          "Medium",
                                                                                          "Large"))

model_urban <- felm(FVC ~ 1 + 
                      WaterDeficit + 
                      WaterDeficit_2 +
                      urban +
                      WaterDeficit:CitySize + 
                      WaterDeficit_2:CitySize +
                      urban:CitySize +
                      WaterDeficit:urban +
                      WaterDeficit_2:urban +
                      WaterDeficit:urban:CitySize +
                      WaterDeficit_2:urban:CitySize +
                      Seasonalityindex + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    | CityID + year | 0 | CityID, 
                    data = merged_dt_WaterDeficit)

summary(model_urban)

output_file <- "FVC_urban_in_out_citysize.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# Compute effects, marginal effects and confidence intervals
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "WaterDeficit",
  "WaterDeficit_2",
  "urbanin",
  "WaterDeficit:CitySizeMedium",
  "WaterDeficit:CitySizeLarge",
  "WaterDeficit_2:CitySizeMedium",
  "WaterDeficit_2:CitySizeLarge",
  "urbanin:CitySizeMedium",
  "urbanin:CitySizeLarge",
  "WaterDeficit:urbanin",
  "WaterDeficit_2:urbanin",
  "WaterDeficit:urbanin:CitySizeMedium",
  "WaterDeficit:urbanin:CitySizeLarge",
  "WaterDeficit_2:urbanin:CitySizeMedium",
  "WaterDeficit_2:urbanin:CitySizeLarge"
)

b <- beta_all[v_names]

citysize_colors <- c(
  "Small" = "#95a5a6",
  "Medium" = "#2ecc71",
  "Large" = "#e74c3c"
)

wd_seq <- seq(0, 2400, length.out = 200)
city_sizes <- c("Small","Medium","Large")

plot_df     <- data.frame()
me_plot_df  <- data.frame()

for(cs in city_sizes){
  for(wd in wd_seq){
    if(cs=="Small"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0
      )
    }else if(cs=="Medium"){
      g_non <- c(
        wd, wd^2, 0,
        wd,0,
        wd^2,0,
        0,0,
        0,0,
        0,0,
        0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        1,0,
        2*wd,0,
        0,0,
        0,0,
        0,0,
        0,0
      )
    }else if(cs=="Large"){
      g_non <- c(
        wd, wd^2, 0,
        0,wd,
        0,wd^2,
        0,0,
        0,0,
        0,0,
        0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,1,
        0,2*wd,
        0,0,
        0,0,
        0,0,
        0,0
      )
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
    me_non      <- as.numeric(crossprod(gme_non, b))
    se_me_non   <- sqrt(as.numeric(t(gme_non) %*% vcov_mat[v_names, v_names] %*% gme_non))
    me_ci_low_non  <- me_non - 1.96*se_me_non
    me_ci_high_non <- me_non + 1.96*se_me_non

    if(cs=="Small"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,
        0,0,
        0,0,
        wd, wd^2,
        0,0,
        0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,
        0,0,
        0,0,
        1, 2*wd,
        0,0,
        0,0
      )
    }else if(cs=="Medium"){
      g_urb <- c(
        wd, wd^2, 1,
        wd,0,
        wd^2,0,
        1,0,
        wd, wd^2,
        wd,0,
        wd^2,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        1,0,
        2*wd,0,
        0,0,
        1, 2*wd,
        1,0,
        2*wd,0
      )
    }else if(cs=="Large"){
      g_urb <- c(
        wd, wd^2, 1,
        0,wd,
        0,wd^2,
        0,1,
        wd, wd^2,
        0,wd,
        0,wd^2
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,1,
        0,2*wd,
        0,0,
        1, 2*wd,
        0,1,
        0,2*wd
      )
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
    me_urb      <- as.numeric(crossprod(gme_urb, b))
    se_me_urb   <- sqrt(as.numeric(t(gme_urb) %*% vcov_mat[v_names, v_names] %*% gme_urb))
    me_ci_low_urb  <- me_urb - 1.96*se_me_urb
    me_ci_high_urb <- me_urb + 1.96*se_me_urb
    
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                CitySize=cs,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                CitySize=cs,
                                LandType="Urban"))
    
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_non,
                                   ci_low=me_ci_low_non,
                                   ci_high=me_ci_high_non,
                                   CitySize=cs,
                                   LandType="Non-urban"))
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_urb,
                                   ci_low=me_ci_low_urb,
                                   ci_high=me_ci_high_urb,
                                   CitySize=cs,
                                   LandType="Urban"))
  }
}

citysize_limits <- tibble(
  CitySize = c("Small","Medium","Large"),
  min_WD = c(
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$CitySize=="Small",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$CitySize=="Medium",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$CitySize=="Large",]$WaterDeficit,na.rm=T)
  ),
  max_WD = c(
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$CitySize=="Small",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$CitySize=="Medium",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$CitySize=="Large",]$WaterDeficit,na.rm=T)
  )
)
print(citysize_limits)

plot_df_clip <- plot_df %>%
  left_join(citysize_limits, by = "CitySize") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

me_plot_df_clip <- me_plot_df %>%
  left_join(citysize_limits, by = "CitySize") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)


# Plot effect curve
plot_df_clip$CitySize <- factor(plot_df_clip$CitySize,
                                levels = c("Small","Medium","Large"))

effect_plot <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC,
                                        color = CitySize, fill=CitySize,
                                        group = interaction(CitySize, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"),  guide = "none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600),
                     labels = scales::number_format(accuracy = 0.01)) +

  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.25, 0.3), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=citysize_colors) +
  scale_fill_manual(values=citysize_colors)

print(effect_plot)

# Plot distribution histogram
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=CitySize), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.011), breaks = seq(0, 0.011, by =0.004)) +
  labs(x = "Water deficit (mm)", y = "") +
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3), 
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 6),
    legend.position = "none"
  ) +
  scale_fill_manual(values = citysize_colors)

print(hist_plot)

effect_citysize_waterdeficit_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_citysize_waterdeficit_FVC)


# Plot marginal effect curve
me_plot_df_clip$CitySize <- factor(me_plot_df_clip$CitySize,
                                   levels = c("Small","Medium","Large"))

marginal_plot <- ggplot(me_plot_df_clip, aes(x = WaterDeficit, y = MarginalEffect,
                                             color = CitySize, fill=CitySize,
                                             group = interaction(CitySize, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = "Marginal effect on FVC",
    color = "", fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(-0.0006, 0.0002),
                     breaks = seq(-0.0006, 0.0002, by = 0.0002),
                     labels = scales::number_format(accuracy = 0.0001)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.75, 0.37), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=citysize_colors)+
  scale_fill_manual(values=citysize_colors)

print(marginal_plot)

marginal_citysize_waterdeficit_FVC <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_citysize_waterdeficit_FVC)





# Development level
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]


merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))

merged_dt_WaterDeficit$DevelopmentLevel <- factor(merged_dt_WaterDeficit$DevelopmentLevel, levels = c("developing", "developed"))


model_urban <- felm(FVC ~ 1 + 
                      WaterDeficit + 
                      WaterDeficit_2 +
                      urban +
                      WaterDeficit:DevelopmentLevel + 
                      WaterDeficit_2:DevelopmentLevel +
                      urban:DevelopmentLevel +
                      WaterDeficit:urban +
                      WaterDeficit_2:urban +
                      WaterDeficit:urban:DevelopmentLevel +
                      WaterDeficit_2:urban:DevelopmentLevel +
                      Seasonalityindex + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    | CityID + year | 0 | CityID, 
                    data = merged_dt_WaterDeficit)


summary(model_urban)

output_file <- "FVC_urban_in_out_level.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# Compute effects, marginal effects and confidence intervals
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "WaterDeficit",
  "WaterDeficit_2",
  "urbanin",
  "WaterDeficit:DevelopmentLeveldeveloped",
  "WaterDeficit_2:DevelopmentLeveldeveloped",
  "urbanin:DevelopmentLeveldeveloped",
  "WaterDeficit:urbanin",
  "WaterDeficit_2:urbanin",
  "WaterDeficit:urbanin:DevelopmentLeveldeveloped",
  "WaterDeficit_2:urbanin:DevelopmentLeveldeveloped"
)

b <- beta_all[v_names]

dev_colors <- c(
  "Global North" = "#0071bc", 
  "Global South" = "#d95218"
)

wd_seq <- seq(0, 2400, length.out = 200)
dev_levels <- c("developing","developed")

plot_df     <- data.frame()
me_plot_df  <- data.frame()

for(dl in dev_levels){
  for(wd in wd_seq){
    if(dl=="developing"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,0,
        0,0,
        0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,0,
        0,0,
        0,0
      )
    }else if(dl=="developed"){
      g_non <- c(
        wd, wd^2, 0,
        wd, wd^2, 0,
        0,0,
        0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        1, 2*wd, 0,
        0,0,
        0,0
      )
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
    me_non      <- as.numeric(crossprod(gme_non, b))
    se_me_non   <- sqrt(as.numeric(t(gme_non) %*% vcov_mat[v_names, v_names] %*% gme_non))
    me_ci_low_non  <- me_non - 1.96*se_me_non
    me_ci_high_non <- me_non + 1.96*se_me_non

    if(dl=="developing"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,
        wd, wd^2,
        0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,0,
        1, 2*wd,
        0,0
      )
    }else if(dl=="developed"){
      g_urb <- c(
        wd, wd^2, 1,
        wd, wd^2, 1,
        wd, wd^2,
        wd, wd^2
      )
      gme_urb <- c(
        1, 2*wd, 0,
        1, 2*wd, 0,
        1, 2*wd,
        1, 2*wd
      )
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
    me_urb      <- as.numeric(crossprod(gme_urb, b))
    se_me_urb   <- sqrt(as.numeric(t(gme_urb) %*% vcov_mat[v_names, v_names] %*% gme_urb))
    me_ci_low_urb  <- me_urb - 1.96*se_me_urb
    me_ci_high_urb <- me_urb + 1.96*se_me_urb
    
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                DevelopmentLevel=dl,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                DevelopmentLevel=dl,
                                LandType="Urban"))
    
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_non,
                                   ci_low=me_ci_low_non,
                                   ci_high=me_ci_high_non,
                                   DevelopmentLevel=dl,
                                   LandType="Non-urban"))
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_urb,
                                   ci_low=me_ci_low_urb,
                                   ci_high=me_ci_high_urb,
                                   DevelopmentLevel=dl,
                                   LandType="Urban"))
  }
}

dev_limits <- tibble(
  DevelopmentLevel = c("developing","developed"),
  min_WD = c(
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$DevelopmentLevel=="developing",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$DevelopmentLevel=="developed",]$WaterDeficit,na.rm=T)
  ),
  max_WD = c(
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$DevelopmentLevel=="developing",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$DevelopmentLevel=="developed",]$WaterDeficit,na.rm=T)
  )
)
print(dev_limits)

plot_df_clip <- plot_df %>%
  left_join(dev_limits, by = "DevelopmentLevel") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

me_plot_df_clip <- me_plot_df %>%
  left_join(dev_limits, by = "DevelopmentLevel") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

plot_df_clip$DevelopmentLevel <- factor(plot_df_clip$DevelopmentLevel, 
                                 levels = c("developing","developed"),
                                 labels = c("Global South","Global North")) 
me_plot_df_clip$DevelopmentLevel <- factor(me_plot_df_clip$DevelopmentLevel,
                                    levels = c("developing","developed"),
                                    labels = c("Global South","Global North"))


# Plot effect curve
effect_plot <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC,
                                        color = DevelopmentLevel, fill=DevelopmentLevel,
                                        group = interaction(DevelopmentLevel, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"),  guide = "none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600),
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
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=dev_colors) +
  scale_fill_manual(values=dev_colors)

print(effect_plot)


# Plot distribution histogram
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=DevelopmentLevel), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.01), breaks = seq(0, 0.01, by =0.004)) +
  labs(x = "Water deficit (mm)", y = "") +
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3), 
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 6),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c( "developed" = "#0071bc", 
                                 "developing" = "#d95218"))

print(hist_plot)

effect_dev_waterdeficit_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_dev_waterdeficit_FVC)


# Plot marginal effect curve
marginal_plot <- ggplot(me_plot_df_clip, aes(x = WaterDeficit, y = MarginalEffect,
                                             color = DevelopmentLevel, fill=DevelopmentLevel,
                                             group = interaction(DevelopmentLevel, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = "Marginal effect on FVC",
    color = "", fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(-0.0004, 0.0002),
                     breaks = seq(-0.0004, 0.0002, by = 0.0002),
                     labels = scales::number_format(accuracy = 0.0001)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.72, 0.25), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=dev_colors)+
  scale_fill_manual(values=dev_colors)

print(marginal_plot)

marginal_dev_waterdeficit_FVC <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_dev_waterdeficit_FVC)




# Country
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]


merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))

merged_dt_WaterDeficit$country <- factor(merged_dt_WaterDeficit$country, 
                                                        levels = c("United States of America",
                                                                   "Europe",
                                                                   "China",
                                                                   "India",
                                                                   "Brazil", "other"))

model_urban <- felm(FVC ~ 1 + 
                      WaterDeficit + 
                      WaterDeficit_2 +
                      urban +
                      WaterDeficit:country + 
                      WaterDeficit_2:country +
                      urban:country +
                      WaterDeficit:urban +
                      WaterDeficit_2:urban +
                      WaterDeficit:urban:country +
                      WaterDeficit_2:urban:country +
                      Seasonalityindex + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    | CityID + year | 0 | CityID, 
                    data = merged_dt_WaterDeficit)
summary(model_urban)

output_file <- "FVC_urban_in_out_country.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# Compute effects, marginal effects and confidence intervals
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "WaterDeficit",
  "WaterDeficit_2",
  "urbanin",
  "WaterDeficit:countryEurope",
  "WaterDeficit:countryChina",
  "WaterDeficit:countryIndia",
  "WaterDeficit:countryBrazil",
  "WaterDeficit:countryother",
  "WaterDeficit_2:countryEurope",
  "WaterDeficit_2:countryChina",
  "WaterDeficit_2:countryIndia",
  "WaterDeficit_2:countryBrazil",
  "WaterDeficit_2:countryother",
  "urbanin:countryEurope",
  "urbanin:countryChina",
  "urbanin:countryIndia",
  "urbanin:countryBrazil",
  "urbanin:countryother",
  "WaterDeficit:urbanin",
  "WaterDeficit_2:urbanin",
  "WaterDeficit:urbanin:countryEurope",
  "WaterDeficit:urbanin:countryChina",
  "WaterDeficit:urbanin:countryIndia",
  "WaterDeficit:urbanin:countryBrazil",
  "WaterDeficit:urbanin:countryother",
  "WaterDeficit_2:urbanin:countryEurope",
  "WaterDeficit_2:urbanin:countryChina",
  "WaterDeficit_2:urbanin:countryIndia",
  "WaterDeficit_2:urbanin:countryBrazil",
  "WaterDeficit_2:urbanin:countryother"
)

b <- beta_all[v_names]

wd_seq <- seq(0, 2400, length.out = 200)
cntry_levels <- c("United States","Europe (excluding Russia)","China","India","Brazil","other")

plot_df     <- data.frame()
me_plot_df  <- data.frame()

for(cnt in cntry_levels){
  for(wd in wd_seq){
    if(cnt=="United States"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(cnt=="Europe (excluding Russia)"){
      g_non <- c(
        wd, wd^2, 0,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        1,0,0,0,0,
        2*wd,0,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(cnt=="China"){
      g_non <- c(
        wd, wd^2, 0,
        0,wd,0,0,0,
        0,wd^2,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,1,0,0,0,
        0,2*wd,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(cnt=="India"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,wd,0,0,
        0,0,wd^2,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,1,0,0,
        0,0,2*wd,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(cnt=="Brazil"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,0,wd,0,
        0,0,0,wd^2,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,0,1,0,
        0,0,0,2*wd,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(cnt=="other"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,0,0,wd,
        0,0,0,0,wd^2,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,0,0,1,
        0,0,0,0,2*wd,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
    me_non      <- as.numeric(crossprod(gme_non, b))
    se_me_non   <- sqrt(as.numeric(t(gme_non) %*% vcov_mat[v_names, v_names] %*% gme_non))
    me_ci_low_non  <- me_non - 1.96*se_me_non
    me_ci_high_non <- me_non + 1.96*se_me_non
    
    if(cnt=="United States"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        wd, wd^2,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        1, 2*wd,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(cnt=="Europe (excluding Russia)"){
      g_urb <- c(
        wd, wd^2, 1,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
        1,0,0,0,0,
        wd, wd^2,
        wd,0,0,0,0,
        wd^2,0,0,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        1,0,0,0,0,
        2*wd,0,0,0,0,
        0,0,0,0,0,
        1, 2*wd,
        1,0,0,0,0,
        2*wd,0,0,0,0
      )
    }else if(cnt=="China"){
      g_urb <- c(
        wd, wd^2, 1,
        0,wd,0,0,0,
        0,wd^2,0,0,0,
        0,1,0,0,0,
        wd, wd^2,
        0,wd,0,0,0,
        0,wd^2,0,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,1,0,0,0,
        0,2*wd,0,0,0,
        0,0,0,0,0,
        1, 2*wd,
        0,1,0,0,0,
        0,2*wd,0,0,0
      )
    }else if(cnt=="India"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,wd,0,0,
        0,0,wd^2,0,0,
        0,0,1,0,0,
        wd, wd^2,
        0,0,wd,0,0,
        0,0,wd^2,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,1,0,0,
        0,0,2*wd,0,0,
        0,0,0,0,0,
        1, 2*wd,
        0,0,1,0,0,
        0,0,2*wd,0,0
      )
    }else if(cnt=="Brazil"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,wd,0,
        0,0,0,wd^2,0,
        0,0,0,1,0,
        wd, wd^2,
        0,0,0,wd,0,
        0,0,0,wd^2,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,0,1,0,
        0,0,0,2*wd,0,
        0,0,0,0,0,
        1, 2*wd,
        0,0,0,1,0,
        0,0,0,2*wd,0
      )
    }else if(cnt=="other"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,0,wd,
        0,0,0,0,wd^2,
        0,0,0,0,1,
        wd, wd^2,
        0,0,0,0,wd,
        0,0,0,0,wd^2
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,0,0,1,
        0,0,0,0,2*wd,
        0,0,0,0,0,
        1, 2*wd,
        0,0,0,0,1,
        0,0,0,0,2*wd
      )
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
    me_urb      <- as.numeric(crossprod(gme_urb, b))
    se_me_urb   <- sqrt(as.numeric(t(gme_urb) %*% vcov_mat[v_names, v_names] %*% gme_urb))
    me_ci_low_urb  <- me_urb - 1.96*se_me_urb
    me_ci_high_urb <- me_urb + 1.96*se_me_urb
    
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                Country=cnt,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                Country=cnt,
                                LandType="Urban"))
    
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_non,
                                   ci_low=me_ci_low_non,
                                   ci_high=me_ci_high_non,
                                   Country=cnt,
                                   LandType="Non-urban"))
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_urb,
                                   ci_low=me_ci_low_urb,
                                   ci_high=me_ci_high_urb,
                                   Country=cnt,
                                   LandType="Urban"))
  }
}

country_limits <- tibble(
  Country = c("United States", "Europe (excluding Russia)", "China", "India", "Brazil", "other"),
  min_WD = c(
    min(merged_dt_WaterDeficit[country == "United States of America"]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "Europe" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "China" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "India" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "Brazil" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "other" ]$WaterDeficit, na.rm = TRUE)
  ),
  max_WD = c(
    max(merged_dt_WaterDeficit[country == "United States of America" ]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[country == "Europe" ]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[country == "China"]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[country == "India"]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[country == "Brazil"]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[country == "other"]$WaterDeficit, na.rm = TRUE)
  )
)

print(country_limits)

plot_df_clip <- plot_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

me_plot_df_clip <- me_plot_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)


plot_df_clip$Country <- factor(plot_df_clip$Country, 
                               levels = c("United States", "Europe (excluding Russia)",  "China", "India", "Brazil" ,"other"),
                               labels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH")) 
me_plot_df_clip$Country <- factor(me_plot_df_clip$Country,
                                  levels = c("United States", "Europe (excluding Russia)",  "China", "India", "Brazil" ,"other"),
                                  labels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH"))

country_colors <- c(
  "US" = "#2ecc71",
  "EURxR" = "#3498db",
  "CHN" = "#CB0505",
  "IND" = "#f39c12",
  "BRA" = "#4A1054",
  "OTH" = "#95a5a6"
)

# Plot effect curve
plot_df_clip$Country <- factor(plot_df_clip$Country, levels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH"))

effect_plot <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC,
                                        color = Country, fill=Country,
                                        group = interaction(Country, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"),  guide = "none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600),
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-1.4, 0.2),
                     breaks = seq(-1.4, 0.2, by = 0.35),
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.14, 0.44),  
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")), 
    legend.key.spacing.x = unit(1, "pt"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 6, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +  
  scale_color_manual(values=country_colors) +
  scale_fill_manual(values=country_colors)

print(effect_plot)


# Plot distribution histogram
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=country), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.025), breaks = seq(0, 0.025, by =0.01)) +
  labs(x = "Water deficit (mm)", y = "") +
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3), 
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 6),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c(
    "United States of America" = "#2ecc71",
    "Europe" = "#3498db",
    "China" = "#CB0505",
    "India" = "#f39c12",
    "Brazil" = "#4A1054",
    "other" = "#95a5a6"
  ))

print(hist_plot)

effect_country_waterdeficit_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_country_waterdeficit_FVC)

# Plot marginal effect curve
me_plot_df_clip$Country <- factor(me_plot_df_clip$Country, levels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH"))

marginal_plot <- ggplot(me_plot_df_clip, aes(x = WaterDeficit, y = MarginalEffect,
                                             color = Country, fill=Country,
                                             group = interaction(Country, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = "Marginal effect on FVC",
    color = "", fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(-0.0012, 0.0010),
                     breaks = seq(-0.0012, 0.0010, by = 0.0006),
                     labels = scales::number_format(accuracy = 0.0001)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.59, 0.95), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.key.spacing.x = unit(1, "pt"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=country_colors)+
  scale_fill_manual(values=country_colors)

print(marginal_plot)

marginal_country_waterdeficit_FVC <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_country_waterdeficit_FVC)




###########################################################################
#                     Part 2: NPP Regression & Visualization
###########################################################################

# Climate zone
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]


merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))
merged_dt_WaterDeficit$ClimateZone <- factor(merged_dt_WaterDeficit$ClimateZone, levels = c("Temperate","Tropical","Cold","Arid"))

model_urban <- felm(NPP ~ 1 + 
                      WaterDeficit +                     
                      WaterDeficit_2 +  
                      urban +  
                      WaterDeficit:ClimateZone + 
                      WaterDeficit_2:ClimateZone + 
                      urban:ClimateZone +  
                      WaterDeficit:urban + 
                      WaterDeficit_2:urban + 
                      WaterDeficit:urban:ClimateZone +  
                      WaterDeficit_2:urban:ClimateZone + 
                      Seasonalityindex + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    | CityID + year | 0 | CityID, 
                    data = merged_dt_WaterDeficit)

summary(model_urban)

output_file <- "NPP_urban_in_out_climate.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# Compute effects, marginal effects and confidence intervals
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "WaterDeficit",
  "WaterDeficit_2",
  "urbanin",
  "WaterDeficit:ClimateZoneTropical",
  "WaterDeficit:ClimateZoneCold",
  "WaterDeficit:ClimateZoneArid",
  "WaterDeficit_2:ClimateZoneTropical",
  "WaterDeficit_2:ClimateZoneCold",
  "WaterDeficit_2:ClimateZoneArid",
  "urbanin:ClimateZoneTropical",
  "urbanin:ClimateZoneCold",
  "urbanin:ClimateZoneArid",
  "WaterDeficit:urbanin",
  "WaterDeficit_2:urbanin",
  "WaterDeficit:urbanin:ClimateZoneTropical",
  "WaterDeficit:urbanin:ClimateZoneCold",
  "WaterDeficit:urbanin:ClimateZoneArid",
  "WaterDeficit_2:urbanin:ClimateZoneTropical",
  "WaterDeficit_2:urbanin:ClimateZoneCold",
  "WaterDeficit_2:urbanin:ClimateZoneArid"
)
b <- beta_all[v_names]

climate_colors <- c(
  "Temperate" = "#2ecc71",
  "Tropical" = "#e74c3c",
  "Arid" = "#95a5a6",
  "Cold" = "#3498db"
)

wd_seq <- seq(0, 2400, length.out = 200)
climates <- c("Temperate","Tropical","Cold","Arid")
plot_df    <- data.frame()
me_plot_df <- data.frame()


for(cz in climates){
  for(wd in wd_seq){
    if(cz=="Temperate"){
      g_non <- c(wd, wd^2, 0,
                 0,0,0,
                 0,0,0,
                 0,0,0,
                 0,0,0,0,0,
                 0,0,0)
      gme_non <- c(1, 2*wd, 0,
                   0,0,0,
                   0,0,0,
                   0,0,0,
                   0,0,0,0,0,
                   0,0,0)
    }else if(cz=="Tropical"){
      g_non <- c(wd, wd^2, 0,
                 wd,0,0,
                 wd^2,0,0,
                 0,0,0,
                 0,0,0,0,0,
                 0,0,0)
      gme_non <- c(1, 2*wd, 0,
                   1,0,0,
                   2*wd,0,0,
                   0,0,0,
                   0,0,0,0,0,
                   0,0,0)
    }else if(cz=="Cold"){
      g_non <- c(wd, wd^2, 0,
                 0,wd,0,
                 0,wd^2,0,
                 0,0,0,
                 0,0,0,0,0,
                 0,0,0)
      gme_non <- c(1, 2*wd, 0,
                   0,1,0,
                   0,2*wd,0,
                   0,0,0,
                   0,0,0,0,0,
                   0,0,0)
    }else if(cz=="Arid"){
      g_non <- c(wd, wd^2, 0,
                 0,0,wd,
                 0,0,wd^2,
                 0,0,0,
                 0,0,0,0,0,
                 0,0,0)
      gme_non <- c(1, 2*wd, 0,
                   0,0,1,
                   0,0,2*wd,
                   0,0,0,
                   0,0,0,0,0,
                   0,0,0)
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
    me_non      <- as.numeric(crossprod(gme_non, b))
    se_me_non   <- sqrt(as.numeric(t(gme_non) %*% vcov_mat[v_names, v_names] %*% gme_non))
    me_ci_low_non  <- me_non - 1.96*se_me_non
    me_ci_high_non <- me_non + 1.96*se_me_non
    

    if(cz=="Temperate"){
      g_urb <- c(wd, wd^2, 1,
                 0,0,0,
                 0,0,0,
                 0,0,0,
                 wd, wd^2,
                 0,0,0,
                 0,0,0)
      gme_urb <- c(1, 2*wd, 0,
                   0,0,0,
                   0,0,0,
                   0,0,0,
                   1, 2*wd,
                   0,0,0,
                   0,0,0)
    }else if(cz=="Tropical"){
      g_urb <- c(wd, wd^2, 1,
                 wd,0,0,
                 wd^2,0,0,
                 1,0,0,
                 wd, wd^2,
                 wd,0,0,
                 wd^2,0,0)
      gme_urb <- c(1, 2*wd, 0,
                   1,0,0,
                   2*wd,0,0,
                   0,0,0,
                   1, 2*wd,
                   1,0,0,
                   2*wd,0,0)
    }else if(cz=="Cold"){
      g_urb <- c(wd, wd^2, 1,
                 0,wd,0,
                 0,wd^2,0,
                 0,1,0,
                 wd, wd^2,
                 0,wd,0,
                 0,wd^2,0)
      gme_urb <- c(1, 2*wd, 0,
                   0,1,0,
                   0,2*wd,0,
                   0,0,0,
                   1, 2*wd,
                   0,1,0,
                   0,2*wd,0)
    }else if(cz=="Arid"){
      g_urb <- c(wd, wd^2, 1,
                 0,0,wd,
                 0,0,wd^2,
                 0,0,1,
                 wd, wd^2,
                 0,0,wd,
                 0,0,wd^2)
      gme_urb <- c(1, 2*wd, 0,
                   0,0,1,
                   0,0,2*wd,
                   0,0,0,
                   1, 2*wd,
                   0,0,1,
                   0,0,2*wd)
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
    me_urb      <- as.numeric(crossprod(gme_urb, b))
    se_me_urb   <- sqrt(as.numeric(t(gme_urb) %*% vcov_mat[v_names, v_names] %*% gme_urb))
    me_ci_low_urb  <- me_urb - 1.96*se_me_urb
    me_ci_high_urb <- me_urb + 1.96*se_me_urb
    
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_NPP=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                ClimateZone=cz,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_NPP=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                ClimateZone=cz,
                                LandType="Urban"))
    
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_non,
                                   ci_low=me_ci_low_non,
                                   ci_high=me_ci_high_non,
                                   ClimateZone=cz,
                                   LandType="Non-urban"))
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_urb,
                                   ci_low=me_ci_low_urb,
                                   ci_high=me_ci_high_urb,
                                   ClimateZone=cz,
                                   LandType="Urban"))
  }
}

climate_limits <- tibble(
  ClimateZone = c("Temperate", "Tropical", "Arid",  "Cold"),
  min_WD = c(
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone=="Temperate",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone=="Tropical",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone=="Arid",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone=="Cold",]$WaterDeficit,na.rm=T)
  ),
  max_WD = c(
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone=="Temperate",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone=="Tropical",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone=="Arid",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$ClimateZone=="Cold",]$WaterDeficit,na.rm=T)
  )
)
print(climate_limits)

plot_df_clip <- plot_df %>%
  left_join(climate_limits, by = "ClimateZone") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

me_plot_df_clip <- me_plot_df %>%
  left_join(climate_limits, by = "ClimateZone") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

# Plot effect curve
plot_df_clip$ClimateZone <- factor(plot_df_clip$ClimateZone, levels = c( "Tropical", "Temperate",  "Cold","Arid"))

effect_plot <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_NPP,
                                        color = ClimateZone, fill=ClimateZone,
                                        group = interaction(ClimateZone, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"),  guide = "none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = expression(Effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "",
    fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.82, 0.85), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=climate_colors) +
  scale_fill_manual(values=climate_colors)

print(effect_plot)

# Plot distribution histogram
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=ClimateZone), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.011), breaks = seq(0, 0.011, by =0.004)) +
  labs(x = "Water deficit (mm)", y = "") +
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3), 
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 6),
    legend.position = "none"
  ) +
  scale_fill_manual(values = climate_colors)

print(hist_plot)

effect_climate_waterdeficit_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_climate_waterdeficit_NPP)

# Plot marginal effect curve
me_plot_df_clip$ClimateZone <- factor(me_plot_df_clip$ClimateZone, levels = c( "Tropical", "Temperate",  "Cold","Arid"))

marginal_plot <- ggplot(me_plot_df_clip, aes(x = WaterDeficit, y = MarginalEffect,
                                             color = ClimateZone, fill=ClimateZone,
                                             group = interaction(ClimateZone, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = expression(Marginal~effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "", fill =""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.75, 0.33), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=climate_colors)+
  scale_fill_manual(values=climate_colors)

print(marginal_plot)

marginal_climate_waterdeficit_NPP <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_climate_waterdeficit_NPP)


# Continent
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]


merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))
merged_dt_WaterDeficit$continent <- factor(merged_dt_WaterDeficit$continent, levels = c( "Asia",
                                                                                         "NorthAmerica",
                                                                                         "Africa",
                                                                                         "Europe",
                                                                                         "Oceania",
                                                                                         "SouthAmerica"))

model_urban <- felm(NPP ~ 1 + 
                      WaterDeficit + WaterDeficit_2 +
                      urban +
                      WaterDeficit:continent + 
                      WaterDeficit_2:continent +
                      urban:continent +
                      WaterDeficit:urban +
                      WaterDeficit_2:urban +
                      WaterDeficit:urban:continent +
                      WaterDeficit_2:urban:continent +
                      Seasonalityindex + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    | CityID + year | 0 | CityID, 
                    data = merged_dt_WaterDeficit)
summary(model_urban)

output_file <- "NPP_urban_in_out_continent.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# Compute effects, marginal effects and confidence intervals
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "WaterDeficit",
  "WaterDeficit_2",
  "urbanin",
  "WaterDeficit:continentNorthAmerica",
  "WaterDeficit:continentAfrica",
  "WaterDeficit:continentEurope",
  "WaterDeficit:continentOceania",
  "WaterDeficit:continentSouthAmerica",
  "WaterDeficit_2:continentNorthAmerica",
  "WaterDeficit_2:continentAfrica",
  "WaterDeficit_2:continentEurope",
  "WaterDeficit_2:continentOceania",
  "WaterDeficit_2:continentSouthAmerica",
  "urbanin:continentNorthAmerica",
  "urbanin:continentAfrica",
  "urbanin:continentEurope",
  "urbanin:continentOceania",
  "urbanin:continentSouthAmerica",
  "WaterDeficit:urbanin",
  "WaterDeficit_2:urbanin",
  "WaterDeficit:urbanin:continentNorthAmerica",
  "WaterDeficit:urbanin:continentAfrica",
  "WaterDeficit:urbanin:continentEurope",
  "WaterDeficit:urbanin:continentOceania",
  "WaterDeficit:urbanin:continentSouthAmerica",
  "WaterDeficit_2:urbanin:continentNorthAmerica",
  "WaterDeficit_2:urbanin:continentAfrica",
  "WaterDeficit_2:urbanin:continentEurope",
  "WaterDeficit_2:urbanin:continentOceania",
  "WaterDeficit_2:urbanin:continentSouthAmerica"
)
b <- beta_all[v_names]

wd_seq <- seq(0, 2400, length.out = 200)
continents <- c("Asia","NorthAmerica","Africa","Europe","Oceania","SouthAmerica")
plot_df    <- data.frame()
me_plot_df <- data.frame()


for(ct in continents){
  for(wd in wd_seq){
    if(ct=="Asia"){
      g_non <- c(wd, wd^2, 0,
                 0,0,0,0,0,
                 0,0,0,0,0,
                 0,0,0,0,0,
                 0,0,
                 0,0,0,0,0,
                 0,0,0,0,0)
      gme_non <- c(1, 2*wd, 0,
                   0,0,0,0,0,
                   0,0,0,0,0,
                   0,0,0,0,0,
                   0,0,
                   0,0,0,0,0,
                   0,0,0,0,0)
    }else if(ct=="NorthAmerica"){
      g_non <- c(wd, wd^2, 0,
                 wd,0,0,0,0,
                 wd^2,0,0,0,0,
                 0,0,0,0,0,
                 0,0,
                 0,0,0,0,0,
                 0,0,0,0,0)
      gme_non <- c(1, 2*wd, 0,
                   1,0,0,0,0,
                   2*wd,0,0,0,0,
                   0,0,0,0,0,
                   0,0,
                   0,0,0,0,0,
                   0,0,0,0,0)
    }else if(ct=="Africa"){
      g_non <- c(wd, wd^2, 0,
                 0,wd,0,0,0,
                 0,wd^2,0,0,0,
                 0,0,0,0,0,
                 0,0,
                 0,0,0,0,0,
                 0,0,0,0,0)
      gme_non <- c(1, 2*wd, 0,
                   0,1,0,0,0,
                   0,2*wd,0,0,0,
                   0,0,0,0,0,
                   0,0,
                   0,0,0,0,0,
                   0,0,0,0,0)
    }else if(ct=="Europe"){
      g_non <- c(wd, wd^2, 0,
                 0,0,wd,0,0,
                 0,0,wd^2,0,0,
                 0,0,0,0,0,
                 0,0,
                 0,0,0,0,0,
                 0,0,0,0,0)
      gme_non <- c(1, 2*wd, 0,
                   0,0,1,0,0,
                   0,0,2*wd,0,0,
                   0,0,0,0,0,
                   0,0,
                   0,0,0,0,0,
                   0,0,0,0,0)
    }else if(ct=="Oceania"){
      g_non <- c(wd, wd^2, 0,
                 0,0,0,wd,0,
                 0,0,0,wd^2,0,
                 0,0,0,0,0,
                 0,0,
                 0,0,0,0,0,
                 0,0,0,0,0)
      gme_non <- c(1, 2*wd, 0,
                   0,0,0,1,0,
                   0,0,0,2*wd,0,
                   0,0,0,0,0,
                   0,0,
                   0,0,0,0,0,
                   0,0,0,0,0)
    }else if(ct=="SouthAmerica"){
      g_non <- c(wd, wd^2, 0,
                 0,0,0,0,wd,
                 0,0,0,0,wd^2,
                 0,0,0,0,0,
                 0,0,
                 0,0,0,0,0,
                 0,0,0,0,0)
      gme_non <- c(1, 2*wd, 0,
                   0,0,0,0,1,
                   0,0,0,0,2*wd,
                   0,0,0,0,0,
                   0,0,
                   0,0,0,0,0,
                   0,0,0,0,0)
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
    me_non      <- as.numeric(crossprod(gme_non, b))
    se_me_non   <- sqrt(as.numeric(t(gme_non) %*% vcov_mat[v_names, v_names] %*% gme_non))
    me_ci_low_non  <- me_non - 1.96*se_me_non
    me_ci_high_non <- me_non + 1.96*se_me_non
    
    if(ct=="Asia"){
      g_urb <- c(wd, wd^2, 1,
                 0,0,0,0,0,
                 0,0,0,0,0,
                 0,0,0,0,0,
                 wd, wd^2,
                 0,0,0,0,0,
                 0,0,0,0,0)
      gme_urb <- c(1, 2*wd, 0,
                   0,0,0,0,0,
                   0,0,0,0,0,
                   0,0,0,0,0,
                   1, 2*wd,
                   0,0,0,0,0,
                   0,0,0,0,0)
    }else if(ct=="NorthAmerica"){
      g_urb <- c(wd, wd^2, 1,
                 wd,0,0,0,0,
                 wd^2,0,0,0,0,
                 1,0,0,0,0,
                 wd, wd^2,
                 wd,0,0,0,0,
                 wd^2,0,0,0,0)
      gme_urb <- c(1, 2*wd, 0,
                   1,0,0,0,0,
                   2*wd,0,0,0,0,
                   0,0,0,0,0,
                   1, 2*wd,
                   1,0,0,0,0,
                   2*wd,0,0,0,0)
    }else if(ct=="Africa"){
      g_urb <- c(wd, wd^2, 1,
                 0,wd,0,0,0,
                 0,wd^2,0,0,0,
                 0,1,0,0,0,
                 wd, wd^2,
                 0,wd,0,0,0,
                 0,wd^2,0,0,0)
      gme_urb <- c(1, 2*wd, 0,
                   0,1,0,0,0,
                   0,2*wd,0,0,0,
                   0,0,0,0,0,
                   1, 2*wd,
                   0,1,0,0,0,
                   0,2*wd,0,0,0)
    }else if(ct=="Europe"){
      g_urb <- c(wd, wd^2, 1,
                 0,0,wd,0,0,
                 0,0,wd^2,0,0,
                 0,0,1,0,0,
                 wd, wd^2,
                 0,0,wd,0,0,
                 0,0,wd^2,0,0)
      gme_urb <- c(1, 2*wd, 0,
                   0,0,1,0,0,
                   0,0,2*wd,0,0,
                   0,0,0,0,0,
                   1, 2*wd,
                   0,0,1,0,0,
                   0,0,2*wd,0,0)
    }else if(ct=="Oceania"){
      g_urb <- c(wd, wd^2, 1,
                 0,0,0,wd,0,
                 0,0,0,wd^2,0,
                 0,0,0,1,0,
                 wd, wd^2,
                 0,0,0,wd,0,
                 0,0,0,wd^2,0)
      gme_urb <- c(1, 2*wd, 0,
                   0,0,0,1,0,
                   0,0,0,2*wd,0,
                   0,0,0,0,0,
                   1, 2*wd,
                   0,0,0,1,0,
                   0,0,0,2*wd,0)
    }else if(ct=="SouthAmerica"){
      g_urb <- c(wd, wd^2, 1,
                 0,0,0,0,wd,
                 0,0,0,0,wd^2,
                 0,0,0,0,1,
                 wd, wd^2,
                 0,0,0,0,wd,
                 0,0,0,0,wd^2)
      gme_urb <- c(1, 2*wd, 0,
                   0,0,0,0,1,
                   0,0,0,0,2*wd,
                   0,0,0,0,0,
                   1, 2*wd,
                   0,0,0,0,1,
                   0,0,0,0,2*wd)
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
    me_urb      <- as.numeric(crossprod(gme_urb, b))
    se_me_urb   <- sqrt(as.numeric(t(gme_urb) %*% vcov_mat[v_names, v_names] %*% gme_urb))
    me_ci_low_urb  <- me_urb - 1.96*se_me_urb
    me_ci_high_urb <- me_urb + 1.96*se_me_urb
    
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_NPP=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                Continent=ct,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_NPP=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                Continent=ct,
                                LandType="Urban"))

    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_non,
                                   ci_low=me_ci_low_non,
                                   ci_high=me_ci_high_non,
                                   Continent=ct,
                                   LandType="Non-urban"))
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_urb,
                                   ci_low=me_ci_low_urb,
                                   ci_high=me_ci_high_urb,
                                   Continent=ct,
                                   LandType="Urban"))
  }
}

continent_limits <- tibble(
  Continent = c("Asia","NorthAmerica","Africa","Europe","Oceania","SouthAmerica"),
  min_WD = c(
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Asia",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="NorthAmerica",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Africa",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Europe",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Oceania",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="SouthAmerica",]$WaterDeficit,na.rm=T)
  ),
  max_WD = c(
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Asia",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="NorthAmerica",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Africa",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Europe",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="Oceania",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$continent=="SouthAmerica",]$WaterDeficit,na.rm=T)
  )
)
print(continent_limits)

plot_df_clip <- plot_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

me_plot_df_clip <- me_plot_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

plot_df_clip$Continent <- factor(plot_df_clip$Continent, 
                                 levels = c("Asia","NorthAmerica","Africa","Europe","Oceania","SouthAmerica"),
                                 labels = c("ASIA","NAM","AFR","EUR","OCE","SAM")) 
me_plot_df_clip$Continent <- factor(me_plot_df_clip$Continent,
                                    levels = c("Asia","NorthAmerica","Africa","Europe","Oceania","SouthAmerica"),
                                    labels = c("ASIA","NAM","AFR","EUR","OCE","SAM"))

continent_colors <- c(
  "ASIA" = "#3498db",
  "AFR" = "#95a5a6",
  "NAM" = "#2ecc71",
  "SAM" = "#4A1054",
  "EUR" = "#e74c3c",
  "OCE" = "#f39c12"
)

# Plot effect curve
plot_df_clip$Continent <- factor(plot_df_clip$Continent,
                                 levels = c("ASIA", "AFR", "NAM", "SAM", "EUR", "OCE"))

effect_plot <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_NPP,
                                        color = Continent, fill=Continent,
                                        group = interaction(Continent, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"),  guide = "none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = expression(Effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "",
    fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.35, 0.25), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.key.spacing.x = unit(1, "pt"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=continent_colors) +
  scale_fill_manual(values=continent_colors)

print(effect_plot)

# Plot distribution histogram
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=continent), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.011), breaks = seq(0, 0.011, by =0.004)) +
  labs(x = "Water deficit (mm)", y = "") +
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3), 
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 6),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("Asia" = "#3498db",
                               "Africa" = "#95a5a6",
                               "NorthAmerica" = "#2ecc71",
                               "SouthAmerica" = "#4A1054",
                               "Europe" = "#e74c3c",
                               "Oceania" = "#f39c12"
  ))

print(hist_plot)

effect_continent_waterdeficit_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_continent_waterdeficit_NPP)

# Plot marginal effect curve
me_plot_df_clip$Continent <- factor(me_plot_df_clip$Continent,
                                    levels = c("ASIA", "AFR", "NAM", "SAM", "EUR", "OCE"))

marginal_plot <- ggplot(me_plot_df_clip, aes(x = WaterDeficit, y = MarginalEffect,
                                             color = Continent, fill=Continent,
                                             group = interaction(Continent, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = expression(Marginal~effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "", fill =""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.45, 0.22), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.key.spacing.x = unit(1, "pt"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=continent_colors)+
  scale_fill_manual(values=continent_colors)

print(marginal_plot)

marginal_continent_waterdeficit_NPP <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_continent_waterdeficit_NPP)




# City size
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]


merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))
merged_dt_WaterDeficit$CitySize <- factor(merged_dt_WaterDeficit$CitySize, levels = c("Small",
                                                                                      "Medium",
                                                                                      "Large"))

model_urban <- felm(NPP ~ 1 + 
                      WaterDeficit + 
                      WaterDeficit_2 +
                      urban +
                      WaterDeficit:CitySize + 
                      WaterDeficit_2:CitySize +
                      urban:CitySize +
                      WaterDeficit:urban +
                      WaterDeficit_2:urban +
                      WaterDeficit:urban:CitySize +
                      WaterDeficit_2:urban:CitySize +
                      Seasonalityindex + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    | CityID + year | 0 | CityID, 
                    data = merged_dt_WaterDeficit)

summary(model_urban)

output_file <- "NPP_urban_in_out_citysize.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)


# Compute effects, marginal effects and confidence intervals
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "WaterDeficit",
  "WaterDeficit_2",
  "urbanin",
  "WaterDeficit:CitySizeMedium",
  "WaterDeficit:CitySizeLarge",
  "WaterDeficit_2:CitySizeMedium",
  "WaterDeficit_2:CitySizeLarge",
  "urbanin:CitySizeMedium",
  "urbanin:CitySizeLarge",
  "WaterDeficit:urbanin",
  "WaterDeficit_2:urbanin",
  "WaterDeficit:urbanin:CitySizeMedium",
  "WaterDeficit:urbanin:CitySizeLarge",
  "WaterDeficit_2:urbanin:CitySizeMedium",
  "WaterDeficit_2:urbanin:CitySizeLarge"
)

b <- beta_all[v_names]

citysize_colors <- c(
  "Small" = "#95a5a6",
  "Medium" = "#2ecc71",
  "Large" = "#e74c3c"
)

wd_seq <- seq(0, 2400, length.out = 200)
city_sizes <- c("Small","Medium","Large")

plot_df     <- data.frame()
me_plot_df  <- data.frame()

for(cs in city_sizes){
  for(wd in wd_seq){
    if(cs=="Small"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0
      )
    }else if(cs=="Medium"){
      g_non <- c(
        wd, wd^2, 0,
        wd,0,
        wd^2,0,
        0,0,
        0,0,
        0,0,
        0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        1,0,
        2*wd,0,
        0,0,
        0,0,
        0,0,
        0,0
      )
    }else if(cs=="Large"){
      g_non <- c(
        wd, wd^2, 0,
        0,wd,
        0,wd^2,
        0,0,
        0,0,
        0,0,
        0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,1,
        0,2*wd,
        0,0,
        0,0,
        0,0,
        0,0
      )
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
    me_non      <- as.numeric(crossprod(gme_non, b))
    se_me_non   <- sqrt(as.numeric(t(gme_non) %*% vcov_mat[v_names, v_names] %*% gme_non))
    me_ci_low_non  <- me_non - 1.96*se_me_non
    me_ci_high_non <- me_non + 1.96*se_me_non

    if(cs=="Small"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,
        0,0,
        0,0,
        wd, wd^2,
        0,0,
        0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,
        0,0,
        0,0,
        1, 2*wd,
        0,0,
        0,0
      )
    }else if(cs=="Medium"){
      g_urb <- c(
        wd, wd^2, 1,
        wd,0,
        wd^2,0,
        1,0,
        wd, wd^2,
        wd,0,
        wd^2,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        1,0,
        2*wd,0,
        0,0,
        1, 2*wd,
        1,0,
        2*wd,0
      )
    }else if(cs=="Large"){
      g_urb <- c(
        wd, wd^2, 1,
        0,wd,
        0,wd^2,
        0,1,
        wd, wd^2,
        0,wd,
        0,wd^2
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,1,
        0,2*wd,
        0,0,
        1, 2*wd,
        0,1,
        0,2*wd
      )
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
    me_urb      <- as.numeric(crossprod(gme_urb, b))
    se_me_urb   <- sqrt(as.numeric(t(gme_urb) %*% vcov_mat[v_names, v_names] %*% gme_urb))
    me_ci_low_urb  <- me_urb - 1.96*se_me_urb
    me_ci_high_urb <- me_urb + 1.96*se_me_urb

    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                CitySize=cs,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                CitySize=cs,
                                LandType="Urban"))
    
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_non,
                                   ci_low=me_ci_low_non,
                                   ci_high=me_ci_high_non,
                                   CitySize=cs,
                                   LandType="Non-urban"))
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_urb,
                                   ci_low=me_ci_low_urb,
                                   ci_high=me_ci_high_urb,
                                   CitySize=cs,
                                   LandType="Urban"))
  }
}

citysize_limits <- tibble(
  CitySize = c("Small","Medium","Large"),
  min_WD = c(
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$CitySize=="Small",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$CitySize=="Medium",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$CitySize=="Large",]$WaterDeficit,na.rm=T)
  ),
  max_WD = c(
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$CitySize=="Small",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$CitySize=="Medium",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$CitySize=="Large",]$WaterDeficit,na.rm=T)
  )
)
print(citysize_limits)

plot_df_clip <- plot_df %>%
  left_join(citysize_limits, by = "CitySize") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

me_plot_df_clip <- me_plot_df %>%
  left_join(citysize_limits, by = "CitySize") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)


# Plot effect curve
plot_df_clip$CitySize <- factor(plot_df_clip$CitySize,
                                levels = c("Small","Medium","Large"))

effect_plot <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC,
                                        color = CitySize, fill=CitySize,
                                        group = interaction(CitySize, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"),  guide = "none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = expression(Effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "",
    fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600),
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
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=citysize_colors) +
  scale_fill_manual(values=citysize_colors)

print(effect_plot)

# Plot distribution histogram
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=CitySize), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.011), breaks = seq(0, 0.011, by =0.004)) +
  labs(x = "Water deficit (mm)", y = "") +
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3), 
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 6),
    legend.position = "none"
  ) +
  scale_fill_manual(values = citysize_colors)

print(hist_plot)

effect_citysize_waterdeficit_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_citysize_waterdeficit_NPP)


# Plot marginal effect curve
me_plot_df_clip$CitySize <- factor(me_plot_df_clip$CitySize,
                                   levels = c("Small","Medium","Large"))

marginal_plot <- ggplot(me_plot_df_clip, aes(x = WaterDeficit, y = MarginalEffect,
                                             color = CitySize, fill=CitySize,
                                             group = interaction(CitySize, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = expression(Marginal~effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "", fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.55, 0.9), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=citysize_colors)+
  scale_fill_manual(values=citysize_colors)

print(marginal_plot)

marginal_citysize_waterdeficit_NPP <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_citysize_waterdeficit_NPP)



# Development level
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]
merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))

merged_dt_WaterDeficit$DevelopmentLevel <- factor(merged_dt_WaterDeficit$DevelopmentLevel, levels = c("developing", "developed"))

model_urban <- felm(NPP ~ 1 + 
                      WaterDeficit + 
                      WaterDeficit_2 +
                      urban +
                      WaterDeficit:DevelopmentLevel + 
                      WaterDeficit_2:DevelopmentLevel +
                      urban:DevelopmentLevel +
                      WaterDeficit:urban +
                      WaterDeficit_2:urban +
                      WaterDeficit:urban:DevelopmentLevel +
                      WaterDeficit_2:urban:DevelopmentLevel +
                      Seasonalityindex + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    | CityID + year | 0 | CityID, 
                    data = merged_dt_WaterDeficit)
summary(model_urban)

output_file <- "NPP_urban_in_out_level.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# Compute effects, marginal effects and confidence intervals
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "WaterDeficit",
  "WaterDeficit_2",
  "urbanin",
  "WaterDeficit:DevelopmentLeveldeveloped",
  "WaterDeficit_2:DevelopmentLeveldeveloped",
  "urbanin:DevelopmentLeveldeveloped",
  "WaterDeficit:urbanin",
  "WaterDeficit_2:urbanin",
  "WaterDeficit:urbanin:DevelopmentLeveldeveloped",
  "WaterDeficit_2:urbanin:DevelopmentLeveldeveloped"
)

b <- beta_all[v_names]

dev_colors <- c(
  "Global North" = "#0071bc", 
  "Global South" = "#d95218"
)

wd_seq <- seq(0, 2400, length.out = 200)
dev_levels <- c("developing","developed")

plot_df     <- data.frame()
me_plot_df  <- data.frame()

for(dl in dev_levels){
  for(wd in wd_seq){
    if(dl=="developing"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,0,
        0,0,
        0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,0,
        0,0,
        0,0
      )
    }else if(dl=="developed"){
      g_non <- c(
        wd, wd^2, 0,
        wd, wd^2, 0,
        0,0,
        0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        1, 2*wd, 0,
        0,0,
        0,0
      )
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
    me_non      <- as.numeric(crossprod(gme_non, b))
    se_me_non   <- sqrt(as.numeric(t(gme_non) %*% vcov_mat[v_names, v_names] %*% gme_non))
    me_ci_low_non  <- me_non - 1.96*se_me_non
    me_ci_high_non <- me_non + 1.96*se_me_non

    if(dl=="developing"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,
        wd, wd^2,
        0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,0,
        1, 2*wd,
        0,0
      )
    }else if(dl=="developed"){
      g_urb <- c(
        wd, wd^2, 1,
        wd, wd^2, 1,
        wd, wd^2,
        wd, wd^2
      )
      gme_urb <- c(
        1, 2*wd, 0,
        1, 2*wd, 0,
        1, 2*wd,
        1, 2*wd
      )
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
    me_urb      <- as.numeric(crossprod(gme_urb, b))
    se_me_urb   <- sqrt(as.numeric(t(gme_urb) %*% vcov_mat[v_names, v_names] %*% gme_urb))
    me_ci_low_urb  <- me_urb - 1.96*se_me_urb
    me_ci_high_urb <- me_urb + 1.96*se_me_urb
    
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                DevelopmentLevel=dl,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                DevelopmentLevel=dl,
                                LandType="Urban"))

    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_non,
                                   ci_low=me_ci_low_non,
                                   ci_high=me_ci_high_non,
                                   DevelopmentLevel=dl,
                                   LandType="Non-urban"))
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_urb,
                                   ci_low=me_ci_low_urb,
                                   ci_high=me_ci_high_urb,
                                   DevelopmentLevel=dl,
                                   LandType="Urban"))
  }
}

dev_limits <- tibble(
  DevelopmentLevel = c("developing","developed"),
  min_WD = c(
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$DevelopmentLevel=="developing",]$WaterDeficit,na.rm=T),
    min(merged_dt_WaterDeficit[merged_dt_WaterDeficit$DevelopmentLevel=="developed",]$WaterDeficit,na.rm=T)
  ),
  max_WD = c(
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$DevelopmentLevel=="developing",]$WaterDeficit,na.rm=T),
    max(merged_dt_WaterDeficit[merged_dt_WaterDeficit$DevelopmentLevel=="developed",]$WaterDeficit,na.rm=T)
  )
)
print(dev_limits)

plot_df_clip <- plot_df %>%
  left_join(dev_limits, by = "DevelopmentLevel") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

me_plot_df_clip <- me_plot_df %>%
  left_join(dev_limits, by = "DevelopmentLevel") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)


plot_df_clip$DevelopmentLevel <- factor(plot_df_clip$DevelopmentLevel, 
                                        levels = c("developing","developed"),
                                        labels = c("Global South","Global North")) 
me_plot_df_clip$DevelopmentLevel <- factor(me_plot_df_clip$DevelopmentLevel,
                                           levels = c("developing","developed"),
                                           labels = c("Global South","Global North"))


# Plot effect curve
effect_plot <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC,
                                        color = DevelopmentLevel, fill=DevelopmentLevel,
                                        group = interaction(DevelopmentLevel, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"),  guide = "none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = expression(Effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "",
    fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600),
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
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=dev_colors) +
  scale_fill_manual(values=dev_colors)

print(effect_plot)


# Plot distribution histogram
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=DevelopmentLevel), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.01), breaks = seq(0, 0.01, by =0.004)) +
  labs(x = "Water deficit (mm)", y = "") +
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3), 
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 6),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c( "developed" = "#0071bc", 
                                "developing" = "#d95218"))

print(hist_plot)

effect_dev_waterdeficit_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_dev_waterdeficit_NPP)

# Plot marginal effect curve
marginal_plot <- ggplot(me_plot_df_clip, aes(x = WaterDeficit, y = MarginalEffect,
                                             color = DevelopmentLevel, fill=DevelopmentLevel,
                                             group = interaction(DevelopmentLevel, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = expression(Marginal~effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "", fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.35, 0.9), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=dev_colors)+
  scale_fill_manual(values=dev_colors)

print(marginal_plot)

marginal_dev_waterdeficit_NPP <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_dev_waterdeficit_NPP)


# Country
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]

merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))

merged_dt_WaterDeficit$country <- factor(merged_dt_WaterDeficit$country, 
                                         levels = c("United States of America",
                                                    "Europe",
                                                    "China",
                                                    "India",
                                                    "Brazil", "other"))

model_urban <- felm(NPP ~ 1 + 
                      WaterDeficit + 
                      WaterDeficit_2 +
                      urban +
                      WaterDeficit:country + 
                      WaterDeficit_2:country +
                      urban:country +
                      WaterDeficit:urban +
                      WaterDeficit_2:urban +
                      WaterDeficit:urban:country +
                      WaterDeficit_2:urban:country +
                      Seasonalityindex + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    | CityID + year | 0 | CityID, 
                    data = merged_dt_WaterDeficit)

summary(model_urban)

output_file <- "NPP_urban_in_out_country.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# Compute effects, marginal effects and confidence intervals
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "WaterDeficit",
  "WaterDeficit_2",
  "urbanin",
  "WaterDeficit:countryEurope",
  "WaterDeficit:countryChina",
  "WaterDeficit:countryIndia",
  "WaterDeficit:countryBrazil",
  "WaterDeficit:countryother",
  "WaterDeficit_2:countryEurope",
  "WaterDeficit_2:countryChina",
  "WaterDeficit_2:countryIndia",
  "WaterDeficit_2:countryBrazil",
  "WaterDeficit_2:countryother",
  "urbanin:countryEurope",
  "urbanin:countryChina",
  "urbanin:countryIndia",
  "urbanin:countryBrazil",
  "urbanin:countryother",
  "WaterDeficit:urbanin",
  "WaterDeficit_2:urbanin",
  "WaterDeficit:urbanin:countryEurope",
  "WaterDeficit:urbanin:countryChina",
  "WaterDeficit:urbanin:countryIndia",
  "WaterDeficit:urbanin:countryBrazil",
  "WaterDeficit:urbanin:countryother",
  "WaterDeficit_2:urbanin:countryEurope",
  "WaterDeficit_2:urbanin:countryChina",
  "WaterDeficit_2:urbanin:countryIndia",
  "WaterDeficit_2:urbanin:countryBrazil",
  "WaterDeficit_2:urbanin:countryother"
)

b <- beta_all[v_names]

wd_seq <- seq(0, 2400, length.out = 200)
cntry_levels <- c("United States","Europe (excluding Russia)","China","India","Brazil","other")

plot_df     <- data.frame()
me_plot_df  <- data.frame()

for(cnt in cntry_levels){
  for(wd in wd_seq){
    if(cnt=="United States"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(cnt=="Europe (excluding Russia)"){
      g_non <- c(
        wd, wd^2, 0,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        1,0,0,0,0,
        2*wd,0,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(cnt=="China"){
      g_non <- c(
        wd, wd^2, 0,
        0,wd,0,0,0,
        0,wd^2,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,1,0,0,0,
        0,2*wd,0,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(cnt=="India"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,wd,0,0,
        0,0,wd^2,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,1,0,0,
        0,0,2*wd,0,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(cnt=="Brazil"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,0,wd,0,
        0,0,0,wd^2,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,0,1,0,
        0,0,0,2*wd,0,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(cnt=="other"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,0,0,wd,
        0,0,0,0,wd^2,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_non <- c(
        1, 2*wd, 0,
        0,0,0,0,1,
        0,0,0,0,2*wd,
        0,0,0,0,0,
        0,0,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
    me_non      <- as.numeric(crossprod(gme_non, b))
    se_me_non   <- sqrt(as.numeric(t(gme_non) %*% vcov_mat[v_names, v_names] %*% gme_non))
    me_ci_low_non  <- me_non - 1.96*se_me_non
    me_ci_high_non <- me_non + 1.96*se_me_non

    if(cnt=="United States"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        wd, wd^2,
        0,0,0,0,0,
        0,0,0,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        1, 2*wd,
        0,0,0,0,0,
        0,0,0,0,0
      )
    }else if(cnt=="Europe (excluding Russia)"){
      g_urb <- c(
        wd, wd^2, 1,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
        1,0,0,0,0,
        wd, wd^2,
        wd,0,0,0,0,
        wd^2,0,0,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        1,0,0,0,0,
        2*wd,0,0,0,0,
        0,0,0,0,0,
        1, 2*wd,
        1,0,0,0,0,
        2*wd,0,0,0,0
      )
    }else if(cnt=="China"){
      g_urb <- c(
        wd, wd^2, 1,
        0,wd,0,0,0,
        0,wd^2,0,0,0,
        0,1,0,0,0,
        wd, wd^2,
        0,wd,0,0,0,
        0,wd^2,0,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,1,0,0,0,
        0,2*wd,0,0,0,
        0,0,0,0,0,
        1, 2*wd,
        0,1,0,0,0,
        0,2*wd,0,0,0
      )
    }else if(cnt=="India"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,wd,0,0,
        0,0,wd^2,0,0,
        0,0,1,0,0,
        wd, wd^2,
        0,0,wd,0,0,
        0,0,wd^2,0,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,1,0,0,
        0,0,2*wd,0,0,
        0,0,0,0,0,
        1, 2*wd,
        0,0,1,0,0,
        0,0,2*wd,0,0
      )
    }else if(cnt=="Brazil"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,wd,0,
        0,0,0,wd^2,0,
        0,0,0,1,0,
        wd, wd^2,
        0,0,0,wd,0,
        0,0,0,wd^2,0
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,0,1,0,
        0,0,0,2*wd,0,
        0,0,0,0,0,
        1, 2*wd,
        0,0,0,1,0,
        0,0,0,2*wd,0
      )
    }else if(cnt=="other"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,0,wd,
        0,0,0,0,wd^2,
        0,0,0,0,1,
        wd, wd^2,
        0,0,0,0,wd,
        0,0,0,0,wd^2
      )
      gme_urb <- c(
        1, 2*wd, 0,
        0,0,0,0,1,
        0,0,0,0,2*wd,
        0,0,0,0,0,
        1, 2*wd,
        0,0,0,0,1,
        0,0,0,0,2*wd
      )
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
    me_urb      <- as.numeric(crossprod(gme_urb, b))
    se_me_urb   <- sqrt(as.numeric(t(gme_urb) %*% vcov_mat[v_names, v_names] %*% gme_urb))
    me_ci_low_urb  <- me_urb - 1.96*se_me_urb
    me_ci_high_urb <- me_urb + 1.96*se_me_urb

    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                Country=cnt,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(WaterDeficit=wd,
                                Effect_on_FVC=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                Country=cnt,
                                LandType="Urban"))

    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_non,
                                   ci_low=me_ci_low_non,
                                   ci_high=me_ci_high_non,
                                   Country=cnt,
                                   LandType="Non-urban"))
    me_plot_df <- rbind(me_plot_df,
                        data.frame(WaterDeficit=wd,
                                   MarginalEffect=me_urb,
                                   ci_low=me_ci_low_urb,
                                   ci_high=me_ci_high_urb,
                                   Country=cnt,
                                   LandType="Urban"))
  }
}

country_limits <- tibble(
  Country = c("United States", "Europe (excluding Russia)", "China", "India", "Brazil", "other"),
  min_WD = c(
    min(merged_dt_WaterDeficit[country == "United States of America"]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "Europe" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "China" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "India" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "Brazil" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "other" ]$WaterDeficit, na.rm = TRUE)
  ),
  max_WD = c(
    max(merged_dt_WaterDeficit[country == "United States of America" ]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[country == "Europe" ]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[country == "China"]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[country == "India"]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[country == "Brazil"]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[country == "other"]$WaterDeficit, na.rm = TRUE)
  )
)

print(country_limits)

plot_df_clip <- plot_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

me_plot_df_clip <- me_plot_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

plot_df_clip$Country <- factor(plot_df_clip$Country, 
                               levels = c("United States", "Europe (excluding Russia)",  "China", "India", "Brazil" ,"other"),
                               labels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH")) 
me_plot_df_clip$Country <- factor(me_plot_df_clip$Country,
                                  levels = c("United States", "Europe (excluding Russia)",  "China", "India", "Brazil" ,"other"),
                                  labels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH"))

country_colors <- c(
  "US" = "#2ecc71",
  "EURxR" = "#3498db",
  "CHN" = "#CB0505",
  "IND" = "#f39c12",
  "BRA" = "#4A1054",
  "OTH" = "#95a5a6"
)

# Plot effect curve
plot_df_clip$Country <- factor(plot_df_clip$Country, levels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH"))

effect_plot <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC,
                                        color = Country, fill=Country,
                                        group = interaction(Country, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"),  guide = "none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = expression(Effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "",
    fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600),
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-2800, 1),
                     breaks = seq(-2800, 1, by = 700)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.49, 0.23), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.key.spacing.x = unit(1, "pt"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=country_colors) +
  scale_fill_manual(values=country_colors)

print(effect_plot)


# Plot distribution histogram
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=country), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.025), breaks = seq(0, 0.025, by =0.01)) +
  labs(x = "Water deficit (mm)", y = "") +
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.3), 
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 6),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c(
    "United States of America" = "#2ecc71",
    "Europe" = "#3498db",
    "China" = "#CB0505",
    "India" = "#f39c12",
    "Brazil" = "#4A1054",
    "other" = "#95a5a6"
  ))

print(hist_plot)

effect_country_waterdeficit_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_country_waterdeficit_NPP)

# Plot marginal effect curve
me_plot_df_clip$Country <- factor(me_plot_df_clip$Country, levels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH"))


marginal_plot <- ggplot(me_plot_df_clip, aes(x = WaterDeficit, y = MarginalEffect,
                                             color = Country, fill=Country,
                                             group = interaction(Country, LandType))) +
  geom_ribbon(aes(ymin=ci_low, ymax=ci_high), alpha=0.1, colour=NA) +
  geom_line(aes(linetype = LandType, linewidth = LandType)) +
  geom_hline(yintercept =0, linetype="dashed", color="gray50") +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() +
  labs(
    x = "",
    y = expression(Marginal~effect~on~NPP~(g~C~m^-2~yr^-1)),
    color = "", fill = ""
  ) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.62, 0.25), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),   
    legend.key.spacing.x = unit(1, "pt"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=country_colors)+
  scale_fill_manual(values=country_colors)

print(marginal_plot)

marginal_country_waterdeficit_NPP <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_country_waterdeficit_NPP)





################################################ Composite figure ###############################################################
crop_bottom <- function(plot, crop_ratio = 0.15) {
  ggdraw() +
    draw_plot(plot, x = 0, y = crop_ratio, width = 1, height = 1 - crop_ratio) +
    theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"))
}

ii_crop  <- crop_bottom(effect_waterdeficit_FVC, crop_ratio = -0.1)
kk_crop  <- crop_bottom(effect_continent_waterdeficit_FVC, crop_ratio = -0.1)
oo_crop  <- crop_bottom(effect_climate_waterdeficit_FVC, crop_ratio = -0.1)

oo5_crop <- crop_bottom(effect_citysize_waterdeficit_FVC, crop_ratio = -0.1)
oo2_crop <- crop_bottom(effect_dev_waterdeficit_FVC, crop_ratio = -0.1)
oo4_crop <- crop_bottom(effect_country_waterdeficit_FVC, crop_ratio = -0.1)

jj_crop  <- crop_bottom(effect_waterdeficit_NPP, crop_ratio = -0.1)
mm_crop  <- crop_bottom(effect_continent_waterdeficit_NPP, crop_ratio = -0.1)
pp_crop  <- crop_bottom(effect_climate_waterdeficit_NPP, crop_ratio = -0.1)

pp5_crop <- crop_bottom(effect_citysize_waterdeficit_NPP, crop_ratio = -0.1)
pp2_crop <- crop_bottom(effect_dev_waterdeficit_NPP, crop_ratio = -0.1)
pp4_crop <- crop_bottom(effect_country_waterdeficit_NPP, crop_ratio = -0.1)

fig_total <- plot_grid(
  ii_crop,  kk_crop,  oo_crop,    
  oo5_crop, oo2_crop, oo4_crop,   
  jj_crop,  mm_crop,  pp_crop,   
  pp5_crop, pp2_crop, pp4_crop,   
  nrow = 4,
  ncol = 3,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1,1,1),
  rel_heights = c(1,1,1,1),
  
  labels = letters[1:12],
  label_size = 8,
  label_fontfamily = "Times",
  label_fontface = "bold",
  label_x = rep(0.11, 12),
  label_y = c(rep(1.03,6), rep(1.06,6)),
  hjust = 0,
  vjust = 1
)

fig_total <- cowplot::ggdraw(fig_total) +
  theme(plot.margin = margin(10, 5, 5, 5, unit = "mm"))

print(fig_total)

leg_grob <- get_legend(effect_plot_legend)
fig_total_with_legend <- ggdraw() +
  draw_plot(fig_total, x = 0, y = 0, width = 1, height = 0.999) +
  draw_grob(
    leg_grob,
    x = 0.07,   
    y = 0.915,    
    width = 0.98,
    height = 0.07 
  )

print(fig_total_with_legend)

ggsave("water_deficit_effect.pdf",
       plot = fig_total_with_legend, width = 160, height = 200, units = "mm", dpi=600)


################################################ Composite figure ###############################################################
ii_crop  <- crop_bottom(marginal_waterdeficit_FVC, crop_ratio = -0.05)
kk_crop  <- crop_bottom(marginal_continent_waterdeficit_FVC, crop_ratio = -0.05)
oo_crop  <- crop_bottom(marginal_climate_waterdeficit_FVC, crop_ratio = -0.05)

oo5_crop <- crop_bottom(marginal_citysize_waterdeficit_FVC, crop_ratio = -0.05)
oo2_crop <- crop_bottom(marginal_dev_waterdeficit_FVC, crop_ratio = -0.05)
oo4_crop <- crop_bottom(marginal_country_waterdeficit_FVC, crop_ratio = -0.05)

jj_crop  <- crop_bottom(marginal_waterdeficit_NPP, crop_ratio = -0.05)
mm_crop  <- crop_bottom(marginal_continent_waterdeficit_NPP, crop_ratio = -0.05)
pp_crop  <- crop_bottom(marginal_climate_waterdeficit_NPP, crop_ratio = -0.05)

pp5_crop <- crop_bottom(marginal_citysize_waterdeficit_NPP, crop_ratio = -0.05)
pp2_crop <- crop_bottom(marginal_dev_waterdeficit_NPP, crop_ratio = -0.05)
pp4_crop <- crop_bottom(marginal_country_waterdeficit_NPP, crop_ratio = -0.05)

fig_total <- plot_grid(
  ii_crop,  kk_crop,  oo_crop,    
  oo5_crop, oo2_crop, oo4_crop,   
  jj_crop,  mm_crop,  pp_crop,    
  pp5_crop, pp2_crop, pp4_crop,   
  nrow = 4,
  ncol = 3,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1,1,1),
  rel_heights = c(1,1,1,1),
  
  labels = letters[1:12],
  label_size = 8,
  label_fontfamily = "Times",
  label_fontface = "bold",
  label_x = rep(0.11,12),
  label_y = c(rep(1.03,6), rep(1.13,6)),
  hjust = 0,
  vjust = 1
)

fig_total <- cowplot::ggdraw(fig_total) +
  theme(plot.margin = margin(10, 5, 5, 5, unit = "mm"))

print(fig_total)


leg_grob <- get_legend(marginal_plot_legend)
fig_total_with_legend <- ggdraw() +
  draw_plot(fig_total, x = 0, y = 0, width = 1, height = 0.999) +
  draw_grob(
    leg_grob,
    x = 0.07,    
    y = 0.915,    
    width = 0.98,
    height = 0.07 
  )

print(fig_total_with_legend)

ggsave("water_deficit_marginal_effect.pdf",
       plot = fig_total_with_legend, width = 160, height = 210, units = "mm", dpi=600)




