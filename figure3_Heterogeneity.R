
# Heterogeneity

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
#                   Three-way interactive quadratic regression
###########################################################################

setwd("XXX") # please set your folder path
merged_dt <- fread("Dataset_matched_buffer10km.csv")

merged_dt$Latitude <- abs(merged_dt$Latitude)
merged_dt$WaterDeficit_2 <- merged_dt$WaterDeficit * merged_dt$WaterDeficit
merged_dt$Seasonalityindex_2 <- merged_dt$Seasonalityindex * merged_dt$Seasonalityindex

###########################################################################
#                   climate zones: FVC ~ Water deficit
###########################################################################
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
    }else if(cz=="Tropical"){
      g_non <- c(
        wd, wd^2, 0,
        wd,0,0,
        wd^2,0,0,
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
    }else if(cz=="Arid"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,wd,
        0,0,wd^2,
        0,0,0,
        0,0,0,0,0,
        0,0,0
      )
    }
    
    pred_non <- as.numeric(crossprod(g_non, b))
    se_non   <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
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
    }
    
    pred_urb <- as.numeric(crossprod(g_urb, b))
    se_urb   <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
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
plot_df_clip$ClimateZone <- factor(plot_df_clip$ClimateZone, levels = c( "Tropical", "Temperate",  "Cold","Arid"))
# Plot effect curve
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=climate_colors) +
  scale_fill_manual(values=climate_colors) 
print(effect_plot)
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=ClimateZone), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.011), breaks = seq(0, 0.011, by =0.004)) +
  labs(x = "", y = "") +
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
effect_climate_waterdeficit_FVC_nolegend <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_climate_waterdeficit_FVC_nolegend)
# Plot effect curve legend
effect_plot_legend <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC,
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
  scale_y_continuous(limits = c(-0.8, 0.8),
                     breaks = seq(-0.8, 0.8, by = 0.2),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), 
    axis.line.x.top = element_blank(), 
    axis.line.y.right = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6),  
    axis.text.x = element_blank(), 
    axis.ticks = element_line(linewidth = 0.3),   
    axis.title.y = element_text(size = 6),  
    legend.position = c(0.2, 0.95),
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),
    legend.key.spacing.x = unit(1, "pt"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, ncol = 2, keywidth = 0.5, keyheight = 0.4)) +
  scale_color_manual(values=climate_colors) +
  scale_fill_manual(values=climate_colors) 
leg_grob <- get_legend(effect_plot_legend)
effect_climate_waterdeficit_FVC <- ggdraw() +
  draw_plot(effect_climate_waterdeficit_FVC_nolegend, x = 0, y = 0, width = 1, height = 0.999) +
  draw_grob(
    leg_grob,
    x = 0.49,    
    y = 0.47,    
    width = 0.98,
    height = 0.07 
  )
print(effect_climate_waterdeficit_FVC)


###########################################################################
#                   continents: FVC ~ Water deficit
###########################################################################
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
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
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
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
plot_df_clip$Continent <- factor(plot_df_clip$Continent, 
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
plot_df_clip$Continent <- factor(plot_df_clip$Continent,
                                 levels = c("ASIA", "AFR", "NAM", "SAM", "EUR", "OCE"))
# Plot effect curve
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
    legend.position = "none",
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
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=continent), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.023), breaks = seq(0, 0.023, by = 0.01)) +
  labs(x = "", y = "") +
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
effect_continent_waterdeficit_FVC_nolegend <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_continent_waterdeficit_FVC_nolegend)
# Plot effect curve legend (continents)
effect_plot_legend <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC,
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
  guides(color = guide_legend(nrow = 3, ncol = 2,keywidth = 0.5, keyheight = 0.4)) +
  scale_color_manual(values=continent_colors) +
  scale_fill_manual(values=continent_colors)
leg_grob <- get_legend(effect_plot_legend)
effect_continent_waterdeficit_FVC_nocity <- ggdraw() +
  draw_plot(effect_continent_waterdeficit_FVC_nolegend, x = 0, y = 0, width = 1, height = 0.999) +
  draw_grob(
    leg_grob,
    x = 0.42,    
    y = 0.50,    
    width = 0.98,
    height = 0.07 
  )
print(effect_continent_waterdeficit_FVC_nocity)
# Plot effect curve legend (urban vs. Non-urban)
effect_plot_legend2 <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC, color = LandType, linetype = LandType, linewidth = LandType)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = ci_low, ymax = ci_high, fill = LandType),
    alpha=0.1, color = NA
  ) +
  scale_linetype_manual(values = c("Non-urban"="dashed", "Urban"="solid"), 
                        labels = c("Non-urban effect curve with 95% CI", 
                                   "Urban effect curve with 95% CI"),
                        guide="none") +
  scale_linewidth_manual(values = c("Non-urban" = 0.25, "Urban" = 0.4), guide="none") +
  theme_classic() + 
  labs(
    x = "Water deficit (mm)",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by = 0.02), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(name = "Effect on FVC", limits = c(-0.2, 0.25),
                     breaks = seq(-0.2, 0.25, by = 0.1),
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
    legend.position = c(0.45, 0.63),
    legend.text = element_text(size = 6, margin = margin(l = 5, unit = "pt")),   
    legend.key.spacing.x = unit(1, "pt"),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = c(
    "Urban" = "gray50",
    "Non-urban" = "gray50"
  )) +
  guides(
    color = guide_legend(
      nrow = 2, ncol = 1,
      keywidth = 0.9,
      keyheight = 0.4,
      override.aes = list(
        linetype = c("Non-urban"="dashed", "Urban"="solid"),
        linewidth = c("Non-urban" = 0.4, "Urban" = 0.4),
        fill = c("Urban" = "gray50",
                 "Non-urban" = "gray50")
      )
    )
  )
leg_grob <- get_legend(effect_plot_legend2)
effect_continent_waterdeficit_FVC <- ggdraw() +
  draw_plot(effect_continent_waterdeficit_FVC_nocity, x = 0, y = 0, width = 1, height = 0.999) +
  draw_grob(
    leg_grob,
    x = 0.12,    
    y = 0.96,    
    width = 0.98,
    height = 0.07 
  )
print(effect_continent_waterdeficit_FVC)


###########################################################################
#                   City sizes: FVC ~ Water deficit
###########################################################################
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
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
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
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
plot_df_clip$CitySize <- factor(plot_df_clip$CitySize,
                                levels = c("Small","Medium","Large"))
# Plot effect curve
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=citysize_colors) +
  scale_fill_manual(values=citysize_colors)
print(effect_plot)
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=CitySize), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.011), breaks = seq(0, 0.011, by =0.004)) +
  labs(x = "", y = "") +
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
effect_citysize_waterdeficit_FVC_nolegend <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_citysize_waterdeficit_FVC_nolegend)
# Plot effect curve legend
effect_plot_legend <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC,
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
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.4)) +
  scale_color_manual(values=citysize_colors) +
  scale_fill_manual(values=citysize_colors)
leg_grob <- get_legend(effect_plot_legend)
effect_citysize_waterdeficit_FVC <- ggdraw() +
  draw_plot(effect_citysize_waterdeficit_FVC_nolegend, x = 0, y = 0, width = 1, height = 0.999) +
  draw_grob(
    leg_grob,
    x = 0.55,    
    y = 0.49,    
    width = 0.98,
    height = 0.07 
  )
print(effect_citysize_waterdeficit_FVC)


###########################################################################
#                   development levels: FVC ~ Water deficit
###########################################################################
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
for(dl in dev_levels){
  for(wd in wd_seq){
    if(dl=="developing"){
      g_non <- c(
        wd, wd^2, 0,
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    if(dl=="developing"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,
        wd, wd^2,
        0,0
      )
    }else if(dl=="developed"){
      g_urb <- c(
        wd, wd^2, 1,
        wd, wd^2, 1,
        wd, wd^2,
        wd, wd^2
      )
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
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
plot_df_clip$DevelopmentLevel <- factor(plot_df_clip$DevelopmentLevel, 
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=dev_colors) +
  scale_fill_manual(values=dev_colors)
print(effect_plot)
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill = DevelopmentLevel), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.01), breaks = seq(0, 0.01, by =0.004)) +
  labs(x = "", y = "") +
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
effect_dev_waterdeficit_FVC_nolegend <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_dev_waterdeficit_FVC_nolegend)
# Plot effect curve legend
effect_plot_legend <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC,
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
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.4)) +
  scale_color_manual(values=dev_colors) +
  scale_fill_manual(values=dev_colors)
leg_grob <- get_legend(effect_plot_legend)
effect_dev_waterdeficit_FVC <- ggdraw() +
  draw_plot(effect_dev_waterdeficit_FVC_nolegend, x = 0, y = 0, width = 1, height = 0.999) +
  draw_grob(
    leg_grob,
    x = 0.40,    
    y = 0.49,    
    width = 0.98,
    height = 0.07 
  )
print(effect_dev_waterdeficit_FVC)


###########################################################################
#                   nations/regions: FVC ~ Water deficit
###########################################################################
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]
merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))
merged_dt_WaterDeficit$country <- factor(merged_dt_WaterDeficit$country, 
                                         levels = c("United States of America",
                                                    "EuropeNR",
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
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "WaterDeficit",
  "WaterDeficit_2",
  "urbanin",
  "WaterDeficit:countryEuropeNR",
  "WaterDeficit:countryChina",
  "WaterDeficit:countryIndia",
  "WaterDeficit:countryBrazil",
  "WaterDeficit:countryother",
  "WaterDeficit_2:countryEuropeNR",
  "WaterDeficit_2:countryChina",
  "WaterDeficit_2:countryIndia",
  "WaterDeficit_2:countryBrazil",
  "WaterDeficit_2:countryother",
  "urbanin:countryEuropeNR",
  "urbanin:countryChina",
  "urbanin:countryIndia",
  "urbanin:countryBrazil",
  "urbanin:countryother",
  "WaterDeficit:urbanin",
  "WaterDeficit_2:urbanin",
  "WaterDeficit:urbanin:countryEuropeNR",
  "WaterDeficit:urbanin:countryChina",
  "WaterDeficit:urbanin:countryIndia",
  "WaterDeficit:urbanin:countryBrazil",
  "WaterDeficit:urbanin:countryother",
  "WaterDeficit_2:urbanin:countryEuropeNR",
  "WaterDeficit_2:urbanin:countryChina",
  "WaterDeficit_2:urbanin:countryIndia",
  "WaterDeficit_2:urbanin:countryBrazil",
  "WaterDeficit_2:urbanin:countryother"
)
b <- beta_all[v_names]
wd_seq <- seq(0, 2400, length.out = 200)
cntry_levels <- c("United States","EuropeNR","China","India","Brazil","other")
plot_df     <- data.frame()
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
    }else if(cnt=="EuropeNR"){
      g_non <- c(
        wd, wd^2, 0,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
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
    }else if(cnt=="EuropeNR"){
      g_urb <- c(
        wd, wd^2, 1,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
        1,0,0,0,0,
        wd, wd^2,
        wd,0,0,0,0,
        wd^2,0,0,0,0
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
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
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
  }
}
country_limits <- tibble(
  Country = c("United States", "EuropeNR", "China", "India", "Brazil", "other"),
  min_WD = c(
    min(merged_dt_WaterDeficit[country == "United States of America"]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "EuropeNR" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "China" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "India" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "Brazil" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "other" ]$WaterDeficit, na.rm = TRUE)
  ),
  max_WD = c(
    max(merged_dt_WaterDeficit[country == "United States of America" ]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[country == "EuropeNR" ]$WaterDeficit, na.rm = TRUE),
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
plot_df_clip$Country <- factor(plot_df_clip$Country, 
                               levels = c("United States", "EuropeNR",  "China", "India", "Brazil" ,"other"),
                               labels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH")) 
country_colors <- c(
  "US" = "#2ecc71",
  "EURxR" = "#3498db",
  "CHN" = "#CB0505",
  "IND" = "#f39c12",
  "BRA" = "#4A1054",
  "OTH" = "#95a5a6"
)
plot_df_clip$Country <- factor(plot_df_clip$Country, levels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH"))
# Plot effect curve
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")), 
    legend.key.spacing.x = unit(1, "pt"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 6, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +  
  scale_color_manual(values=country_colors) +
  scale_fill_manual(values=country_colors)
print(effect_plot)
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
    "EuropeNR" = "#3498db",
    "China" = "#CB0505",
    "India" = "#f39c12",
    "Brazil" = "#4A1054",
    "other" = "#95a5a6"
  ))
print(hist_plot)
effect_country_waterdeficit_FVC_nolegend <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_country_waterdeficit_FVC_nolegend)
# Plot effect curve legend
effect_plot_legend <- ggplot(plot_df_clip, aes(x = WaterDeficit, y = Effect_on_FVC,
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
  guides(color = guide_legend(nrow = 3, ncol = 2, keywidth = 0.5, keyheight = 0.4)) +
  scale_color_manual(values=country_colors) +
  scale_fill_manual(values=country_colors)
leg_grob <- get_legend(effect_plot_legend)
effect_country_waterdeficit_FVC <- ggdraw() +
  draw_plot(effect_country_waterdeficit_FVC_nolegend, x = 0, y = 0, width = 1, height = 0.999) +
  draw_grob(
    leg_grob,
    x = 0.55,    
    y = 0.48,    
    width = 0.98,
    height = 0.07 
  )
print(effect_country_waterdeficit_FVC)





###########################################################################
#                   climate zones: NPP ~ Water deficit
###########################################################################
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
for(cz in climates){
  for(wd in wd_seq){
    if(cz=="Temperate"){
      g_non <- c(wd, wd^2, 0,
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
    }else if(cz=="Cold"){
      g_non <- c(wd, wd^2, 0,
                 0,wd,0,
                 0,wd^2,0,
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    if(cz=="Temperate"){
      g_urb <- c(wd, wd^2, 1,
                 0,0,0,
                 0,0,0,
                 0,0,0,
                 wd, wd^2,
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
    }else if(cz=="Cold"){
      g_urb <- c(wd, wd^2, 1,
                 0,wd,0,
                 0,wd^2,0,
                 0,1,0,
                 wd, wd^2,
                 0,wd,0,
                 0,wd^2,0)
    }else if(cz=="Arid"){
      g_urb <- c(wd, wd^2, 1,
                 0,0,wd,
                 0,0,wd^2,
                 0,0,1,
                 wd, wd^2,
                 0,0,wd,
                 0,0,wd^2)
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
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
plot_df_clip$ClimateZone <- factor(plot_df_clip$ClimateZone, levels = c( "Tropical", "Temperate",  "Cold","Arid"))
# Plot effect curve
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
    y = "Effect on NPP",
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),   
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=climate_colors) +
  scale_fill_manual(values=climate_colors)
print(effect_plot)
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=ClimateZone), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.011), breaks = seq(0, 0.011, by =0.004)) +
  labs(x = "", y = "") +
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


###########################################################################
#                   continents: NPP ~ Water deficit
###########################################################################
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
    }else if(ct=="NorthAmerica"){
      g_non <- c(wd, wd^2, 0,
                 wd,0,0,0,0,
                 wd^2,0,0,0,0,
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
    }else if(ct=="Europe"){
      g_non <- c(wd, wd^2, 0,
                 0,0,wd,0,0,
                 0,0,wd^2,0,0,
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
    }else if(ct=="SouthAmerica"){
      g_non <- c(wd, wd^2, 0,
                 0,0,0,0,wd,
                 0,0,0,0,wd^2,
                 0,0,0,0,0,
                 0,0,
                 0,0,0,0,0,
                 0,0,0,0,0)
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    if(ct=="Asia"){
      g_urb <- c(wd, wd^2, 1,
                 0,0,0,0,0,
                 0,0,0,0,0,
                 0,0,0,0,0,
                 wd, wd^2,
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
    }else if(ct=="Africa"){
      g_urb <- c(wd, wd^2, 1,
                 0,wd,0,0,0,
                 0,wd^2,0,0,0,
                 0,1,0,0,0,
                 wd, wd^2,
                 0,wd,0,0,0,
                 0,wd^2,0,0,0)
    }else if(ct=="Europe"){
      g_urb <- c(wd, wd^2, 1,
                 0,0,wd,0,0,
                 0,0,wd^2,0,0,
                 0,0,1,0,0,
                 wd, wd^2,
                 0,0,wd,0,0,
                 0,0,wd^2,0,0)
    }else if(ct=="Oceania"){
      g_urb <- c(wd, wd^2, 1,
                 0,0,0,wd,0,
                 0,0,0,wd^2,0,
                 0,0,0,1,0,
                 wd, wd^2,
                 0,0,0,wd,0,
                 0,0,0,wd^2,0)
    }else if(ct=="SouthAmerica"){
      g_urb <- c(wd, wd^2, 1,
                 0,0,0,0,wd,
                 0,0,0,0,wd^2,
                 0,0,0,0,1,
                 wd, wd^2,
                 0,0,0,0,wd,
                 0,0,0,0,wd^2)
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
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
plot_df_clip$Continent <- factor(plot_df_clip$Continent, 
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
plot_df_clip$Continent <- factor(plot_df_clip$Continent,
                                 levels = c("ASIA", "AFR", "NAM", "SAM", "EUR", "OCE"))
# Plot effect curve
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
    y = "Effect on NPP",
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
    legend.position = "none",
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
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=continent), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.011), breaks = seq(0, 0.011, by =0.004)) +
  labs(x = "", y = "") +
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




###########################################################################
#                   city sizes: NPP ~ Water deficit
###########################################################################
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
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
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
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
plot_df_clip$CitySize <- factor(plot_df_clip$CitySize,
                                levels = c("Small","Medium","Large"))
# Plot effect curve
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
    y = "Effect on NPP",
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=citysize_colors) +
  scale_fill_manual(values=citysize_colors)
print(effect_plot)
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=CitySize), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.011), breaks = seq(0, 0.011, by =0.004)) +
  labs(x = "", y = "") +
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
###########################################################################
#                   development levels: NPP ~ Water deficit
###########################################################################
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
for(dl in dev_levels){
  for(wd in wd_seq){
    if(dl=="developing"){
      g_non <- c(
        wd, wd^2, 0,
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    if(dl=="developing"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,
        wd, wd^2,
        0,0
      )
    }else if(dl=="developed"){
      g_urb <- c(
        wd, wd^2, 1,
        wd, wd^2, 1,
        wd, wd^2,
        wd, wd^2
      )
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
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
plot_df_clip$DevelopmentLevel <- factor(plot_df_clip$DevelopmentLevel, 
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
    y = "Effect on NPP",
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=dev_colors) +
  scale_fill_manual(values=dev_colors)
print(effect_plot)
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill = DevelopmentLevel), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0, 0.01), breaks = seq(0, 0.01, by =0.004)) +
  labs(x = "", y = "") +
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


###########################################################################
#                   nations/regions: NPP ~ Water deficit
###########################################################################
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]
merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))
merged_dt_WaterDeficit$country <- factor(merged_dt_WaterDeficit$country, 
                                         levels = c("United States of America",
                                                    "EuropeNR",
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
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "WaterDeficit",
  "WaterDeficit_2",
  "urbanin",
  "WaterDeficit:countryEuropeNR",
  "WaterDeficit:countryChina",
  "WaterDeficit:countryIndia",
  "WaterDeficit:countryBrazil",
  "WaterDeficit:countryother",
  "WaterDeficit_2:countryEuropeNR",
  "WaterDeficit_2:countryChina",
  "WaterDeficit_2:countryIndia",
  "WaterDeficit_2:countryBrazil",
  "WaterDeficit_2:countryother",
  "urbanin:countryEuropeNR",
  "urbanin:countryChina",
  "urbanin:countryIndia",
  "urbanin:countryBrazil",
  "urbanin:countryother",
  "WaterDeficit:urbanin",
  "WaterDeficit_2:urbanin",
  "WaterDeficit:urbanin:countryEuropeNR",
  "WaterDeficit:urbanin:countryChina",
  "WaterDeficit:urbanin:countryIndia",
  "WaterDeficit:urbanin:countryBrazil",
  "WaterDeficit:urbanin:countryother",
  "WaterDeficit_2:urbanin:countryEuropeNR",
  "WaterDeficit_2:urbanin:countryChina",
  "WaterDeficit_2:urbanin:countryIndia",
  "WaterDeficit_2:urbanin:countryBrazil",
  "WaterDeficit_2:urbanin:countryother"
)
b <- beta_all[v_names]
wd_seq <- seq(0, 2400, length.out = 200)
cntry_levels <- c("United States","EuropeNR","China","India","Brazil","other")
plot_df     <- data.frame()
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
    }else if(cnt=="EuropeNR"){
      g_non <- c(
        wd, wd^2, 0,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
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
    }else if(cnt=="EuropeNR"){
      g_urb <- c(
        wd, wd^2, 1,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
        1,0,0,0,0,
        wd, wd^2,
        wd,0,0,0,0,
        wd^2,0,0,0,0
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
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
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
  }
}
country_limits <- tibble(
  Country = c("United States", "EuropeNR", "China", "India", "Brazil", "other"),
  min_WD = c(
    min(merged_dt_WaterDeficit[country == "United States of America"]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "EuropeNR" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "China" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "India" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "Brazil" ]$WaterDeficit, na.rm = TRUE),
    min(merged_dt_WaterDeficit[country == "other" ]$WaterDeficit, na.rm = TRUE)
  ),
  max_WD = c(
    max(merged_dt_WaterDeficit[country == "United States of America" ]$WaterDeficit, na.rm = TRUE),
    max(merged_dt_WaterDeficit[country == "EuropeNR" ]$WaterDeficit, na.rm = TRUE),
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
plot_df_clip$Country <- factor(plot_df_clip$Country, 
                               levels = c("United States", "EuropeNR",  "China", "India", "Brazil" ,"other"),
                               labels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH")) 
country_colors <- c(
  "US" = "#2ecc71",
  "EURxR" = "#3498db",
  "CHN" = "#CB0505",
  "IND" = "#f39c12",
  "BRA" = "#4A1054",
  "OTH" = "#95a5a6"
)
plot_df_clip$Country <- factor(plot_df_clip$Country, levels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH"))
# Plot effect curve
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
    y = "Effect on NPP",
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.key.spacing.x = unit(1, "pt"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=country_colors) +
  scale_fill_manual(values=country_colors)
print(effect_plot)
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
    "EuropeNR" = "#3498db",
    "China" = "#CB0505",
    "India" = "#f39c12",
    "Brazil" = "#4A1054",
    "other" = "#95a5a6"
  ))
print(hist_plot)
effect_country_waterdeficit_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_country_waterdeficit_NPP)



###########################################################################
#                   climate zones: FVC ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                         .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))
merged_dt_Seasonalityindex$ClimateZone <- factor(merged_dt_Seasonalityindex$ClimateZone, levels = c("Temperate","Tropical","Cold","Arid"))
model_urban <- felm(FVC ~ 1 + 
                      Seasonalityindex +                        
                      Seasonalityindex_2 +  
                      urban +  
                      Seasonalityindex:ClimateZone + 
                      Seasonalityindex_2:ClimateZone + 
                      urban:ClimateZone +  
                      Seasonalityindex:urban + 
                      Seasonalityindex_2:urban + 
                      Seasonalityindex:urban:ClimateZone +  
                      Seasonalityindex_2:urban:ClimateZone + 
                      WaterDeficit + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    |  CityID + year | 0 | CityID, 
                    data = merged_dt_Seasonalityindex)
summary(model_urban)
output_file <- "FVC_urban_in_out_climate.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "Seasonalityindex",
  "Seasonalityindex_2",
  "urbanin",
  "Seasonalityindex:ClimateZoneTropical",
  "Seasonalityindex:ClimateZoneCold",
  "Seasonalityindex:ClimateZoneArid",
  "Seasonalityindex_2:ClimateZoneTropical",
  "Seasonalityindex_2:ClimateZoneCold",
  "Seasonalityindex_2:ClimateZoneArid",
  "urbanin:ClimateZoneTropical",
  "urbanin:ClimateZoneCold",
  "urbanin:ClimateZoneArid",
  "Seasonalityindex:urbanin",
  "Seasonalityindex_2:urbanin",
  "Seasonalityindex:urbanin:ClimateZoneTropical",
  "Seasonalityindex:urbanin:ClimateZoneCold",
  "Seasonalityindex:urbanin:ClimateZoneArid",
  "Seasonalityindex_2:urbanin:ClimateZoneTropical",
  "Seasonalityindex_2:urbanin:ClimateZoneCold",
  "Seasonalityindex_2:urbanin:ClimateZoneArid"
)
b <- beta_all[v_names]
climate_colors <- c(
  "Temperate" = "#2ecc71",
  "Tropical"  = "#e74c3c",
  "Arid"      = "#95a5a6",
  "Cold"      = "#3498db"
)
wd_seq <- seq(0, 0.1, length.out = 200)
climates <- c("Temperate","Tropical","Cold","Arid")
plot_df <- data.frame()
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
    }else if(cz=="Tropical"){
      g_non <- c(
        wd, wd^2, 0,
        wd,0,0,
        wd^2,0,0,
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
    }else if(cz=="Arid"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,wd,
        0,0,wd^2,
        0,0,0,
        0,0,0,0,0,
        0,0,0
      )
    }
    
    pred_non <- as.numeric(crossprod(g_non, b))
    se_non   <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
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
    }
    
    pred_urb <- as.numeric(crossprod(g_urb, b))
    se_urb   <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_FVC=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                ClimateZone=cz,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_FVC=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                ClimateZone=cz,
                                LandType="Urban"))
  }
}
climate_limits <- tibble(
  ClimateZone = c("Temperate", "Tropical", "Arid",  "Cold"),
  min_WD = c(
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Temperate", ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Tropical", ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Arid", ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Cold", ]$Seasonalityindex, na.rm = TRUE)
  ),
  max_WD = c(
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Temperate", ]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Tropical", ]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Arid", ]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Cold", ]$Seasonalityindex, na.rm = TRUE)
  )
)
print(climate_limits)
plot_df_clip <- plot_df %>%
  left_join(climate_limits, by = "ClimateZone") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
plot_df_clip$ClimateZone <- factor(plot_df_clip$ClimateZone, levels = c( "Tropical", "Temperate",  "Cold","Arid"))
# Plot effect curve
effect_plot <- ggplot(plot_df_clip, aes(x = Seasonalityindex, y = Effect_on_FVC,
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=climate_colors) +
  scale_fill_manual(values=climate_colors) 
print(effect_plot)
hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill=ClimateZone), binwidth = 0.0001, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.03)) +
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by =200)) +
  labs(x = "", y = "") +
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
effect_climate_Seasonalityindex_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_climate_Seasonalityindex_FVC)


###########################################################################
#                   continents: FVC ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                         .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))
merged_dt_Seasonalityindex$continent <- factor(merged_dt_Seasonalityindex$continent, levels = c( "Asia",
                                                                                                 "NorthAmerica",
                                                                                                 "Africa",
                                                                                                 "Europe",
                                                                                                 "Oceania",
                                                                                                 "SouthAmerica"))
model_urban <- felm(FVC ~ 1 + 
                      Seasonalityindex + Seasonalityindex_2 +
                      urban +
                      Seasonalityindex:continent + 
                      Seasonalityindex_2:continent +
                      urban:continent +
                      Seasonalityindex:urban +
                      Seasonalityindex_2:urban +
                      Seasonalityindex:urban:continent +
                      Seasonalityindex_2:urban:continent +
                      WaterDeficit + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    |  CityID + year | 0 | CityID, 
                    data = merged_dt_Seasonalityindex)
summary(model_urban)
output_file <- "FVC_urban_in_out_continent.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "Seasonalityindex",
  "Seasonalityindex_2",
  "urbanin",
  "Seasonalityindex:continentNorthAmerica",
  "Seasonalityindex:continentAfrica",
  "Seasonalityindex:continentEurope",
  "Seasonalityindex:continentOceania",
  "Seasonalityindex:continentSouthAmerica",
  "Seasonalityindex_2:continentNorthAmerica",
  "Seasonalityindex_2:continentAfrica",
  "Seasonalityindex_2:continentEurope",
  "Seasonalityindex_2:continentOceania",
  "Seasonalityindex_2:continentSouthAmerica",
  "urbanin:continentNorthAmerica",
  "urbanin:continentAfrica",
  "urbanin:continentEurope",
  "urbanin:continentOceania",
  "urbanin:continentSouthAmerica",
  "Seasonalityindex:urbanin",
  "Seasonalityindex_2:urbanin",
  "Seasonalityindex:urbanin:continentNorthAmerica",
  "Seasonalityindex:urbanin:continentAfrica",
  "Seasonalityindex:urbanin:continentEurope",
  "Seasonalityindex:urbanin:continentOceania",
  "Seasonalityindex:urbanin:continentSouthAmerica",
  "Seasonalityindex_2:urbanin:continentNorthAmerica",
  "Seasonalityindex_2:urbanin:continentAfrica",
  "Seasonalityindex_2:urbanin:continentEurope",
  "Seasonalityindex_2:urbanin:continentOceania",
  "Seasonalityindex_2:urbanin:continentSouthAmerica"
)
b <- beta_all[v_names]
wd_seq <- seq(0, 0.1, length.out = 200)
continents <- c("Asia","NorthAmerica","Africa","Europe","Oceania","SouthAmerica")
plot_df     <- data.frame()
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
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
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_FVC=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                Continent=ct,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_FVC=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                Continent=ct,
                                LandType="Urban"))
  }
}
continent_limits <- tibble(
  Continent = c("Asia","NorthAmerica","Africa","Europe","Oceania","SouthAmerica"),
  min_WD = c(
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Asia",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="NorthAmerica",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Africa",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Europe",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Oceania",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="SouthAmerica",]$Seasonalityindex,na.rm=T)
  ),
  max_WD = c(
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Asia",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="NorthAmerica",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Africa",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Europe",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Oceania",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="SouthAmerica",]$Seasonalityindex,na.rm=T)
  )
)
print(continent_limits)
plot_df_clip <- plot_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
plot_df_clip$Continent <- factor(plot_df_clip$Continent, 
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
plot_df_clip$Continent <- factor(plot_df_clip$Continent,
                                 levels = c("ASIA", "AFR", "NAM", "SAM", "EUR", "OCE"))
# Plot effect curve
effect_plot <- ggplot(plot_df_clip, aes(x = Seasonalityindex, y = Effect_on_FVC,
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
    legend.position = "none",
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
hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill=continent), binwidth = 0.0001, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.03)) +
  scale_y_continuous(limits = c(0, 1100), breaks = seq(0, 1100, by = 300)) +
  labs(x = "", y = "") +
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
effect_continent_Seasonalityindex_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_continent_Seasonalityindex_FVC)



###########################################################################
#                   city sizes: FVC ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                         .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))
merged_dt_Seasonalityindex$CitySize <- factor(merged_dt_Seasonalityindex$CitySize, levels = c("Small",
                                                                                              "Medium",
                                                                                              "Large"))
model_urban <- felm(FVC ~ 1 + 
                      Seasonalityindex + 
                      Seasonalityindex_2 +
                      urban +
                      Seasonalityindex:CitySize + 
                      Seasonalityindex_2:CitySize +
                      urban:CitySize +
                      Seasonalityindex:urban +
                      Seasonalityindex_2:urban +
                      Seasonalityindex:urban:CitySize +
                      Seasonalityindex_2:urban:CitySize +
                      WaterDeficit + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    |  CityID + year | 0 | CityID, 
                    data = merged_dt_Seasonalityindex)
summary(model_urban)
output_file <- "FVC_urban_in_out_citysize.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "Seasonalityindex",
  "Seasonalityindex_2",
  "urbanin",
  "Seasonalityindex:CitySizeMedium",
  "Seasonalityindex:CitySizeLarge",
  "Seasonalityindex_2:CitySizeMedium",
  "Seasonalityindex_2:CitySizeLarge",
  "urbanin:CitySizeMedium",
  "urbanin:CitySizeLarge",
  "Seasonalityindex:urbanin",
  "Seasonalityindex_2:urbanin",
  "Seasonalityindex:urbanin:CitySizeMedium",
  "Seasonalityindex:urbanin:CitySizeLarge",
  "Seasonalityindex_2:urbanin:CitySizeMedium",
  "Seasonalityindex_2:urbanin:CitySizeLarge"
)
b <- beta_all[v_names]
citysize_colors <- c(
  "Small" = "#95a5a6",
  "Medium" = "#2ecc71",
  "Large" = "#e74c3c"
)
wd_seq <- seq(0, 0.1, length.out = 200)
city_sizes <- c("Small","Medium","Large")
plot_df     <- data.frame()
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
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
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_FVC=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                CitySize=cs,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_FVC=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                CitySize=cs,
                                LandType="Urban"))
  }
}
citysize_limits <- tibble(
  CitySize = c("Small","Medium","Large"),
  min_WD = c(
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$CitySize=="Small",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$CitySize=="Medium",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$CitySize=="Large",]$Seasonalityindex,na.rm=T)
  ),
  max_WD = c(
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$CitySize=="Small",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$CitySize=="Medium",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$CitySize=="Large",]$Seasonalityindex,na.rm=T)
  )
)
print(citysize_limits)
plot_df_clip <- plot_df %>%
  left_join(citysize_limits, by = "CitySize") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
plot_df_clip$CitySize <- factor(plot_df_clip$CitySize,
                                levels = c("Small","Medium","Large"))
# Plot effect curve
effect_plot <- ggplot(plot_df_clip, aes(x = Seasonalityindex, y = Effect_on_FVC,
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=citysize_colors) +
  scale_fill_manual(values=citysize_colors)
print(effect_plot)
hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill=CitySize), binwidth = 0.0002, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.03)) +
  labs(x = "", y = "") +
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
effect_citysize_Seasonalityindex_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_citysize_Seasonalityindex_FVC)


###########################################################################
#                   development levels: FVC ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                         .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))
merged_dt_Seasonalityindex$DevelopmentLevel <- factor(merged_dt_Seasonalityindex$DevelopmentLevel, levels = c("developing", "developed"))
model_urban <- felm(FVC ~ 1 + 
                      Seasonalityindex + 
                      Seasonalityindex_2 +
                      urban +
                      Seasonalityindex:DevelopmentLevel + 
                      Seasonalityindex_2:DevelopmentLevel +
                      urban:DevelopmentLevel +
                      Seasonalityindex:urban +
                      Seasonalityindex_2:urban +
                      Seasonalityindex:urban:DevelopmentLevel +
                      Seasonalityindex_2:urban:DevelopmentLevel +
                      WaterDeficit + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    |  CityID + year | 0 | CityID, 
                    data = merged_dt_Seasonalityindex)
summary(model_urban)
output_file <- "FVC_urban_in_out_level.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "Seasonalityindex",
  "Seasonalityindex_2",
  "urbanin",
  "Seasonalityindex:DevelopmentLeveldeveloped",
  "Seasonalityindex_2:DevelopmentLeveldeveloped",
  "urbanin:DevelopmentLeveldeveloped",
  "Seasonalityindex:urbanin",
  "Seasonalityindex_2:urbanin",
  "Seasonalityindex:urbanin:DevelopmentLeveldeveloped",
  "Seasonalityindex_2:urbanin:DevelopmentLeveldeveloped"
)
b <- beta_all[v_names]
dev_colors <- c(
  "Global North" = "#0071bc", 
  "Global South" = "#d95218"
)
wd_seq <- seq(0, 0.1, length.out = 200)
dev_levels <- c("developing","developed")
plot_df     <- data.frame()
for(dl in dev_levels){
  for(wd in wd_seq){
    if(dl=="developing"){
      g_non <- c(
        wd, wd^2, 0,
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    if(dl=="developing"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,
        wd, wd^2,
        0,0
      )
    }else if(dl=="developed"){
      g_urb <- c(
        wd, wd^2, 1,
        wd, wd^2, 1,
        wd, wd^2,
        wd, wd^2
      )
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_FVC=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                DevelopmentLevel=dl,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_FVC=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                DevelopmentLevel=dl,
                                LandType="Urban"))
  }
}
dev_limits <- tibble(
  DevelopmentLevel = c("developing","developed"),
  min_WD = c(
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$DevelopmentLevel=="developing",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$DevelopmentLevel=="developed",]$Seasonalityindex,na.rm=T)
  ),
  max_WD = c(
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$DevelopmentLevel=="developing",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$DevelopmentLevel=="developed",]$Seasonalityindex,na.rm=T)
  )
)
print(dev_limits)
plot_df_clip <- plot_df %>%
  left_join(dev_limits, by = "DevelopmentLevel") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
plot_df_clip$DevelopmentLevel <- factor(plot_df_clip$DevelopmentLevel, 
                                        levels = c("developing","developed"),
                                        labels = c("Global South","Global North")) 
# Plot effect curve
effect_plot <- ggplot(plot_df_clip, aes(x = Seasonalityindex, y = Effect_on_FVC,
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=dev_colors) +
  scale_fill_manual(values=dev_colors)
print(effect_plot)
hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill = DevelopmentLevel), binwidth =  0.0001,  alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.03)) +
  scale_y_continuous(limits = c(0, 360), breaks = seq(0, 360, by =100)) +
  labs(x = "", y = "") +
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
effect_dev_Seasonalityindex_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_dev_Seasonalityindex_FVC)


###########################################################################
#                   nations/regions: FVC ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                         .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))
merged_dt_Seasonalityindex$country <- factor(merged_dt_Seasonalityindex$country, 
                                             levels = c("United States of America",
                                                        "EuropeNR",
                                                        "China",
                                                        "India",
                                                        "Brazil", "other"))
model_urban <- felm(FVC ~ 1 + 
                      Seasonalityindex + 
                      Seasonalityindex_2 +
                      urban +
                      Seasonalityindex:country + 
                      Seasonalityindex_2:country +
                      urban:country +
                      Seasonalityindex:urban +
                      Seasonalityindex_2:urban +
                      Seasonalityindex:urban:country +
                      Seasonalityindex_2:urban:country +
                      WaterDeficit + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    |  CityID + year | 0 | CityID, 
                    data = merged_dt_Seasonalityindex)
summary(model_urban)
output_file <- "FVC_urban_in_out_country.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "Seasonalityindex",
  "Seasonalityindex_2",
  "urbanin",
  "Seasonalityindex:countryEuropeNR",
  "Seasonalityindex:countryChina",
  "Seasonalityindex:countryIndia",
  "Seasonalityindex:countryBrazil",
  "Seasonalityindex:countryother",
  "Seasonalityindex_2:countryEuropeNR",
  "Seasonalityindex_2:countryChina",
  "Seasonalityindex_2:countryIndia",
  "Seasonalityindex_2:countryBrazil",
  "Seasonalityindex_2:countryother",
  "urbanin:countryEuropeNR",
  "urbanin:countryChina",
  "urbanin:countryIndia",
  "urbanin:countryBrazil",
  "urbanin:countryother",
  "Seasonalityindex:urbanin",
  "Seasonalityindex_2:urbanin",
  "Seasonalityindex:urbanin:countryEuropeNR",
  "Seasonalityindex:urbanin:countryChina",
  "Seasonalityindex:urbanin:countryIndia",
  "Seasonalityindex:urbanin:countryBrazil",
  "Seasonalityindex:urbanin:countryother",
  "Seasonalityindex_2:urbanin:countryEuropeNR",
  "Seasonalityindex_2:urbanin:countryChina",
  "Seasonalityindex_2:urbanin:countryIndia",
  "Seasonalityindex_2:urbanin:countryBrazil",
  "Seasonalityindex_2:urbanin:countryother"
)
b <- beta_all[v_names]
wd_seq <- seq(0, 0.1, length.out = 200)
cntry_levels <- c("United States","EuropeNR","China","India","Brazil","other")
plot_df     <- data.frame()
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
    }else if(cnt=="EuropeNR"){
      g_non <- c(
        wd, wd^2, 0,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
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
    }else if(cnt=="EuropeNR"){
      g_urb <- c(
        wd, wd^2, 1,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
        1,0,0,0,0,
        wd, wd^2,
        wd,0,0,0,0,
        wd^2,0,0,0,0
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
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_FVC=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                Country=cnt,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_FVC=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                Country=cnt,
                                LandType="Urban"))
  }
}
country_limits <- tibble(
  Country = c("United States", "EuropeNR", "China", "India", "Brazil", "other"),
  min_WD = c(
    min(merged_dt_Seasonalityindex[country == "United States of America"]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[country == "EuropeNR" ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[country == "China" ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[country == "India" ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[country == "Brazil" ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[country == "other" ]$Seasonalityindex, na.rm = TRUE)
  ),
  max_WD = c(
    max(merged_dt_Seasonalityindex[country == "United States of America" ]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[country == "EuropeNR" ]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[country == "China"]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[country == "India"]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[country == "Brazil"]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[country == "other"]$Seasonalityindex, na.rm = TRUE)
  )
)
print(country_limits)
plot_df_clip <- plot_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
plot_df_clip$Country <- factor(plot_df_clip$Country, 
                               levels = c("United States", "EuropeNR",  "China", "India", "Brazil" ,"other"),
                               labels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH")) 
country_colors <- c(
  "US" = "#2ecc71",
  "EURxR" = "#3498db",
  "CHN" = "#CB0505",
  "IND" = "#f39c12",
  "BRA" = "#4A1054",
  "OTH" = "#95a5a6"
)
plot_df_clip$Country <- factor(plot_df_clip$Country, levels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH"))
# Plot effect curve
effect_plot <- ggplot(plot_df_clip, aes(x = Seasonalityindex, y = Effect_on_FVC,
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=country_colors) +
  scale_fill_manual(values=country_colors)
print(effect_plot)
hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill=country), binwidth = 0.0002, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.03)) +
  scale_y_continuous(limits = c(0, 550), breaks = seq(0, 550, by =200)) +
  labs(x = "Rainfall seasonality", y = "") +
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
    "EuropeNR" = "#3498db",
    "China" = "#CB0505",
    "India" = "#f39c12",
    "Brazil" = "#4A1054",
    "other" = "#95a5a6"
  ))
print(hist_plot)
effect_country_Seasonalityindex_FVC <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_country_Seasonalityindex_FVC)


###########################################################################
#                   climate zones: NPP ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude',
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population',
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)),
                                                         .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))
merged_dt_Seasonalityindex$ClimateZone <- factor(merged_dt_Seasonalityindex$ClimateZone, levels = c("Temperate","Tropical","Cold","Arid"))
model_urban <- felm(NPP ~ 1 +
                      Seasonalityindex +
                      Seasonalityindex_2 +
                      urban +
                      Seasonalityindex:ClimateZone +
                      Seasonalityindex_2:ClimateZone +
                      urban:ClimateZone +
                      Seasonalityindex:urban +
                      Seasonalityindex_2:urban +
                      Seasonalityindex:urban:ClimateZone +
                      Seasonalityindex_2:urban:ClimateZone +
                      WaterDeficit + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    |  CityID + year | 0 | CityID,
                    data = merged_dt_Seasonalityindex)
summary(model_urban)
output_file <- "NPP_urban_in_out_climate.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "Seasonalityindex",
  "Seasonalityindex_2",
  "urbanin",
  "Seasonalityindex:ClimateZoneTropical",
  "Seasonalityindex:ClimateZoneCold",
  "Seasonalityindex:ClimateZoneArid",
  "Seasonalityindex_2:ClimateZoneTropical",
  "Seasonalityindex_2:ClimateZoneCold",
  "Seasonalityindex_2:ClimateZoneArid",
  "urbanin:ClimateZoneTropical",
  "urbanin:ClimateZoneCold",
  "urbanin:ClimateZoneArid",
  "Seasonalityindex:urbanin",
  "Seasonalityindex_2:urbanin",
  "Seasonalityindex:urbanin:ClimateZoneTropical",
  "Seasonalityindex:urbanin:ClimateZoneCold",
  "Seasonalityindex:urbanin:ClimateZoneArid",
  "Seasonalityindex_2:urbanin:ClimateZoneTropical",
  "Seasonalityindex_2:urbanin:ClimateZoneCold",
  "Seasonalityindex_2:urbanin:ClimateZoneArid"
)
b <- beta_all[v_names]
climate_colors <- c(
  "Temperate" = "#2ecc71",
  "Tropical"  = "#e74c3c",
  "Arid"      = "#95a5a6",
  "Cold"      = "#3498db"
)
wd_seq <- seq(0, 0.1, length.out = 200)
climates <- c("Temperate","Tropical","Cold","Arid")
plot_df <- data.frame()
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
    }else if(cz=="Tropical"){
      g_non <- c(
        wd, wd^2, 0,
        wd,0,0,
        wd^2,0,0,
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
    }else if(cz=="Arid"){
      g_non <- c(
        wd, wd^2, 0,
        0,0,wd,
        0,0,wd^2,
        0,0,0,
        0,0,0,0,0,
        0,0,0
      )
    }
    
    pred_non <- as.numeric(crossprod(g_non, b))
    se_non   <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
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
    }
    
    pred_urb <- as.numeric(crossprod(g_urb, b))
    se_urb   <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_NPP=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                ClimateZone=cz,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_NPP=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                ClimateZone=cz,
                                LandType="Urban"))
  }
}
climate_limits <- tibble(
  ClimateZone = c("Temperate", "Tropical", "Arid",  "Cold"),
  min_WD = c(
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Temperate", ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Tropical", ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Arid", ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Cold", ]$Seasonalityindex, na.rm = TRUE)
  ),
  max_WD = c(
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Temperate", ]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Tropical", ]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Arid", ]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$ClimateZone == "Cold", ]$Seasonalityindex, na.rm = TRUE)
  )
)
print(climate_limits)
plot_df_clip <- plot_df %>%
  left_join(climate_limits, by = "ClimateZone") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
plot_df_clip$ClimateZone <- factor(plot_df_clip$ClimateZone, levels = c( "Tropical", "Temperate",  "Cold","Arid"))

effect_plot <- ggplot(plot_df_clip, aes(x = Seasonalityindex, y = Effect_on_NPP,
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
    y = "Effect on NPP",
    color = "",
    fill = ""
  ) +
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
    legend.position = "none",
    legend.direction = "vertical",
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=climate_colors) +
  scale_fill_manual(values=climate_colors)
print(effect_plot)

hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill=ClimateZone), binwidth = 0.0001, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.03)) +
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by =200)) +
  labs(x = "", y = "") +
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
effect_climate_Seasonalityindex_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_climate_Seasonalityindex_NPP)


###########################################################################
#                   continents: NPP ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude',
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population',
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)),
                                                         .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))
merged_dt_Seasonalityindex$continent <- factor(merged_dt_Seasonalityindex$continent, levels = c( "Asia",
                                                                                                 "NorthAmerica",
                                                                                                 "Africa",
                                                                                                 "Europe",
                                                                                                 "Oceania",
                                                                                                 "SouthAmerica"))
model_urban <- felm(NPP ~ 1 +
                      Seasonalityindex + Seasonalityindex_2 +
                      urban +
                      Seasonalityindex:continent +
                      Seasonalityindex_2:continent +
                      urban:continent +
                      Seasonalityindex:urban +
                      Seasonalityindex_2:urban +
                      Seasonalityindex:urban:continent +
                      Seasonalityindex_2:urban:continent +
                      WaterDeficit + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    |  CityID + year | 0 | CityID,
                    data = merged_dt_Seasonalityindex)
summary(model_urban)
output_file <- "NPP_urban_in_out_continent.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "Seasonalityindex",
  "Seasonalityindex_2",
  "urbanin",
  "Seasonalityindex:continentNorthAmerica",
  "Seasonalityindex:continentAfrica",
  "Seasonalityindex:continentEurope",
  "Seasonalityindex:continentOceania",
  "Seasonalityindex:continentSouthAmerica",
  "Seasonalityindex_2:continentNorthAmerica",
  "Seasonalityindex_2:continentAfrica",
  "Seasonalityindex_2:continentEurope",
  "Seasonalityindex_2:continentOceania",
  "Seasonalityindex_2:continentSouthAmerica",
  "urbanin:continentNorthAmerica",
  "urbanin:continentAfrica",
  "urbanin:continentEurope",
  "urbanin:continentOceania",
  "urbanin:continentSouthAmerica",
  "Seasonalityindex:urbanin",
  "Seasonalityindex_2:urbanin",
  "Seasonalityindex:urbanin:continentNorthAmerica",
  "Seasonalityindex:urbanin:continentAfrica",
  "Seasonalityindex:urbanin:continentEurope",
  "Seasonalityindex:urbanin:continentOceania",
  "Seasonalityindex:urbanin:continentSouthAmerica",
  "Seasonalityindex_2:urbanin:continentNorthAmerica",
  "Seasonalityindex_2:urbanin:continentAfrica",
  "Seasonalityindex_2:urbanin:continentEurope",
  "Seasonalityindex_2:urbanin:continentOceania",
  "Seasonalityindex_2:urbanin:continentSouthAmerica"
)
b <- beta_all[v_names]
wd_seq <- seq(0, 0.1, length.out = 200)
continents <- c("Asia","NorthAmerica","Africa","Europe","Oceania","SouthAmerica")
plot_df     <- data.frame()
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
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
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_NPP=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                Continent=ct,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_NPP=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                Continent=ct,
                                LandType="Urban"))
  }
}
continent_limits <- tibble(
  Continent = c("Asia","NorthAmerica","Africa","Europe","Oceania","SouthAmerica"),
  min_WD = c(
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Asia",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="NorthAmerica",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Africa",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Europe",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Oceania",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="SouthAmerica",]$Seasonalityindex,na.rm=T)
  ),
  max_WD = c(
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Asia",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="NorthAmerica",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Africa",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Europe",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="Oceania",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$continent=="SouthAmerica",]$Seasonalityindex,na.rm=T)
  )
)
print(continent_limits)
plot_df_clip <- plot_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
plot_df_clip$Continent <- factor(plot_df_clip$Continent,
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
plot_df_clip$Continent <- factor(plot_df_clip$Continent,
                                 levels = c("ASIA", "AFR", "NAM", "SAM", "EUR", "OCE"))

effect_plot <- ggplot(plot_df_clip, aes(x = Seasonalityindex, y = Effect_on_NPP,
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
    y = "Effect on NPP",
    color = "",
    fill = ""
  ) +
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
    legend.position = "none",
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

hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill=continent), binwidth = 0.0001, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.03)) +
  scale_y_continuous(limits = c(0, 1100), breaks = seq(0, 1100, by = 300)) +
  labs(x = "", y = "") +
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
effect_continent_Seasonalityindex_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_continent_Seasonalityindex_NPP)


###########################################################################
#                   city sizes: NPP ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                         .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))
merged_dt_Seasonalityindex$CitySize <- factor(merged_dt_Seasonalityindex$CitySize, levels = c("Small",
                                                                                              "Medium",
                                                                                              "Large"))
model_urban <- felm(NPP ~ 1 + 
                      Seasonalityindex + 
                      Seasonalityindex_2 +
                      urban +
                      Seasonalityindex:CitySize + 
                      Seasonalityindex_2:CitySize +
                      urban:CitySize +
                      Seasonalityindex:urban +
                      Seasonalityindex_2:urban +
                      Seasonalityindex:urban:CitySize +
                      Seasonalityindex_2:urban:CitySize +
                      WaterDeficit + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    |  CityID + year | 0 | CityID, 
                    data = merged_dt_Seasonalityindex)
summary(model_urban)
output_file <- "NPP_urban_in_out_citysize.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "Seasonalityindex",
  "Seasonalityindex_2",
  "urbanin",
  "Seasonalityindex:CitySizeMedium",
  "Seasonalityindex:CitySizeLarge",
  "Seasonalityindex_2:CitySizeMedium",
  "Seasonalityindex_2:CitySizeLarge",
  "urbanin:CitySizeMedium",
  "urbanin:CitySizeLarge",
  "Seasonalityindex:urbanin",
  "Seasonalityindex_2:urbanin",
  "Seasonalityindex:urbanin:CitySizeMedium",
  "Seasonalityindex:urbanin:CitySizeLarge",
  "Seasonalityindex_2:urbanin:CitySizeMedium",
  "Seasonalityindex_2:urbanin:CitySizeLarge"
)
b <- beta_all[v_names]
citysize_colors <- c(
  "Small" = "#95a5a6",
  "Medium" = "#2ecc71",
  "Large" = "#e74c3c"
)
wd_seq <- seq(0, 0.1, length.out = 200)
city_sizes <- c("Small","Medium","Large")
plot_df     <- data.frame()

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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
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
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_NPP=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                CitySize=cs,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_NPP=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                CitySize=cs,
                                LandType="Urban"))
  }
}

citysize_limits <- tibble(
  CitySize = c("Small","Medium","Large"),
  min_WD = c(
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$CitySize=="Small",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$CitySize=="Medium",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$CitySize=="Large",]$Seasonalityindex,na.rm=T)
  ),
  max_WD = c(
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$CitySize=="Small",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$CitySize=="Medium",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$CitySize=="Large",]$Seasonalityindex,na.rm=T)
  )
)
print(citysize_limits)
plot_df_clip <- plot_df %>%
  left_join(citysize_limits, by = "CitySize") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
plot_df_clip$CitySize <- factor(plot_df_clip$CitySize,
                                levels = c("Small","Medium","Large"))
# Plot effect curve
effect_plot <- ggplot(plot_df_clip, aes(x = Seasonalityindex, y = Effect_on_NPP,
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
    y = "Effect on NPP",
    color = "",
    fill = ""
  ) +
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=citysize_colors) +
  scale_fill_manual(values=citysize_colors)
print(effect_plot)
hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill=CitySize), binwidth = 0.0002, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.03)) +
  labs(x = "", y = "") +
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
effect_citysize_Seasonalityindex_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_citysize_Seasonalityindex_NPP)

###########################################################################
#                   development levels: NPP ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                         .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))
merged_dt_Seasonalityindex$DevelopmentLevel <- factor(merged_dt_Seasonalityindex$DevelopmentLevel, levels = c("developing", "developed"))
model_urban <- felm(NPP ~ 1 + 
                      Seasonalityindex + 
                      Seasonalityindex_2 +
                      urban +
                      Seasonalityindex:DevelopmentLevel + 
                      Seasonalityindex_2:DevelopmentLevel +
                      urban:DevelopmentLevel +
                      Seasonalityindex:urban +
                      Seasonalityindex_2:urban +
                      Seasonalityindex:urban:DevelopmentLevel +
                      Seasonalityindex_2:urban:DevelopmentLevel +
                      WaterDeficit + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    |  CityID + year | 0 | CityID, 
                    data = merged_dt_Seasonalityindex)
summary(model_urban)
output_file <- "NPP_urban_in_out_level.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "Seasonalityindex",
  "Seasonalityindex_2",
  "urbanin",
  "Seasonalityindex:DevelopmentLeveldeveloped",
  "Seasonalityindex_2:DevelopmentLeveldeveloped",
  "urbanin:DevelopmentLeveldeveloped",
  "Seasonalityindex:urbanin",
  "Seasonalityindex_2:urbanin",
  "Seasonalityindex:urbanin:DevelopmentLeveldeveloped",
  "Seasonalityindex_2:urbanin:DevelopmentLeveldeveloped"
)
b <- beta_all[v_names]
dev_colors <- c(
  "Global North" = "#0071bc", 
  "Global South" = "#d95218"
)
wd_seq <- seq(0, 0.1, length.out = 200)
dev_levels <- c("developing","developed")
plot_df     <- data.frame()

for(dl in dev_levels){
  for(wd in wd_seq){
    if(dl=="developing"){
      g_non <- c(
        wd, wd^2, 0,
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
    if(dl=="developing"){
      g_urb <- c(
        wd, wd^2, 1,
        0,0,0,
        wd, wd^2,
        0,0
      )
    }else if(dl=="developed"){
      g_urb <- c(
        wd, wd^2, 1,
        wd, wd^2, 1,
        wd, wd^2,
        wd, wd^2
      )
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_NPP=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                DevelopmentLevel=dl,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_NPP=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                DevelopmentLevel=dl,
                                LandType="Urban"))
  }
}

dev_limits <- tibble(
  DevelopmentLevel = c("developing","developed"),
  min_WD = c(
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$DevelopmentLevel=="developing",]$Seasonalityindex,na.rm=T),
    min(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$DevelopmentLevel=="developed",]$Seasonalityindex,na.rm=T)
  ),
  max_WD = c(
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$DevelopmentLevel=="developing",]$Seasonalityindex,na.rm=T),
    max(merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$DevelopmentLevel=="developed",]$Seasonalityindex,na.rm=T)
  )
)
print(dev_limits)
plot_df_clip <- plot_df %>%
  left_join(dev_limits, by = "DevelopmentLevel") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
plot_df_clip$DevelopmentLevel <- factor(plot_df_clip$DevelopmentLevel, 
                                        levels = c("developing","developed"),
                                        labels = c("Global South","Global North")) 
# Plot effect curve
effect_plot <- ggplot(plot_df_clip, aes(x = Seasonalityindex, y = Effect_on_NPP,
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
    y = "Effect on NPP",
    color = "",
    fill = ""
  ) +
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=dev_colors) +
  scale_fill_manual(values=dev_colors)
print(effect_plot)
hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill = DevelopmentLevel), binwidth =  0.0001,  alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.03)) +
  scale_y_continuous(limits = c(0, 360), breaks = seq(0, 360, by =100)) +
  labs(x = "", y = "") +
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
effect_dev_Seasonalityindex_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_dev_Seasonalityindex_NPP)

###########################################################################
#                   nations/regions: NPP ~ Rainfall seasonality
###########################################################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                         .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))
merged_dt_Seasonalityindex$country <- factor(merged_dt_Seasonalityindex$country, 
                                             levels = c("United States of America",
                                                        "EuropeNR",
                                                        "China",
                                                        "India",
                                                        "Brazil", "other"))
model_urban <- felm(NPP ~ 1 + 
                      Seasonalityindex + 
                      Seasonalityindex_2 +
                      urban +
                      Seasonalityindex:country + 
                      Seasonalityindex_2:country +
                      urban:country +
                      Seasonalityindex:urban +
                      Seasonalityindex_2:urban +
                      Seasonalityindex:urban:country +
                      Seasonalityindex_2:urban:country +
                      WaterDeficit + AverageTemperature +
                      TemperatureRange + WindSpeed + Elevation +
                      Latitude + Precipitation + SoilMoisture +
                      SoilPH + GDP_per_capita_PPP + HDI +
                      Population + ImperviousSurface +
                      HumanSettlement + CityArea
                    |  CityID + year | 0 | CityID, 
                    data = merged_dt_Seasonalityindex)
summary(model_urban)
output_file <- "NPP_urban_in_out_country.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)
beta_all <- coef(model_urban)
vcov_mat <- vcov(model_urban)
v_names <- c(
  "Seasonalityindex",
  "Seasonalityindex_2",
  "urbanin",
  "Seasonalityindex:countryEuropeNR",
  "Seasonalityindex:countryChina",
  "Seasonalityindex:countryIndia",
  "Seasonalityindex:countryBrazil",
  "Seasonalityindex:countryother",
  "Seasonalityindex_2:countryEuropeNR",
  "Seasonalityindex_2:countryChina",
  "Seasonalityindex_2:countryIndia",
  "Seasonalityindex_2:countryBrazil",
  "Seasonalityindex_2:countryother",
  "urbanin:countryEuropeNR",
  "urbanin:countryChina",
  "urbanin:countryIndia",
  "urbanin:countryBrazil",
  "urbanin:countryother",
  "Seasonalityindex:urbanin",
  "Seasonalityindex_2:urbanin",
  "Seasonalityindex:urbanin:countryEuropeNR",
  "Seasonalityindex:urbanin:countryChina",
  "Seasonalityindex:urbanin:countryIndia",
  "Seasonalityindex:urbanin:countryBrazil",
  "Seasonalityindex:urbanin:countryother",
  "Seasonalityindex_2:urbanin:countryEuropeNR",
  "Seasonalityindex_2:urbanin:countryChina",
  "Seasonalityindex_2:urbanin:countryIndia",
  "Seasonalityindex_2:urbanin:countryBrazil",
  "Seasonalityindex_2:urbanin:countryother"
)
b <- beta_all[v_names]
wd_seq <- seq(0, 0.1, length.out = 200)
cntry_levels <- c("United States","EuropeNR","China","India","Brazil","other")
plot_df     <- data.frame()

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
    }else if(cnt=="EuropeNR"){
      g_non <- c(
        wd, wd^2, 0,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
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
    }
    
    pred_non  <- as.numeric(crossprod(g_non, b))
    se_non    <- sqrt(as.numeric(t(g_non) %*% vcov_mat[v_names, v_names] %*% g_non))
    ci_low_non  <- pred_non - 1.96*se_non
    ci_high_non <- pred_non + 1.96*se_non
    
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
    }else if(cnt=="EuropeNR"){
      g_urb <- c(
        wd, wd^2, 1,
        wd,0,0,0,0,
        wd^2,0,0,0,0,
        1,0,0,0,0,
        wd, wd^2,
        wd,0,0,0,0,
        wd^2,0,0,0,0
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
    }
    
    pred_urb  <- as.numeric(crossprod(g_urb, b))
    se_urb    <- sqrt(as.numeric(t(g_urb) %*% vcov_mat[v_names, v_names] %*% g_urb))
    ci_low_urb  <- pred_urb - 1.96*se_urb
    ci_high_urb <- pred_urb + 1.96*se_urb
    
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_NPP=pred_non,
                                ci_low=ci_low_non,
                                ci_high=ci_high_non,
                                Country=cnt,
                                LandType="Non-urban"))
    plot_df <- rbind(plot_df,
                     data.frame(Seasonalityindex=wd,
                                Effect_on_NPP=pred_urb,
                                ci_low=ci_low_urb,
                                ci_high=ci_high_urb,
                                Country=cnt,
                                LandType="Urban"))
  }
}

country_limits <- tibble(
  Country = c("United States", "EuropeNR", "China", "India", "Brazil", "other"),
  min_WD = c(
    min(merged_dt_Seasonalityindex[country == "United States of America"]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[country == "EuropeNR" ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[country == "China" ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[country == "India" ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[country == "Brazil" ]$Seasonalityindex, na.rm = TRUE),
    min(merged_dt_Seasonalityindex[country == "other" ]$Seasonalityindex, na.rm = TRUE)
  ),
  max_WD = c(
    max(merged_dt_Seasonalityindex[country == "United States of America" ]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[country == "EuropeNR" ]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[country == "China"]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[country == "India"]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[country == "Brazil"]$Seasonalityindex, na.rm = TRUE),
    max(merged_dt_Seasonalityindex[country == "other"]$Seasonalityindex, na.rm = TRUE)
  )
)
print(country_limits)
plot_df_clip <- plot_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
plot_df_clip$Country <- factor(plot_df_clip$Country, 
                               levels = c("United States", "EuropeNR",  "China", "India", "Brazil" ,"other"),
                               labels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH")) 
country_colors <- c(
  "US" = "#2ecc71",
  "EURxR" = "#3498db",
  "CHN" = "#CB0505",
  "IND" = "#f39c12",
  "BRA" = "#4A1054",
  "OTH" = "#95a5a6"
)
plot_df_clip$Country <- factor(plot_df_clip$Country, levels = c("US","EURxR","CHN", "IND", "BRA" ,"OTH"))
# Plot effect curve
effect_plot <- ggplot(plot_df_clip, aes(x = Seasonalityindex, y = Effect_on_NPP,
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
    y = "Effect on NPP",
    color = "",
    fill = ""
  ) +
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
    legend.position = "none",
    legend.direction = "vertical", 
    legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "pt")),  
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  scale_color_manual(values=country_colors) +
  scale_fill_manual(values=country_colors)
print(effect_plot)
hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill=country), binwidth = 0.0002, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.03)) +
  scale_y_continuous(limits = c(0, 550), breaks = seq(0, 550, by =200)) +
  labs(x = "Rainfall seasonality", y = "") +
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
    "EuropeNR" = "#3498db",
    "China" = "#CB0505",
    "India" = "#f39c12",
    "Brazil" = "#4A1054",
    "other" = "#95a5a6"
  ))
print(hist_plot)
effect_country_Seasonalityindex_NPP <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_country_Seasonalityindex_NPP)



crop_bottom <- function(plot, crop_ratio = 0.15) {
  ggdraw() +
    draw_plot(plot, x = 0, y = crop_ratio, width = 1, height = 1 - crop_ratio) +
    theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"))
}

crop_bottom_cm <- function(plot, crop_cm = 0.3) {
  ggdraw() +
    draw_plot(plot, x = 0, y = crop_cm/2, width = 1, height = 1 - crop_cm/2) +
    theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"))
}

col1_1 <- crop_bottom(effect_continent_waterdeficit_FVC, crop_ratio = -0.2)  # continent
col1_2 <- crop_bottom(effect_climate_waterdeficit_FVC,   crop_ratio = -0.2)  # climate
col1_3 <- crop_bottom(effect_citysize_waterdeficit_FVC,  crop_ratio = -0.2)  # city size
col1_4 <- crop_bottom(effect_dev_waterdeficit_FVC,       crop_ratio = -0.2)  # development level
col1_5 <- crop_bottom(effect_country_waterdeficit_FVC,   crop_ratio = -0.2)  # country/region

col2_1 <- crop_bottom(effect_continent_waterdeficit_NPP, crop_ratio = -0.2)
col2_2 <- crop_bottom(effect_climate_waterdeficit_NPP,   crop_ratio = -0.2)
col2_3 <- crop_bottom(effect_citysize_waterdeficit_NPP,  crop_ratio = -0.2)
col2_4 <- crop_bottom(effect_dev_waterdeficit_NPP,       crop_ratio = -0.2)
col2_5 <- crop_bottom(effect_country_waterdeficit_NPP,   crop_ratio = -0.2)

col3_1 <- crop_bottom(effect_continent_Seasonalityindex_FVC, crop_ratio = -0.2)
col3_2 <- crop_bottom(effect_climate_Seasonalityindex_FVC,   crop_ratio = -0.2)
col3_3 <- crop_bottom(effect_citysize_Seasonalityindex_FVC,  crop_ratio = -0.2)
col3_4 <- crop_bottom(effect_dev_Seasonalityindex_FVC,       crop_ratio = -0.2)
col3_5 <- crop_bottom(effect_country_Seasonalityindex_FVC,   crop_ratio = -0.2)

col4_1 <- crop_bottom(effect_continent_Seasonalityindex_NPP, crop_ratio = -0.2)
col4_2 <- crop_bottom(effect_climate_Seasonalityindex_NPP,   crop_ratio = -0.2)
col4_3 <- crop_bottom(effect_citysize_Seasonalityindex_NPP,  crop_ratio = -0.2)
col4_4 <- crop_bottom(effect_dev_Seasonalityindex_NPP,       crop_ratio = -0.2)
col4_5 <- crop_bottom(effect_country_Seasonalityindex_NPP,   crop_ratio = -0.2)

fig_total <- plot_grid(
  col1_1, col2_1, col3_1, col4_1,
  col1_2, col2_2, col3_2, col4_2,
  col1_3, col2_3, col3_3, col4_3,
  col1_4, col2_4, col3_4, col4_4,
  col1_5, col2_5, col3_5, col4_5,
  nrow = 5,
  ncol = 4,
  align = "hv",
  axis = "tblr",
  rel_widths  = c(1, 1, 1, 1),
  rel_heights = c(1, 1, 1, 1, 1),
  labels = letters[1:20],
  label_size = 8,
  label_fontfamily = "Times",
  label_fontface = "bold",
  label_x = rep(0.11, 20),
  label_y = c(rep(1.03, 10), rep(1.06, 10)),
  hjust = 0,
  vjust = 1
)
fig_total <- cowplot::ggdraw(fig_total) +
  theme(plot.margin = margin(3, 3, 3, 3, unit = "mm"))

print(fig_total)

ggsave("./Fig3_5x4.pdf",
       plot = fig_total, width = 160, height = 170, units = "mm", dpi = 300)


