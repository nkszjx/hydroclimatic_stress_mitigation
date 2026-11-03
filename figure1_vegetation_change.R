

# Vegetation change

rm(list=ls())  
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(scales)
library(sf)
library(rcolors)
library(magrittr)
library(grid)
library(scattermore)



###########################################################################
#       Continental-averaged temporal variations of urban FVC
###########################################################################

setwd("XXX") # please set your folder path
merged_dt <- fread("Dataset_matched_buffer10km.csv")


urban_FVC_in <- merged_dt[urban == "in"]
urban_FVC_in_mean <- urban_FVC_in[, .(
  N = .N,
  FVC_mean = mean(FVC, na.rm = TRUE)
), by = .(year, continent)]

urban_FVC_in_mean[continent == "NorthAmerica", continent := "North America"]
urban_FVC_in_mean[continent == "SouthAmerica", continent := "South America"]

continent_colors <- c(
  "Asia" = "#3498db",
  "Africa" = "#95a5a6",
  "North America" = "#2ecc71",
  "South America" = "#9b59b6",
  "Europe" = "#e74c3c",
  "Oceania" = "#f39c12"
)

order_conti <- names(continent_colors)
conti_label <- c("Asia", "Africa", "North\nAmerica", "South\nAmerica", "Europe", "Oceania")


urban_FVC_in_mean$continent <- factor(urban_FVC_in_mean$continent, levels = order_conti)

p_a <- ggplot(urban_FVC_in_mean, aes(x = year, y = FVC_mean, color = continent, group = continent)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, stroke = 0, shape = 19) +
  scale_x_continuous(
    name = "Year",
    breaks = seq(2000, 2020, 5),
    expand = c(0.05, 0.05)
  ) +
  scale_y_continuous(
    name = expression(FVC[Urban]),
    limits = c(0.2, 0.53),
    breaks = seq(0.2, 0.53, 0.1),
    expand = c(0.02, 0.02)
  ) +
  scale_color_manual(
    values = continent_colors,
    breaks = order_conti,
    labels = conti_label,
    name = "Continent"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
    # axis.text.x = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black", margin = margin(t = 5)),
    # axis.title.y = element_text(size = 12, family = "serif", color = "black", margin = margin(r = 5)),
    # axis.ticks.length = unit(0.05, "cm"),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = c(0.50, 0.85),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = "black", lineheight = 0.7),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 3, ncol = 2, byrow = TRUE))

print(p_a)



###########################################################################
#       Continental-averaged temporal variations of non-urban FVC
###########################################################################

nonurban_FVC_out <- merged_dt[urban == "out"]
nonurban_FVC_out_mean <- nonurban_FVC_out[, .(
  N = .N,
  FVC_mean = mean(FVC, na.rm = TRUE)
), by = .(year, continent)]

nonurban_FVC_out_mean[continent == "NorthAmerica", continent := "North America"]
nonurban_FVC_out_mean[continent == "SouthAmerica", continent := "South America"]

continent_colors <- c(
  "Asia" = "#3498db",
  "Africa" = "#95a5a6",
  "North America" = "#2ecc71",
  "South America" = "#9b59b6",
  "Europe" = "#e74c3c",
  "Oceania" = "#f39c12"
)

order_conti <- names(continent_colors)
conti_label <- c("Asia", "Africa", "North\nAmerica", "South\nAmerica", "Europe", "Oceania")
nonurban_FVC_out_mean$continent <- factor(nonurban_FVC_out_mean$continent, levels = order_conti)

p_b <- ggplot(nonurban_FVC_out_mean, aes(x = year, y = FVC_mean, color = continent, group = continent)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, stroke = 0, shape = 19) +
  scale_x_continuous(
    name = "Year",
    breaks = seq(2000, 2020, 5),
    expand = c(0.05, 0.05)
  ) +
  scale_y_continuous(
    name = expression(FVC[Non-urban]),
    limits = c(0.2, 0.55),
    breaks = seq(0.2, 0.55, 0.1),
    expand = c(0.02, 0.02)
  ) +
  scale_color_manual(
    values = continent_colors,
    breaks = order_conti,
    labels = conti_label,
    name = "Continent"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
    # axis.text.x = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black", margin = margin(t = 5)),
    # axis.title.y = element_text(size = 12, family = "serif", color = "black", margin = margin(r = 5)),
    # axis.ticks.length = unit(0.05, "cm"),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = "black", lineheight = 0.7),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 3, ncol = 2, byrow = TRUE))

print(p_b)




###########################################################################
#       Continental-averaged temporal variations of urban NPP
###########################################################################

urban_NPP_in <- merged_dt[urban == "in"]

urban_NPP_in_mean <- urban_NPP_in[, .(
  N = .N,
  NPP_mean = mean(NPP, na.rm = TRUE)
), by = .(year, continent)]

urban_NPP_in_mean[continent == "NorthAmerica", continent := "North America"]
urban_NPP_in_mean[continent == "SouthAmerica", continent := "South America"]

continent_colors <- c(
  "Asia" = "#3498db",
  "Africa" = "#95a5a6",
  "North America" = "#2ecc71",
  "South America" = "#9b59b6",
  "Europe" = "#e74c3c",
  "Oceania" = "#f39c12"
)

order_conti <- names(continent_colors)

conti_label <- c("Asia", "Africa", "North\nAmerica", "South\nAmerica", "Europe", "Oceania")
urban_NPP_in_mean$continent <- factor(urban_NPP_in_mean$continent, levels = order_conti)

p_c <- ggplot(urban_NPP_in_mean, aes(x = year, y = NPP_mean, color = continent, group = continent)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, stroke = 0, shape = 19) +
  scale_x_continuous(
    name = "Year",
    breaks = seq(2000, 2020, 5),
    expand = c(0.05, 0.05)
  ) +
  scale_y_continuous(
    name = expression(NPP[Urban]~(g~C~m^-2~yr^-1)),
    limits = c(100, 600),
    breaks = seq(100, 600, 200),
    expand = c(0.02, 0.02)
  ) +
  scale_color_manual(
    values = continent_colors,
    breaks = order_conti,
    labels = conti_label,
    name = "Continent"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
    # axis.text.x = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black", margin = margin(t = 5)),
    # axis.title.y = element_text(size = 12, family = "serif", color = "black", margin = margin(r = 5)),
    # axis.ticks.length = unit(0.05, "cm"),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = "black", lineheight = 0.7),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 3, ncol = 2, byrow = TRUE))

print(p_c)




###########################################################################
#       Continental-averaged temporal variations of non-urban NPP
###########################################################################

nonurban_NPP_out <- merged_dt[urban == "out"]

nonurban_NPP_out_mean <- nonurban_NPP_out[, .(
  N = .N,
  NPP_mean = mean(NPP, na.rm = TRUE)
), by = .(year, continent)]

nonurban_NPP_out_mean[continent == "NorthAmerica", continent := "North America"]
nonurban_NPP_out_mean[continent == "SouthAmerica", continent := "South America"]

continent_colors <- c(
  "Asia" = "#3498db",
  "Africa" = "#95a5a6",
  "North America" = "#2ecc71",
  "South America" = "#9b59b6",
  "Europe" = "#e74c3c",
  "Oceania" = "#f39c12"
)

order_conti <- names(continent_colors)
conti_label <- c("Asia", "Africa", "North\nAmerica", "South\nAmerica", "Europe", "Oceania")

nonurban_NPP_out_mean$continent <- factor(nonurban_NPP_out_mean$continent, levels = order_conti)

p_d <- ggplot(nonurban_NPP_out_mean, aes(x = year, y = NPP_mean, color = continent, group = continent)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, stroke = 0, shape = 19) +
  scale_x_continuous(
    name = "Year",
    breaks = seq(2000, 2020, 5),
    expand = c(0.05, 0.26)
  ) +
  scale_y_continuous(
    name = expression(NPP[Non-urban]~(g~C~m^-2~yr^-1)),
    limits = c(200, 900),
    breaks = seq(200, 900, 200),
    expand = c(0.02, 0.02)
  ) +
  scale_color_manual(
    values = continent_colors,
    breaks = order_conti,
    labels = conti_label,
    name = "Continent"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 12, color = "black"),
    # axis.text.y = element_text(size = 11, family = "serif", color = "black"),
    axis.title = element_text(size = 13, color = "black", margin = margin(t = 5)),
    # axis.title.y = element_text(size = 12, family = "serif", color = "black", margin = margin(r = 5)),
    # axis.ticks.length = unit(0.05, "cm"),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = "black", lineheight = 0.7),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 3, ncol = 2, byrow = TRUE))

print(p_d)




###########################################################################
#                Temporal dynamics of ∆FVC over 2000–2020
###########################################################################

out_dt <- merged_dt[
  urban == "out" & year %in% c(2000,2005,2010,2015,2020),
  .(CityID, FVC, NPP, year)
]
city_year_out <- out_dt[,
                        .(
                          FVC_out = mean(FVC, na.rm = TRUE),
                          NPP_out = mean(NPP, na.rm = TRUE)
                        ),
                        by = .(CityID, year)
]


in_dt <- merged_dt[
  urban == "in" & year %in% c(2000,2005,2010,2015,2020),
  .(CityID, FVC, NPP, year, Latitude, Longitude)
]
city_year_in <- in_dt[,
                      .(
                        FVC_in = mean(FVC, na.rm = TRUE),
                        NPP_in = mean(NPP, na.rm = TRUE),
                        Latitude = mean(Latitude, na.rm = TRUE),
                        Longitude = mean(Longitude, na.rm = TRUE)
                      ),
                      by = .(CityID, year)
]


city_final_dt_year <- merge(
  city_year_in,
  city_year_out,
  by = c("CityID","year"),
  all = FALSE
)
city_final_dt_year[, `:=`(
  delta_FVC = FVC_in - FVC_out,
  delta_NPP = NPP_in - NPP_out
)]


year_summary <- city_final_dt_year[,
                                   .(
                                     dFVC_mean = mean(delta_FVC, na.rm = TRUE),
                                     dFVC_sd05 = 0.5 * sd(delta_FVC, na.rm = TRUE),
                                     dNPP_mean = mean(delta_NPP, na.rm = TRUE),
                                     dNPP_sd05 = 0.5 * sd(delta_NPP, na.rm = TRUE),
                                     n_city = .N
                                   ),
                                   by = .(year)
]


lm_fvc <- lm(delta_FVC ~ year, data = city_final_dt_year)
slope_fvc <- coef(lm_fvc)[["year"]]
p_fvc <- summary(lm_fvc)$coefficients["year","Pr(>|t|)"]
lab_fvc <- sprintf("Slope=%.3e", slope_fvc)
ptext_fvc <- ifelse(p_fvc<0.001, "p<0.001", sprintf("p=%.3f",p_fvc))

lm_npp <- lm(delta_NPP ~ year, data = city_final_dt_year)
slope_npp <- coef(lm_npp)[["year"]]
p_npp <- summary(lm_npp)$coefficients["year","Pr(>|t|)"]
lab_npp <- sprintf("Slope=%.3e", slope_npp)
ptext_npp <- ifelse(p_npp<0.001, "p<0.001", sprintf("p=%.3f",p_npp))


color_fvc <- "#1f77b4"   # 蓝色
color_npp <- "#ff7f0e"   # 橙色


p_e <- ggplot() +
  geom_segment(data = year_summary,
               aes(x = year, xend = year,
                   y = dFVC_mean - dFVC_sd05,
                   yend = dFVC_mean + dFVC_sd05),
               color = color_fvc) +
  geom_point(data = year_summary,
             aes(x = year, y = dFVC_mean),
             shape = 21, fill = color_fvc, size = 3, stroke = 1) +
  geom_smooth(data = city_final_dt_year,
              aes(x = year, y = delta_FVC),
              method = "lm", se = FALSE, color = "black") +
  scale_x_continuous(breaks = c(2000,2005,2010,2015,2020), limits = c(1999,2021)) +
  labs(x = "Year", y = expression(Delta~FVC)) +
  annotate("text", x = 2010, y = Inf, label = lab_fvc,
           hjust = 0.5, vjust = 2, size = 4.5) +
  annotate("text", x = 2010, y = -Inf, label = ptext_fvc,
           hjust = 0.5, vjust = -1.2, size = 5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.title = element_text(size = 13, colour = "black"),
        axis.text  = element_text(size = 12, colour = "black"))

print(p_e)


###########################################################################
#                Temporal dynamics of ∆NPP over 2000–2020
###########################################################################

p_i <- ggplot() +
  geom_segment(data = year_summary,
               aes(x = year, xend = year,
                   y = dNPP_mean - dNPP_sd05,
                   yend = dNPP_mean + dNPP_sd05),
               color = color_npp) +
  geom_point(data = year_summary,
             aes(x = year, y = dNPP_mean),
             shape = 21, fill = color_npp, size = 3, stroke = 1) +
  geom_smooth(data = city_final_dt_year,
              aes(x = year, y = delta_NPP),
              method = "lm", se = FALSE, color = "black") +
  scale_x_continuous(breaks = c(2000,2005,2010,2015,2020), limits = c(1999,2021)) +
  labs(x = "Year", y = expression(Delta * NPP ~ "(g C" ~ m^-2 ~ yr^-1 * ")")) +
  annotate("text", x = 2010, y = Inf, label = lab_npp,
           hjust = 0.5, vjust = 2, size = 4.5) +
  annotate("text", x = 2010, y = -Inf, label = ptext_npp,
           hjust = 0.5, vjust = -1.2, size = 5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.title = element_text(size = 13, colour = "black"),
        axis.text  = element_text(size = 12, colour = "black"))

print(p_i)





###########################################################################
#   City-level comparison of mean FVC between urban and non-urban areas
###########################################################################

out_dt <- merged_dt[
  urban == "out" & year %in% c(2000,2005,2010,2015,2020),
  .(CityID, FVC, NPP, year)
]

city_year_dt <- out_dt[,
                       .(
                         FVC_year_mean = mean(FVC, na.rm = TRUE),
                         NPP_year_mean = mean(NPP, na.rm = TRUE)
                       ),
                       by = .(CityID, year)
]


city_out_final <- city_year_dt[,
                               .(
                                 FVC_out = mean(FVC_year_mean, na.rm = TRUE),
                                 NPP_out = mean(NPP_year_mean, na.rm = TRUE)
                               ),
                               by = .(CityID)
]


in_dt <- merged_dt[
  urban == "in" & year %in% c(2000,2005,2010,2015,2020),
  .(CityID, FVC, NPP, year, Latitude, Longitude)
]

city_in_final <- in_dt[,
                       .(
                         FVC_in = mean(FVC, na.rm = TRUE),
                         NPP_in = mean(NPP, na.rm = TRUE),
                         Latitude = mean(Latitude, na.rm = TRUE),
                         Longitude = mean(Longitude, na.rm = TRUE)
                       ),
                       by = .(CityID)
]


city_final_dt <- merge(
  city_in_final,
  city_out_final,
  by = "CityID",
  all = FALSE
)

city_final_dt$cha_FVC <- city_final_dt$FVC_in-city_final_dt$FVC_out
city_final_dt$cha_NPP <- city_final_dt$NPP_in-city_final_dt$NPP_out


city_final_dt[, region_cat_fvc := fcase(
  FVC_out < FVC_in, "below_diag",   
  FVC_out > FVC_in, "above_diag",   
  default = "on_diag"
)]

cor_test_res <- cor.test(city_final_dt$FVC_in, city_final_dt$FVC_out)
r_val <- round(cor_test_res$estimate, 2)
n_total <- nrow(city_final_dt)
n_below <- sum(city_final_dt$region_cat_fvc == "below_diag")
n_above <- sum(city_final_dt$region_cat_fvc == "above_diag")
pct_below <- round(n_below / n_total * 100, 0)
pct_above <- round(n_above / n_total * 100, 0)

text_r <- sprintf("r = %.2f***", r_val)
y_pos_ratio <- max(city_final_dt$FVC_out)*1
x_start <- 0.04

p_f <- ggplot(city_final_dt, aes(x = FVC_in, y = FVC_out)) +
  # geom_point(aes(color = region_cat_fvc), shape = 1, size = 2.2, alpha = 0.6) +
  geom_scattermore(aes(color = region_cat_fvc),
                   pointsize = 4,
                   pixels = c(800, 800),
                   interpolate = TRUE) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color="black") +
  scale_color_manual(values = c("below_diag"="#dd3333", "above_diag"="#3377bb", "on_diag"="gray50")) +
  annotate("text", x = x_start, y = max(city_final_dt$FVC_out)*1, label = text_r, hjust=0, size=4.5, colour="black") +
  annotate("text", x = x_start+0.460, y = y_pos_ratio, label = paste0(pct_above,"%"), hjust=0, size=4.5, colour="#3377bb") +
  annotate("text", x = x_start+0.595, y = y_pos_ratio, label = " : ", hjust=0, size=4.5, colour="black") +
  annotate("text", x = x_start+0.645, y = y_pos_ratio, label = paste0(pct_below,"%"), hjust=0, size=4.5, colour="#dd3333") +
  labs(
    x = expression(FVC[Urban]),
    y = expression(FVC[Non-urban]),
    subtitle = ""
  ) +
  scale_x_continuous(limits = c(0, 0.9), 
                     breaks = seq(0, 0.9, by =0.3),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 0.9),
                     breaks = seq(0, 0.9, by = 0.3),
                     labels = scales::number_format(accuracy = 0.1),
                     expand = c(0.01, 0.01)) +
  theme_bw() +
  theme(
    plot.subtitle = element_text(hjust=0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = "none",
    axis.title = element_text(size=13, colour = "black"),
    axis.text = element_text(size=12, colour = "black")
  )
print(p_f)


###########################################################################
#   City-level comparison of mean NPP between urban and non-urban areas
###########################################################################

city_final_dt[, region_cat_npp := fcase(
  NPP_out < NPP_in, "below_diag",
  NPP_out > NPP_in, "above_diag",
  default = "on_diag"
)]
cor_npp <- cor.test(city_final_dt$NPP_in, city_final_dt$NPP_out)
r_npp <- round(cor_npp$estimate,2)
n_below_npp <- sum(city_final_dt$region_cat_npp=="below_diag")
n_total_npp <- nrow(city_final_dt)
pct_below_npp <- round(n_below_npp/n_total_npp*100,0)
pct_above_npp <- 100-pct_below_npp

text_r_npp <- sprintf("r = %.2f***", r_npp)
y_pos_ratio_npp <- max(city_final_dt$NPP_out)*1
x_start_npp <- 72

p_j <- ggplot(city_final_dt, aes(x=NPP_in, y=NPP_out))+
  # geom_point(aes(color=region_cat_npp), shape=1, size=2.2, alpha=0.6)+
  geom_scattermore(aes(color=region_cat_npp),
                   pointsize = 4,
                   pixels = c(800,800),
                   interpolate = TRUE) +
  geom_smooth(method="lm", se=FALSE, color="black")+
  geom_abline(slope=1, intercept=0, linetype="dashed", color="black")+
  scale_color_manual(values=c("below_diag"="#dd3333","above_diag"="#3377bb","on_diag"="gray50"))+
  annotate("text", x = x_start_npp, y = max(city_final_dt$NPP_out)*1, label = text_r_npp, hjust=0, size=4.5, colour="black") +
  annotate("text", x = x_start_npp+1250, y = y_pos_ratio_npp, label = paste0(pct_below_npp,"%"), hjust=0, size=4.5, colour="#dd3333") +
  annotate("text", x = x_start_npp+1150, y = y_pos_ratio_npp, label = " : ", hjust=0, size=4.5, colour="black") +
  annotate("text", x = x_start_npp+900, y = y_pos_ratio_npp, label = paste0(pct_above_npp,"%"), hjust=0, size=4.5, colour="#3377bb") +
  labs(x=expression(NPP[Urban]~~(g~C~m^{-2}~yr^{-1})), y=expression(NPP[Non-urban]~~(g~C~m^{-2}~yr^{-1})), subtitle="")+
  scale_x_continuous(limits = c(0, 1700), 
                     breaks = seq(0, 1700, by =400),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1700),
                     breaks = seq(0, 1700, by = 400),
                     expand = c(0.01, 0.01)) +
  theme_bw()+
  theme(legend.position="none", 
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
        axis.ticks = element_line(color = "black", size = 0.5),
        plot.subtitle=element_text(hjust=0.5),
        axis.title = element_text(size=13, colour = "black"),
        axis.text = element_text(size=12, colour = "black"))
print(p_j)




###########################################################################
#                     Spatial distribution of mean ∆FVC
###########################################################################

# ---------------------- 1 读取shp数据 ----------------------
world_gdf <- st_read("word_map.shp")
city_df  <- city_final_dt

world_gdf <- st_transform(world_gdf, 4326)
city_gdf  <- st_as_sf(city_df, coords = c("Longitude", "Latitude"), crs = 4326)

draw_theme <- theme_bw() +
  theme(
    panel.grid      = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title      = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    plot.margin     = margin(2, 2, 2, 2),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, size = 0.75)
  )


make_rwg_cmap <- function(vmin, vmax, n = 256, red = "#8b0000", white = "white", green = "#006400") {
  zero_pos <- (0 - vmin) / (vmax - vmin)
  cols <- colorRampPalette(c(red, white, green))(n)
  n1 <- max(1, round(zero_pos * n))
  n2 <- n - n1
  c(colorRampPalette(c(red, white))(n1),
    colorRampPalette(c(white, green))(n2))
}


city_gdf$x <- st_coordinates(city_gdf)[, 1]
city_gdf$y <- st_coordinates(city_gdf)[, 2]


vmin_fvc <- -0.67
vmax_fvc <-  0.73
cols_fvc <- make_rwg_cmap(vmin_fvc, vmax_fvc)

p_g <- ggplot() +
  draw_theme +
  geom_sf(data = world_gdf, fill = "white", color = "#999999", linewidth = 0.3) +
  geom_scattermore(
    data = city_gdf,
    aes(x = x, y = y, color = cha_FVC),
    pointsize = 4,
    pixels = c(1600, 800),
    interpolate = TRUE
  ) +
  scale_color_gradientn(
    colors = cols_fvc,
    limits = c(vmin_fvc, vmax_fvc),
    oob    = scales::squish,
    # name   = expression(Delta * FVC == FVC[urban] - FVC[non-urban]),
    name   = expression(Delta * FVC),
    guide  = guide_colorbar(
      direction       = "horizontal",
      barwidth        = unit(35, "mm"),
      barheight       = unit(2, "mm"),
      title.position  = "top",
      title.hjust     = 0.5,
      ticks.colour    = "black",
      frame.colour    = "black",
      frame.linewidth = 0.4
    )
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 75), expand = FALSE) +
  theme(
    aspect.ratio = 0.4,
    legend.position   = c(0.02, 0.16),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.title      = element_text(size = 12),
    legend.text       = element_text(size = 12),
    legend.margin     = margin(1, 1, 1, 1)
  )

print(p_g)



###########################################################################
#                     Spatial distribution of mean ∆NPP
###########################################################################
vmin_npp <- -1135.485
vmax_npp <-  402.416
zero_norm_pos <- (0 - vmin_npp) / (vmax_npp - vmin_npp)
message(sprintf("0值对应的色板归一化位置：%.2f", zero_norm_pos))

cols_npp <- make_rwg_cmap(vmin_npp, vmax_npp)

p_k <- ggplot() +
  draw_theme +
  geom_sf(data = world_gdf, fill = "white", color = "#999999", linewidth = 0.3) +
  # geom_sf(data = city_gdf, aes(color = cha_NPP),
  #         size = 1.6, shape = 16, stroke = 0) +
  geom_scattermore(
    data = city_gdf,
    aes(x = x, y = y, color = cha_NPP),
    pointsize = 4,
    pixels = c(1600, 800),
    interpolate = TRUE
  ) +
  scale_color_gradientn(
    colors = cols_npp,
    limits = c(vmin_npp, vmax_npp),
    oob    = scales::squish,
    name   = expression(atop(Delta * NPP, "(g C" ~ m^-2 ~ yr^-1 * ")")), 
    guide  = guide_colorbar(
      direction       = "horizontal",
      barwidth        = unit(35, "mm"),
      barheight       = unit(2, "mm"),
      title.position  = "top",
      title.hjust     = 0.5,
      ticks.colour    = "black",
      frame.colour    = "black",
      frame.linewidth = 0.4
    )
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 75), expand = FALSE) +
  theme(
    aspect.ratio = 0.4,
    legend.position      = c(0.02, 0.16),
    legend.justification = c(0, 0),
    legend.background    = element_blank(),
    legend.title         = element_text(size = 12),
    legend.text          = element_text(size = 12),
    legend.margin        = margin(1, 1, 1, 1)
  )

print(p_k)



###########################################################################
#                     Latitudinal variation of mean ∆FVC
###########################################################################
df <- as.data.frame(city_final_dt)

bin_width <- 2
df <- df %>%
  mutate(lat_bin = cut(Latitude,
                       breaks = seq(-60, 76, by = bin_width), 
                       include.lowest = TRUE,
                       labels = FALSE))

lat_profile <- df %>%
  group_by(lat_bin) %>%
  summarise(
    lat_mean = mean(Latitude, na.rm = TRUE),
    dFVC_mean = mean(cha_FVC, na.rm = TRUE),
    n = n(),
    se = sd(cha_FVC, na.rm = TRUE)/sqrt(n),
    ci95_low = dFVC_mean - 1.96*se,
    ci95_high = dFVC_mean + 1.96*se
  ) %>%
  filter(n >= 5) 


p_h <- ggplot(lat_profile) +
  geom_ribbon(aes(xmin = ci95_low, xmax = ci95_high, y = lat_mean),
              fill = "gray50", alpha = 0.4) +
  geom_path(aes(x = dFVC_mean, y = lat_mean), color = "black", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "black", linetype = 2) +
  scale_y_continuous(
    limits = c(-60,75),
    breaks = seq(-60,75,25),
    name = expression(Latitude~(degree*N))
  ) +
  scale_x_continuous(name = expression(Delta~FVC)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.title = element_text(size=13, colour = "black"),
    axis.text = element_text(size=12, colour = "black")
  )

print(p_h)


###########################################################################
#                     Latitudinal variation of mean ∆NPP
###########################################################################

df <- as.data.frame(city_final_dt)

bin_width <- 2
df <- df %>%
  mutate(lat_bin = cut(Latitude,
                       breaks = seq(-60, 76, by = bin_width), 
                       include.lowest = TRUE,
                       labels = FALSE))

lat_profile <- df %>%
  group_by(lat_bin) %>%
  summarise(
    lat_mean = mean(Latitude, na.rm = TRUE),
    dNPP_mean = mean(cha_NPP, na.rm = TRUE),
    n = n(),
    se = sd(cha_NPP, na.rm = TRUE)/sqrt(n),
    ci95_low = dNPP_mean - 1.96*se,
    ci95_high = dNPP_mean + 1.96*se
  ) %>%
  filter(n >= 5) 


p_l <- ggplot(lat_profile) +
  geom_ribbon(aes(xmin = ci95_low, xmax = ci95_high, y = lat_mean),
              fill = "gray50", alpha = 0.4) +
  geom_path(aes(x = dNPP_mean, y = lat_mean), color = "black", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "black", linetype = 2) +
  scale_y_continuous(
    limits = c(-60,75),
    breaks = seq(-60,75,25),
    name = expression(Latitude~(degree*N))
  ) +
  # scale_x_continuous(name = expression(Delta~NPP)) +
  scale_x_continuous(
    name = expression(Delta * NPP ~ "(g C" ~ m^-2 ~ yr^-1 * ")")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.title = element_text(size=13, colour = "black"),
    axis.text = element_text(size=12, colour = "black")
  )
print(p_l)



strip_margin <- function(p) {
  p + theme(plot.margin = margin(2, 4, 2, 4, unit = "pt"))
}

p_1 <- strip_margin(p_a)
p_2 <- strip_margin(p_b)
p_3 <- strip_margin(p_c)
p_4 <- strip_margin(p_d)
p_5 <- strip_margin(p_f)
p_6 <- strip_margin(p_g)
p_7 <- strip_margin(p_h)
p_8 <- strip_margin(p_e)
p_9 <- strip_margin(p_j)
p_10 <- strip_margin(p_k)
p_11 <- strip_margin(p_l)
p_12 <- strip_margin(p_i)

aligned_59 <- align_plots(
  p_5, p_9,
  align = "v",
  axis  = "tb"
)
p_5_aligned <- aligned_59[[1]]
p_9_aligned <- aligned_59[[2]]

p_6 <- p_6 + theme(aspect.ratio = 0.425)
p_10 <- p_10 + theme(aspect.ratio = 0.419)

make_row <- function(row_list, rel_widths) {
  plot_grid(
    plotlist   = row_list,
    nrow       = 1, ncol = 4,
    align      = "h",
    axis       = "tb",
    rel_widths = rel_widths
  )
}

col_w <- c(2.2, 2.2, 4.2, 1.6)

row0_plot <- make_row(list(p_1, p_2, p_3, p_4), rel_widths = c(1, 1, 1, 1))
row1_plot <- make_row(list(p_8, p_5_aligned, p_6, p_7), rel_widths = col_w)
row2_plot <- make_row(list(p_12, p_9_aligned, p_10, p_11), rel_widths = col_w)


fig_3x4 <- plot_grid(
  row0_plot, row1_plot, row2_plot,
  ncol        = 1,
  align       = "v",
  axis        = "lr",
  rel_heights = c(1, 1, 1)
)

top_pad <- 0.03

fig_3x4_padded <- ggdraw() +
  draw_plot(fig_3x4, x = 0, y = 0, width = 1, height = 1 - top_pad)

row_h   <- (1 - top_pad) / 3
y_row1  <- 1 - top_pad              
y_row2  <- 1 - top_pad - row_h      
y_row3  <- 1 - top_pad - 2 * row_h  

y_tag_row0 <- rep(y_row1 + 0.020, 4)
y_tag_row1 <- rep(y_row2 + 0.003, 4)
y_tag_row2 <- rep(y_row3 + 0.003, 4)

x_tag_row0 <- c(0.005, 0.255, 0.505, 0.755)  
x_tag_row1 <- c(0.005, 0.240, 0.435, 0.850)  
x_tag_row2 <- c(0.005, 0.240, 0.435, 0.850)  

fig_3x4_padded <- fig_3x4_padded +
  draw_plot_label(
    label    = c("a","b","c","d",
                 "e","f","g","h",
                 "i","j","k","l"),
    x        = c(x_tag_row0,
                 x_tag_row1,
                 x_tag_row2),
    y        = c(y_tag_row0,
                 y_tag_row1,
                 y_tag_row2),
    size     = 14, fontface = "bold", family = "Times",
    hjust    = 0, vjust = 1
  )

panel_w_mm <- 62
panel_h_mm <- 60
inflate_w  <- 1.35
inflate_h  <- 1.15

fig_w_mm <- panel_w_mm * 4 * inflate_w + 22
fig_h_mm <- (panel_h_mm * 3 * inflate_h + 20) / (1 - top_pad) 

ggsave("./Fig1_3x4.pdf",
       plot   = fig_3x4_padded,
       width  = fig_w_mm,
       height = fig_h_mm,
       units  = "mm",
       dpi    = 300)



