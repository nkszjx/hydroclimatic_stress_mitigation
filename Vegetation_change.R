

# Changes in FVC and NPP from 2000 to 2020

library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(scales)



###########################################################################
#       Continental-averaged temporal variations of urban FVC
###########################################################################

setwd("XXX") # please set your folder path
merged_dt <- fread("Dataset_matched_buffer10km.csv")


urban_FVC_in <- urban_FVC[urban == "in"]
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

p1 <- ggplot(urban_FVC_in_mean, aes(x = year, y = FVC_mean, color = continent, group = continent)) +
  geom_line(linewidth = 0.5, alpha = 1) +
  geom_point(size = 1.5, stroke = 0, shape = 19) +
  scale_x_continuous(
    name = "",
    breaks = seq(2000, 2020, 5),
    expand = c(0.05, 0.05)
  ) +
  scale_y_continuous(
    name = "Urban FVC",
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
    panel.border = element_rect(linewidth = 0.6, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 11, family = "serif", color = "black"),
    axis.title.x = element_text(size = 12, family = "serif", color = "black", margin = margin(t = 5)),
    axis.title.y = element_text(size = 12, family = "serif", color = "black", margin = margin(r = 5)),
    axis.ticks.length = unit(0.05, "cm"),
    axis.ticks = element_line(color = "black", size = 0.4),
    legend.position = c(0.50, 0.85),
    legend.title = element_blank(),
    legend.text = element_text(size = 11, family = "serif", color = "black", lineheight = 0.7),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 3, ncol = 2, byrow = TRUE))

p1

ggsave("Continental_averaged_urban_FVC.pdf", 
       plot = p1, 
       width = 60, 
       height = 60, 
       units = "mm",
       dpi = 1000)





###########################################################################
#       Continental-averaged temporal variations of non-urban FVC
###########################################################################

merged_dt <- fread("Dataset_match_buffer10.csv")

nonurban_FVC_out <- nonurban_FVC[urban == "out"]
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

p2 <- ggplot(nonurban_FVC_out_mean, aes(x = year, y = FVC_mean, color = continent, group = continent)) +
  geom_line(linewidth = 0.5, alpha = 1) +
  geom_point(size = 1.5, stroke = 0, shape = 19) +
  scale_x_continuous(
    name = "",
    breaks = seq(2000, 2020, 5),
    expand = c(0.05, 0.05)
  ) +
  scale_y_continuous(
    name = "Non-urban FVC",
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
    panel.border = element_rect(linewidth = 0.6, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 11, family = "serif", color = "black"),
    axis.title.x = element_text(size = 12, family = "serif", color = "black", margin = margin(t = 5)),
    axis.title.y = element_text(size = 12, family = "serif", color = "black", margin = margin(r = 5)),
    axis.ticks.length = unit(0.05, "cm"),
    axis.ticks = element_line(color = "black", size = 0.4),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, family = "serif", color = "black", lineheight = 0.7),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 3, ncol = 2, byrow = TRUE))

p2

ggsave("Continental_averaged_nonurban_FVC.pdf", 
       plot = p2, 
       width = 60, 
       height = 60, 
       units = "mm",
       dpi = 1000)





###########################################################################
#       Continental-averaged temporal variations of urban NPP
###########################################################################

merged_dt <- fread("Dataset_match_buffer10.csv")

urban_NPP_in <- urban_NPP[urban == "in"]

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

p3 <- ggplot(urban_NPP_in_mean, aes(x = year, y = NPP_mean, color = continent, group = continent)) +
  geom_line(linewidth = 0.5, alpha = 1) +
  geom_point(size = 1.5, stroke = 0, shape = 19) +
  scale_x_continuous(
    name = "",
    breaks = seq(2000, 2020, 5),
    expand = c(0.05, 0.05)
  ) +
  scale_y_continuous(
    name = expression(Urban~NPP~(g~C~m^-2~yr^-1)),
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
    panel.border = element_rect(linewidth = 0.6, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 11, family = "serif", color = "black"),
    axis.title.x = element_text(size = 12, family = "serif", color = "black", margin = margin(t = 5)),
    axis.title.y = element_text(size = 12, family = "serif", color = "black", margin = margin(r = 5)),
    axis.ticks.length = unit(0.05, "cm"),
    axis.ticks = element_line(color = "black", size = 0.4),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, family = "serif", color = "black", lineheight = 0.7),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 3, ncol = 2, byrow = TRUE))

p3

ggsave("Continental_averaged_urban_NPP.pdf", 
       plot = p3, 
       width = 63, 
       height = 60, 
       units = "mm",
       dpi = 1000)



###########################################################################
#       Continental-averaged temporal variations of nonurban NPP
###########################################################################

merged_dt <- fread("Dataset_match_buffer10.csv")

nonurban_NPP_out <- nonurban_NPP[urban == "out"]

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

p2 <- ggplot(nonurban_NPP_out_mean, aes(x = year, y = NPP_mean, color = continent, group = continent)) +
  geom_line(linewidth = 0.5, alpha = 1) +
  geom_point(size = 1.5, stroke = 0, shape = 19) +
  scale_x_continuous(
    name = "",
    breaks = seq(2000, 2020, 5),
    expand = c(0.05, 0.05)
  ) +
  scale_y_continuous(
    name = expression(paste("Non-urban")~NPP~(g~C~m^-2~yr^-1)),
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
    panel.border = element_rect(linewidth = 0.6, color = "black"),
    axis.text.x = element_text(size = 11, family = "serif", color = "black"),
    axis.text.y = element_text(size = 11, family = "serif", color = "black"),
    axis.title.x = element_text(size = 12, family = "serif", color = "black", margin = margin(t = 5)),
    axis.title.y = element_text(size = 12, family = "serif", color = "black", margin = margin(r = 5)),
    axis.ticks.length = unit(0.05, "cm"),
    axis.ticks = element_line(color = "black", size = 0.4),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, family = "serif", color = "black", lineheight = 0.7),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 3, ncol = 2, byrow = TRUE))

p4

ggsave("Continental_averaged_nonurban_NPP.pdf", 
       plot = p4, 
       width = 63, 
       height = 62.5, 
       units = "mm",
       dpi = 1000)


