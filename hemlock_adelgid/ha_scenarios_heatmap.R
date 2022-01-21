library(tidyverse)
library(ggplot2)

# import data from csv files produced by ha_group_1_scenarios.m
# and ha_group_2_scenarios.m, saved as tibble, with
# column names added and 0 values replaces with NaN
# (indicating that 0.1 proportion tips alive was not
# reached within 15 year simulation)
group_1 <- read_csv("ha_group_1_deadtime.csv", col_names = FALSE) %>%
  rename(H = X1, A = X2, Years = X3) %>%
  mutate(Years = recode(Years, `0` = NaN))
  
group_2 <- read_csv("ha_group_2_deadtime.csv", col_names = FALSE) %>%
  rename(H = X1, A = X2, Years = X3) %>%
  mutate(Years = recode(Years, `0` = NaN))

# heat map for group 1
heat1 <- ggplot(group_1, aes(A, H)) +
  geom_tile(aes(fill = Years), color = "grey70", size = 0.1) +
  scale_fill_viridis_c(na.value = "white") +
  labs(y = "Proportion tips alive", x = "Adelgid density (per cm)") +
  theme_bw(base_size = 14) +
  theme(
    axis.ticks = element_line(size = 0.5),
    plot.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-.3, "cm"))
heat1
ggsave("./ha_group_1_heatmap.tiff", heat1,
       width = 6, height = 5, dpi = 600, bg = "white")

# heat map for group 2
heat2 <- ggplot(group_2, aes(A, H)) +
  geom_tile(aes(fill = Years), color = "grey70", size = 0.1) +
  scale_fill_viridis_c(na.value = "white") +
  labs(y = "Proportion tips alive", x = "Adelgid density (per cm)") +
  theme_bw(base_size = 14) +
  theme(
    axis.ticks = element_line(size = 0.5),
    plot.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-.3, "cm"))
heat2
ggsave("./ha_group_2_heatmap.tiff", heat2,
       width = 6, height = 5, dpi = 600, bg = "white")
