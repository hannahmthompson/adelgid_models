library(tidyverse)
#library(ggplot2)

# data from excel files, saved in data frames
group1_raw <- read_csv("HA_time_to_likely_death_group_1.csv") 

group1 <- pivot_longer(group1_raw, cols= colnames(group1_raw)[-c(1)],
               names_to = "H", values_to = "years")
  
colnames(group1)[1] <- c("A")
  
  
group2_raw <- read_csv("HA_time_to_likely_death_group_2.csv")

group1 <- pivot_longer(group1_raw, cols= colnames(group1_raw)[-c(1)], names_to = "H", values_to = "years")

colnames(group1)[1] <- c("A")
group1["years"]["years" > 15 | "years" == 0] <- NaN


# replace values that indicate that proportion tips alive does not reach 0.15
# within 15 year test with NaN
group1_raw[group1_raw > 15 | group1_raw==0] <- NaN
group2_raw[group2_raw > 15 | group2_raw==0] <- NaN

#reformqte
group1_melt <- melt(group1_raw)
group2_melt <- melt(group2_raw)
colnames(group1_melt) <- c("A","H","Years")
colnames(group2_melt) <- c("A","H","Years")

heat1 <- ggplot(group1_melt, aes(A, H)) + 
  geom_tile(aes(fill=Years), color="grey90", size=0.1) + 
  scale_fill_viridis_c(na.value="white", limits = c(0,15)) +  
  labs(y = "Proportion tips alive", x= "Adelgid density (per cm)") +
  theme_bw(base_size=14) +
  theme(
  axis.ticks=element_line(size=0.5),
  plot.background=element_blank(),
  panel.border=element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.ticks.length = unit(-.3, "cm"))
heat1

heat1b <- ggplot(group1_melt, aes(A, H)) + 
  geom_tile(aes(fill=Years), color="grey30", size=0.1) + 
  scale_fill_viridis_c(na.value="white") +  
  labs(y = "Proportion tips alive", x= "Adelgid density (per cm)") +
  theme_bw(base_size=14) +
  theme(
    axis.ticks=element_line(size=0.5),
    plot.background=element_blank(),
    panel.border=element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-.3, "cm"))
heat1b

heat2 <- ggplot(group2_melt, aes(A, H)) + 
  geom_tile(aes(fill=Years), color="black", size=0.1) + 
  scale_fill_viridis_c(na.value="white", limits = c(0,15)) +  
  labs(y = "Proportion tips alive", x= "Adelgid density (per cm)") +
  theme_bw(base_size=14) +
  theme(
    axis.ticks=element_line(size=0.5),
    plot.background=element_blank(),
    panel.border=element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-.3, "cm"))
heat2


