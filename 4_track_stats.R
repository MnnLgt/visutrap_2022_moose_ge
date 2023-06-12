#------------------------------
# UVP6 MOOSE-GE track data
#------------------------------
# Track analysis
# 08-02-2023
# Manon Laget
#------------------------------


library(tidyverse)
library(cluster)
library(languageR)
library(ggExtra)
library(gridExtra)
library(cowplot)

path <- "C:\\Users\\Manon Laget\\Desktop\\UVP6_data\\projects\\paper_methodo\\analyses\\"
setwd(path)


# 1. Load tracks ----

# Read track files
tracks <- list.files(path = 'data\\track_data\\', pattern = '*.tsv',
                     full.names = TRUE) %>%
lapply(read.csv, sep = '\t') %>%
  bind_rows() %>%
  mutate(deployment = paste(Cycle, Depth, sep = "_")) %>%
  mutate(vertical_speed = 
           ifelse(angle_mean > 180, vertical_speed, -vertical_speed))

# Add path to vignettes
tracks$path_to_vig <- file.path(
  paste0(path, "data\\track_data\\img_to_show\\", tracks$vig_name))

# Read EcoTaxa file and remove living organisms
ecotaxa <- read.csv('data\\track_data\\ecotaxa\\ecotaxa_export_moose.tsv', 
                    sep = '\t') %>%
  select(track_id = 1, everything())

tracks <- left_join(tracks, ecotaxa, by = 'track_id') %>%
  filter(object_annotation_category != "false"
         & object_annotation_category != "dubious") %>%
  filter(esd_um >= 600)

counts_living <- tracks %>% filter(object_annotation_category != "true",
                                   object_annotation_category != "t001") %>%
  group_by(deployment) %>% dplyr::summarise(count = n())

tracks <- tracks %>% filter(object_annotation_category != "living<"
                            & object_annotation_category != "copepod" 
                            & object_annotation_category != "Rhizaria"
                            & object_annotation_category != "Ctenophora<Metazoa"
                            & object_annotation_category != "t010")

counts_particles <- tracks %>% group_by(deployment) %>% dplyr::summarise(n = n())


# 2. Create track categories ----

# Assign category depending on vertical speed
tracks <- tracks %>% 
  mutate(type = case_when(
    orientation == 'desc' ~ 'sinking',
    orientation == 'mix' ~ 'suspended',
    TRUE ~ 'ascending'
  ))

# Save this new dataframe which will be used for further analyses
write_tsv(tracks, 'results\\track_data.tsv')

# Mean trap depth
depths <- read.csv('data\\trap_data\\trap_depth_moose22_c1_500.csv') %>%
  mutate(deployment = 'c1_500') %>%
  bind_rows(read.csv('data\\trap_data\\trap_depth_moose22_c2_200.csv') %>%
              mutate(deployment = 'c2_200')) %>%
  bind_rows(read.csv('data\\trap_data\\trap_depth_moose22_c2_500.csv') %>%
              mutate(deployment = 'c2_500')) %>%
  # Calculate speed at each timestep
  group_by(deployment) %>% 
  dplyr::summarize(mean_trap_depth = mean(averaged_depth, na.rm = T))


# 3. Track statistics ----

# Basic statistics
mean(tracks$length)
sd(tracks$length)
max(tracks$length)
median(tracks$length)
max(tracks$esd_um)

# Number of tracks for each deployment
tracks %>% group_by(deployment) %>%
  dplyr::summarise(count = n())


# 4. Plots ----

cols <- c("#fc8d62", "#8da0cb", "#66c2a5")

# Size and speed spectra

(p1 <- ggplot(tracks, aes(x = esd_um, y = vertical_speed, color = type)) +
    geom_point(alpha = 0.6) + 
    scale_color_manual('Track type', values = cols, 
                       limits = c("ascending", "suspended", "sinking"),
                       labels = c("ascending", "suspended", "sinking")) +
    scale_x_log10(minor_breaks = NULL) +
    annotation_logticks(sides="b") +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = c(0.85, 0.85),
      legend.background =element_rect(size=0.5, linetype="solid", 
                                      colour ="black")) +
    xlab('Equivalent Spherical Diameter (µm)') +
    ylab('Vertical velocity (m/d)'))

(p1 <- ggMarginal(p1, type = "histogram", groupColour = T, groupFill = T))


# Proportion of each type for each deployment
tracks$type <- factor(tracks$type, levels = c('ascending', 'suspended', 'sinking'))
(p2 <- ggplot(tracks,
       aes(x = deployment, y = ..count../sum(..count..), 
           color = type, fill = type)) +
  geom_bar(position = "fill", width = 0.6, alpha = 0.5) +
  scale_fill_manual('Track type', values = cols, 
                    limits = c("ascending", "suspended", "sinking"),
                    labels = c("ascending", "suspended", "sinking")) +
  scale_color_manual('Track type', values = cols, 
                     limits = c("ascending", "suspended", "sinking"),
                     labels = c("ascending", "suspended", "sinking")) +
  scale_x_discrete(
    breaks = c("moose22_c1_500", "moose22_c2_200", "moose22_c2_500"),
    labels = c("C1 500 m", "C2 200 m", "C2 500 m")) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)) +
  xlab('Deployment') +
  ylab('Proportion of track types'))

# Mean angle
(p3 <- ggplot(tracks, aes(x = angle_mean, y = abs(vertical_speed), color = type)) +
  geom_point(alpha = 0.6) + 
  scale_color_manual('Track type', values = cols, 
                     limits = c("ascending", "suspended", "sinking"),
                     labels = c("ascending", "suspended", "sinking")) +
  coord_polar(start = -90 * pi / 180, direction = -1) +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)) +
  xlab('Mean angle (°)') +
  ylab('Absolute vertical velocity (m/d)'))

(all <- ggdraw() +
  draw_plot(p1, 0, 0, .5, 1) +
  draw_plot(p3, .5, 0, .5, .6) +
  draw_plot(p2, .53, .62, .45, .35) +
  draw_plot_label(c("(a)", "(b)", "(c)"), c(0, 0.5, 0.5), c(1, 1, 0.6), size = 15))

ggsave('figures\\plots_stats.pdf', all, height = 6, width = 12)
ggsave('figures\\plots_stats.png', all, height = 6, width = 12)
