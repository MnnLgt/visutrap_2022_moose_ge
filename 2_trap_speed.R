#------------------------------
# UVP6 MOOSE-GE track data
#------------------------------
# Trap speeds
# 23-01-2023
# Manon Laget
#------------------------------


library(tidyverse)
library(plyr)
library(lubridate)
library(cowplot)

path <- "C:\\Users\\Manon Laget\\Desktop\\UVP6_data\\projects\\paper_methodo\\analyses\\"
setwd(path)


# 1. Read depth files ----

depths <- read.csv('data\\trap_data\\trap_depth_moose22_c1_500.csv') %>%
  mutate(deployment = 'c1_500') %>%
  bind_rows(read.csv('data\\trap_data\\trap_depth_moose22_c2_200.csv') %>%
              mutate(deployment = 'c2_200')) %>%
  bind_rows(read.csv('data\\trap_data\\trap_depth_moose22_c2_500.csv') %>%
              mutate(deployment = 'c2_500')) %>%
  mutate(datetime = ymd_hms(datetime, tz = "UTC")) %>%
  # Calculate speed at each timestep
  group_by(deployment, sequence) %>% 
  mutate(dt = as.numeric(datetime - lag(datetime, 1)),
         ddepth = averaged_depth - lag(averaged_depth, 1)) %>% 
  ungroup() %>%
  mutate(speed = (ddepth/dt)*100)

# Calculate mean speed and depth
summ_depths <- depths %>% group_by(deployment, sequence) %>% 
  dplyr::summarize(trap_depth_sd = sd(averaged_depth, na.rm = T),
            trap_speed_mean = mean(abs(speed), na.rm = T),
            trap_speed_sd = sd(abs(speed), na.rm = T))


# 2. Read track files ----

tracks <- read.csv('data\\trap_data\\trap_speeds_moose22_c1_500.csv') %>%
  mutate(deployment = 'c1_500') %>%
  bind_rows(read.csv('data\\trap_data\\trap_speeds_moose22_c2_200.csv') %>%
              mutate(deployment = 'c2_200')) %>%
  bind_rows(read.csv('data\\trap_data\\trap_speeds_moose22_c2_500.csv') %>%
              mutate(deployment = 'c2_500')) %>%
  mutate(object_id = track_id)

# Read EcoTaxa files
ecotaxa <- read_tsv('data\\track_data\\ecotaxa\\ecotaxa_export_moose.tsv') %>%
  select(object_id = 1, everything())

# Remove living organisms
tracks <- left_join(tracks, ecotaxa, by = 'object_id') %>%
  filter(object_annotation_category != "false"
         & object_annotation_category != "dubious" 
         & object_annotation_category != "living<" 
         & object_annotation_category!="copepod" 
         & object_annotation_category!="Rhizaria"
         & object_annotation_category!="Ctenophora<Metazoa"
         & object_annotation_category!="t010") %>%
  filter(object_esd >= 600)


# 3. Sd depth and mean speed per sequence ----

ggplot(summ_depths, aes(x = sequence, y = trap_depth_sd, color = deployment)) +
  geom_point() +
  xlab('Sequence') +
  ylab("Sd trap depth per sequence (m)") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggplot(summ_depths, aes(x = sequence, y = trap_speed_mean, color = deployment)) +
  geom_point() +
  xlab('Sequence') +
  ylab("Mean trap speed per sequence (cm/s)") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# 4. Fig examples of oscillation ----

## 4.1. Mean oscillation

oscillations <- depths %>% 
  dplyr::summarize(mean_osc = mean(abs(ddepth), na.rm = T),
                   sd_osc = sd(abs(ddepth), na.rm = T))

amplitudes <- depths %>% group_by(deployment, sequence) %>%
  dplyr::summarize(
    amplitude = max(averaged_depth, na.rm = T) - min(averaged_depth, na.rm = T))
mean(amplitudes$amplitude)
sd(amplitudes$amplitude)
max(amplitudes$amplitude)
min(amplitudes$amplitude)

## 4.2. Fig examples of oscillation

p1 <- depths %>% 
  filter(deployment == 'c1_500', sequence == '20220909-200000') %>%
  ggplot(aes(x = datetime, y = averaged_depth)) +
  geom_point(size = 1, alpha = 0.7) +
  xlab('Time (UTC)') +
  ylab('Trap depth (m)') +
  ggtitle('C1 - 500 m - Sequence 20220909-200000') +
  theme_bw()

p2 <- depths %>% 
  filter(deployment == 'c2_200', sequence == '20220913-160000') %>%
  ggplot(aes(x = datetime, y = averaged_depth)) +
  geom_point(size = 1, alpha = 0.7) +
  xlab('Time (UTC)') +
  ylab('Trap depth (m)') +
  ggtitle('C2 - 200 m - Sequence 20220913-160000') +
  theme_bw()
  
(poscill <- plot_grid(p1, p2, labels = c('(a)', '(b)'), ncol = 1))
ggsave("figures\\trap_oscillations.pdf",
       poscill, width = 10, height = 5)


# 5. Correlation between trap speed and track speed ---- 

## 5.1. On sequence-binned data (absolute)

tracks %>% group_by(deployment, sequence) %>% 
  dplyr::summarise(mean_track_speed = mean(abs(track_speed))) %>%
  full_join(
    (depths %>% group_by(deployment, sequence) %>%
       dplyr::summarise(mean_trap_speed = mean(abs(speed), na.rm = T)))) %>%
  ggplot(aes(x = mean_trap_speed, y = mean_track_speed, color = deployment)) +
  geom_point() +
  geom_smooth() +
  xlab("Mean absolute trap speed per sequence (cm/s)") +
  ylab("Mean absolute track speed per sequence (m/d)") +
  facet_wrap(~deployment, scales = "free")

tracks %>% group_by(deployment, sequence) %>% 
  dplyr::summarise(mean_track_speed = mean(abs(track_speed))) %>%
  full_join((depths %>% group_by(deployment, sequence) %>%
               dplyr::summarise(
                 sd_trap_depth = sd(abs(averaged_depth), na.rm = T)))) %>%
  ggplot(aes(x = sd_trap_depth, y = mean_track_speed, color = deployment)) +
  geom_point() +
  geom_smooth() +
  xlab("Sd trap depth per sequence (m)") +
  ylab("Mean absolute track speed per sequence (m/d)") +
  facet_wrap(~deployment, scales = "free")

tracks %>% group_by(deployment, sequence) %>% 
  dplyr::summarise(mean_track_speed = mean(abs(track_speed))) %>%
  full_join((depths %>% group_by(deployment, sequence) %>%
               dplyr::summarise(mean_trap_speed = mean(abs(speed), na.rm = T)))) %>%
  ddply(., .(deployment), summarise,
        corr=(cor.test(mean_trap_speed, mean_track_speed,
                       method = "kendall")), name=names(corr))

## 5.2. Speed of the trap during the course of the track

(p3 <- ggplot(tracks, aes(x = trap_speed, y = track_speed)) +
  geom_point() +
  xlab("Mean trap speed during the course of the track (cm/s)") +
  ylab("Mean track speed (m/d)") +
  theme_bw() +
  facet_wrap(~deployment, scales = "free"))

tracks %>%
  ddply(., .(deployment), summarise,
        corr=(cor.test(trap_speed, track_speed,
                       method = "kendall", exact = FALSE)), name=names(corr))

## 5.3 Speed of the trap 30 s before the track

df30 <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df30) <- c('deployment', 'datetime_ini', 'mean_track_speed', 
                    'mean_trap_speed', 'sd_depth')

for (i in 1:dim(tracks)[1]){
  dt2 <- ymd_hms(tracks$datetime_ini[i], tz = "UTC")
  dt1 <- dt2 - seconds(30)
  sub <- depths[depths$datetime>=dt1 & depths$datetime<=dt2,]
  mean_speed <- mean(abs(sub$speed))
  sd_depth <- sd(sub$averaged_depth)
  
  df30[nrow(df30) + 1,] = c(tracks$deployment[i], dt2, 
                            tracks$track_speed[i],
                            mean_speed, sd_depth)
}

(p4 <- ggplot(df30, 
       aes(x = as.numeric(mean_trap_speed), 
           y = as.numeric(mean_track_speed))) +
  geom_point() +
  xlab("Mean absolute trap speed 30 s before the track (cm/s)") +
  ylab("Mean track speed (m/d)") +
  theme_bw() +
  facet_wrap(~deployment, scales = "free"))

df30 %>%
  mutate(mean_trap_speed = as.numeric(mean_trap_speed),
         mean_track_speed = as.numeric(mean_track_speed)) %>%
  ddply(., .(deployment), summarise,
      corr=(cor.test(mean_trap_speed, mean_track_speed,
                     method="kendall", exact = FALSE)), name=names(corr))

(prelations <- plot_grid(p3, p4, labels =  c('(a)', '(b)'), ncol = 1))
ggsave("figures\\trap_speeds_relationships.pdf",
       prelations, width = 10, height = 5)

