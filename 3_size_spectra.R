#------------------------------
# UVP6 MOOSE-GE track data
#------------------------------
# Size spectra comparison
# 08-02-2023
# Manon Laget
#------------------------------


library(tidyverse)
library(lubridate)
library(nbssr)
library(lsmeans)
library(latex2exp)

path <- "C:\\Users\\Manon Laget\\Desktop\\UVP6_data\\projects\\paper_methodo\\analyses\\"
setwd(path)


# 1. Loading UVP5 data ----

# Select stations associated to cycles
casts_lion <- c('008', '009', '010', '011', '012', '013', '014', '015', '016',
                '017', '018', '019')
casts_dyf <- c('028', '029', '030', '031', '032', '033', '034', '035', '036', 
               '037')

# UVP5 file with all stations
uvp5 <- read.csv('data\\particle_data\\export_detailed_uvp5.tsv', sep = '\t',
                 fileEncoding = 'latin1', check.names = F) %>%
  # Select leg 1
  filter(substr(Profile, 1, 4) == "leg1") %>%
  # Remove columns we don't need
  # Seems to have a problem with class 128-161... we'll start from here
  select(-(`Sampled volume [L]`:`LPM (128-161 µm) [# l-1]`)) %>%
  # Remove size classes > 2 mm and columns beyond
  select(-(`LPM (2.05-2.58 mm) [# l-1]`:`extrames20`)) %>%
  filter(substr(Profile, 7, 9) %in% casts_lion |
           substr(Profile, 7, 9) %in% casts_dyf)

max(uvp5$`Depth [m]`)

# Create a vector of values representing size classes
classes <- substr(colnames(uvp5)[-c(1:5)], 6, 15)
classes <- str_split_fixed(classes, " ", 2)[,1]
classes <- data.frame(str_split_fixed(classes, "-", 2))
classes <- (as.numeric(classes$X1) + as.numeric(classes$X2)) /2
classes[8:11] <- classes[8:11] * 1000
  

# 2. VisuTrap data -- only data acquired at low frequency are loaded ----

uvp6vt <- read.csv('data\\particle_data\\export_detailed_uvp6_vt.tsv', sep = '\t',
                 fileEncoding = 'latin1', check.names = F) %>%
  # Remove columns we don't need
  select(-(`Sampled volume [L]`:`LPM (128-161 µm) [# l-1]`)) %>%
  # Remove size classes >2 mm and columns beyond
  select(-(`LPM (2.05-2.58 mm) [# l-1]`:`extrames20`)) %>%
  group_by(Profile) %>%
  select(Profile, `Depth [m]`, starts_with("LPM")) %>%
  summarize_if(is.numeric, mean)


# 3. Comparison of size spectra for the duration of deployments ----

# Create NBSS for each UVP6 deployment

df_nbss <- bind_rows(
  (nbss(x = classes, w = unlist(uvp6vt[1,-c(1,2)]), type = "abundance") %>%
     mutate(uvp = 'uvp6_vt', deployment = 'c1_500')),
  (nbss(x = classes, w = unlist(uvp6vt[2,-c(1,2)]), type = "abundance") %>%
     mutate(uvp = 'uvp6_vt', deployment = 'c2_200')),
  (nbss(x = classes, w = unlist(uvp6vt[3,-c(1,2)]), type = "abundance") %>%
     mutate(uvp = 'uvp6_vt', deployment = 'c2_500')))

# Create NBSS for each UVP5 deployment

nbss_tmp <- uvp5 %>%
  # Take UVP5 casts corresponding to cycle LION (c1)
  filter(substr(Profile, 7, 9) %in% casts_lion) %>%
  # Select depths +/- 10 m of trap depth
  filter(`Depth [m]` > mean(
    (uvp6vt %>% filter(Profile == "c1_500_all"))$`Depth [m]`) - 10,
    `Depth [m]` < mean(
      (uvp6vt %>% filter(Profile == "c1_500_all"))$`Depth [m]`) + 10) %>%
  # Group by datetime and mean
  select(starts_with("LPM")) %>%
  summarize_all(mean)

df_nbss <- bind_rows(
  df_nbss,
  (nbss(x = classes, w = unlist(nbss_tmp[1,]), type = "abundance") %>%
     mutate(uvp = 'uvp5', deployment = 'c1_500')))

nbss_tmp <- uvp5 %>%
  # Take UVP5 casts corresponding to cycle DYF (c2)
  filter(substr(Profile, 7, 9) %in% casts_dyf) %>%
  # Select depths +/- 10 m of trap depth
  filter(`Depth [m]` > mean(
    (uvp6vt %>% filter(Profile == "c2_500_all"))$`Depth [m]`) - 10,
    `Depth [m]` < mean(
      (uvp6vt %>% filter(Profile == "c2_500_all"))$`Depth [m]`) + 10) %>%
  # Group by datetime and mean
  select(starts_with("LPM")) %>%
  summarize_all(mean)

df_nbss <- bind_rows(
  df_nbss,
  (nbss(x = classes, w = unlist(nbss_tmp[1,]), type = "abundance") %>%
     mutate(uvp = 'uvp5', deployment = 'c2_500')))

nbss_tmp <- uvp5 %>%
  # Take UVP5 casts corresponding to cycle DYF (c2)
  filter(substr(Profile, 7, 9) %in% casts_dyf) %>%
  # Select depths +/- 10 m of trap depth
  filter(`Depth [m]` > mean(
    (uvp6vt %>% filter(Profile == "c2_200_all"))$`Depth [m]`) - 10,
    `Depth [m]` < mean(
      (uvp6vt %>% filter(Profile == "c2_200_all"))$`Depth [m]`) + 10) %>%
  # Group by datetime and mean
  select(starts_with("LPM")) %>%
  summarize_all(mean)

df_nbss <- bind_rows(
  df_nbss,
  (nbss(x = classes, w = unlist(nbss_tmp[1,]), type = "abundance") %>%
     mutate(uvp = 'uvp5', deployment = 'c2_200'))) %>%
  mutate(norm_y = norm_y*1000)


# 4. Plot NBSS ---- 

(p <- ggplot(df_nbss, aes(
  x = bin, y = norm_y, color = uvp)) +
  geom_point() +
  geom_path() +
  scale_y_log10() +
  scale_x_log10(minor_breaks = NULL) +
  annotation_logticks(sides="bl") +
  scale_color_discrete('UVP') +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)) +
  xlab("ESD (µm)") +
  ylab(TeX("Normalized abundance (µm µm$^{-1}$ m$^{-3}$)")) +
   facet_wrap(~deployment))

ggsave('figures\\nass.pdf', p, height = 4, width = 12)


# 5. Test for differences in slope with ANOVA

mod2 <- df_nbss %>% filter(deployment == 'c1_500') %>% 
  lm(log10(norm_y) ~ log10(bin) * uvp, data = .)
anova(mod2)

mod2 <- df_nbss %>% filter(deployment == 'c2_200', norm_y != 0) %>% 
  lm(log10(norm_y) ~ log10(bin) * uvp, data = .)
anova(mod2)

mod2 <- df_nbss %>% filter(deployment == 'c2_500', norm_y != 0) %>% 
  lm(log10(norm_y) ~ log10(bin) * uvp, data = .)
anova(mod2)


# 6. Plot concentrations ----

uvp5 %>% filter(`Depth [m]` < 2000) %>%
  mutate(deployment = case_when(
    substr(Profile, 7, 9) %in% casts_dyf ~ 'DYFAMED',
    substr(Profile, 7, 9) %in% casts_lion ~ 'LION')) %>%
  group_by(deployment, `Depth [m]`) %>%
  summarize_if(is.numeric, median) %>%
  ungroup() %>%
  mutate(conc_tot = rowSums(.[3:13])) %>% 
  ggplot(aes(x = `Depth [m]`, y = conc_tot, color = deployment)) +
  geom_point() +
  geom_path()
  

