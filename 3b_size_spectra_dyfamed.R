#------------------------------
# UVP6 MOOSE-GE track data
#------------------------------
# Size spectra comparison
# 19-04-2023
# Manon Laget
#------------------------------


library(tidyverse)
library(lubridate)
library(nbssr)
library(lsmeans)
library(latex2exp)

path <- "C:\\Users\\Manon Laget\\Desktop\\UVP6_data\\projects\\paper_methodo\\analyses\\"
setwd(path)


# 1. UVP6 data -- VT (only data acquired at low frequency are loaded) and stock ----

uvp6vt <- read.csv('data\\particle_data\\export_detailed_dyfamed.tsv', sep = '\t',
                   fileEncoding = 'latin1', check.names = F) %>%
  # Remove columns we don't need
  select(-(`Sampled volume [L]`:`LPM (64-128 µm) [# l-1]`)) %>%
  # Remove size classes > 2 mm and columns beyond
  select(-(`LPM (2.05-4.1 mm) [# l-1]`:`extrames20`)) %>%
  group_by(Profile) %>%
  select(Profile, `Depth [m]`, starts_with("LPM")) %>%
  summarize_if(is.numeric, mean) %>%
  filter(Profile != 'dyfamed_visutrap_01')

classes <- substr(colnames(uvp6vt)[-c(1:2)], 6, 15)
classes <- str_split_fixed(classes, " ", 2)[,1]
classes <- data.frame(str_split_fixed(classes, "-", 2))
classes <- (as.numeric(classes$X1) + as.numeric(classes$X2)) /2
classes[3:4] <- classes[3:4] * 1000


# 2. Comparison of size spectra for the duration of deployments ----

df_nbss <- bind_rows(
  (nbss(x = classes, w = unlist(uvp6vt[1,-c(1,2)]), type = "abundance") %>%
     mutate(uvp = 'uvp6_vt', deployment = 'February')),
  (nbss(x = classes, w = unlist(uvp6vt[2,-c(1,2)]), type = "abundance") %>%
     mutate(uvp = 'uvp6_vt', deployment = 'April')),
  (nbss(x = classes, w = unlist(uvp6vt[3,-c(1,2)]), type = "abundance") %>%
     mutate(uvp = 'uvp6_stock', deployment = 'February')),
  (nbss(x = classes, w = unlist(uvp6vt[4,-c(1,2)]), type = "abundance") %>%
     mutate(uvp = 'uvp6_stock', deployment = 'April')))

(p <- df_nbss %>% 
    ggplot(aes( x = bin, y = norm_y, color = uvp)) +
    geom_point() +
    geom_path() +
    scale_y_log10() +
    scale_x_log10(minor_breaks = NULL) +
    annotation_logticks(sides="bl") +
    scale_color_discrete('UVP') +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)) +
    xlab("ESD (µm)") +
    ylab(TeX("Normalized abundance (µm µm$^{-1}$ m$^{-3}$)")) +
  facet_wrap(~factor(deployment, c('February', 'April'))))

ggsave('figures\\nass_dyfamed.pdf', p, height = 4, width = 12)


# 3. Test for differences in slope with ANOVA

mod2 <- df_nbss %>% filter(deployment == 'February') %>% 
  lm(log10(norm_y) ~ log10(bin) * uvp, data = .)
anova(mod2)

mod2 <- df_nbss %>% filter(deployment == 'April', norm_y != 0) %>% 
  lm(log10(norm_y) ~ log10(bin) * uvp, data = .)
anova(mod2)

