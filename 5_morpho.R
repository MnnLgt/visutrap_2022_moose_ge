#------------------------------
# UVP6 MOOSE-GE track data
#------------------------------
# Particle morphological space
# 14-02-2023
# Manon Laget
#------------------------------


# Load packages
library(tidyverse)
library(patchwork)
library(morphr)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(bestNormalize)
library(languageR)
library(cowplot)
library(latex2exp)
library(ggpubr)
library(rstatix)


# Load tracks ----

# Set working directory
path <- "C:\\Users\\Manon Laget\\Desktop\\UVP6_data\\projects\\paper_methodo\\analyses\\"
setwd(path)

# Read track file
tracks <- read_tsv('results\\track_data.tsv')


# 1. Morphological space of particles ----

# Using morphological descriptors, particles are projected in a low-dimensional
# space. This is performed with the *morphor* package.

## 1.1. Create morphological space

# We use *morphospace* function to perform a PCA on selected descriptors. Then,
# some images of the dataset are displayed at their position in the
# morphological space. Here we show PC1 vs PC2 and PC3 vs PC4.

# Select features

X <- tracks %>%
  mutate(fractal = 2*log(perim_px)/log(area_px)) %>%
  dplyr::select(track_id, area_convex, area_um, circularity,
                elongation, extent, fractal, major_px, 
                meangrey, mediangrey, minor_px, perim_px, 
                perimmajor, rangegrey, skewgrey, solidity, stdgrey)

row.names(X) <- X$track_id
X <- X %>% dplyr::select(-track_id)
pairscor.fnc(X)

# Transform data
X <- sapply(X,\(x)yeojohnson(x)$x.t)
w <- rep(1, dim(X)[1])

# Perform PCA with morphospace function
morpho <- morphospace(X, weights = w)

# Plot images in the morphological space and variable contribution to PCs

vig <- tracks$path_to_vig
cols <- c("#fc8d62", "#8da0cb", "#66c2a5")

(p1 <- ggmorph_tile(morpho, vig, dimensions = c(1, 2),
                    steps = 13, n_imgs = 12, scale = 0.007))

(p2 <- ggmorph_tile(morpho, vig, dimensions = c(2, 3),
                    steps = 13, n_imgs = 12, scale = 0.007))

(p3 <- fviz_pca_biplot(morpho, axes = c(1, 2), repel = T, label = "var", 
                       pointsize = 1, labelsize = 4,
                       habillage = as.factor(tracks$type),
                       addEllipses = TRUE, ellipse.level = 0.5, title = "",
                       palette = cols,
                       col.var = gray(0.2), legend.title = "Type"))

(p4 <- fviz_pca_biplot(morpho, axes = c(2, 3), repel = T, label = "var", 
                       pointsize = 1, labelsize = 4,
                       habillage = as.factor(tracks$type),
                       addEllipses = TRUE, ellipse.level = 0.5, title = "",
                       palette = cols,
                       col.var = gray(0.2), legend.title = "Type"))

(p <- plot_grid(p1, p2, p3, p4, 
                labels = c('(a)', '(b)', '(c)', '(d)'), ncol = 2))
ggsave('figures\\pca_ind_var.pdf', p, height = 10, width = 12)
ggsave('figures\\pca_ind_var.png', p, height = 10, width = 12)

## 1.2. Variable contribution to PCs

pcont1 <- fviz_contrib(morpho, choice = "var", axes = 1, top = 10) +
  theme(axis.text = element_text(size = 14))
pcont2 <- fviz_contrib(morpho, choice = "var", axes = 2, top = 10) +
  theme(axis.text = element_text(size = 14))
pcont3 <- fviz_contrib(morpho, choice = "var", axes = 3, top = 10) +
  theme(axis.text = element_text(size = 14))

(pcont <- plot_grid(pcont1, pcont2, pcont3, 
                    labels = c('(a)', '(b)', '(c)'), ncol = 3))
ggsave('figures\\pca_var_contributions.pdf', pcont, height = 4, width = 12)


# 2. Relationships PCs/track features ----

scores_ind <- as.data.frame(morpho$ind$coord) %>%
  bind_cols(tracks$type) %>%
  bind_cols(tracks$vertical_speed) %>%
  select(type = 6, vertical_speed = 7, everything())

scores_ind$type <- factor(scores_ind$type, 
                          levels = c('ascending', 'suspended', 'sinking'))

## 2.1. Differences between groups

aov_pc1 <- aov(Dim.1 ~ type, data = scores_ind)
summary(aov_pc1)
TukeyHSD(aov_pc1)

aov_pc2 <- aov(Dim.2 ~ type, data = scores_ind)
summary(aov_pc2)
TukeyHSD(aov_pc2)

aov_pc3 <- aov(Dim.3 ~ type, data = scores_ind)
summary(aov_pc3)
TukeyHSD(aov_pc3)

## 2.2. Correlation between vertical speed and PCs

c1 <- cor.test(scores_ind[,1], scores_ind$vertical_speed)
c2 <- cor.test(scores_ind[,2], scores_ind$vertical_speed)
c3 <- cor.test(scores_ind[,3], scores_ind$vertical_speed)
p.adjust(c(c1$p.value, c2$p.value, c3$p.value), method = "holm")

## 2.3. Figs

# Distribution of PCs according to track type

stat.test <- tukey_hsd(scores_ind, formula = Dim.1 ~ type) %>% 
  add_xy_position(x = "type")
(bp1 <- ggplot(scores_ind, aes(x = type, y = Dim.1)) + 
    geom_boxplot(aes(fill = type, color = type), alpha = 0.6, outlier.alpha =  0) +
    stat_pvalue_manual(stat.test, label="p.adj.signif", y.position = c(7.2, 8.4, 9.6),
                       size = 4) +
    scale_fill_manual(values = cols, 
                      limits = c("ascending", "suspended", "sinking"),
                      labels = c("ascending", "suspended", "sinking")) +
    scale_color_manual('Track type', values = cols, 
                       limits = c("ascending", "suspended", "sinking"),
                       labels = c("ascending", "suspended", "sinking")) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      legend.position="none") +
    xlab("") +
    ylab(TeX("PC1 (particle density)")))

stat.test <- tukey_hsd(scores_ind, formula = Dim.2 ~ type) %>% 
  add_xy_position(x = "type")
(bp2 <- ggplot(scores_ind, aes(x = type, y = Dim.2)) + 
    geom_boxplot(aes(fill = type, color = type), alpha = 0.6, outlier.alpha =  0) +
    stat_pvalue_manual(stat.test, label="p.adj.signif", y.position = c(7.2, 8.4, 9.6),
                       size = 4) +
    scale_fill_manual(values = cols, 
                      c("suspended", "sinking", "ascending")) +
    scale_color_manual(values = cols, 
                       c("suspended", "sinking", "ascending")) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      legend.position="none") +
    xlab("") +
    ylab(TeX("PC2 (particle size)")))

stat.test <- tukey_hsd(scores_ind, formula = Dim.3 ~ type) %>% 
  add_xy_position(x = "type")
(bp3 <- ggplot(scores_ind, aes(x = type, y = Dim.3)) + 
    geom_boxplot(aes(fill = type, color = type), alpha = 0.6, outlier.alpha =  0) +
    scale_fill_manual(values = cols, 
                      c("suspended", "sinking", "ascending")) +
    scale_color_manual(values = cols, 
                       c("suspended", "sinking", "ascending")) +
    stat_pvalue_manual(stat.test, label="p.adj.signif", y.position = c(4.2, 5.4, 6.6),
                       size = 4) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      legend.position="none") +
    xlab("") +
    ylab(TeX("PC3 (particle roundness)")))

(allbp <- plot_grid(bp1, bp2, bp3,
                    labels = c('(a)', '(b)', '(c)'), 
                    ncol = 3, hjust = 0))
ggsave('figures\\pca_boxplots.pdf', allbp, height = 4, width = 12)
ggsave('figures\\pca_boxplots.png', allbp, height = 4, width = 12)


