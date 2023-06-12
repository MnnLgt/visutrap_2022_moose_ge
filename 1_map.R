#------------------------------
# UVP6 MOOSE-GE track data
#------------------------------
# Map of deployments
# 23-01-2023
# Manon Laget
#------------------------------


# Load packages
library(readxl)
library(marmap)
library(oce)

path <- "C:\\Users\\Manon Laget\\Desktop\\UVP6_data\\projects\\paper_methodo\\analyses\\"
setwd(path)


# Load sediment trap line coordinates
coords_lion <- read_excel("data\\trap_data\\MOOSE_GE_2022_LION_iridium.xls")
coords_dyf <- read_excel("data\\trap_data\\MOOSE_GE_2022_DYF_iridium.xls")

# Coordinates of LION and DYFAMED sites
special_points = data.frame(name = c("LION", "DYFAMED"),
                            lon = c(4.716780235, 7.888608),
                            lat = c(41.984222, 43.421357))

# Get bathymetry data from NOAA
med <- getNOAA.bathy(lon1 = -3, lon2 = 15,
                        lat1 = 35, lat2 = 47, resolution = 1)

# Creating color palettes
greys <- c(grey(0.6), grey(0.93), grey(0.99))
col <- oceColorsGebco(120)[1:100]

png("figures\\map.png", width = 10, height = 6, unit = 'in', res = 480)  

# Drawing palette
drawPalette(zlim = c(min(med), 0),
            zlab = "Depth (m)",
            col = col, 
            at = seq(-3000, 0, 1000), 
            cex = 0.8,
            pos = 4, 
            las = 0,
            drawTriangles = FALSE)

# Drawing map
plot(med, image = TRUE, land = TRUE, xlim=c(0, 12), ylim=c(38.5, 45.5),
     xlab = "Longitude (°)",
     ylab = "Latitude (°)",
     bpal = list(c(0, max(med), greys), c(min(med), 0 , col)),
     deep = c(-3000, -1000, -200),
     shallow = c(-3000, -1000, -200),
     step = c(0, 0, 0),
     lwd = c(0.1, 0.1, 0.1), lty = c(1, 1, 1),
     col = c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2"),
     drawlabel = c(T, T, T))
axis(side = 1, c(-2, 14), labels = F, lwd.ticks = 0)
axis(side = 2, c(38, 46), labels = F, lwd.ticks = 0)
axis(side = 3, c(-2, 14), labels = F, lwd.ticks = 0)
axis(side = 4, c(38, 46), labels = F, lwd.ticks = 0)

# Add sediment trap tracks

# LION
lines(coords_lion$Longitude, coords_lion$Latitude, lwd = 3)
text(special_points[1, 'lon'] - 0.3, y = special_points[1,'lat'] - 0.3, 
     labels = special_points[1,'name'])

# DYFAMED
lines(coords_dyf$Longitude, coords_dyf$Latitude, lwd = 3)
text(special_points[2, 'lon'] - 0.3, y = special_points[2,'lat'] - 0.3, 
     labels = special_points[2,'name'])

dev.off()
