## Water PCB concentrations mapping.
# Data were obtained from EPA and contractors from PCB Superfund
# sites in USA

# Install packages
install.packages("ggplot2")
install.packages("devtools")
install.packages("dplyr")
install.packages("stringr")
install.packages("maps")
install.packages("mapdata")
install.packages("ggmap")
install.packages("usethis")
install.packages("GISTools")
install.packages("rgeos")
install.packages("ggsn")
install.packages("ggrepel")
install.packages("ggpp")

# Load libraries
library(dplyr)
library(usethis)
library(devtools)
library(ggmap) # function map_data
library(maps)
library(ggplot2)
library(leaflet)
library(raster)
library(GISTools)
library(rgeos)
library(ggsn)
library(ggrepel)
library(reshape2)
library(ggpmisc)

# Data in pg/L
wdc.0 <- read.csv("WaterDataCongenerAroclor062322.csv")

# Map locations -----------------------------------------------------------
us <- map_data("usa")
states <- map_data("state")

# Find number of samples per state to be included as table in maps
wdc.1 <- wdc.0 %>%
  group_by(StateSampled) %>%
  summarise(n = n())
colnames(wdc.1) <- c("State", "# samples")

# Creates map of US with locations from congener data
ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = "lightblue") +
  coord_fixed(1.3) +
  theme_nothing() +
  geom_path(data = states, aes(x = long, y = lat, group = group),
             colour = "white") +
  geom_polygon(color = "black", fill = NA) +
  geom_point(data = wdc.0, aes(x = Longitude, y = Latitude), color = "black",
             size = 1.2, shape = 20) +
  annotate(geom = 'table', x = -57, y = 32,
           label = list(wdc.1), size = 2.5) # add table with info

# Map total PCB concentrations --------------------------------------------
# Determine a time range (e.g., 2010 - 2019) and few locations
# Prepare data
# (i) Remove samples (rows) with total PCBs  = 0
wdc.2 <- wdc.0[!(rowSums(wdc.0[, c(12:115)], na.rm = TRUE)==0),]
# Get tPCB and coordinates
wdc.tPCB <- data.frame(cbind(wdc.2$Latitude, wdc.2$Longitude,
                            rowSums(wdc.2[, c(12:115)], na.rm = TRUE)))
# Name the columns
colnames(wdc.tPCB) <- c("Latitude", "Longitude", "tPCB")
# Average tPCB per site
wdc.tPCB.mean <- aggregate(tPCB ~ Latitude + Longitude, data = wdc.tPCB, mean)

# Plot map + tPCB
# Cannot include legend
ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = "lightblue") +
  coord_fixed(1.3) +
  theme_nothing() +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "white") +
  geom_polygon(color = "black", fill = NA) +
  geom_point(data = wdc.tPCB.mean, aes(x = Longitude, y = Latitude,
                                       size = tPCB), color = "red")

# Portland maps
# Select only OR
wdc.OR <- subset(wdc.0, wdc.0$StateSampled == "OR")
# Select only from Portland
wdc.PO <- subset(wdc.OR, SiteSampled!="NA")

# Create map
PO.box <- make_bbox(lon = wdc.PO$Longitude, lat = wdc.PO$Latitude, f = 0.8)
PO.map <- get_stamenmap(bbox = PO.box, zoom = 10)

# Plot map + sampling locations
ggmap(PO.map) +
  geom_point(data = wdc.PO, aes(x = Longitude, y = Latitude), shape = 21,
             color = "red",
             fill = "white", size = 1.75, stroke = 0.75)

# Plot map + mean tPCB
# Remove samples (rows) with total PCBs  = 0
wdc.PO.1 <- wdc.PO[!(rowSums(wdc.PO[, c(12:115)],
                                   na.rm = TRUE)==0),] # sum of PCB1 to PCB209
# Get tPCB and coordinates
tPCB.PO <- data.frame(cbind(wdc.PO.1$Latitude, wdc.PO.1$Longitude,
                         rowSums(wdc.PO.1[, c(12:115)], na.rm = TRUE)))
# Name the columns
colnames(tPCB.PO) <- c("Latitude", "Longitude", "tPCB")
# Average tPCB per site
tPCB.mean <- aggregate(tPCB ~ Latitude + Longitude, data = tPCB.PO, mean)

# (3) Plot map + tPCB
ggmap(PO.map) +
  geom_point(data = tPCB.mean, aes(x = Longitude, y = Latitude,
                                   size = tPCB), alpha = 1, color  = "red") +
  scale_size_area(breaks = c(100, 125, 150, 175, 200),
                  name = "Ave. PCBs \n2018-2019 (pg/L)") +
  xlab("Longitude") +
  ylab("Latitude")
 


# Map with ggmap
WI.box <- make_bbox(lon = w.WI$Long, lat = w.WI$Lat, f = 0.6)
wi.map <- get_stamenmap(bbox = WI.box, zoom = 11)


# Remove samples (rows) with total PCBs  = 0
w.WI.t <- w.WI[!(rowSums(w.WI[,
                              c(12:115)],
                         na.rm = TRUE)==0),] # sum of PCB1 to PCB209
site.sampled <- w.WI.t$SiteSampled
w.WI.t <- subset(w.WI.t, select = -c(ID:AroclorCongener))
w.WI.t <- subset(w.WI.t, select = -c(AroclorA1016:AroclorA1260))
# Get mean congener per site, excluding zeros
tPCB <- rowSums(w.WI.t, na.rm = TRUE)
tPCB <- data.frame(cbind(site.sampled, tPCB))
tPCB$tPCB <- as.numeric(as.character(tPCB$tPCB))
tPCB.mean <- aggregate(tPCB ~ site.sampled, data = tPCB, mean)
# add coordinates
tPCB.mean <- data.frame(c(tPCB.mean, wi.coord))

# (3) Plot map + tPCB
ggmap(wi.map) +
  geom_point(data = tPCB.mean, aes(x = Long, y = Lat,
                                   size = tPCB), alpha = 0.5) +
  scale_size_area(breaks = c(250, 500, 750, 1000, 1500),
                  labels = c(250, 500, 750, 1000, 1500),
                  name = "PCBs ng/L") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(title = "Fox River PCBs water concentration (pg/L) 2010-2018")



# select only WI
w.WI <- subset(w, w$StateSampled == "WI")

# Map with ggmap
WI.box <- make_bbox(lon = w.WI$Long, lat = w.WI$Lat, f = 0.6)
wi.map <- get_stamenmap(bbox = WI.box, zoom = 11)

# (1) Plot map
ggmap(wi.map)

# (2) Plot map + sampling locations
ggmap(wi.map) +
  geom_point(data = w.WI, aes(x = Long, y = Lat), shape = 21,
  color = "red",
  fill = "white", size = 1.75, stroke = 0.75)

# Prepare congener data for plotting
# Get coordinates per site
LW <- subset(w.WI, w.WI$SiteSampled == 'LakeWinnebago')
LW <- data.frame(c(LW[1,6], LW[1,7]))
OU1 <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit1')
OU1 <- data.frame(c(OU1[1,6], OU1[1,7]))
OU2A <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit2A')
OU2A <- data.frame(c(OU2A[1,6], OU2A[1,7]))
OU2B <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit2B')
OU2B <- data.frame(c(OU2B[1,6], OU2B[1,7]))
OU2C <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit2C')
OU2C <- data.frame(c(OU2C[1,6], OU2C[1,7]))
OU3 <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit3')
OU3 <- data.frame(c(OU3[1,6], OU3[1,7]))
wi.coord <- rbind(LW, OU1, OU2A, OU2B, OU2C, OU3)

# Total PCBs
# # remove samples (rows) with total PCBs  = 0
w.WI.t <- w.WI[!(rowSums(w.WI[,
                           c(12:115)],
                         na.rm = TRUE)==0),] # sum of PCB1 to PCB209
site.sampled <- w.WI.t$SiteSampled
w.WI.t <- subset(w.WI.t, select = -c(ID:AroclorCongener))
w.WI.t <- subset(w.WI.t, select = -c(AroclorA1016:AroclorA1260))
# Get mean congener per site, excluding zeros
tPCB <- rowSums(w.WI.t, na.rm = TRUE)
tPCB <- data.frame(cbind(site.sampled, tPCB))
tPCB$tPCB <- as.numeric(as.character(tPCB$tPCB))
tPCB.mean <- aggregate(tPCB ~ site.sampled, data = tPCB, mean)
# add coordinates
tPCB.mean <- data.frame(c(tPCB.mean, wi.coord))

# (3) Plot map + tPCB
ggmap(wi.map) +
  geom_point(data = tPCB.mean, aes(x = Long, y = Lat,
                              size = tPCB), alpha = 0.5) +
  scale_size_area(breaks = c(250, 500, 750, 1000, 1500),
                  labels = c(250, 500, 750, 1000, 1500),
                  name = "PCBs ng/L") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(title = "Fox River PCBs water concentration (pg/L) 2010-2018")

# Congener maps
# Select congener and remove samples with = 0 and NA for selected congener
w.WI.2 <- subset(w.WI, w.WI$PCB1 != 0 & w.WI$PCB1 != "NA")
# Get mean congener per site, excluding zeros
PCB1 <- aggregate(PCB1 ~ SiteSampled, data = w.WI.2, mean)
PCB1 <- data.frame(c(PCB1, wi.coord))

# (4) Plot map + congener
ggmap(wi.map) +
  geom_point(data = PCB1, aes(x = Long, y = Lat,
                              size = PCB1), alpha = 0.5) +
  scale_size_area(breaks = c(0.1, 1, 2, 4, 6),
                 labels = c(0.1, 1, 2, 4, 6),
                 name = "PCB 1 ng/L") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(title = "Fox River PCB 1 water concentration (pg/L) 2010-2018")
  #geom_label_repel(data = PCB1, aes(x = Long, y = Lat, label = SiteSampled),
  #                 fill = "white", box.padding = unit(0.3, "lines"),
  #                 label.padding = unit(0.15, "lines"),
  #                 segment.color = "black", segment.size = 1)
                   
                   
