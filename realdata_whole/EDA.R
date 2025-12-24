# ===================================
# 1. Load libraries
# ===================================
rm(list = ls())
library(ggmap)
library(sf)
library(ggplot2)
library(viridis)
library(dplyr)
library(maps)
library(cowplot)
library(ggrepel)

library(moments)
library(nortest)  


# ===================================
# 2. Stadia Maps API
# ===================================
register_stadiamaps("868f9432-846b-4444-9e4e-da68b56ab6b2")

# ===================================
# 3. Disable S2
# ===================================
sf_use_s2(FALSE)

# ===================================
# 4. Load data
# ===================================
frm  <- read.csv("Gold_clean.csv")
cmaq <- read.csv("CMAQ_clean.csv")

names(frm)  <- c("lat", "lon", "AQI")
names(cmaq) <- c("lat", "lon", "PM25")

frm <- frm %>% filter(AQI <= 70)

frm_sf  <- st_as_sf(frm,  coords = c("lon","lat"), crs = 4326)
cmaq_sf <- st_as_sf(cmaq, coords = c("lon","lat"), crs = 4326)

# ===================================
# 5. US boundary polygon
# ===================================
usa_states_raw <- map("state", plot = FALSE, fill = TRUE)
usa_states <- st_as_sf(usa_states_raw)
st_crs(usa_states) <- 4326
usa_states <- st_make_valid(usa_states)
usa_union  <- st_union(usa_states)

# ===================================
# 6. Remove ocean points
# ===================================
frm_land_sf  <- frm_sf[lengths(st_intersects(frm_sf,  usa_states)) > 0, ]
cmaq_land_sf <- cmaq_sf[lengths(st_intersects(cmaq_sf, usa_states)) > 0, ]

frm_land <- frm_land_sf %>%
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2]) %>% st_drop_geometry()

cmaq_land <- cmaq_land_sf %>%
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2]) %>% st_drop_geometry()

# ===================================
# 7. Background map
# ===================================
bbox_vec <- c(left=-125, bottom=24, right=-66, top=50)
map_bg <- get_stadiamap(bbox_vec, zoom=5, maptype="stamen_terrain_background")

# ===================================
# 8. Cities labels
# ===================================
cities <- data.frame(
  name = c("New York","Washington DC","Chicago",
           "Los Angeles","San Francisco",
           "Houston","Dallas","Miami",
           "Atlanta","Phoenix","Seattle","Boston"),
  lon = c(-74.0060,-77.0369,-87.6298,-118.2437,-122.4194,
          -95.3698,-96.7970,-80.1918,-84.5880,-112.0740,
          -122.3321,-71.0589),
  lat = c(40.4128,38.9072,41.8781,34.0522,37.7749,
          29.7604,32.7767,25.7617,33.2490,33.4484,
          47.6062,42.3601)
)

# ============================================
# 9. Define CA / NE / SE regions + polygons
# ============================================
region_states <- list(
  CA = c("california"),
  NE = c("new york", "massachusetts", "pennsylvania", "new jersey"),
  SE = c("georgia", "florida", "alabama", "south carolina")
)

region_polys <- list()

for(reg in names(region_states)){
  poly <- usa_states %>%
    filter(ID %in% region_states[[reg]]) %>%
    st_union() %>% st_make_valid()
  
  region_polys[[reg]] <- poly
}

region_labels <- data.frame(
  region = c("CA", "NE", "SE"),
  lon = c(-117, -79.5, -86.5),
  lat = c( 40.5,  45,   36.5)
)


# ============================================
# 10. FRM Map (with region highlights)
# ============================================
p_frm <- ggmap(map_bg) +
  geom_sf(data = usa_states, fill = NA, color = "white", size = 0.4,
          inherit.aes = FALSE) +
  geom_point(data = frm_land,
             aes(x = lon, y = lat, color = AQI),
             size = 2.3, alpha = 0.9) +
  scale_color_viridis(option = "plasma", name = "AQI",
                      guide = guide_colorbar(
                        label.theme = element_text(size = 12),   
                        title.theme = element_text(size = 12)    
                      )) +
  coord_sf(crs = 4326, datum = NA) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text  = element_blank())

# ⭐ Add region highlights (yellow borders)
for(reg in names(region_polys)){
  p_frm <- p_frm +
    geom_sf(data = region_polys[[reg]],
            fill = NA, color = "yellow", size = 1.4,
            inherit.aes = FALSE)
}

# ⭐ Region labels
p_frm <- p_frm +
  geom_text(data = region_labels,
            aes(x = lon, y = lat, label = region),
            color = "yellow", fontface = "bold", size = 6)

# ⭐⭐⭐ Add city labels LAST (so they appear ON TOP)
p_frm <- p_frm +
  geom_text_repel(data = cities,
                  aes(x = lon, y = lat, label = name),
                  size = 3.2, fontface = "bold",
                  color = "black",
                  box.padding = 0.3,
                  bg.color = "white")


# ============================================
# 11. CMAQ Grid + Map (with region highlights)
# ============================================
grid <- st_make_grid(cmaq_land_sf, n = c(160,100), what = "polygons") %>% st_sf()
nearest_idx <- st_nearest_feature(grid, cmaq_land_sf)
grid$PM25 <- cmaq_land_sf$PM25[nearest_idx]
grid_clipped <- st_intersection(grid, usa_union)

p_cmaq <- ggmap(map_bg) +
  geom_sf(data = grid_clipped,
          aes(fill = PM25),
          color = NA, alpha = 0.9,
          inherit.aes = FALSE) +
  geom_sf(data = usa_states,
          fill = NA, color = "white", size = 0.4,
          inherit.aes = FALSE) +
  scale_fill_viridis(option = "plasma", name = "µg/m³",
                     guide = guide_colorbar(
                       label.theme = element_text(size = 12),   
                       title.theme = element_text(size = 12)   
                     )) +
  coord_sf(crs = 4326, datum = NA) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text  = element_blank())

# ⭐ highlight regions
for(reg in names(region_polys)){
  p_cmaq <- p_cmaq +
    geom_sf(data = region_polys[[reg]],
            fill = NA, color = "yellow", size = 1.4,
            inherit.aes = FALSE)
}

# ⭐ region labels
p_cmaq <- p_cmaq +
  geom_text(data = region_labels,
            aes(x = lon, y = lat, label = region),
            color = "yellow", fontface = "bold", size = 6)

# ⭐⭐⭐ city labels LAST (top layer)
p_cmaq <- p_cmaq +
  geom_text_repel(data = cities,
                  aes(x = lon, y = lat, label = name),
                  size = 3.2, fontface = "bold",
                  color = "black",
                  box.padding = 0.3,
                  bg.color = "white")


# ============================================
# 12. Save plots
# ============================================
ggsave("US_FRM_highlight.png", p_frm, width = 7, height = 3, dpi = 320)
ggsave("US_CMAQ_highlight.png", p_cmaq, width = 7, height = 3, dpi = 320)

cat("✔ CA / NE / SE highlighted and labeled successfully! Cities are on TOP.\n")








# ================================
# 1. Skewness & Kurtosis
# ================================
cat("FRM skewness =", skewness(frm$AQI), "\n")
cat("FRM kurtosis =", kurtosis(frm$AQI), "\n\n")

cat("CMAQ skewness =", skewness(cmaq$PM25), "\n")
cat("CMAQ kurtosis =", kurtosis(cmaq$PM25), "\n")


# ================================
# 2. Normality tests
# Anderson-Darling: good for medium-large n
# ================================
cat("\n=== Normality Tests ===\n")

cat("\nFRM: Anderson-Darling test\n")
print(ad.test(frm$AQI))

cat("\nCMAQ: Anderson-Darling test\n")
print(ad.test(cmaq$PM25))




# ================================
# 1. FRM histogram
# ================================
h1 <- ggplot(frm, aes(x = AQI)) +
  geom_histogram(aes(y = ..density..), bins = 20, 
                 fill = "skyblue", color = "black") +
  geom_density(color = "red", size = 1.2) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(size = 17),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 17),
    axis.title   = element_text(size = 17),
    axis.text    = element_text(size = 17)
  ) 
  #ggtitle("FRM")


# ================================
# 2. CMAQ histogram
# ================================
h2 <- ggplot(cmaq, aes(x = PM25)) +   
  geom_histogram(aes(y = ..density..), bins = 30, 
                 fill = "lightgreen", color = "black") +
  geom_density(color = "red", size = 1.2) +
  labs(
    x = expression(PM[2.5]),   # ← 下标 2.5
    y = "density"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(size = 17),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 17),
    axis.title   = element_text(size = 17),
    axis.text    = element_text(size = 17)
  )



# ================================
# 3. Combine side-by-side
# ================================
combined_hist <- plot_grid(h1, h2, ncol = 2, align = "h")
combined_hist


# ================================
# 4. Save figure
# ================================
ggsave("hist_combined.png",
       combined_hist, width = 7, height = 2.5)


# ================================
# 1. FRM Q-Q plot
# ================================
qq1 <- ggplot(frm, aes(sample = AQI)) +
  stat_qq(size = 1) +
  stat_qq_line(color = "red", linewidth = 1) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(size = 17),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 17),
    axis.title   = element_text(size = 17),
    axis.text    = element_text(size = 17)
  ) 
  #ggtitle("FRM")


# ================================
# 2. CMAQ Q-Q plot
# ================================
qq2 <- ggplot(cmaq, aes(sample = PM25)) +   # 注意列名
  stat_qq(size = 1) +
  stat_qq_line(color = "red", linewidth = 1) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(size = 17),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 17),
    axis.title   = element_text(size = 17),
    axis.text    = element_text(size = 17)
  ) 
  #ggtitle("CMAQ")


# ================================
# 3. Combine side by side
# ================================
combined_qq <- plot_grid(qq1, qq2, ncol = 2, align = "h")
combined_qq


# ================================
# 4. Save
# ================================
ggsave("qq_combined.png",
       combined_qq, width = 7, height = 2.5)


