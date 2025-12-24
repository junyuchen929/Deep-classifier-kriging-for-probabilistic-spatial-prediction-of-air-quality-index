# ===================================
# Final DK Map Plot System (CA/NE/SE)
# ===================================

rm(list = ls())
setwd("/Users/beauniverse/Documents/GitHub/Deep-classifier-kriging-for-probabilistic-spatial-prediction-of-air-quality-index/realdata_subset")

library(sf)
library(ggplot2)
library(viridis)
library(dplyr)
library(cowplot)
library(ggrepel)
library(maps)
sf_use_s2(FALSE)

regions <- c("CA", "NE", "SE")

usa_states_raw <- map("state", plot = FALSE, fill = TRUE)
usa_states <- st_as_sf(usa_states_raw)
st_crs(usa_states) <- 4326
usa_states <- st_make_valid(usa_states)

region_states <- list(
  CA = c("california"),
  NE = c("new york", "massachusetts", "pennsylvania", "new jersey"),
  SE = c("georgia", "florida", "alabama", "south carolina")
)

for(reg in regions){
  message("\n============================")
  message("Plotting Region: ", reg)
  message("============================")
  
  df_pred <- read.csv(paste0("predicted_results_heatmap_", reg, ".csv"))
  df_frm  <- read.csv(paste0("FRM_", reg, "_subset.csv"))
  
  pred_sf <- st_as_sf(df_pred, coords=c("longitude","latitude"), crs=4326)
  frm_sf  <- st_as_sf(df_frm,  coords=c("lon","lat"), crs=4326)
  
  # compute quantiles
  samples <- grep("^sample_", names(df_pred), value=TRUE)
  df_pred$Q025 <- apply(df_pred[,samples],1,quantile,0.025)
  df_pred$Q50  <- apply(df_pred[,samples],1,quantile,0.50)
  df_pred$Q975 <- apply(df_pred[,samples],1,quantile,0.975)
  
  # region polygon
  region_poly <- usa_states %>%
    filter(ID %in% region_states[[reg]]) %>%
    st_union() %>% st_make_valid()
  
  # grid only inside region
  bbox <- st_bbox(region_poly)
  grid_reg <- st_make_grid(region_poly, cellsize = (bbox[4]-bbox[2])/200,
                           what="polygons") %>% st_sf() %>%
    st_intersection(region_poly)
  
  nearest_idx <- st_nearest_feature(grid_reg, pred_sf)
  grid_reg$Q025 <- df_pred$Q025[nearest_idx]
  grid_reg$Q50  <- df_pred$Q50[nearest_idx]
  grid_reg$Q975 <- df_pred$Q975[nearest_idx]
  
  get_limits <- function(reg){
    if (reg == "CA") return(c(0, 70))
    if (reg == "NE") return(c(0, 45))
    return(c(5, 40))   # SE or others
  }
  
  plot_fun <- function(value, title){
    
    lims <- get_limits(reg)
    
    ggplot() +
      geom_sf(data = grid_reg, aes(fill = .data[[value]]), color = NA) +
      geom_sf(data = region_poly, fill = NA, color = "white", size=0.3) +
      scale_fill_viridis(option="plasma",
                         name="AQI",
                         limits = lims,
                         guide = guide_colorbar(
                           label.theme = element_text(size = 17),  
                           title.theme = element_text(size = 17)    
                         )) +
      coord_sf(crs=4326, datum=NA) +
      #ggtitle(title) +
      theme_minimal() +
      theme(axis.text=element_blank(),
            axis.title=element_blank(),
            panel.grid=element_blank(),
            plot.title = element_text(size=10))
  }
  
  
  p_frm <- ggplot() +
    geom_sf(data = region_poly, fill="gray95", color="white") +
    geom_sf(data = frm_sf, aes(color = AQI), size=2) +
    scale_color_viridis(option="plasma",
                        name="AQI",
                        limits = get_limits(reg),
                          guide = guide_colorbar(
                            label.theme = element_text(size = 17),   
                            title.theme = element_text(size = 17)    
                          )) +
    coord_sf(crs=4326, datum=NA) +
    #ggtitle(paste0("Observed FRM - ", reg)) +
    theme_minimal() +
    theme(axis.text=element_blank(),
          axis.title=element_blank(),
          panel.grid=element_blank(),
          plot.title = element_text(size=10))
  
  
  # ---- predicted quantiles ----
  p_q50  <- plot_fun("Q50","Pred Q50")
  p_q025 <- plot_fun("Q025","Pred Q2.5")
  p_q975 <- plot_fun("Q975","Pred Q97.5")
  
  # ---- combine figure ----
  final <- plot_grid(p_frm,p_q50,p_q025,p_q975,
                     ncol=2, align="hv")
  
  out_file <- paste0("DK_FRM_Quantile_", reg, ".png")
  ggsave(out_file, final, width=7, height=6, dpi=400)
  message("✔ Saved: ", out_file)
}

message("\n🎯 All regions finished! Maps ready for publication!\n")







# =====================================================
# Extreme AQI Probability Map: P(AQI > 40)
# =====================================================

for(reg in regions){
  message("\n============================")
  message("Plotting Extreme AQI Probability: ", reg)
  message("============================")
  
  # -----------------------------
  # Load data
  # -----------------------------
  df_pred <- read.csv(paste0("predicted_results_heatmap_", reg, ".csv"))
  pred_sf <- st_as_sf(df_pred, coords=c("longitude","latitude"), crs=4326)
  
  # -----------------------------
  # Identify posterior samples
  # -----------------------------
  sample_cols <- grep("^sample_", names(df_pred), value = TRUE)
  
  # -----------------------------
  # Compute P(AQI > 40)
  # -----------------------------
  prob_extreme <- apply(
    df_pred[, sample_cols],
    1,
    function(x) mean(x > 40)
  )
  
  df_pred$Prob_GT40 <- prob_extreme
  
  # -----------------------------
  # Region polygon
  # -----------------------------
  region_poly <- usa_states %>%
    filter(ID %in% region_states[[reg]]) %>%
    st_union() %>%
    st_make_valid()
  
  # -----------------------------
  # Grid within region
  # -----------------------------
  bbox <- st_bbox(region_poly)
  grid_reg <- st_make_grid(
    region_poly,
    cellsize = (bbox[4] - bbox[2]) / 200,
    what = "polygons"
  ) %>%
    st_sf() %>%
    st_intersection(region_poly)
  
  # -----------------------------
  # Map probability to grid
  # -----------------------------
  nearest_idx <- st_nearest_feature(grid_reg, pred_sf)
  grid_reg$Prob_GT40 <- df_pred$Prob_GT40[nearest_idx]
  
  # -----------------------------
  # Plot probability map
  # -----------------------------
  p_prob <- ggplot() +
    geom_sf(data = grid_reg, aes(fill = Prob_GT40), color = NA) +
    geom_sf(data = region_poly, fill = NA, color = "white", size = 0.3) +
    scale_fill_viridis(
      option = "plasma",
      name = expression(P(AQI > 40)),
      # limits = c(0, 1),
      # breaks = c(0, 0.25, 0.5, 0.75, 1)
      # breaks = pretty(grid_reg$Prob_GT40, n = 4),
      oob = scales::squish,
      guide = guide_colorbar(
        label.theme = element_text(size = 14),   
        title.theme = element_text(size = 14)    
      )
    ) +
    coord_sf(crs = 4326, datum = NA) +
    #ggtitle(paste0("Probability of Extreme AQI (>", 40, ") - ", reg)) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      # plot.title = element_text(size = 14)
    )
  
  # -----------------------------
  # Save figure
  # -----------------------------
  out_file <- paste0("DK_Prob_AQI_GT40_", reg, ".png")
  ggsave(out_file, p_prob, width = 5, height = 3, dpi = 400)
  message("✔ Saved: ", out_file)
}

message("\n✅ Extreme AQI probability maps completed!\n")



