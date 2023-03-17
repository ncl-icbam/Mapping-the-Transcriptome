#' This script will generate the Spatial Autocorrelation figure
#' 

# ---------------------------------------------------------------------------- #
## A. Prepare for SH calculations ----
# ---------------------------------------------------------------------------- #
### Get Seurat-normalised data ----
data <- SA_counts_Nor

### Select a set of genes to test ----
## The list of 2000 genes from Seurat comes after variance stabilising and 
## selecting the most variable genes. We want to narrow this further down and 
## get a pool of genes that are variable and expressed in more than 2/3 of the 
## tissue area.
#### Filter for number of locations expressed ----
SH_genes <-  data %>%
  t() %>%
  as.data.frame() %>%
  filter(rowSums(. != 0) > 1000) %>% 
  rownames(.)

### Add geometries and convert to sp format ----
dt_SH <- data %>%
  dplyr::select(all_of(SH_genes)) %>%
  rownames_to_column(var = "Barcode") %>%            # barcodes from row names to column
  left_join(polygons[,c("geom_pol", "Barcode")]) %>% # merge with geometries
  column_to_rownames(var = "Barcode") %>%            # barcodes back to row names
  st_as_sf(., sf_column_name = "geom_pol") %>%       # data frame to sf
  as(., "Spatial")                                   # sf to sp - GWmodel still works with sp objects

### Set plot output prefix ----
prfx <- "countNorm"

## Housekeeping
rm(data)

# ---------------------------------------------------------------------------- #
## B. Run Multiple Linear Regression ----
# ---------------------------------------------------------------------------- #
### Create a set of combinations to model ----
model_vars <- expand.grid(gene_1 = SH_genes, 
                          gene_2 = SH_genes) %>% 
  mutate(model_vars = paste(gene_1, gene_2, sep = "~")) %>% 
  filter(!(gene_1 == gene_2))

r_list <- vector(mode = "list", length = nrow(model_vars))
names(r_list) <- model_vars$model_vars

### Run linear regression for multiple pairs ----
for (mv in model_vars$model_vars) {
  message("Working on model: ", mv)
  m = lm(mv, data = SA_counts_Nor)
  m_summ <- summary(m)
  
  r_list[mv] <- m_summ$adj.r.squared
  
  rm(mv, m, m_summ)
}

lm_Rsq <- data.frame(r_list) %>% 
  t() %>% 
  as.data.frame()
colnames(lm_Rsq) <- "R_sqrd"
lm_Rsq <- lm_Rsq %>%
  arrange(desc(R_sqrd))

rm(r_list)
# ---------------------------------------------------------------------------- #
## C. Run GWR ----
# ---------------------------------------------------------------------------- #
### Select response and predictor variables with R^2 above 0.5 ----
model_vars <- lm_Rsq %>%
  filter(R_sqrd > 0.40) %>%
  rownames_to_column(var = "model") %>% 
  mutate(model = gsub("\\.", "~", model))
mdl = "ENSG00000091513~ENSG00000197971"
for (mdl in model_vars$model) {
  ### Determine the kernel bandwidth ----
  bw <- bw.gwr(mdl,
               approach = "AIC",
               adaptive = T,
               data = dt_SH)
  
  ## The bandwidth is an adaptive bandwidth (kernel size) indicating that its size 
  ## varies
  message("Bandwidth (bw): ", bw)
  
  ### Run GWR ----
  gwr <- gwr.basic(mdl, 
                   adaptive = T,
                   data = dt_SH,
                   bw = bw)
  
  gwr.tab <- apply(gwr$SDF@data[, 1:7], 2, summary)
  gwr.tab <- round(gwr.tab, 1)
  print(t(gwr.tab[,1:2]))
  
  ### Prepare the data for plotting ----
  predictorVar <- gsub(".+?~", "", mdl)
  
  gwr_sf <- st_as_sf(gwr$SDF) %>% 
    mutate(signif0 = factor(if_else(abs(.data[["Intercept_TV"]]) > 1.96, 
                                    "Signif.", 
                                    "Not-signif.")),
           signif1 = factor(if_else(abs(.data[[paste0(predictorVar, "_TV")]]) > 1.96, 
                                   "Signif.", 
                                   "Not-signif.")),
           b0_geom = if_else(signif0 == "Signif.", 
                             geometry, 
                             NA),
           b1_geom = if_else(signif1 == "Signif.", 
                            geometry, 
                            NA))
  
  for (i in c(0,1)) {
    tab = round(summary(gwr$lm)$coefficients, 3)
    colnames(tab)[3:4] = c("t-value", "p-value")
    rownames(tab) = c("Intercept", "Covariate")
    
    tab = data.frame(tab)
    tab1 = tab[1,]
    tab2 = tab[2,]
    
    if (i == 0) {
      fill <-  "Intercept"
      b_geom <- "b0_geom"
      tit = expression(""*beta[0]*"")
      table = tab1
    } else {
      fill  <- predictorVar
      b_geom <- "b1_geom"
      tit = expression(""*beta[1]*"")
      table = tab2
    }
    
    ### Plot: Map Coefficients that flip ----
    (ggplot(gwr_sf) + 
      geom_sf(aes(geometry = geometry, 
                  fill = get(fill)),
              colour = NA) +
      geom_sf(aes(geometry = get(b_geom)),
              fill = alpha("white", 0),
              colour = "black",
              linewidth = 0.5) +
      scale_fill_gradient2(high = "#B2182B",
                           mid = "#F7F7F7",
                           low = "#2166AC",
                           midpoint = 0,
                           n.breaks = 7) +
      # scale_fill_binned_c4a_div(palette = "brewer.rd_bu", 
      #                           mid = 0, 
      #                           show.limits = TRUE, 
      #                           n.breaks = 11) +
      labs(x = "",
           y = "",
           fill = tit) +
      annotation_custom(tableGrob(table), 
                        xmin = 150, xmax = 450, 
                        ymin = -570, ymax = -500) +
      my_theme_s) %>%
      print()
    
    ggsave(paste0("./data/graphics_out/", 
                  prfx, "_SHFlip_", mdl, "beta", i, ".svg"),
           device = "svg",
           width = grDevices::dev.size(units = "in")[1],
           height = grDevices::dev.size(units = "in")[2],
           units = "in",
           dpi = 300)
  }
  
  
  
  rm(mdl, predictorVar, gwr, gwr.tab, gwr_sf, bw, tit, i, fill)
}

# high = "#67001F"
# low = "#053061"

