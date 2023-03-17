#' This script will generate the Spatial Autocorrelation figure
#' 

# ---------------------------------------------------------------------------- #
## A. Prepare for SA calculations ----
# ---------------------------------------------------------------------------- #
### Get spot names ----
nb_names <- polygons$Barcode

### Get counts for the 2000 most variable genes ----
### The counts are normalised 
### It is important at the last step to order the rows as the vector of names above
pattern <- grep("obs_Nor|obs$", names(.GlobalEnv), value = TRUE)

for (p in pattern) {
  dt <- get(p)
  
  out_name <- paste0("SA_counts", gsub("obs", "", p))
  
  dt_SA <- dt %>% 
    filter(rownames(.) %in% var_genes) %>%
    t() %>%
    as.data.frame() %>%
    .[nb_names,]
  
  assign(out_name, as.data.frame(dt_SA))
  
  rm(p, dt, out_name, dt_SA)
}

### Find neighbours ----
neighbours_knn <- knn2nb(knearneigh(polygons$geom_cntd, k = 6), 
                         row.names = nb_names)

### Calculate neighbour weights ----
neighbours_w_exp <- nb2listwdist(neighbours_knn, polygons$geom_cntd,
                                 type = "idw", style = "raw", alpha = 2,
                                 zero.policy = TRUE)

### Set plot output prefix ----
prfx <- "countNorm"

### Set active dataset ----
data <- SA_counts_Nor

# ---------------------------------------------------------------------------- #
## B. Calculate GLOBAL Moran's I for all genes ----
# ---------------------------------------------------------------------------- #
### Prepare the data ----
moran.n <- length(neighbours_w_exp$neighbours)
moran.S0 <- Szero(neighbours_w_exp)
### Calculate Moran's I ----
moran.multi <- apply(data, 2, function(x) moran(x, 
                                                  neighbours_w_exp,
                                                  n = moran.n,
                                                  S0 = moran.S0)) %>% 
  purrr::reduce(bind_rows) %>% 
  as.data.frame()

rownames(moran.multi) <- colnames(data)
### Order output ----
moran.multi <- moran.multi[order(moran.multi$I, decreasing = TRUE),]
### Select genes ----
## Genes with either high positive or close to zero Moran's I
SA_genes_pos <- rownames(moran.multi)[1:10]
SA_genes_zero <- tail(rownames(moran.multi), n = 10)
SA_genes <-  cbind(SA_genes_pos, SA_genes_zero)

# ---------------------------------------------------------------------------- #
## C. Calculate LOCAL Moran's I ----
# ---------------------------------------------------------------------------- #
## For the genes with highest, positive Moran's I
### Prepare the data ----
gCounts <- data[,SA_genes_pos]
### Calculate Moran's I ----
moran.multi.local <- lapply(gCounts, 
                            function(x) {
                              names(x) <- rownames(gCounts)
                              localmoran(x,
                                         listw = neighbours_w_exp, 
                                         spChk = TRUE)
                            })

rm(gCounts)

# ---------------------------------------------------------------------------- #
## D. GLOBAL Moran's I plots and tables ----
# ---------------------------------------------------------------------------- #
### Plot: Moran plots ----
## These plots are for the top 10 genes with positive SA and the 10 genes with 
## the lowest Moran's I (which means there is no SA or is alsmost no SA)
for (g in SA_genes) {
  message("Working on: ", g)
  
  if (g %in% SA_genes_pos) {
    SA.type <- "SA.pos"
  } else if (g %in% SA_genes_zero) {
    SA.type <-  "SA.zero"
  }
  
  ## Get gene expression
  gene.exp <- data[,g] %>%
    setNames(rownames(data))
  
  ## Create Moran Plot (mp)
  mp <- moran.plot(x = gene.exp, 
                   listw = neighbours_w_exp,
                   labels = FALSE,
                   zero.policy = TRUE,
                   spChk = TRUE)
  
  ## Plot mp
  ggplot(mp, aes(x = x, y = wx)) + geom_point(shape = 1) + 
    geom_smooth(formula = y ~ x, method = "lm") + 
    geom_hline(yintercept = mean(mp$wx), lty = 2) + 
    geom_vline(xintercept = mean(mp$x), lty = 2) + theme_minimal() + 
    geom_point(data = mp[mp$is_inf,], aes(x = x, y = wx), shape = 9) +
    labs(x = "Gene expression",
         y = "Spatially lagged gene expression") + 
    theme_bw() + 
    theme(panel.border = element_rect(colour = "black", linewidth = rel(4))) + 
    my_theme
  
  ## Save mp
  ggsave(paste0("./data/graphics_out/", prfx, "_MoranPlot_", SA.type, "_", g, ".svg"),
         device = "svg",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)
  
  rm(g, SA.type, mp, gene.exp)
}

### Tables: Moran stats ----
for (g in SA_genes) {
  message("Working on: ", g)
  
  if (g %in% SA_genes_pos) {
    SA.type <- "SA.pos"
  } else if (g %in% SA_genes_zero) {
    SA.type <-  "SA.zero"
  }
  
  ## Calculate Moran's I tests with MC perm. and Z-score
  ## Z-score
  moran.stat.Z <- moran.test.ste(x = data[,g], 
                             listw = neighbours_w_exp) %>% 
    unlist() %>% 
    t() %>%
    as.data.frame() %>% 
    select(starts_with("estimate"), p.value) %>% 
    dplyr::rename("Moran's I" = "estimate.Moran I statistic",
           "P-value" = "p.value", 
           "Expected" = "estimate.Expectation", 
           "Variance" = "estimate.Variance") %>%
    mutate(`Test type` = "Z-score",
           `Moran's I` = round(as.numeric(`Moran's I`), 4),
           Expected = round(as.numeric(Expected), 6),
           Variance = round(as.numeric(Variance), 6),) %>% 
    relocate(`Test type`, `Moran's I`, `P-value`, .before = 1)
  
  ## Monte Carlo permutations
  moran.stat.MC <- moran.mc(x = data[,g], 
                            listw = neighbours_w_exp, 
                            nsim = 999) %>% 
    unlist() %>% 
    t() %>%
    as.data.frame() %>% 
    select(statistic.statistic, p.value) %>% 
    dplyr::rename("Moran's I" = "statistic.statistic",
                  "P-value" = "p.value") %>%
    mutate(`Test type` = "MC perm",
           "Expected" = "N/A", 
           "Variance" = "N/A",
           `Moran's I` = round(as.numeric(`Moran's I`), 4)) %>% 
    relocate(`Test type`, `Moran's I`, `P-value`, .before = 1)
  
  tb <- rbind(moran.stat.Z, moran.stat.MC)
  
  ## Create table
  tb <-  tb %>%
    flextable::flextable() %>%
    flextable::theme_vanilla() %>%
    flextable::fontsize(size = 12) %>%
    flextable::fontsize(size = 14, part = "header") %>%
    flextable::align_text_col(align = "center") %>%
    flextable::set_caption(caption = "SA statistics (Moran's I)")
  
  ## Save table
  flextable::save_as_image(tb, path = paste0("./data/graphics_out/", 
                                             prfx, "_MoranStatsTbl_", 
                                             SA.type, "_", g, ".png"))
  
  rm(tb, SA.type, moran.stat.Z, moran.stat.MC, g)
}


### Plot: Map Expression ---- 
## Plot from Seurat-normalised counts and visualise them as log2-transformed
dt_type <- "Nor"
dt_fill <- "Normalised "
  
for (g in SA_genes) {
  message("Working on: ", g)
  
  if (g %in% SA_genes_pos) {
    SA.type <- "SA.pos"
  } else if (g %in% SA_genes_zero) {
    SA.type <-  "SA.zero"
  }
  
  ## Create a test counts df to map the gene expression
  test.counts <- data %>%
    dplyr::select(.data[[g]]) %>% 
    rownames_to_column(var = "Barcode") %>%
    left_join(polygons["Barcode"])
  
  ## Map gene expression
  ggplot(data = test.counts) + 
    geom_sf(aes(geometry = geom_pol,
                fill = log2(get(g) + 1))) + 
    scale_fill_distiller(palette = "YlGnBu") + 
    labs(x = "",
         y = "",
         fill = paste0(dt_fill, "\ngene expression")) + 
    my_theme_s
  
  ## Save gene expression
  ggsave(paste0("./data/graphics_out/", 
                prfx, 
                "_MapGeneExpr_", 
                dt_type, "_",
                SA.type, "_", g, ".svg"),
         device = "svg",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)
  
  rm(g, SA.type, test.counts)
  
}

rm(dt_type, dt_fill)



# ---------------------------------------------------------------------------- #
## E. LOCAL Moran's I plots and tables ----
# ---------------------------------------------------------------------------- #
pattern <- names(moran.multi.local)

for (p in pattern) {
  message("Working on: ", p)
### Prepare the data ----
  SA.type <- "SA.pos"
  
  quadr <- attr(moran.multi.local[[p]], "quadr")
  
  moran.local.map <- moran.multi.local[[p]] %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Barcode") %>% 
    dplyr::rename("p.value" = `Pr(z != E(Ii))`) %>% 
    left_join(polygons[,"Barcode"]) %>%
    mutate(Quadrant = quadr$pysal,
           colours = Quadrant,
           colours = as.factor(case_when(p.value > 0.05 ~ "#ffffff",
                                         Quadrant == "Low-Low" ~ "#0066ff",
                                         Quadrant == "High-Low" ~ "#ffb3b3",
                                         Quadrant == "Low-High" ~ "#99c2ff",
                                         Quadrant == "High-High" ~ "#ff0000")),
           label = as.factor(case_when(p.value > 0.05 ~ "Not signif.",
                                       Quadrant == "Low-Low" ~ "Low-Low",
                                       Quadrant == "High-Low" ~ "High-Low",
                                       Quadrant == "Low-High" ~ "Low-High",
                                       Quadrant == "High-High" ~ "High-High")),
           signif = as.factor(case_when(p.value < 0.05 ~ "Signif.",
                              p.value >= 0.05 ~ "Not Signif.")),
           b_geom = if_else(signif == "Signif.", 
                            geom_pol, 
                            NA)) 
  
### Plot: Map Local Moran ----
  ggplot(data = moran.local.map) +
    geom_sf(aes(geometry = geom_pol, fill = Ii), 
            colour = NA) + 
    scale_fill_viridis_c(option = "inferno") + 
    geom_sf(aes(geometry = b_geom),
            fill = NA,
            colour = "#E0E0E0",
            linewidth = 0.25) +
    labs(x = "",
         y = "",
         fill = "Local Moran\nIi") +
    my_theme_s
  
  ## Save
  ggsave(paste0("./data/graphics_out/", 
                prfx, 
                "_MapMoranLoc_",
                SA.type, "_", p, ".svg"),
         device = "svg",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)
  
### Plot: Map LISA clusters ----
  ggplot(data = moran.local.map) +
    geom_sf(aes(geometry = geom_pol, fill = colours), 
            colour = "grey30") + 
    scale_fill_identity(guide = "legend",
                        labels = moran.local.map$label,
                        breaks = moran.local.map$colours) +
    labs(x = "",
         y = "",
         fill = "Local Moran\nIi Clusters") +
    my_theme_s
  
  ## Save
  ggsave(paste0("./data/graphics_out/", 
                prfx, 
                "_MapMoranLocClusts_",
                SA.type, "_", p, ".svg"),
         device = "svg",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)
  
  rm(p, SA.type, quadr, moran.local.map)
}

# ---------------------------------------------------------------------------- #
## F. Negative control ----
# ---------------------------------------------------------------------------- #
### Prepare the data ----
g = "ENSG00000197971"

order <- sample(1:nrow(data), size = nrow(data), replace = FALSE)

N_control <- data %>%
  dplyr::select(all_of(g)) %>%
  rownames_to_column(var = "Barcode") %>% 
  mutate(Barcode = Barcode[order]) %>%
  left_join(polygons[, c("Barcode", "geom_cntd")])

rownames(N_control) <- N_control$Barcode

neighbours_knn_N_control <- knn2nb(knearneigh(N_control$geom_cntd, k = 6), 
                                   row.names = rownames(N_control))

neighbours_N_Control <- nb2listwdist(neighbours_knn_N_control, N_control$geom_cntd,
                                     type = "idw", style = "raw", alpha = 2,
                                     zero.policy = TRUE)

dt_fill <- "Normalised"
dt_type <- "Nor"
SA.type <- "SA.pos"

### Plot: Map Expression ----
## Map gene expression
ggplot(data = N_control) + 
  geom_sf(aes(geometry = geom_pol,
              fill = log2(get("ENSG00000197971") + 1))) + 
  scale_fill_distiller(palette = "YlGnBu") + 
  labs(x = "",
       y = "",
       fill = paste0(dt_fill, "\ngene expression")) + 
  my_theme_s

## Save gene expression
ggsave(paste0("./data/graphics_out/", 
              prfx, 
              "_MapGeneExpr_", 
              dt_type, "_",
              SA.type, "_", g, "_NegCntrl.svg"),
       device = "svg",
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 300)

### Tables: Moran stats ----
## Calculate Moran's I tests with MC perm. and Z-score
## Z-score
moran.stat.Z <- moran.test.ste(x = N_control[,g], 
                               listw = neighbours_N_Control) %>% 
  unlist() %>% 
  t() %>%
  as.data.frame() %>% 
  select(starts_with("estimate"), p.value) %>% 
  dplyr::rename("Moran's I" = "estimate.Moran I statistic",
                "P-value" = "p.value", 
                "Expected" = "estimate.Expectation", 
                "Variance" = "estimate.Variance") %>%
  mutate(`Test type` = "Z-score",
         `Moran's I` = round(as.numeric(`Moran's I`), 4),
         Expected = round(as.numeric(Expected), 6),
         Variance = round(as.numeric(Variance), 6),) %>% 
  relocate(`Test type`, `Moran's I`, `P-value`, .before = 1)

## Monte Carlo permutations
moran.stat.MC <- moran.mc(x = N_control[,g], 
                          listw = neighbours_N_Control, 
                          nsim = 999) %>% 
  unlist() %>% 
  t() %>%
  as.data.frame() %>% 
  select(statistic.statistic, p.value) %>% 
  dplyr::rename("Moran's I" = "statistic.statistic",
                "P-value" = "p.value") %>%
  mutate(`Test type` = "MC perm",
         "Expected" = "N/A", 
         "Variance" = "N/A",
         `Moran's I` = round(as.numeric(`Moran's I`), 4)) %>% 
  relocate(`Test type`, `Moran's I`, `P-value`, .before = 1)

tb <- rbind(moran.stat.Z, moran.stat.MC)

## Create table
tb <-  tb %>%
  flextable::flextable() %>%
  flextable::theme_vanilla() %>%
  flextable::fontsize(size = 12) %>%
  flextable::fontsize(size = 14, part = "header") %>%
  flextable::align_text_col(align = "center") %>%
  flextable::set_caption(caption = "SA statistics (Moran's I)")

## Save table
flextable::save_as_image(tb, path = paste0("./data/graphics_out/", 
                                           prfx, "_MoranStatsTbl_", 
                                           SA.type, "_", g, "_NegCntrl.png"))

### Plot: Map Local Moran ----
## For the selected gene
SA_genes <- "ENSG00000197971"
#### Prepare the data ----
gCounts <- N_control[,SA_genes]
#### Calculate Moran's I ----
names(gCounts) <- N_control$Barcode
moran.multi.local <- localmoran(gCounts,
                                listw = neighbours_N_Control, 
                                spChk = TRUE)

rm(gCounts)

moran.local.map <- moran.multi.local %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Barcode") %>% 
  dplyr::rename("p.value" = `Pr(z != E(Ii))`) %>% 
  left_join(polygons[,"Barcode"]) %>% 
  mutate(signif = as.factor(case_when(p.value < 0.05 ~ "Signif.",
                                      p.value >= 0.05 ~ "Not Signif.")),
         b_geom = if_else(signif == "Signif.", 
                          geom_pol, 
                          NA))

#### Plot: Map Local Moran ----
ggplot(data = moran.local.map) +
  geom_sf(aes(geometry = geom_pol, fill = Ii), 
          colour = "grey30") + 
  scale_fill_viridis_c(option = "inferno", limits = c(-0.5, 1)) +
  geom_sf(aes(geometry = b_geom),
          fill = NA,
          colour = "#E0E0E0",
          linewidth = 0.25) +
  labs(x = "",
       y = "",
       fill = "Local Moran\nIi") +
  my_theme_s

## Save local Moran
ggsave(paste0("./data/graphics_out/", 
              prfx, 
              "_MapMoranLoc_",
              "_NegCntrl", "_", g, ".svg", ".svg"),
       device = "svg",
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 300)

rm(g, dt_fill, dt_type, SA.type, order,
   moran.stat.Z, moran.stat.MC, tb, neighbours_N_Control, 
   neighbours_knn_N_control, N_control)

# ---------------------------------------------------------------------------- #
## G. No SA example ----
# ---------------------------------------------------------------------------- #
## For an example of no SA we further filtered the genes and selected genes that
## were expressed in more than 2000 spots. This was done because from the 2000
## most variable genes that were selected initially, the genes with no SA are 
## genes that are not present in more that 100 locations on the tissue. This 
## means that are very sparse and no meaningful demonstration of no SA is 
## available.

### Prepare the data ----
noSA_dt <- data %>%
  t() %>% 
  as.data.frame() %>% 
  filter(rowSums(. != 0) > 2000) %>% 
  t() %>% 
  as.data.frame()

SA.type <- "SA.zero"

### Calculate Moran's I ----
## Global Moran's I
moran.multi <- apply(noSA_dt, 2, function(x) moran(x, 
                                                neighbours_w_exp,
                                                n = moran.n,
                                                S0 = moran.S0)) %>% 
  purrr::reduce(bind_rows) %>% 
  as.data.frame()

rownames(moran.multi) <- colnames(noSA_dt)
### Order output ----
moran.multi <- moran.multi[order(moran.multi$I, decreasing = TRUE),]
### Select genes ----
## Genes with close to zero Moran's I --> select the last one
noSA_g <- tail(rownames(moran.multi), n = 1)


### Tables: Moran stats ----
## Calculate Moran's I tests with MC perm. and Z-score
## Z-score
moran.stat.Z <- moran.test.ste(x = noSA_dt[,noSA_g], 
                               listw = neighbours_w_exp) %>% 
  unlist() %>% 
  t() %>%
  as.data.frame() %>% 
  select(starts_with("estimate"), p.value) %>% 
  dplyr::rename("Moran's I" = "estimate.Moran I statistic",
                "P-value" = "p.value", 
                "Expected" = "estimate.Expectation", 
                "Variance" = "estimate.Variance") %>%
  mutate(`Test type` = "Z-score",
         `Moran's I` = round(as.numeric(`Moran's I`), 4),
         Expected = round(as.numeric(Expected), 6),
         Variance = round(as.numeric(Variance), 6),) %>% 
  relocate(`Test type`, `Moran's I`, `P-value`, .before = 1)

## Monte Carlo permutations
moran.stat.MC <- moran.mc(x = noSA_dt[,noSA_g], 
                          listw = neighbours_w_exp, 
                          nsim = 999) %>% 
  unlist() %>% 
  t() %>%
  as.data.frame() %>% 
  select(statistic.statistic, p.value) %>% 
  dplyr::rename("Moran's I" = "statistic.statistic",
                "P-value" = "p.value") %>%
  mutate(`Test type` = "MC perm",
         "Expected" = "N/A", 
         "Variance" = "N/A",
         `Moran's I` = round(as.numeric(`Moran's I`), 4)) %>% 
  relocate(`Test type`, `Moran's I`, `P-value`, .before = 1)

tb <- rbind(moran.stat.Z, moran.stat.MC)

## Create table
tb <-  tb %>%
  flextable::flextable() %>%
  flextable::theme_vanilla() %>%
  flextable::fontsize(size = 12) %>%
  flextable::fontsize(size = 14, part = "header") %>%
  flextable::align_text_col(align = "center") %>%
  flextable::set_caption(caption = "SA statistics (Moran's I)")

## Save table
flextable::save_as_image(tb, path = paste0("./data/graphics_out/", 
                                           prfx, "_MoranStatsTbl_", 
                                           SA.type, "_", noSA_g, ".png"))

rm(tb, moran.stat.Z, moran.stat.MC, moran.multi)

### Plot: Map Expression ---- 
## Plot from Seurat-normalised counts and visualise them as log2-transformed
dt_fill <- "Normalised"
dt_type <- "Nor"

## Create a test counts df to map the gene expression
test.counts <- noSA_dt %>%
  dplyr::select(all_of(noSA_g)) %>% 
  rownames_to_column(var = "Barcode") %>%
  left_join(polygons["Barcode"])

## Map gene expression
ggplot(data = test.counts) + 
  geom_sf(aes(geometry = geom_pol,
              fill = log2(get(noSA_g) + 1))) + 
  scale_fill_distiller(palette = "YlGnBu") + 
  labs(x = "",
       y = "",
       fill = paste0(dt_fill, "\ngene expression")) + 
  my_theme_s

## Save gene expression
ggsave(paste0("./data/graphics_out/", 
              prfx, 
              "_MapGeneExpr_", 
              dt_type, "_",
              SA.type, "_", noSA_g, ".svg"),
       device = "svg",
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 300)

rm(noSA_g, noSA_dt, SA.type, test.counts, dt_type, dt_fill)


# ---------------------------------------------------------------------------- #
## H. Clean up the environment ----
# ---------------------------------------------------------------------------- #
rm(data)
