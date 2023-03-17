#' This script will generate the Modifiable Aerial Unit Problem figure
#' 

# ---------------------------------------------------------------------------- #
## A. Prepare for MAUP calculations ----
# ---------------------------------------------------------------------------- #
### Get counts to be used ----
# Select either the Normalised or the Normalised-and-Scaled gene counts
select <- SA_counts_Nor %>% 
  t() %>%
  as.data.frame() %>%
  filter(rowSums(. != 0) > 1000) %>% 
  rownames(.)
  
cor.input <- SA_counts_Nor %>%
  dplyr::select(all_of(select))

### Set plot output prefix ----
prfx <- "countNorm"

### Get a set of colours ----
## Colourblind-friendly colours using cols4all
## c4a palette: c4a.bu_br_bivs
cols_grid <- c("#4992C8", "#AECBEA")

## c4a palette: misc.okabe
cols_lay <- c("#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", 
              "#CC79A7")

# ---------------------------------------------------------------------------- #
## B. Calculate multiple correlations ----
# ---------------------------------------------------------------------------- #
## Spearman correlations between all genes
cor.multi <- cor(cor.input, cor.input, method = "spearman")
diag(cor.multi) <- 0 #set diagonal to zero

g1 <- seq_len(nrow(cor.multi))
g2 <- max.col(abs(cor.multi))

cor.multi.Rs <- cor.multi %>%
  as.data.frame() %>%
  mutate(gene_X = rownames(.), 
         gene_Y = colnames(.)[g2],
         R = .[cbind(seq_along(g1), g2)]) %>% 
  arrange(desc(abs(R))) %>% 
  select(gene_X, gene_Y, R) %>% 
  filter(R < 0.5 & R > 0.3)

sampleInts <- sample(nrow(cor.multi.Rs), size = 10, replace = FALSE)

for (i in sampleInts) {
  (ggscatter(cor.input, 
             x = cor.multi.Rs[i, 1], 
             y = cor.multi.Rs[i, 2], 
             add = "reg.line", conf.int = TRUE, 
             cor.coef = TRUE, cor.method = "spearman", 
             xlab = "", 
             ylab = "", 
             ggtheme = my_theme)) %>%
  print()
  
  ggsave(paste0("./data/graphics_out/",
                prfx,
                "_MAUPCorls_",
                cor.multi.Rs[i, 1], "_", cor.multi.Rs[i, 2], ".svg"),
         device = "svg",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)
  
  rm(i)
}

rm(g1, g2)

## randomly selecting cor.multi.Rs[45,] as suitable for visualisation

# ---------------------------------------------------------------------------- #
## C. Format the different zonations ----
# ---------------------------------------------------------------------------- #
### Prepare the data ----
## Get the range of the locations
rangeX <- range(meta_D[meta_D[, "Section"] == 1, "Spot_X"])
rangeY <- range(meta_D[meta_D[, "Section"] == 1, "Spot_Y"])


## We need a 361-groups grid. sqrt(361) = 19 so we need roughly 19 groups on the
## X and 19 on the Y axis. The X has a range of 115 and the Y of 76. Therefore
## 115/19 = 6.05 (so we need to break every 6 spots) and 76/19 = 4 (so we need
## to break every 4 spots).

## Following the same logic we divide the tissue into horizontal layers and store 
## the information in columns denoting groups.

dt_MAUP <- polygons %>% 
  left_join(meta_D[,c("Barcode", "Spot_X", "Spot_Y")], by = "Barcode") %>%
  mutate(group_X = as.numeric(cut(.data[["Spot_X"]], 
                                  breaks = seq(rangeX[1], rangeX[2], by = 6), 
                                  include.lowest = TRUE)),
         group_Y = as.numeric(cut(.data[["Spot_Y"]], 
                                  breaks = seq(rangeY[1], rangeY[2], by = 4), 
                                  include.lowest = TRUE)),
         group_361 = factor(paste(group_X, group_Y, sep = "_"))) %>%    # 361 groups
  mutate(group_Y = as.numeric(cut(.data[["Spot_Y"]], 
                                  breaks = seq(rangeY[1], rangeY[2], by = 2), 
                                  include.lowest = TRUE)),
         group_38 = factor(group_Y)) %>%                        # 38 groups
  mutate(group_Y = as.numeric(cut(.data[["Spot_Y"]], 
                                  breaks = seq(rangeY[1], rangeY[2], by = 4), 
                                  include.lowest = TRUE)),
         group_19 = factor(group_Y))                            # 19 groups



## The grid-like grouping will have a two-colour scheme
cols <- rep(cols_grid, 181)

colours <- expand.grid(group_X = seq(1:19),
                       group_Y = seq(1:19)) %>%
  mutate(group_361 = factor(paste(group_X, group_Y, sep = "_"))) %>%
  mutate(cols = cols[1:nrow(.)]) %>%
  dplyr::select(cols, group_361)

dt_MAUP <- dt_MAUP %>%
  left_join(colours)

### Plot: Map zones ----
groups <- c("group_361", "group_38", "group_19", "layer")

for (g in groups) {
  
  if (g == "group_361") {
    g <- "cols"
    cols <- dt_MAUP[[g]]
  } else if (g == "layer") {
    cols <- cols_lay
  } else {
    cols <- rep(cols_grid, length(unique(dt_MAUP[[g]]))) 
  }
  
  (ggplot(data = dt_MAUP) +
    geom_sf(aes(geometry = geom_pol, fill = .data[[g]]), 
            colour = "black") + 
    scale_fill_manual(values = cols) +
    labs(x = "",
         y = "") + 
    my_theme_s + 
    theme(legend.position = "none")) %>%
    print()
  
  ggsave(paste0("./data/graphics_out/",
                prfx,
                "_MAUPMaps_", g, ".svg"),
         device = "svg",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)
  
  rm(g, cols)
}

rm(rangeX, rangeY)
# ---------------------------------------------------------------------------- #
## D. Add zones to the correlations input ----
# ---------------------------------------------------------------------------- #
### Prepare the data ----
## Add zones to the correlation input df
cor.input <- cor.input %>%
  rownames_to_column(var = "Barcode") %>%
  left_join(dt_MAUP[,c("Barcode", "group_361", "group_38", "group_19", "layer")]) %>%
  st_drop_geometry() %>% 
  dplyr::select(matches("[^geom_pol]")) %>%
  column_to_rownames(var = "Barcode") 

## Because the tissue is not a proper square the number of expected groups might
## not match the number of groups we get back. Let's have a look then:
groups <- c("group_361", "group_38", "group_19")

for (g in groups) {
  expected <- gsub("group_", "", g)
  
  gs <- length(unique(cor.input[,g]))
  
  message("Expected ", expected, " groups. Got ", gs, " instead.")
  
  # rm(expected, gs, g)
}

groups <- c(groups, "layer")
## Generate the grouped correlation inputs
for (g in groups) {
  message("Working on: ", g)
  
  cIn <- cor.input %>%                                          # get the correlation input
    dplyr::select(matches(paste0("ENSG|", g))) %>%              # select one group type
    group_by_at(vars(g)) %>%                                    # group by group type and gene
    summarise(across(matches("ENSG"), tibble::lst(sum, mean)))  # summarise per group per column
  
  for (ends in c("_sum", "_mean")) {
    out <- cIn %>%
      dplyr::select(ends_with(ends)) %>% 
      ungroup() %>%
      select_all(~gsub(ends, "", .))
    
    assign(paste0("cor.input_", g, ends), out)
  }
  
  rm(g, cIn, ends, out)
}

### Calculate correlations ----
pattern <- grep("_mean|_sum", names(.GlobalEnv), value = TRUE)

for (p in pattern) {
  message("Working on: ", p)
  input <- get(p)
  
  for (i in sampleInts) {
    sName <- gsub("cor.input_", "", p)
    
    (ggscatter(input, 
               x = cor.multi.Rs[45, 1], 
               y = cor.multi.Rs[45, 2], 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman", 
               xlab = "", 
               ylab = "", 
               ggtheme = my_theme)) %>%
      print()
    
    ggsave(paste0("./data/graphics_out/", 
                  prfx, 
                  "_MAUPCorls_",  
                  cor.multi.Rs[45, 1], "_", cor.multi.Rs[45, 2], 
                  "_", sName, ".svg"),
           device = "svg",
           width = grDevices::dev.size(units = "in")[1],
           height = grDevices::dev.size(units = "in")[2],
           units = "in",
           dpi = 300)
    
    rm(i, sName)
  }
  
  rm(p, input)
}

rm(groups, cor.multi.Rs)

# ---------------------------------------------------------------------------- #
# E. Empty tissue map ----
# ---------------------------------------------------------------------------- #

ggplot(data = dt_MAUP) +
  geom_sf(aes(geometry = geom_pol), 
          colour = "black") + 
  scale_fill_manual(values = "grey36") +
  labs(x = "",
       y = "") + 
  my_theme_s + 
  theme(legend.position = "none")
