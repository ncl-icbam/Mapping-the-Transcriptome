#' This script is to load the input data
#' 
#' The data are based on the DLPFC data set (see below info about the authors).
#' 

## Load gene count data ----
## Raw counts data as outputted by spaceranger.
obs <- readRDS("./data/rObjects/DLPFC_151673_inputD.rds")
colnames(obs) <- gsub("-", ".", colnames(obs))

## Load location metadata ----
## Location information of spots (X-Y coordinates in image pixels is included).
meta_D <- readRDS("./data/rObjects/DLPFC_151673_inputMD.rds") %>% 
  mutate(Barcode = gsub("-", ".", Barcode))

## Load data as Seurat object ----
obs_seu <- readRDS("./data/rObjects/DLPFC_151673_SeuratObject.rds")

## Load ground truth (whole data set) ----
## The DLPFC dataset comes with a per spot manual annotation by the original 
## authors:
## Maynard et al., 2020 (https://www.nature.com/articles/s41593-020-00787-0)
g_truth_all <- readRDS("./data/rObjects/DLPFC_151673_gTruth.rds")

## Load the sf object for locations ----
## An sf data frame that gathers different types of metadata for each spot and 
## also contains two geo-columns that enclose spatial information. The dataset 
## loaded here is pre-filtered to match the ground truth by keeping only spots
## that were also manually annotated by the authors of the DLPFC data set.
polygons <- readRDS("./data/rObjects/DLPFC_151673_polygons_gTruth-filtered.rds")

## Load my custom ggplot theme (for maps) ----
my_theme_s <- readRDS("./data/rObjects/my_ggplot_theme_spatial.rds")

## Load my custom ggplot theme (for other plots) ----
my_theme <- readRDS("./data/rObjects/my_ggplot_theme.rds")

## Work with 1/3 of locations ----
polygons_part <- polygons %>% 
  filter(polygons$pixel_x < 250)

obs_part <- obs[,polygons_part$Barcode]

meta_D_part <- meta_D %>% 
  right_join(polygons_part[,c("Barcode")]) %>% 
  select(-geom_pol)

g_truth_part <- g_truth_all %>% 
  right_join(polygons_part[,c("Barcode")]) %>% 
  select(-geom_pol)
