#' This script is to pre-process the data that are going to be used to generate 
#' figures.
#' 
#' 

## -------------------------------------------------------------------------- ##
## Keep only locations that are present in the ground truth dataset too ----
## -------------------------------------------------------------------------- ##
obs_filt <- obs %>% 
  select(g_truth_all$Barcode) 

## -------------------------------------------------------------------------- ##
## Remove genes with low counts and in low number of locations ----
## -------------------------------------------------------------------------- ##
obs_filt <-  obs_filt %>% 
  mutate(rSums = rowSums(.),
         nLoc = rowSums(. != 0),
         average = rSums / nLoc) %>% 
  filter(rSums != 0) %>%
  filter(!(nLoc < 20 & average < 1.5)) %>% 
  dplyr::select(-rSums, -nLoc, -average)

## -------------------------------------------------------------------------- ##
## Normalise counts  using the standard Seurat method----
## -------------------------------------------------------------------------- ##
obs_seu_Nor <- obs_seu %>%
  .[rownames(obs_filt),] %>%
  ## normalise for cell to cell differences
  Seurat::NormalizeData(normalization.method = "LogNormalize",
                        scale.factor = 10000,
                        verbose = TRUE)

obs_seu_NorVst <- obs_seu_Nor %>%
  ## find genes to use for clustering
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

obs_seu_NorVstSc <- obs_seu_NorVst %>%
  ## scale genes so level of expression is the similar for all genes
  ScaleData(verbose = TRUE,
            vars.to.regress = c("detected",
                                "total"))

## -------------------------------------------------------------------------- ##
## Extract the normalised and variance stabilised counts from Seurat object ----
## -------------------------------------------------------------------------- ##
pattern <- grep("obs_seu_", names(.GlobalEnv), value = TRUE)

for (p in pattern) {
  dt <- get(p)
  
  out_name <- paste0("obs_", gsub("obs_seu_", "", p))
  
  if (grepl("Sc", out_name, fixed = TRUE)) {
    dt_counts <- dt@assays[["originalexp"]]@scale.data %>% 
      as.matrix()
  } else {
    dt_counts <- dt@assays[["originalexp"]]@data %>% 
      as.matrix()
  }
  
  spotNames <- dt$Barcode %>% 
    gsub("-", ".", .)
  
  colnames(dt_counts) <- spotNames[match(names(spotNames), 
                                            colnames(dt_counts))]
  
  assign(out_name, as.data.frame(dt_counts))
  
  rm(p, dt, out_name, dt_counts, spotNames)
}

## -------------------------------------------------------------------------- ##
## Extract the 2000 most variable genes from the Seurat object ----
## -------------------------------------------------------------------------- ##
var_genes <- obs_seu_NorVst@assays[["originalexp"]]@var.features

## -------------------------------------------------------------------------- ##
## Clean up the environment ----
## -------------------------------------------------------------------------- ##
## Remove unneeded Seurat objects
rm(list = pattern)

## Remove the Noramlised and Vst dataframe because the gene counts are the same
rm(obs_NorVst)

## Remove the filtered observations data frame that is not needed anymore
rm(obs_filt)
