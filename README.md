# Mapping-the-Transcriptome

Code and data for the paper: "Mapping the Transcriptome - realising the full potential of spatial transcriptomics data analysis"

***

## Contents

This repository contains all the required files to recreate the figures from the publication. Below is a detailed list of the repository's contents.

1. **data**: Contains the data used for the analysis and the graphical outputs.
   1. **graphics_out**: graphical outputs
      1. **illustrator_exports**: the final figures submitted.
      2. **illustrator_files**: the .ai files where we put together the figures.
      3. **selected_svg_MAUP**: the .svg files that where selected to be used in the MAUP figure, as outputed from RStudio.
      4. **selected_svg_SA**: the .svg files that where selected to be used in the SA figure as outputed from RStudio.
      5. **selected_svg_SH**: the .svg files that where selected to be used in the SH figure as outputed from RStudio.
   2. **rObjects**: data used for the analysis. A more detailed description of each file below can be found inside script: *load_data.R* located inside the *R* folder.
      1. DLPFC_151673_SeuratObject.rds
      2. DLPFC_151673_gTruth.rds
      3. DLPFC_151673_inputD.rds
      4. DLPFC_151673_inputMD.rds
      5. DLPFC_151673_most_variable_genes.rds
      6. DLPFC_151673_polygons_gTruth-filtered.rds
      7. my_ggplot_theme.rds
      8. my_ggplot_theme_spatial.rds
2. **R**: COntains the scripts used to run the analysis, explore the results and produce teh graphs that were subsequently used in the figures.
   1. **figure_MAUP.R**: code for generating graphs used in the MAUP figure.
   2. **figure_SA.R**: code for generating graphs used in the SA figure.
   3. **figure_SH.R**: code for generating graphs used in the SH figure.
   4. **load_data.R**: code to load the required data.
   5. **load_packages.R**: code to load the required packages.
   6. **moran.test.ste.R**: a modification of the *moran.test* function from the *GWmodel* package that does not round really low *p*-values to zero.
   7. **preprocess_data.R**: code that preprocesses the data before the analysis.
3. **README.md**: the readme file describing the repository's contents.
4. **figures.Rmd**: An *RMarkdown* file that sums up the code that produced the graphs that where finally selected and comments on the process.
5. **figures.Proj**: The R project file.
6. **figures.html**: the HTML as rendered from the *RMarkdown* file (*figures.Rmd*).

***

## Guidelines of use

If you are looking to recreate the figures we suggest you start by running the *figures.Rmd* file and then divie into the individual figure scripts (*figure_MAUP.R*, *figure_SA.R*, *figure_SH.R*) to further explore.
