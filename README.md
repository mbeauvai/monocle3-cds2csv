# monocle3-cds2csv
A lightweight tool to convert Monocle3 CellDataSet (CDS) objects from RDS format to CSV files for easy data analysis and sharing.

## Installation

To use this script, ensure you have the following R libraries installed:

- `dplyr` (MIT License)
- `igraph` (GPL-2 | GPL-3 License)
- `SingleCellExperiment` (Artistic-2.0 License)
- `dynwrap` (MIT License)

### Install Required Libraries
Run the following commands in your R console to install the required libraries:

```R
install.packages("dplyr")  # Install from CRAN
install.packages("igraph") # Install from CRAN
BiocManager::install("SingleCellExperiment") # Install from Bioconductor
install.packages("dynwrap") # Install from CRAN
