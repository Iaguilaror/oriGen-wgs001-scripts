## load libraries
library("dplyr")
library("ggplot2")
library("cowplot")

##get output file as the last element fro the args vecto
output_file <- "all_admixture.svg"

# Get all .rds file paths
rds_files <- list.files( path = ".",
                         pattern = "\\.rds$",
                         full.names = TRUE )

rds_files %>% 
  sort( decreasing = TRUE )

# Read all .rds files into a list
my_plot_list <- lapply( rds_files, readRDS )

## make a grid with every K
grid.p <- plot_grid( plotlist = my_plot_list,
                     ncol = 1 )

##save plot as svg
ggsave(filename = output_file,
       plot = grid.p,
       device = "svg",
       width = 14.4,
       height = 28.8,
       units = "cm", dpi = 300)