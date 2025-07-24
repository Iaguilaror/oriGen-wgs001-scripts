# Load pkgs
pacman::p_load( "vroom", "dplyr", "tidyr", "ggplot2", "pheatmap", "tibble" )

### prepare Functions
# Define the function for creating the heatmap
create_heatmap <- function( matrix_data, 
                            main_title = "Heatmap", 
                            x_axis_title = "Samples", 
                            y_axis_title = "Genes", 
                            legend_title = "Legend", 
                            palette_colors = c( "lightblue", "white", "lightpink" ),
                            num_colors = 10 ) {
  
  # Create a pastel color-blind-friendly palette
  mycolors.v <- colorRampPalette( palette_colors )( num_colors )
  
  # Generate the heatmap
  heatmap <- pheatmap( matrix_data,
                       scale = "none",
                       cluster_rows = TRUE,
                       cluster_cols = TRUE,
                       color = mycolors.v,
                       main = main_title,
                       show_rownames = FALSE, # Remove gene names
                       show_colnames = FALSE, # Remove sample names
                       legend = TRUE )        # Enable the legend
  
  # Add axis titles (manually with grid graphics)
  grid::grid.text( y_axis_title, x = 0.02, y = 0.5,
                   rot = 90, gp = grid::gpar( fontsize = 14 ) ) # Y-axis title
  
  # Add legend title (custom using grid graphics)
  grid::grid.text( legend_title, x = 0.98, y = 0.4,
                   rot = 90, gp = grid::gpar( fontsize = 12, fontface = "bold" ) )
  
  # Return the heatmap object
  return( heatmap )
}

# load all affected genes
all_data.df <- vroom( file = "allsamples_affected_genes_dataframe.tsv" )

# Si no hay evidencia de DUP, vamos a llenar el valor con 1 ( asumimos que la profundidad fue 1 en referencia (ni mas ni menos) )
# dup df
dup.df <- all_data.df %>% 
  pivot_wider( id_cols = Gene,
               names_from = "sample",
               values_from = "dup" )

dup.df[ is.na( dup.df ) ] <- 1

# Si no hay evidencia de DUP, vamos a llenar el valor con 1 ( asumimos que la profundidad fue 1 en referencia (ni mas ni menos) )
# del df
del.df <- all_data.df %>% 
  pivot_wider( id_cols = Gene,
               names_from = "sample",
               values_from = "del" )

del.df[ is.na( del.df ) ] <- 1

del.df

# do heatmap
# Prepare the data for pheatmap
del.matrix <- del.df %>%
  column_to_rownames("Gene") %>% # Set the 'Gene' column as row names
  as.matrix()                   # Convert the tibble to a matrix

# Create Del matrix
# Open a PNG device to save the plot
png( filename = "DEL_heatmap.png", width = 800, height = 800 )

create_heatmap( matrix_data = del.matrix,
                main_title = "DEL matrix",
                x_axis_title = "Samples",
                y_axis_title = "Genes",
                legend_title = "CNVpytor rd",
                palette_colors = c( "lightblue", "white", "lightpink" ), num_colors = 10 )


dev.off( )

# Create DUP matrix
# Prepare the data for pheatmap
dup.matrix <- dup.df %>%
  column_to_rownames("Gene") %>% # Set the 'Gene' column as row names
  as.matrix()                   # Convert the tibble to a matrix

# Create Del matrix
# Open a PNG device to save the plot
png( filename = "DUP_heatmap.png", width = 800, height = 800 )

create_heatmap( matrix_data = dup.matrix,
                main_title = "DUP matrix",
                x_axis_title = "Samples",
                y_axis_title = "Genes",
                legend_title = "CNVpytor rd",
                palette_colors = c( "lightblue", "white", "lightpink" ), num_colors = 10 )


dev.off( )