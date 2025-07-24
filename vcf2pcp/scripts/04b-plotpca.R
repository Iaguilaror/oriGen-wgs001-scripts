## Load libraries
pacman::p_load( "vroom", "dplyr", "tidyr", "ggplot2", "svglite", "cowplot", "scales", "stringr" )

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## For debugging only
# args[1] <- "allwgs.evec.gz" ## evec.gz file
# args[2] <- "allwgs.eval" ## eval file
# args[3] <- "allwgs.tracy_widom_statistics" ## tracy_widom_statistics file
# args[4] <- "smartpca.stdout" ## smartpca.stdout file
# args[5] <- "parallel_plot.svg" ## .parallel_plot.svg file
# args[6] <- "OriGen_Phase1_1Kg_and_AiGA.tsv" ## smartpca.stdout file

## Passing args to named objects
evec_file <- args[1]
eval_file <- args[2]
tw_file <- args[3]
smartpca_stdout_file <- args[4]
svg_file <- args[5]
reference_file <- args[6]

## Load data
rawdata.df <- vroom( file = evec_file, skip = 1, col_names = FALSE )

# get ncols
lastcol_n <- ncol( rawdata.df )

rawdata.df <- rawdata.df %>% 
  rename(sample = 1, lastcol = lastcol_n ) %>%
  rename_at( vars( 2:(lastcol_n -1) ), ~ paste0("PC", seq_along(.) ) )

## make a screeplot
##Load eval data data
raw_evaldata.df <- read.table(file = eval_file,
                              header = F,
                              sep = "\t", stringsAsFactors = F)

## add a column with numbers
raw_evaldata.df$component_number <- as.numeric(rownames(raw_evaldata.df))

#filter evalues of 0 or less
raw_evaldata.df <- raw_evaldata.df %>% filter(V1 > 0)

## identify last significant PC
## Should be able to load the table title to avoid missing the pvalue column due
# to changes in smartpca version
tw.df <- read.table(file = tw_file,
                    header = F,
                    sep = "\t", stringsAsFactors = F)

# Extract number of last significant PC
last_significant_pc <- tw.df %>%
  filter(V5 < 0.01) %>%
  select(V1) %>% max()

last_significant_pvalue <- tw.df %>%
  filter(V5 < 0.01) %>%
  select(V5) %>% max()

## create message to inform last significant PC and p-value
scree_title <- paste0("Scree plot - last significant PC is PC",
                      last_significant_pc,
                      "\nwith p-value =",
                      last_significant_pvalue)

## plot scree
scree.p <- raw_evaldata.df %>% 
  head( n = last_significant_pc ) %>% 
  ggplot( data = .,
          mapping = aes(x = component_number, y = V1, group = 1)) +
  geom_vline(xintercept = last_significant_pc,
             lty = "dashed", color = "orange4") +
  geom_line( color = "blue" , size = 0.5) +
  geom_point( shape=19, color = "red4", size = 1) +
  ggtitle(label = scree_title) +
  scale_y_continuous(name = "eigenvalues") +
  scale_x_continuous(breaks = 1:max(raw_evaldata.df$component_number),
                     labels = 1:max(raw_evaldata.df$component_number)) +
  theme_classic() #+
#  theme(axis.text.x = element_text(angle = 90, size = 5),
#        plot.title = element_text(size = 10))

## calculate % of variance for each PC
significant_pc.df <- tw.df %>%
  filter(V5 < 0.01) %>%
  select(V1, V2)

significant_pc.df$V1 <- paste0("PC", significant_pc.df$V1)
significant_pc.df$explained_variance_proportion <- significant_pc.df$V2 / sum(significant_pc.df$V2)
# calculate as percentage
significant_pc.df$variance_proportion_percent <- percent(significant_pc.df$explained_variance_proportion)
significant_pc.df$cumproportion <- cumsum(significant_pc.df$explained_variance_proportion)

## plot variance proportions
variance.p <- ggplot(data = significant_pc.df,
                     aes(x = V1,
                         y = explained_variance_proportion,
                         label = variance_proportion_percent) ) +
  geom_bar(stat = "identity", color = "black", fill = NA) +
  geom_text(size = 3,
            position = position_stack(vjust = 0.5)) +
  geom_line(aes(y = cumproportion, group = 1), color = "red4", alpha = 0.5) +
  geom_point(aes(y = cumproportion), color = "red4") +
  scale_y_continuous(limits = c(0,1)) +
  ggtitle(label = "Explained variance for each significant PC") +
  theme_bw() +
  theme(axis.title.x = element_blank())

## save significant pc dataframe
o_file <- gsub(pattern = ".svg", replacement = ".significant_pc.tsv", svg_file)
write.table(x = significant_pc.df, file = o_file,
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)

## process data for parallel coordinate plot
# remove first line, and remove last column
plotable.df <- rawdata.df

## create a vector with the PC number required
pc_series <- paste0("PC", 1:(ncol(plotable.df)-2) )

##rename cols
colnames(plotable.df) <- c("sample", pc_series, "tag")

# retag data by supercontinental groups
the_reference <- vroom( file = reference_file )

# 
plotable.df <-  plotable.df %>% 
  separate( col = sample,
            into = c( "project", "sample" ),
            sep = "-|:", remove = TRUE ) %>% 
  select( -project ) %>% 
  left_join( x =.,
             y = the_reference,
             by = "sample" ) %>% 
  mutate( tag = ifelse( test = `Population code` == "ORI",
                        yes = `Population code`,
                        no = `Superpopulation code` ) ) %>% 
  select( -`Superpopulation code`, -`Population code` )

## from wide to long format
long_plotable.df <- gather(plotable.df,
                           component_number, value,
                           pc_series, factor_key=F)

##transform PC to factor and order for correct plotting
long_plotable.df$component_number <- factor(long_plotable.df$component_number,
                                            levels = pc_series)

plot_title <- "Parallel Coordinate Plot"

## plot PCP
PCP.p <- ggplot(data = long_plotable.df,
                aes(x = component_number,
                    y = value, group = sample, color = tag)) +
  geom_line( alpha = 0.1 ) +
  ggtitle(label = plot_title) +
  theme(legend.position = "none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

## arrange grid
grid_below <- plot_grid(scree.p, variance.p,
                        nrow = 1, rel_widths = c(0.3,0.7))
grid.p <- plot_grid(PCP.p, grid_below, ncol = 1)

## save plot as svg
ggsave(filename = svg_file,
       plot = grid.p,
       device = "svg",
       width = 10, height = 7 , units = "in",
       dpi = 300)

## save plot as png
png_file <- gsub(pattern = ".svg", replacement = ".png", svg_file)
ggsave(filename = png_file,
       plot = grid.p,
       device = "png",
       width = 10, height = 7 , units = "in",
       dpi = 300)

##save the long dataframe
## create filename
o_file <- gsub(pattern = ".svg", replacement = ".tsv", svg_file)
write.table(x = long_plotable.df, file = o_file,
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)

##save the wide dataframe
## create filename
o_file <- gsub(pattern = ".svg", replacement = ".PCA_df.tsv", svg_file)
write.table(x = plotable.df, file = o_file,
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)

### Calculate PCP summarized
max_evec <- max(long_plotable.df$value)
min_evec <- min(long_plotable.df$value)

### Calculate mean of each region
means_per_region.df <- long_plotable.df %>%
  group_by(component_number, tag) %>%
  summarise( mean_per_region = mean(value),
             max_per_region = max(value),
             min_per_region = min(value))

means.p <- ggplot( data = means_per_region.df,
                   aes( x = component_number, y = mean_per_region,
                        group = tag, color = tag) ) +
  geom_line() +
  scale_y_continuous(limits = c(min_evec, max_evec)) +
  ggtitle(label = "Parallel Coordinate Plot") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

means_maxmin.p <- ggplot( data = means_per_region.df,
                          aes( x = component_number, y = mean_per_region,
                               group = tag, color = tag) ) +
  geom_line() +
  geom_line( aes( y = max_per_region),
             lty = "dashed", alpha = 0.3 ) +
  geom_line( aes( y = min_per_region),
             lty = "dashed", alpha = 0.3) +
  scale_y_continuous(limits = c(min_evec, max_evec)) +
  ggtitle(label = "Parallel Coordinate Plot") +
  theme(axis.text.x = element_text(angle = 90))

## save PCP plots in long format
# put into grido
grid_pcp.p <- plot_grid(means.p, means_maxmin.p, ncol = 1)

## save plot
ggsave(filename = svg_file,
       plot = grid_pcp.p,
       device = "svg", bg = "white",
       width = 10, height = 14 , units = "in",
       dpi = 300)

### Calculate many PC perspectives
## generate PC 1 vs 2
PCA.p <- ggplot(data = plotable.df, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = tag),
             size = 1, shape = 21) +
  theme_bw()

# as PC1 vs every other PC
# define main function
my_pca_function <- function( pca_data, pc_x, pc_y, x_name, y_name, the_tag) {
  ##pc_x and pc_y must be a string with the colname to plot
  tmp.df <- data.frame( X1 = pc_x,
                        X2 = pc_y,
                        X3 = the_tag )
  
  all_points.p <- ggplot(data = tmp.df,
                      mapping = aes( x = tmp.df[[1]],
                                     y = tmp.df[[2]],
                                     fill = tmp.df[[3]] ) ) +
    geom_hline( yintercept = 0, lty = "dashed" ) +
    geom_vline( xintercept = 0, lty = "dashed" ) +
    geom_point( alpha = 0.5, shape = 21 ) +
#    scale_x_continuous(name = x_name, limits = c(min_evec,max_evec)) +
#    scale_y_continuous(name = y_name, limits = c(min_evec,max_evec)) +
    scale_x_continuous(name = x_name) +
    scale_y_continuous(name = y_name) +
    theme_linedraw( )
  
  faceted.p <- all_points.p +  facet_wrap( ~ tag )
  
  panel_pca_.p <- plot_grid( all_points.p + theme( legend.position = "none" ), 
                            faceted.p, nrow = 1, rel_widths = c( 0.3, 0.7 ) )
  
}

##looping trough every PC after 1
# define starting PC col number
start_pc_coln <- 3
# define ending PC col number
end_pc_coln <- ncol(plotable.df) -1

# get a vector of posible PCs
the_PCs <- colnames( plotable.df[ , c(-1, -lastcol_n) ] )

# Create all possible combinations of two different elements
combinations <- combn(the_PCs, 2)

# Convert the combinations to a dataframe
combinations_df <- as.data.frame(t(combinations))

# Rename the columns
for (i in 1:nrow( combinations_df ) ) {
  
  PC_x_name <- combinations_df[i, 1]
  PC_y_name <- combinations_df[i, 2]
  message( paste( "[debug-R] plotting", PC_x_name, "vs", PC_y_name ) )
  
  ## generate plot
  my_plot.p <- my_pca_function( pc_x = plotable.df[ ,PC_x_name],
                                pc_y = plotable.df[ ,PC_y_name],
                                x_name = PC_x_name,
                                y_name = PC_y_name,
                                the_tag = plotable.df[ ,"tag"] )
  
  ## generate dynamic out file name
  #generate extension
  new_ext <- paste0(".",PC_x_name,"vs",PC_y_name,".svg")
  
  out_file <- gsub(pattern = "parallel_plot.svg",
                   replacement = paste0( "bipca", new_ext ),
                   svg_file)
  
  ## save plot
  ggsave(filename = out_file,
         plot = my_plot.p,
         device = "svg",
         width = 14, height = 7 , units = "in",
         dpi = 300)
}
