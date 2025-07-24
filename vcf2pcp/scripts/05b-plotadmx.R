## load libraries
pacman::p_load( "vroom", "dplyr", "tidyr", "stringr",
                "ggplot2", "svglite", "cowplot", "ggsci", "RColorBrewer" )

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
# # Comment for production mode only
# args[1] <- "allwgs.5.Q" ## .Q file
# args[2] <- "allwgs.popinfo.txt" ## popinfo file
# args[3] <- "OriGen_Phase1_1Kg_and_AiGA.tsv"

## Passing args to named objects
q_file <- args[1]
pop_file <- args[2]
sample_list <- args[3]

# read proportions
q.df <- vroom( file = q_file, col_names = FALSE )

# get the k
k_value <- ncol( q.df )

# read the popinfo
pop.df <- vroom( file = pop_file, col_names = FALSE )

full_tag.df <- vroom( file = sample_list, col_names = TRUE )

retagged.df <- left_join( x = pop.df,
                          y = full_tag.df,
                          by = c( "X2" = "sample" ) ) %>% 
  mutate( X3 = ifelse( test = `Population code` == "ORI",
                       yes = "ORI",
                       no = X3 ) ) %>% 
  select( X1:X3 )

data.df <- bind_cols( retagged.df, q.df )

## generate the colname vector
#dynamically set of numbered groups
vector_of_groups <- paste0("group", 4:ncol(data.df) - 3)
names_for_cols <- c("plink_code","sample","region",vector_of_groups)
colnames(data.df) <- names_for_cols

# Save the tagged Q files
str_replace( string = q_file,
             pattern = ".Q",
             replacement = ".admixture_proportion.tsv" ) %>% 
write.table( x = data.df,
             file = .,
             append = FALSE, quote = FALSE,
             sep = "\t", row.names = FALSE,
             col.names = TRUE )

## pass from wide to long data
long.df <- data.df %>% pivot_longer(cols = 4:ncol(data.df),
                                    names_to = "group_k",
                                    values_to = "proportion")

# Create a vector of 10 pastel harmonic colors
# pastel_colors <- hue_pal()(10)
# pastel_colors <- alpha(pastel_colors, 0.5)  # Adjust for pastel effect



# create a manual pallete
# mycolors.v <- c( "red", "green", "blue", "yellow", "purple",
#                  "brown", "orange", "skyblue", "black", "gray" )
# mycolors.v <- pastel_colors
mycolors.v <- brewer.pal(10, "Set3")

names( mycolors.v ) <-  paste0( "group", 1:10 ) 

# transparency in colors
mycolors.v <- alpha(mycolors.v, 0.9) 

## Define a custom function for plotting
myplotting_function <- function(data, kval, region_name) {
  
  # count number of samples
  nsamples.tmp <- data %>% 
    pull( sample ) %>% 
    unique( ) %>% 
    length( )
  
  # find the most abundant group by cumm proportion
  rank_proportions.tmp <- data %>% 
    group_by( group_k ) %>% 
    summarise( sum_proportion = sum( proportion ) ) %>% 
    ungroup( ) %>% 
    arrange( -sum_proportion )
  
  # get the order of groups
  groups_ordered.v <- rank_proportions.tmp %>% 
    pull( group_k )
  
  # reorder samples by the ranked proportions
  wide.tmp <- data %>% 
    select( sample, group_k, proportion ) %>% 
    mutate( group_k = factor( group_k, levels = groups_ordered.v ) )
  
  ordered.tmp <- wide.tmp %>% 
    arrange( group_k, proportion ) %>% 
    pull( sample ) %>% 
    unique( )
  
  data %>% 
    mutate( group_k = factor( group_k, levels = rev( groups_ordered.v ) ) ) %>%
    ggplot( ., aes(x = sample, y = proportion)) +
    geom_col(aes(fill = group_k)) +
    scale_y_continuous(name = "Proportion", 
                       breaks = seq(0,1,by = 0.2),
                       labels = paste0(seq(0,1,by = 0.2) * 100, "%"),
                       expand = c(0, 0)) +
    scale_x_discrete( limits = ordered.tmp ) +
    ggtitle(label = region_name,
            subtitle = paste( "n =", nsamples.tmp ) ) +
    # scale_fill_npg( name = paste( "K=", kval ) ) +
    scale_fill_manual( name = paste("K=", kval),
                       limits = groups_ordered.v,
                       values = mycolors.v ) +
    theme_bw() +
    theme(
      text = element_text(size = 5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 5),
      axis.title.x = element_blank(),
      panel.background = element_blank(),
      legend.key.size = unit(0.3,"cm"),
      legend.box.margin=margin(-10,-10,-10,-10))
  
}

## define a vector with all available regions
allregions <- unique(long.df$region)

# reorder if AMI is present or not
if ("AMI" %in% allregions) {
  
  ordered_allregions.v <- c( "AFR", "SAS", "EAS", "AMI", "AMR", "ORI", "EUR" )
  
} else {
  
  ordered_allregions.v <- c( "AFR", "SAS", "EAS", "AMR", "ORI_REF", "ORI_problem", "EUR" )
  
}

# keep only regions pressent in the allregions
allregions <- ordered_allregions.v[ ordered_allregions.v %in% allregions ]

##create an empty list
my_plot_list <- list()

## plot every region
for (i in 1:length(allregions)) {
  region_in_turn <- allregions[i]
  # Debug message
  message(paste("[..] plotting", region_in_turn))
  #filter data for required region
  region.df <- long.df %>% 
    filter(region == region_in_turn)
  
  myplot.p <- myplotting_function( data = region.df, kval = k_value, region_name = region_in_turn)
  
  ## modify every panel but the last
  if (i != length(allregions)) {
    myplot.p <- myplot.p + theme(legend.position = "none",
                                 plot.margin = unit(c(0.1,0,0,0), "cm"))
  }
  
  ## ## modify every panel but the first
  if (i != 1) {
    myplot.p <- myplot.p + theme(axis.text.y = element_blank(),
                                 axis.title.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 plot.margin = unit(c(0.1,0,0,0), "cm")
    )
  }
  
  ## modify the last panel
  if (i == length(allregions)) {
    myplot.p <- myplot.p + theme(plot.margin = unit(c(0.1,0.5,0,0), "cm")
    )
  }
  
  ## pass the plot to plotlist
  my_plot_list[[i]] <- myplot.p
}

# calculate the number of samples
samples_by_number.df <- data.df %>% 
  filter( !is.na(region) ) %>% 
  group_by( region ) %>% 
  summarise( n = n( ) ) %>% 
  ungroup(  ) %>% 
  mutate( region = factor( region, levels = ordered_allregions.v ) ) %>% 
  arrange( region ) %>% 
  mutate( fraction = n / sum(n) ) %>% 
  mutate( rank = row_number( ) ) %>% 
  mutate( adj_width = ifelse( test = rank == max(rank),
                              yes = fraction + 0.08,
                              no = fraction - ( 0.08 / (max(rank) -1) ) ) )

## plot every plot together
allregions.p <- plot_grid(plotlist = my_plot_list, nrow = 1,
                          rel_widths = samples_by_number.df$adj_width )

##save plot as svg
str_replace( string = q_file,
             pattern = ".Q",
             replacement = ".admixture_plot.svg" ) %>% 
  ggsave(filename = .,
         plot = allregions.p,
         device = "svg",
         width = 14.4,
         height = 7.2,
         units = "cm", dpi = 300)

## save plot object for use in another R session
str_replace( string = q_file,
             pattern = ".Q",
             replacement = ".admixture_plot.rds" ) %>% 
  saveRDS( allregions.p, file = . )

## Diagnostic
# a plot showing all groups
##create an empty list
my_plot_list_diagnostic <- list()

## plot every region
for (i in 1:length(allregions)) {
  
  region_in_turn <- allregions[i]
  
  # Debug message
  message(paste("[..] plotting", region_in_turn))
  
  #filter data for required region
  region.df <- long.df %>% 
    filter(region == region_in_turn)
  
  myplot.p <- myplotting_function( data = region.df, kval = k_value, region_name = region_in_turn)
  
  myplot.p <- myplot.p +
    theme( legend.position = "bottom",
           legend.box = "horizontal" ) +
    guides( fill = guide_legend( nrow = 2 ) )
  
  ## pass the plot to plotlist
  my_plot_list_diagnostic[[i]] <- myplot.p
}

## plot every plot together
allregions_diagnostic.p <- plot_grid( plotlist = my_plot_list_diagnostic, nrow = 1 )

##save plot as svg
str_replace( string = q_file,
             pattern = ".Q",
             replacement = ".admixture_plot.diagnostic.svg" ) %>% 
  ggsave(filename = .,
         plot = allregions_diagnostic.p,
         device = "svg",
         width = 28,
         height = 10,
         units = "cm", dpi = 600)
