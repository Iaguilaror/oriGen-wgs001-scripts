# load packs
pacman::p_load( "ggplot2", "patchwork", "dplyr", "cowplot" )

# load base plots ----

afr.p <- readRDS( file = "gnomad_AF_afr_scatter.rds" )

amr.p <- readRDS( file = "gnomad_AF_amr_scatter.rds" )

eas.p <- readRDS( file = "gnomad_AF_eas_scatter.rds" )

nfe.p <- readRDS( file = "gnomad_AF_nfe_scatter.rds" )

sas.p <- readRDS( file = "gnomad_AF_sas_scatter.rds" )

mcps_raw.p <- readRDS( file = "mcps_AF_RAW_scatter.rds" )

# crate themes ----
toprow_theme <- theme( plot.title = element_blank( ),
                       axis.title = element_blank( ),
                       plot.subtitle = element_blank( ),
                       plot.caption = element_blank( ),
                       legend.position = "none",
                       axis.title.x = element_blank( ),
                       # axis.ticks.x = element_blank( ),
                       axis.text.x = element_blank( ),
                       axis.text = element_text( size = 30 ) ) 

bottomrow_theme <- theme( plot.title = element_blank( ),
                          axis.title = element_blank( ),
                          plot.subtitle = element_blank( ),
                          plot.caption = element_blank( ),
                          legend.position = "none",
                          axis.title.x = element_blank( ),
                          axis.text = element_text( size = 30, hjust = 0.5 ) 
                          # axis.ticks.x = element_blank( ),
                          # axis.text.x = element_blank( )
) 

common_theme <- theme( panel.grid = element_blank( ) )
  
side_theme <- theme( axis.text.y = element_blank( )  )

# params -----
pop_text_size <- 15
pop_x_pos <- 1

# Build top row first ----

# 
toprow.p <- ( eas.p + 
                # annotate( "text", label = "EAS",
                #                 x = pop_x_pos, y = 0.05,
                #                 hjust = 1, size = pop_text_size ) +
                toprow_theme + common_theme +
                afr.p + 
                # annotate( "text", label = "AFR",
                #                   x = pop_x_pos, y = 0.05,
                #                   hjust = 1, size = pop_text_size ) +
                toprow_theme + side_theme + common_theme +
                nfe.p + 
                # annotate( "text", label = "NF-EUR",
                #                   x = pop_x_pos, y = 0.05,
                #                   hjust = 1, size = pop_text_size ) +
                toprow_theme + side_theme + common_theme )

bottomrow.p <- ( sas.p  + 
                   # annotate( "text", label = "SAS",
                   #                  x = pop_x_pos, y = 0.05,
                   #                  hjust = 1, size = pop_text_size ) +
                   bottomrow_theme + common_theme +
                   amr.p + 
                   # annotate( "text", label = "AMR",
                   #                   x = pop_x_pos, y = 0.05,
                   #                   hjust = 1, size = pop_text_size ) +
                   bottomrow_theme + side_theme + common_theme +
                   mcps_raw.p + 
                   # annotate( "text", label = "MCPS",
                   #                        x = pop_x_pos, y = 0.05,
                   #                        hjust = 1, size = pop_text_size ) +
                   bottomrow_theme + side_theme + common_theme )

full_panel.p <- toprow.p / bottomrow.p

# full_panel.p

ggsave( filename = "full_panel.svg", plot = full_panel.p, height = 5, width = 10 )

#
# Extract only the legend
legend.p <- get_legend( mcps_raw.p +
                          theme( legend.background = element_rect( fill = NA,
                                                                   color = NA ) ) ) %>% 
  plot_grid( )

ggsave(  filename = "legend.svg", plot = legend.p, height = 5, width = 5  )
