pacman::p_load( "patchwork", "ggplot2", "cowplot", "dplyr" )

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

side_theme <- theme( axis.text.y = element_blank( )  )

# params -----
pop_text_size <- 15
pop_x_pos <- 1

# Build top row first ----

eas.p <- readRDS( file = "test/results/06-b-af_plot/gnomad_AF_eas_scatter.rds" )

afr.p <- readRDS( file = "test/results/06-b-af_plot/gnomad_AF_afr_scatter.rds" )

nfe.p <- readRDS( file = "test/results/06-b-af_plot/gnomad_AF_nfe_scatter.rds" )

# 
toprow.p <- ( eas.p + 
                # annotate( "text", label = "EAS",
                #                 x = pop_x_pos, y = 0.05,
                #                 hjust = 1, size = pop_text_size ) +
                toprow_theme +
                afr.p + 
                # annotate( "text", label = "AFR",
                #                   x = pop_x_pos, y = 0.05,
                #                   hjust = 1, size = pop_text_size ) +
                toprow_theme + side_theme +
                nfe.p + 
                # annotate( "text", label = "NF-EUR",
                #                   x = pop_x_pos, y = 0.05,
                #                   hjust = 1, size = pop_text_size ) +
                toprow_theme + side_theme )

ggsave( filename = "toprow.png", plot = toprow.p, height = 7, width = 21 )
ggsave( filename = "toprow.svg", plot = toprow.p, height = 7, width = 21 )

rm( list = c( "eas.p", "afr.p", "nfe.p", "toprow.p", "toprow_theme" ) )
gc( )

# Bottom ----

amr.p <- readRDS( file = "test/results/06-b-af_plot/gnomad_AF_amr_scatter.rds" )

sas.p <- readRDS( file = "test/results/06-b-af_plot/gnomad_AF_sas_scatter.rds" )

mcps_raw.p <- readRDS( file = "test/results/06-b-af_plot/mcps_AF_RAW_scatter.rds" )

bottomrow.p <- ( sas.p  + 
                   # annotate( "text", label = "SAS",
                   #                  x = pop_x_pos, y = 0.05,
                   #                  hjust = 1, size = pop_text_size ) +
                   bottomrow_theme +
                   amr.p + 
                   # annotate( "text", label = "AMR",
                   #                   x = pop_x_pos, y = 0.05,
                   #                   hjust = 1, size = pop_text_size ) +
                   bottomrow_theme + side_theme +
                   mcps_raw.p + 
                   # annotate( "text", label = "MCPS",
                   #                        x = pop_x_pos, y = 0.05,
                   #                        hjust = 1, size = pop_text_size ) +
                   bottomrow_theme + side_theme )

ggsave( filename = "bottomrow.png", plot = bottomrow.p, height = 7, width = 21 )
ggsave( filename = "bottomrow.svg", plot = bottomrow.p, height = 7, width = 21 )

# Extract only the legend
legend.p <- get_legend( mcps_raw.p +
                          theme( legend.background = element_rect( fill = NA,
                                                                   color = NA ) ) ) %>% 
  plot_grid( )

ggsave(  filename = "legend.svg", plot = legend.p, height = 5, width = 5  )

# to test only, because we run out of RAM

# full_panel.p <- toprow.p / bottomrow.p
# 
# ggsave( filename = "full_panel.png", plot = full_panel.p, height = 14, width = 24 )
