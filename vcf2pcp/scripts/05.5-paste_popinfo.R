pacman::p_load( "vroom", "dplyr", "tidyr" )

id_order <- vroom( file = "prob_and_ref.fam", delim = "\t", col_names = FALSE )

popinfo <- vroom( file = "allwgs.popinfo.txt", delim = "\t", col_names = FALSE )

popinfo_ordered <- popinfo %>%
  # filter( X2 %in% id_order$X1 ) %>%
  slice( match( id_order$X1, X2 ) )

write.table( x = popinfo_ordered, file = "pasted_allwgs.popinfo.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = F )