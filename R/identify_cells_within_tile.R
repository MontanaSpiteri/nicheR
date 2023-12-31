#' Identify Cells Within Tile
#'
#'
identify_cells_within_tile <- function(start_tile, df, tile_height, x_coord, y_coord, anno){

    tile <- rbind(start_tile, start_tile+tile_height)

    df %>% dplyr::filter(.data[[x_coord]] >= tile[1,1] & .data[[x_coord]] < tile[2,1] &
                             .data[[y_coord]] >= tile[1,2] & .data[[y_coord]] < tile[2,2]) %>%
        dplyr::count(.data[[anno]], .drop=FALSE)

}
