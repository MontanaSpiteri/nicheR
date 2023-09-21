#' Identify Cells Within Tile
#'
#' @param start_tile
#' @param df_data
#' @param tile_height
#' @param x_coord
#' @param y_coord
#'
#' @return
#' @export
#'
#' @examples
identify_cells_within_tile <- function(start_tile, df, tile_height, x_coord, y_coord){

    tile <- rbind(start_tile, start_tile+tile_height)

    df %>% dplyr::filter(!!sym(x_coord) >= tile[1,1] & !!sym(x_coord) < tile[2,1] &
                             !!sym(y_coord) >= tile[1,2] & !!sym(y_coord) < tile[2,2]) %>%
        dplyr::count(Anno, .drop=FALSE)
}


