#' Identify Cells Within Tile
#'
#' @param start_tile
#' @param df
#' @param tile_height
#'
#' @return
#' @export
#'
#' @examples
identify_cells_within_tile <- function(start_tile, df, tile_height){

    tile <- rbind(start_tile, start_tile+tile_height)

    df %>% dplyr::filter(x_local_px>=tile[1,1] & x_local_px<tile[2,1] &
                             y_local_px>=tile[1,2] & y_local_px<tile[2,2]) %>%
        dplyr::count(Anno, .drop=FALSE)

}
