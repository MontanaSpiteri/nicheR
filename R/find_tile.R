#' Find Tiles
#'
#' @param x
#' @param df_final_1
#' @param tile_shift
#'
#' @return
#' @export
#'
#' @examples
find_tile <- function(x, df_final_1, tile_shift) {

    which(df_final_1[,c("x")]<x[1] & x[1]<df_final_1[,c("x")]+tile_shift &
              x[2]> df_final_1[,c("y")] & x[2]<df_final_1[,c("y")]+tile_shift)

}
