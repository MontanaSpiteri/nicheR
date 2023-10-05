#' Identify Spatial Niches
#'
#' Find niches of commonly co-occurring cell types in spatial transcriptomics or proteomics data.
#'
#' @param spe A SpatialExperiment object. Polygon/x and y spatial co-ordinates, annotation and cell barcode information stored in metadata.
#' @param fov A character string specifying the name of the fields of view column in spe (Default: 'fov').
#' @param anno A character string specifying the name of the annotation column in spe (Default: 'Anno').
#' @param x_coord A character string specifying the name of the column containing x spatial co-ordinates (Default: 'x_local_px').
#' @param y_coord A character string specifying the name of the column containing y spatial co-ordinates (Default: 'y_local_px').
#' @param tile_height A vector specifying tile height (Default: 1000).
#' @param tile_shift A vector specifying tile shift (Default: 100).
#' @param image_size Vectors specifying image size (Default: c(4400, 4400)).
#' @param num_clusters A vector specifying the number of clusters (Default: 9).
#' @param cores A vector specifying the number of cores (Default: 10).
#' @param prob A logical value indicating if the user in interested in calculating probability (Default: FALSE).
#'
#' @import mclust
#' @import parallel
#' @import magrittr
#' @import dplyr
#' @importFrom rDist parDist
#'
#' @return Dataframe identifying spatial niches assignment or when `prob=TRUE` the probability of each niche for each tile.
#' @export
#'
#' @details `find_niches` uses the cell type positions and provided cell type annotation as input to identify niches of
#' commonly co-occurring cell types. `find_niches` works using a window approach. First each FOV, if applicable, or the entire region
#' are divided into overlapping windows of size `tile_height` times `tile_height`. These windows shifted by `tile_shift` in x and y
#' direction resulting in overlaps. In each window, we count the number of different cell types.
#' Using a Bray-Curtis dissimilarity implemented a hierarchical clustering tree is built, which can be used to determine `num_clusters`
#' niches. For each `tile_shift` by `tile_shift` tile of each region is overlapped by multiple windows with
#' an associated niche label used to determine the probability of each tile to be a member of a certain niche or using
#' a majority voting strategy to identify a niche label. Note that edges of FOVs or entire regions are expanded to avoid edge effects. Hence, for
#' FOVs that are adjacent treat them like one continuous region.
#'
#' @examples
find_niches <- function(spe, fov="fov", anno="Anno", x_coord="x_local_px", y_coord="y_local_px", tile_height=1000, tile_shift=100,
                        image_size=c(4400, 4400), num_clusters=9, cores=10, prob=FALSE) {

    df_polygons <- as.data.frame(colData(spe)[, c(fov, x_coord, y_coord, anno)])

    start=0-(tile_height-tile_shift)
    image_size=image_size+tile_height-tile_shift

    all_tiles <- expand.grid(seq(start,image_size[1],tile_shift),
                             seq(start,image_size[2],tile_shift))
    all_tiles_list <- as.list(as.data.frame(t(all_tiles)))

    index <- split(1:nrow(df_polygons), df_polygons[,fov])
    df_polygons_split <- lapply(index, function(x) df_polygons[x, ])

    cells_in_tiles <- lapply(df_polygons_split, function(df) mclapply(all_tiles_list, function(x)
        identify_cells_within_tile(x, df, tile_height, x_coord, y_coord, anno)[,2], mc.cores=cores))
    # per fov find all cells within the defined tile (starting point defined in all_tiles)
    cells_in_tiles <- unlist(cells_in_tiles, recursive = FALSE)
    cells_in_tiles <- do.call(rbind, cells_in_tiles)

    colnames(cells_in_tiles) <- levels(df_polygons[,anno])

    all_tiles <- data.frame(x=rep(all_tiles[,1], length(index)), y=rep(all_tiles[,2], length(index)),
                            fov=rep(names(index), each=nrow(all_tiles)))
    index_remove <- rowSums(cells_in_tiles) < 3
    all_tiles <- all_tiles[!index_remove, ]
    cells_in_tiles <- cells_in_tiles[!index_remove, ]
    rownames(cells_in_tiles) <- paste0(all_tiles[,1], "_", all_tiles[,2], "_", all_tiles[,3])
    # find name for each tile

    dist_tiles <- parDist(cells_in_tiles, method = 'bray', threads = cores)
    # calculate distance between tiles
    clusters <- hclust(dist_tiles, method="ward.D2")
    clusters <- cutree(clusters, k = num_clusters)
    # cluster tiles
    clusters_df <- as.data.frame(all_tiles)
    clusters_df$cluster <- as.factor(clusters)
    colnames(clusters_df) <- c("x", "y", "fov", "cluster")
    # make data frame for each tile starting point and cluster

    all_tiles$x <- as.numeric(all_tiles$x)
    all_tiles$y <- as.numeric(all_tiles$y)
    index <- split(1:nrow(all_tiles), all_tiles$fov)
    all_tiles_list <- lapply(index, function(x) all_tiles[x, ])
    all_tiles_list <- lapply(all_tiles_list, function(x) as.list(as.data.frame(t(x))))
    # split all tiles by fov
    clusters_df <- lapply(index, function(x) clusters_df [x, ])
    # split dataframe for each tile starting point and cluster by fov

    cluster_in_shift <- mapply(X=all_tiles_list, Y=clusters_df, function(X,Y) mclapply(X, function(x)
        average_clusters_within_shift(x, Y, tile_shift, tile_height)[,2], mc.cores=cores))
    # find all tiles that overlap starting point and collect clusters do so by fov separately
    cluster_in_shift_x <- lapply(cluster_in_shift, function(x) do.call(rbind, x))

    if(prob){

        final_cluster <- lapply(cluster_in_shift_x, function(y)
            t(apply(y, 1, function(x) x/sum(x))))
        # identify most common cluster for the tile do so by fov separately
        all_tiles <- lapply(all_tiles_list, function(x) do.call(rbind, x))
        df_final <- mapply(X=all_tiles, Y=final_cluster, function(X,Y)
            cbind(data.frame(x=as.numeric(X[,1]),
                             y=as.numeric(X[,2]),
                             fov=X[,3], id=rownames(Y)), Y), SIMPLIFY = FALSE)
        df_final <- do.call(rbind, df_final)


    } else {

        final_cluster <- lapply(cluster_in_shift_x, function(y)
            apply(y, 1, function(x) c(1:num_clusters)[which.max(x)]))
        # identify most common cluster for the tile do so by fov separately
        all_tiles <- lapply(all_tiles_list, function(x) do.call(rbind, x))
        df_final <- mapply(X=all_tiles, Y=final_cluster, function(X,Y)
            data.frame(x=as.numeric(X[,1]),
                       y=as.numeric(X[,2]),
                       fov=X[,3],
                       Cluster=as.factor(Y), id=names(Y)), SIMPLIFY = FALSE)
        df_final <- do.call(rbind, df_final)

    }

    return(df_final)

}
