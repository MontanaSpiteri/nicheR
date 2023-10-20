#' Identify Spatial Niches
#'
#' Find niches of commonly co-occurring cell types in spatial transcriptomics or proteomics data.
#'
#' @param spe A SpatialExperiment object. Polygon centroid coordinates (when applicable), global and/or local x/y spatial co-ordinates, annotation and cell barcode information stored in column data.
#' @param fov Optional (Default: NULL). If applicable, a character string specifying the name of the fields of view column in spe.
#' @param anno A character string specifying the name of the annotation column in spe (Default: 'Anno').
#' @param x_coord A character string specifying the name of the column containing x spatial co-ordinates (Default: 'x'). If `fov=NULL`, these will be global co-ordinates, else, local co-ordinates.
#' @param y_coord A character string specifying the name of the column containing y spatial co-ordinates (Default: 'y'). If `fov=NULL`, these will be global co-ordinates, else, local co-ordinates.
#' @param tile_height A vector specifying tile height (Default: 1000).
#' @param tile_shift A vector specifying tile shift (Default: 100).
#' @param num_clusters A vector specifying the number of clusters (Default: 9).
#' @param cores A vector specifying the number of cores (Default: 10).
#' @param prob A logical value indicating if the user in interested in calculating probability (Default: FALSE).
#' @param n_permutations A vector specifying the number of permutations to use for null hypothesis testing (Default: 100).
#'
#' @import mclust
#' @import SpatialExperiment
#' @import parallel
#' @import dplyr
#' @importFrom parallelDist parDist
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
find_niches <- function(spe, fov=NULL, anno="Anno", x_coord="x", y_coord="y", tile_height=1000, tile_shift=100,
                        num_clusters=9, cores=10, prob=FALSE, n_permutations=100) {

    colData(spe)[,anno] <- as.factor(colData(spe)[, anno])
    # ensure anno is of class factor
    if (!is.null(fov)) {
        df_polygons <- as.data.frame(colData(spe)[, c(fov, x_coord, y_coord, anno)])
    } else {
        df_polygons <- as.data.frame(colData(spe)[, c(x_coord, y_coord, anno)])
    }

    x_range <- range(df_polygons[, x_coord])
    y_range <- range(df_polygons[, y_coord])
    width <- x_range[2] - x_range[1]
    height <- y_range[2] - y_range[1]
    image_size <- c(width, height)
    #estimate approximate image size

    start=0-(tile_height-tile_shift)
    image_size=image_size+tile_height-tile_shift

    all_tiles <- expand.grid(seq(start,image_size[1],tile_shift),
                             seq(start,image_size[2],tile_shift))
    all_tiles_list <- as.list(as.data.frame(t(all_tiles)))

    if (!is.null(fov)) {
        index <- split(1:nrow(df_polygons), df_polygons[,fov])
        df_polygons_split <- lapply(index, function(x) df_polygons[x, ])
    } else {
        df_polygons_split <- list(df_polygons)
    }
    # adjust how data is split depending on 'fov' availability

    cells_in_tiles <- lapply(df_polygons_split, function(df) mclapply(all_tiles_list, function(x)
        identify_cells_within_tile(x, df, tile_height, x_coord, y_coord, anno)[,2], mc.cores=cores))
    # per fov find all cells within the defined tile (starting point defined in all_tiles)
    cells_in_tiles <- unlist(cells_in_tiles, recursive = FALSE)
    cells_in_tiles <- do.call(rbind, cells_in_tiles)

    colnames(cells_in_tiles) <- levels(df_polygons[,anno])

    if (!is.null(fov)) {
        all_tiles <- data.frame(x=rep(all_tiles[,1], length(index)), y=rep(all_tiles[,2], length(index)),
                                fov=rep(names(index), each=nrow(all_tiles)))
    } else {
        all_tiles <- data.frame(x=all_tiles[,1], y=all_tiles[,2]) # 'fov' is absent, so we skip it
    }
    # construct the 'all_tiles' dataframe differently based on 'fov' availability

    index_remove <- rowSums(cells_in_tiles) < 3
    all_tiles <- all_tiles[!index_remove, ]
    cells_in_tiles <- cells_in_tiles[!index_remove, ]
    if (nrow(cells_in_tiles) == 0) {
        stop("No cells remain after filtering tiles with fewer than 3 cells. Please adjust your parameters or input data.")
    }

    if (!is.null(fov)) {
        rownames(cells_in_tiles) <- paste0(all_tiles[,1], "_", all_tiles[,2], "_", all_tiles[,3])
    } else {
        rownames(cells_in_tiles) <- paste0(all_tiles[,1], "_", all_tiles[,2]) # excludes 'fov'
    }
    # find name for each tile

    dist_tiles <- parDist(cells_in_tiles, method = 'bray', threads = cores)
    dist_matrix <- as.matrix(dist_tiles)
    # calculate distance between tiles
    clusters <- hclust(dist_tiles, method="ward.D2")
    clusters <- cutree(clusters, k = num_clusters)
    # cluster tiles
    clusters_df <- as.data.frame(all_tiles)
    clusters_df$cluster <- as.factor(clusters)
    if (!is.null(fov)) {
        colnames(clusters_df) <- c("x", "y", "fov", "cluster")
    } else {
        colnames(clusters_df) <- c("x", "y", "cluster") # excludes 'fov'
    }
    # make data frame for each tile starting point and cluster

    current_dist_matrix <- as.matrix(dist_tiles)
    p_value <- compute_pvalue(clusters_df, current_dist_matrix, n_permutations)
    #compute pvalues

    all_tiles$x <- as.numeric(all_tiles$x)
    all_tiles$y <- as.numeric(all_tiles$y)

    if (!is.null(fov)) {
        index <- split(1:nrow(all_tiles), all_tiles$fov)
        all_tiles_list <- lapply(index, function(x) all_tiles[x, ])
        all_tiles_list <- lapply(all_tiles_list, function(x) as.list(as.data.frame(t(x))))
        clusters_df <- lapply(index, function(x) clusters_df[x, ])
    } else {
        all_tiles_list <- list(all_tiles)
        clusters_df <- list(clusters_df)
    }
    # process 'all_tiles_list' and 'clusters_df' based on 'fov' availability

    cluster_in_shift <- mapply(X=all_tiles_list, Y=clusters_df, function(X,Y) mclapply(X, function(x)
        average_clusters_within_shift(x, Y, tile_shift, tile_height)[,2], mc.cores=cores))
    # find all tiles that overlap starting point and collect clusters do so by 'fov' separately - if applicable
    cluster_in_shift_x <- lapply(cluster_in_shift, function(x) do.call(rbind, x))

    if(prob){

        final_cluster <- lapply(cluster_in_shift_x, function(y)
            t(apply(y, 1, function(x) x/sum(x))))
        # identify most common cluster for the tile do so by 'fov' separately - if applicable
        all_tiles <- lapply(all_tiles_list, function(x) do.call(rbind, x))
        df_final <- mapply(X=all_tiles, Y=final_cluster, function(X,Y)
            cbind(data.frame(x=as.numeric(X[,1]),
                             y=as.numeric(X[,2]),
                             fov=X[,3], id=rownames(Y)), Y), SIMPLIFY = FALSE)
        df_final <- do.call(rbind, df_final)


    } else {

        final_cluster <- lapply(cluster_in_shift_x, function(y)
            apply(y, 1, function(x) c(1:num_clusters)[which.max(x)]))
        # identify most common cluster for the tile do so by 'fov' separately - if applicable
        all_tiles <- lapply(all_tiles_list, function(x) do.call(rbind, x))
        df_final <- mapply(X=all_tiles, Y=final_cluster, function(X,Y)
            data.frame(x=as.numeric(X[,1]),
                       y=as.numeric(X[,2]),
                       fov=X[,3],
                       Cluster=as.factor(Y), id=names(Y)), SIMPLIFY = FALSE)
        df_final <- do.call(rbind, df_final)

    }
    df_final$p_value <- p_value
    if (is.null(fov) && "fov" %in% names(df_final)) {
        df_final <- df_final[, !names(df_final) %in% "fov"]
    }
    # check if 'fov' exists in the final dataframe and remove if not present

    return(df_final)

}
