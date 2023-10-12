#' Find Uniquely Expressed Genes per Niche
#'
#' @param spe A SpatialExperiment object with gene expression data stored in assays(spe). Ensure cell IDs are in colnames(spe).
#' @param niches_output A dataframe object from the find_niches output. Includes 'Cluster' column indicating niches.
#' @param assay A character vector specifying the name of the gene expression assay in spe (Default: 'counts').
#'
#' @import limma
#' @import SpatialExperiment
#'
#' @return List of top differentially expressed genes per niche (Clusters)
#' @export
#'
#' @examples
niche_expression <- function(spe, niches_output, assay="counts") {

    if (!(assay %in% names(assays(spe)))) {
        stop(paste("The specified assay:", assay, "is not found in the provided SpatialExperiment object."))
    }
    if (!"Cluster" %in% colnames(niches_output)) {
        stop("'Cluster' column not found in the provided niches_output.")
    }

    expression_data <- assays(spe)[[assay]]
    niches_output <- niches_output[!is.na(niches_output$id), ]
    id_order <- match(colnames(spe), niches_output$id)
    niches_output <- niches_output[!is.na(id_order), ]

    differential_results <- list()
    # for each cluster found in the find_niches function
    for (cluster in unique(niches_output$Cluster)) {

        # identify which cells are in the Cluster vs not (NonCluster)
        cluster_cells <- which(niches_output$Cluster == cluster)
        non_cluster_cells <- which(niches_output$Cluster != cluster)

        all_cells <- c(cluster_cells, non_cluster_cells)
        design <- model.matrix(~ 0 + factor(c(rep(1, length(cluster_cells)), rep(0, length(non_cluster_cells)))))
        colnames(design) <- c("Cluster", "NonCluster")

        combined_data <- expression_data[, all_cells]

        #fit linear model
        fit <- lmFit(combined_data, design)

        contrast_matrix <- makeContrasts(Cluster-NonCluster, levels=design)
        fit2 <- contrasts.fit(fit, contrast_matrix)
        fit2 <- eBayes(fit2)

        top_genes <- topTable(fit2, number=Inf)
        differential_results[[paste0("Cluster_", cluster)]] <- top_genes
    }

    return(differential_results)
}
