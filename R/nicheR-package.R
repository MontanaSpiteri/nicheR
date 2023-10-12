#' @importFrom parallelDist parDist
#' @import mclust
#' @import SpatialExperiment
#' @import mclust
#' @import dplyr
#' @import scater
#' @import scran
#' @import parallel
#' @import infotheo
NULL

#' Method to identify niches of co-occuring cell types from spatial
#' transcriptomics and proteomics data with single cell resolution.
#'
#' `nicheR` implements a convolutional network type approach to identify niches
#'  in spatial omics data with single cell resolution. `nicheR` works with the
#'  \link{SpatialExperiment::SpatialExperiment} class and requires cell type
#'  annotations as inputs.
#'
#' Key neighborhood analysis functions include \code{\link{find_niches}}, \code{\link{cell_to_tile}}
#'
#'
#'
#'
#' @author  Montana Spiteri \email{spiteri.m@wehi.edu.au}
#' @name nicheR
#' @docType package
#' @aliases nicheR nicheR-package
#' @keywords internal
#'
#'
#'
NULL
