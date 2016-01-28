#' Cluster Variables in a Spectra Object Composed of NMR Data
#' 
#' This function carries out clustering of variables around latent
#' variables.  Briefly, as used here the function finds NMR
#' frequencies that are related to each other.  Peaks may be
#' related by being part of the same spin system/compound,
#' or depending upon the experimental design, the peaks
#' identified may be related because they arise from
#' separate compounds which co-vary because they are on
#' related metabolic pathways.  The latent variables which drive the
#' clustering are the principal components.
#'
#' @param spectra An object of class \code{Spectra}.
#' Should be shrinkwrapped, see the examples.
#' 
#' @param mode Character. One of \code{c("clv", "kmeans")}
#'
#' @param method Character. One of \code{c("directional", "local")}.
#' \code{"directional"} identifies variables that are both positively
#' and negatively correlated with the associated latent variable.
#' \code{"local"} identifies variables based upon a positive
#' correlation with the latent variable.
#'
#' @param clust Only applies to \code{mode = "kmeans"}. Integer.  
#' The desired number of clusters.
#'
#' @param nmax Only applies to \code{mode = "clv"}. Integer.  
#' The maximum number of partitions to be considered.
#'
#' @return An object of class \code{"clv"}.  The contents of this
#' object depend upon \code{mode}.  For \code{mode = "clv"}
#' a number of partitions are returned, while for 
#' \code{mode = "kmeans"} only a single partition is returned.
#' 
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @references E. Vigneau and E. Qannari. "Clustering of variables around latent components."
#' \emph{Communications in Statistics - Simulation and Computation}, 32(4):1131-1150, 2003. 
#' For a good introduction, see E. Vigneau, M. Chen and E. Qannari
#' "ClustVarLV: An R Package for the Clustering of Variables Around Latent 
#' Variables" \emph{The R Journal} Vol. 7/2, December 2015 
#'
#' @seealso \code{\link{plotSpectraCLV}} or \code{\link{extractCLVclusters}}
#' for ways to plot the results.
#' 
#' @section Note:
#' \code{mode = "clv"} will be much slower than \code{mode = "kmeans"}.
#' For larger data sets it might be wise to run it in batch mode.
#'
#' @section Details:
#' \code{method = "local"} will identify variables that are positively
#' correlated with the associated latent variable.  In practice, this
#' means frequencies in a cluster all co-vary positively.  This means
#' they are either peaks from the same compound, or peaks from different
#' compounds that rise and fall in concentration with the experimental
#' treatment.  This method is best for identifying the compounds.
#' For \code{method = "directional"}, positive and negative
#' correlations are permitted.  The interpretation is as for \code{"local"}
#' but in addition, frequencies (peaks) will be included that arise
#' from compounds that go down in concentration when other compounds go
#' up in concentration, due to the experimental treatment.  This
#' method is more useful for fleshing out specific metabolic changes.
#'
#' @keywords multivariate
#' @export
#' @importFrom ClustVarLV CLV
#' @importFrom ClustVarLV CLV_kmeans
#' @importFrom ChemoSpec chkSpectra
#'
#' @examples
#'
#' require("ClustVarLV")
#' require("ChemoSpec")
#' data(metMUD2)
#'
#' # Shrinkwrap the spectra
#' myt <- "shrinkwrapSpectra Demo"
#' sw <- shrinkwrapSpectra(spectra = metMUD2, sigma = 0.01,
#'	thres = 4.85, xlim = c(1.25, 1.74), main = myt, plot = FALSE)
#'
#' # CLV analysis
#' clv1 <- clvSpectra(sw$Spectra)
#' plot(clv1, type = "delta") # suggests 4 clusters
#' clv2 <- clvSpectra(sw$Spectra, mode = "kmeans", clust = 4)
#' \dontrun{
#' # To see the differences between modes, compare the following:
#' names(clv1)
#' names(clv2)
#' str(clv1)
#' str(clv2)
#' }

clvSpectra <- function(spectra, mode = "clv", method = "local",
	clust = NULL, nmax = 20) {

	if (!requireNamespace("ClustVarLV", quietly = TRUE)) {
		stop("You need to install package ClustVarLV to use this function")
		}	
	if (missing(spectra)) stop("No spectral data provided")
	chkSpectra(spectra)
	
	if (mode == "clv") {
		res <- CLV(spectra$data, method = method, nmax = nmax)
		}
		
	if (mode == "kmeans") {
		if (is.null(clust)) stop("You must supply a value for clust.")
		res <- CLV_kmeans(spectra$data, method = method, clust = clust)
		}
	
	return(res)
	}

