#' 
#' Extract the Clusters Identified in a CLV Analysis of a Spectra Object Composed of NMR Data
#' 
#' This function extracts the clusters of frequencies identified in a CLV analysis
#' of a \code{Spectra} object containing NMR data.  The clusters of frequencies
#' are optionally presented as a dotplot for visual comparison.
#' Clusters of frequences (peaks) may be
#' related by being part of the same spin system/compound,
#' or depending upon the experimental design, the peaks
#' identified may be related because they arise from
#' separate compounds which co-vary because they are on
#' related metabolic pathways.
#'
#' @param spectraORIG An object of class \code{Spectra}.  The data to be plotted.
#' Should not be shrink-wrapped.
#' 
#' @param spectraSW An object of class \code{Spectra}. The source of the frequencies
#' of interest.  Should be shrink-wrapped.  See \code{\link{shrinkwrapSpectra}}.
#' 
#' @param clv An object of class \code{"clv"}.  Created by
#' \code{\link{CLV}} or \code{\link{CLV_kmeans}} in package \code{ClustVarLV}.
#' In this case, created via \code{\link{clvSpectra}}.
#'
#' @param whichP Integer.  The desired partition. See Details.
#'
#' @param plot Logical.  Should a plot be made?
#'
#'
#' @param ... Additional arguments to be passed to the plotting routines.
#' \code{xlim} would be a typical example, since you usually need to zoom
#' in to see what's going on.
#'
#' @return A list with elements equal to the number of clusters. Each element
#' contains the frequencies associated with that cluster.  Side effect is a plot.
#'
#' @seealso \code{\link{plotSpectraCLV}}
#'
#' @section Details:
#' Objects of class \code{"clv"} contain different
#' entries if they were created by function \code{CLV} or \code{CLV_kmeans}.
#' For \code{mode = "clv"}
#' a number of partitions are returned, while for 
#' \code{mode = "kmeans"} only a single partition is returned.
#' In the latter case argument
#' \code{whichP} is ignored.
#' 
#' 
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @references E. Vigneau and E. Qannari. "Clustering of variables around latent components."
#' \emph{Communications in Statistics - Simulation and Computation}, 32(4):1131-1150, 2003. 
#' For a really good introduction, see E. Vigneau, M. Chen and E. Qannari
#' "ClustVarLV: An R Package for the Clustering of Variables Around Latent 
#' Variables" \emph{The R Journal} Vol. 7/2, December 2015 
#'
#' @keywords multivariate
#' @export
#' @importFrom reshape2 melt
#' @importFrom lattice dotplot panel.dotplot panel.abline
#' @importFrom ChemoSpec chkSpectra
#'
#' @examples
#'
#' require("ClustVarLV")
#' require("ChemoSpec")
#' data(metMUD2)
#'
#' # Shrinkwrap the spectra
#' sw <- shrinkwrapSpectra(spectra = metMUD2, sigma = 0.01,
#' thres = 4.85, xlim = c(1.25, 1.74), plot = FALSE)
#'
#' # Now CLV analysis
#' clv1 <- clvSpectra(sw$Spectra)
#' plot(clv1, type = "delta") # suggests 4 clusters
#' freqs <- extractCLVclusters(metMUD2, sw$Spectra, clv = clv1, whichP = 4,
#' main = "Peak Clusters in metMUD2")
#'

extractCLVclusters <- function(spectraORIG, spectraSW, clv = NULL, whichP = 1, plot = TRUE, ...) {

	if (!requireNamespace("ClustVarLV", quietly = TRUE)) {
		stop("You need to install package ClustVarLV to use this function")
		}	
	if (missing(spectraORIG)) stop("No spectral data provided")
	if (missing(spectraSW)) stop("No shrinkwrapped spectral data provided")
	if (missing(clv)) stop("No clv object provided")
	chkSpectra(spectraORIG)
	chkSpectra(spectraSW)
	
	# Map the freqs of interest back onto the original spectrum
	part <- get_partition(clv, K = whichP)
	nc <- length(unique(part))
	res <- vector("list", nc)
	names(res) <- paste("Cluster", 1:nc, sep = " ")
	for (i in 1:nc) {
		keep <- which(part == i)
		res[[i]] <- unique(round(spectraSW$freq[keep], digits = 2))
		}
	
	if (plot) {
		# Take the list and fill in w/NA and make a data frame
		ml <- max(lengths(res))
		M <- matrix(NA_real_, nrow = nc, ncol = ml)
		for (i in 1:nc) {
			cl <- length(res[[i]])
			M[i,] <- c(res[[i]], rep(NA_real_, ml - cl))
			}
		DF <- as.data.frame(M)
		DF$cluster <- names(res)
		DFm <- melt(DF, id.vars = "cluster")
		vlines <- seq(from = min(DFm$value, na.rm = TRUE), to = max(DFm$value, na.rm = TRUE), by = 0.1)
		p <- dotplot(cluster ~ value, xlab = "ppm", cex = 0.5, data = DFm, ...,
			panel = function(...) {
				panel.abline(v = vlines, col = "gray90", ...)
				panel.dotplot(...)
				})
		plot(p)
		}
	
	return(res)
	}
