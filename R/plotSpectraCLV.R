#' 
#' Plot the Results of a CLV Analysis of a Spectra Object Composed of NMR Data
#' 
#' This function plots the frequencies identified in a CLV analysis
#' of a Spectra Object containing NMR data.  The frequencies corresponding
#' to a particular partition and cluster are highlighted in red.
#' 
#'
#' @param spectraORIG An object of class \code{Spectra}.  The data to be plotted.
#' Should not be shrink-wrapped.
#' 
#' @param spectraSW An object of class \code{Spectra}. The source of the frequencies
#' of interest.  Should be shrink-wrapped.  See \code{\link{shrinkwrapSpectra}}.
#' 
#' @param clv An object of class \code{"clv"}.  Created by
#' \code{\link{CLV}} or \code{\link{CLV_kmeans}} in package \code{ClustVarLV}.
#' In this case, most easily created via \code{\link{clvSpectra}}.
#'
#' @param whichP Integer.  The desired partition. See Details.
#'
#' @param whichC Integer.  The desired cluster.
#'
#' @param whichS Integer.  The spectrum in \code{spectraORIG} to be plotted.
#' Can also be \code{"colSums"} in which case the sum of the spectra in
#' \code{spectraORIG} is plotted.
#'
#' @param ... Additional arguments to be passed to the plotting routines.
#' \code{xlim} would be a typical example, since you usually need to zoom
#' in to see what's going on.
#'
#' @return None. Side effect is a plot.
#'
#' @section Details:
#' Objects of class \code{"clv"} contain different
#' entries if they were created by function \code{CLV} or \code{CLV_kmeans}.
#' If \code{clv} was created by \code{CLV_kmeans} then argument
#' \code{whichP} is ignored, as there is only one partition.
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
#' @importFrom ClustVarLV get_partition
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
#'	thres = 4.85, xlim = c(1.25, 1.74), plot = FALSE)
#'
#' # CLV analysis
#' clv1 <- clvSpectra(sw$Spectra)
#' plot(clv1, type = "delta") # suggests 4 clusters
#'
#' # Compare these plots (other regions are interesting too)
#' plotSpectraCLV(metMUD2, sw$Spectra, clv = clv1, whichP = 4, whichC = 1, xlim = c(0.8, 2.1),
#' main = "Peak Cluster 1 from metMUD2")
#' plotSpectraCLV(metMUD2, sw$Spectra, clv = clv1, whichP = 4, whichC = 2, xlim = c(0.8, 2.1),
#' main = "Peak Cluster 2 from metMUD2")
#' plotSpectraCLV(metMUD2, sw$Spectra, clv = clv1, whichP = 4, whichC = 3, xlim = c(0.8, 2.1),
#' main = "Peak Cluster 3 from metMUD2")
#' plotSpectraCLV(metMUD2, sw$Spectra, clv = clv1, whichP = 4, whichC = 4, xlim = c(0.8, 2.1),
#' main = "Peak Cluster 4 from metMUD2")
#'

plotSpectraCLV <- function(spectraORIG, spectraSW, clv = NULL, whichP = 1, whichC = 1, whichS = 1, ...) {

	if (!requireNamespace("ClustVarLV", quietly = TRUE)) {
		stop("You need to install package ClustVarLV to use this function")
		}	
	if (missing(spectraORIG)) stop("No spectral data provided")
	if (missing(spectraSW)) stop("No shrinkwrapped spectral data provided")
	if (missing(clv)) stop("No clv object provided")
	chkSpectra(spectraORIG)
	chkSpectra(spectraSW)
	
	if (whichC > whichP) stop("whichC must be less than whichP.")
	
	# Map the freqs of interest back onto the original spectrum
	part <- get_partition(clv, K = whichP)
	keep <- which(part == whichC)
	keepf <- spectraSW$freq[keep]
	myc <- spectraORIG$freq %in% keepf
	myc <- ifelse(myc == TRUE, "red", "black")

	# set up the plot
	np <- length(spectraORIG$freq)
	ind1 <- 1:(np-1)
	ind2 <- 2:np
	
	if (whichS == "colSums") spec <- colSums(spectraORIG$data)
	if (!whichS == "colSums") spec <- spectraORIG$data[whichS,]

	plot(spectraORIG$freq, spec, type = "n",
		xlab = spectraORIG$unit[1], ylab = spectraORIG$unit[2], ...)
	segments(spectraORIG$freq[ind1], spec[ind1], spectraORIG$freq[ind2], spec[ind2],  col = myc)
	}

