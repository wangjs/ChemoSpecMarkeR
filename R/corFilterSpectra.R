#' @title Filter a Spectra Object Composed of NMR Data
#' 
#' @description
#' This function computes the correlation matrix of
#' the intensities in the \code{Spectra} object.  Any 
#' pair of peaks with a correlation below the specified threshold
#' (argument \code{lvl}) 
#' are removed from the \code{Spectra} object, giving
#' a smaller object which is composed of only highly
#' positively correlated peaks. 
#'
#' @param spectra An object of class \code{Spectra}.
#' Should be shrinkwrapped, see the examples.
#' 
#' @param lvl Numeric. The correlation level.  One value in (0...1).
#'
#' @return An object of class \code{Spectra}.
#' 
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @references E. Vigneau and E. Qannari. "Clustering of variables around latent components."
#' \emph{Communications in Statistics - Simulation and Computation}, 32(4):1131-1150, 2003. 
#' For a good introduction, see E. Vigneau, M. Chen and E. Qannari
#' "ClustVarLV: An R Package for the Clustering of Variables Around Latent 
#' Variables" \emph{The R Journal} Vol. 7/2, December 2015 
#'
#' @keywords multivariate
#' @export
#' @importFrom ChemoSpec chkSpectra
#' @importFrom stats cor
#'
#' @examples
#'
#' require("ClustVarLV")
#' require("ChemoSpec")
#' data(metMUD2)
#'
#'sw <- shrinkwrapSpectra(spectra = metMUD2, sigma = 0.02,
#'	thres = 4.5, xlim = c(1.25, 1.74), xlab = "ppm", plot = FALSE,
#'	main = expression(Shrinkwrap~of~metMUD2~sigma == 0.02))
#'
#'tst <- corFilterSpectra(sw$Spectra)
#'length(sw$Spectra$freq) # original no. of frequencies
#'length(tst$freq) # after filtering
#'
corFilterSpectra <- function(spectra, lvl = 0.9) {
	
	chkSpectra(spectra)
	if (length(lvl) != 1) stop("lvl must be a single value")
	if ((lvl < 0) | (lvl > 1)) stop("lvl should be in (0...1)")
	
	# Compute cor matrix and then apply threshold
	roc <- cor(spectra$data)
	roc <- apply(roc, 1, function(x) {ifelse(x > lvl, x, NA)})
	diag(roc) <- NA
	
	# Determine which freqs to keep
	keep <- rep(FALSE, nrow(roc))
	for (i in 1:nrow(roc)) {
		if (is.na(any(roc[i,] > lvl))) next	# entire row is NA
		if (any(roc[i,] > lvl)) keep[i] <- TRUE		
		}

	# Update the spectra object
	spectra$freq <- spectra$freq[keep]
	spectra$data <- spectra$data[,keep]
	
	chkSpectra(spectra)
	return(spectra)
	}
