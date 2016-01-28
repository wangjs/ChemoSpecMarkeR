#' Shrinkwrap a Spectra Object Composed of NMR Data
#' 
#' This function identifies the central part of an NMR peak using
#' a modification of the convolution product / second
#' derivative method described in the
#' reference.  The rest of the frequencies are then discarded,
#' resulting in a much smaller data set.
#' To ensure that all peaks in a data set are
#' represented, the spectra are added together as a first step.
#' In this context, a peak means any peak, ignoring any larger structure.
#' For instance, all three
#' peaks in a triplet will be identified as separate peaks.  This
#' facilitates database searching in a subsequent step (and
#' in a crowded spectrum there is no knowing what goes with what
#' anyway).
#'
#' @param spectra An object of class \code{Spectra}.
#' 
#' @param thres Numeric.  Any peaks with a maximum below this value
#' will be discarded.  Determine by inspection (see examples).
#'
#' @param plot Logical.  Option for plot showing bucketed peaks
#' on the original spectrum.
#'
#' @param method One of \code{c("theoretical", "empirical")}, giving the method for
#' computing the second derivative.
#'
#' @param sigma Only applies to \code{method = "theoretical"}.  Numeric.
#' Peak width at half maximum of the curve (FWHM).  Equivalent to
#' 2 * gamma, the half-width at half-maximum (HWHM).
#' In simpler terms, \code{sigma} should correspond to typical peak widths in
#' the spectrum. The units of \code{sigma}
#' are the units of the frequency values in the \code{Spectra Object}.
#' The default value is 0.0005 ppm, as recommended in the reference.
#' You probably need to adjust this parameter.  See Details and the examples.
#' 
#' @param W  Only applies to \code{method = "theoretical"}.  Integer.
#' The number of data points in the window used to
#' compute the second derivative.  Due to the shape of a Lorentzian
#' curve, the second derivative declines dramatically beyond a
#' certain point and time can be saved by not computing values
#' outside a certain range.  The reference recommends a value of 2000
#' for a typical NMR data set.
#'
#' @param ... Additional arguments to be passed to the plotting routines.
#' Most likely you will want to adjust \code{xlim} so you can see 
#' details in the plots.
#'
#' @return A list consisting of a modified \code{Spectra} object, and
#' the indices of the peak boundaries.  Indices correspond to the
#' original \code{Spectra} object, not the one returned.
#'
#' @section Details:
#' The theoretical method is based upon
#' equation 2 in the reference, and requires an assumption about the value
#' of sigma.  However, this function does not extend the zero crossings
#' as described in the reference because the same effect can be accomplished
#' by adjusting the value of \code{sigma}.  The results of the second
#' derivative calculation are very sensitive to the value of \code{sigma}.
#' It is recommended that you experiment to find a suitable value.
#' You can do so using this function and examining a few regions
#' of interest to get \code{thres} and \code{sigma} the way you want them.
#' The empirical method uses the
#' Savitzky-Golay second derivative and makes no assumptions, so 
#' it is more robust.  This method also typically identifies many
#' more peaks than the theoretical method.
#' 
#' 
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @references Daniel Jacob, Catherine Deborde, and Annick Moing.
#' "An Efficient Spectra Processing Method for Metabolite Identification
#' from 1H-NMR Metabolomics Data."
#' Analytical and Bioanalytical Chemistry vol. 405 (2013) pgs. 5049-5061.
#' DOI: 10.1007/s00216-013-6852-y
#'
#' @keywords utilities
#' @export
#' @importFrom EMD extrema
#' @importFrom ChemoSpec chkSpectra
#' @importFrom graphics abline lines plot points rect segments
#' @importFrom stats na.omit
#'
#' @examples
#' 
#' require("ChemoSpec")
#' data(metMUD2)
#'
#' # Inspect to determine suitable threshold
#' s <- colSums(metMUD2$data)
#' x <- metMUD2$freq
#' plot(x, s, type = "l", ylim = c(4, 8),
#' 	xlab = "ppm", ylab = "colSums", main = "Checking the Threshold")
#' th <- 4.5
#' abline(h = th, col = "green")
#'
#' # Shrinkwrap using theoretical method
#' # Note the effect of sigma on the results
#'
#' sw <- shrinkwrapSpectra(spectra = metMUD2, sigma = 0.01,
#' thres = th, xlim = c(1.25, 1.74), xlab = "ppm",
#' main = expression(Shrinkwrap~of~metMUD2~sigma == 0.01))
#'
#' sw <- shrinkwrapSpectra(spectra = metMUD2, sigma = 0.02,
#' thres = th, xlim = c(1.25, 1.74), xlab = "ppm",
#' main = expression(Shrinkwrap~of~metMUD2~sigma == 0.02))
#'
#' sw <- shrinkwrapSpectra(spectra = metMUD2, sigma = 0.01,
#' thres = th, xlim = c(0.85, 1.1), xlab = "ppm",
#' main = expression(Shrinkwrap~of~metMUD2~sigma == 0.01))
#'
#' sw <- shrinkwrapSpectra(spectra = metMUD2, sigma = 0.02,
#' thres = th, xlim = c(0.85, 1.1), xlab = "ppm",
#' main = expression(Shrinkwrap~of~metMUD2~sigma == 0.02))
#'
#' # Compare to empirical method (finds more peaks at a given thres)
#'
#' sw <- shrinkwrapSpectra(spectra = metMUD2, thres = 5,
#' xlim = c(0.85, 1.1), method = "empirical", xlab = "ppm",
#' main = "Shrinkwrap of metMUD2, method = empirical")
#'
#' \dontrun{
#'
#' # A much larger data set
#' data(SrE.NMR)
#'
#' # Inspect and set the threshold
#' s <- colSums(SrE.NMR$data)
#' x <- SrE.NMR$freq
#' plot(x, s, type = "l", ylim = c(-1e6, 1e7),
#' 	xlab = "ppm", ylab = "colSums", main = "Checking the Threshold")
#' th <- 5e5
#' abline(h = th, col = "green")
#'
#' sw <- shrinkwrapSpectra(spectra = SrE.NMR, sigma = 0.01,
#' thres = th, xlim = c(1.2, 1.4), xlab = "ppm",
#' main = expression(Shrinkwrap~of~SrE.NMR~sigma == 0.01))
#'
#' sw <- shrinkwrapSpectra(spectra = SrE.NMR, sigma = 0.02,
#' thres = th, xlim = c(1.2, 1.4), xlab = "ppm",
#' main = expression(Shrinkwrap~of~SrE.NMR~sigma == 0.02))
#'
#' sw <- shrinkwrapSpectra(spectra = SrE.NMR,
#' thres = th, xlim = c(1.2, 1.4), xlab = "ppm", method = "empirical",
#' main = "Shrinkwrap of SrE.NMR, method = empirical")
#' }

shrinkwrapSpectra <- function(spectra, thres = 2.0, plot = TRUE,
	method = "theoretical", sigma = 0.0005, W = 2000, ...) {

	if (!requireNamespace("EMD", quietly = TRUE)) {
		stop("You need to install package EMD to use this function")
		}	

	# Helper function from github.com/bryanhanson/HandyStuff
	vectorizeByRow <- function(IN) {
		OUT <- rep(NA_real_, length(IN))
		nc <- ncol(IN)
		nr <- nrow(IN)
		a <- seq(1, length(IN), nc)
		b <- a + nc - 1
		for (n in 1:length(a)) {
			OUT[a[n]:b[n]] <- IN[n,]
			}
		OUT
	}
	
	if (missing(spectra)) stop("No spectral data provided")
	chkSpectra(spectra)
	mydots <- list(...)

	x <- spectra$freq
	s <- colSums(spectra$data)

	# Get zero crossings which define the edges of the peaks	
	cp <- CP(S = s, X = x, sigma = sigma, W = W, method = method)
	zeros <- extrema(cp)
	zeros <- zeros$cross # a matrix with 2 columns
	
	# Clean up the zero crossings
	# The zero crossings in noisy regions are challenging.
	# From the front of the data, remove crossings up to the point
	# where we see the first downward crossing in the cp data, this
	# would be the first peak, even if it is just noise.
	# Then from the end of the data, remove back to the last 
	# upward crossing, this would be the last peak, even if just noise
	
	zeros <- vectorizeByRow(zeros)
	neg <- cp[zeros] < 0.0
	pos <- cp[zeros] > 0.0
	dropStart <- dropEnd <- NA_integer_
	for (i in 2:length(pos)) { # check the start of the sequence
		if ((neg[i]) & (pos[i-1])) {
			dropStart <- i - 2
			}
		if (!is.na(dropStart)) break # found one, we're done
		}

	for (i in length(pos):2) { # check the end of the sequence
		if ((neg[i-1]) & (pos[i])) {
			dropEnd <- i + 1
			}
		if (!is.na(dropEnd)) break # found one, we're done
		}
	
	if (dropStart > 1) zeros <- zeros[-c(1:dropStart)]
	if (dropEnd < length(pos)) zeros <- zeros[-c(dropEnd:length(pos))]
		
	# # An NMR peak has 2 zero crossings, and 4 associated indices
	# # Keep the outer most set of indices
	
	zeros <- matrix(zeros, ncol = 2, byrow = TRUE)
	keep <- as.logical(1:nrow(zeros) %% 2)
	zeros[keep, 2] <- NA
	zeros[!keep, 1] <- NA
		
	zeros <- na.omit(vectorizeByRow(zeros))
	
	# This next line seems to be necesary; the code above
	# sometimes leaves half a 'peak' in the results otherwise
	if (!(length(zeros) %% 2 == 0)) zeros <- zeros[-length(zeros)]
	
	zeros <- matrix(zeros, ncol = 2, byrow = TRUE)
	
	# Discard bins that are below the threshold
	keep <- NA_integer_ # a vector of zero crossing indices
	dropRow <- NA_integer_ # a vector or rows of zeros that will need to be dropped
	for (i in 1:nrow(zeros)) {
		z <- zeros[i,1]:zeros[i,2]
		if (max(s[z]) < thres) {
			dropRow <- c(dropRow, i)
			next
			}
		keep <- c(keep, z)
		}
	keep <- keep[-1]
	dropRow <- dropRow[-1]
	message("You discarded ", length(dropRow), " peaks which were below the threshold.")
	zeros <- zeros[-dropRow,]

	if (nrow(zeros) == 0L) stop("You discarded all your peaks. Set the threshold lower.")
	message("A total of ", nrow(zeros), " peaks were found.")

	pb <- vectorizeByRow(zeros)
	
	if (plot) {
		x <- spectra$freq
		s <- colSums(spectra$data)
		
		if (!"ylim" %in% names(mydots)) yl <- range(s)
		if ("ylim" %in% names(mydots)) yl <- mydots$ylim
		
 		plot(x, s, type = "n", ylab = "colSums", ...)
    
		yb <- yl[1] # compute rectangles to shade the peaks
		yt <- yl[2]
		re <- as.logical((1:length(pb)) %% 2) # right edge logical
		xl <- x[pb[!re]]
		xr <- x[pb[re]]
		rect(xl, yb, xr, yt, col = "azure")

		lines(x, s)
		points(x, s, cex = 0.5)
		abline(h = thres, col = "green")	
		}

	# Subset the spectra object to include only the peaks found
	spectra$freq <- spectra$freq[keep]
	spectra$data <- spectra$data[,keep]
	
	chkSpectra(spectra)
		
	# Note: PeakBounds are indices of zero crossings in the ORIGINAL spectra object
	return(list(Spectra = spectra, PeakBounds = pb))
	}

