#' Compute the Second Derivative of a Lorentzian Function
#'
#' Computes the second derivative of a Lorentzian curve using
#' the either the theoretical form, or the Savitzky-Golay
#' empirical method.
#' 
#' @param x A numeric vector of x values from a Lorentzian curve.
#'
#' @param x0 Numeric.  Only applies to \code{method = "theoretical"}.
#' Center frequency of the peak; also called the location parameter.
#' 
#'
#' @param sigma Only applies to \code{method = "theoretical"}.  Numeric.
#' Peak width at half maximum of the curve (FWHM).  Equivalent to
#' 2 * gamma, the half-width at half-maximum (HWHM). The units of sigma
#' are the units of the data passed (\code{x, x0}).  The default value
#' is 0.0005 ppm, as recommended in the reference.
#' 
#'
#' @param y Only applies to \code{method = "empirical"}.
#' A numeric vector of y values from a Lorentzian curve.  
#'
#' @param method One of \code{c("theoretical", "empirical")}, giving the method for
#' computing the second derivative.
#'
#' @return A numeric vector containing the values of the second derivative.
#'
#' @section Details:
#' The theoretical method is based upon
#' equation 2 in the reference, and requires an assumption about the value
#' of sigma.  The results of the second
#' derivative calculation are very sensitive to the value of sigma.
#' It is recommended that you experiment to find a suitable value.
#' The empirical method uses the
#' Savitzky-Golay second derivative and makes no assumptions.
#' 
#' @author Brian K. Saulnier, Bryan A. Hanson,  DePauw University. \email{hanson@@depauw.edu}
#'
#' @references Daniel Jacob, Catherine Deborde, and Annick Moing.
#' "An Efficient Spectra Processing Method for Metabolite Identification
#' from 1H-NMR Metabolomics Data."
#' Analytical and Bioanalytical Chemistry vol. 405 (2013) pgs. 5049-5061.
#' DOI: 10.1007/s00216-013-6852-y
#'
#' @keywords utilities
#' @export
#' @importFrom signal sgolayfilt
#'
#' @examples
#' require("SpecHelpers")
#' require("signal")
#'
#' # Create a basic Lorentzian curve with a single peak
#' loren1 <- data.frame(x0 = 0, area = 1, gamma = 0.5)
#' lorentz1 <- makeSpec(loren1, plot = FALSE, type = "lorentz",
#' dd = 100, x.range = c(-10, 10))
#'
#' # Figure out the value of sigma on the scale of the data.  0.0005 ppm is the default,
#' # but in this case, the data are given in Hertz and gamma = 0.5.
#' # Sigma is by definition half of gamma, 0.25.
#' # NOTE: sigma has a big effect on the shape of the 2nd derivative
#' # when using the theoretical method!
#'
#' # Compute the second derivative using each method
#' x <- lorentz1[1,] # Frequency values
#' y <- lorentz1[2,] # Intensity values
#' d2Theo <- SDL(x = x, x0 = x[1000], sigma = 0.25)
#' d2Emp <- SDL(y = y, method = "empirical")
#'
#' # Now plot
#' ylabel <- "data (black), theoretical 2nd deriv (red), empirical 2nd deriv (blue)"
#' myt <- "Second Derivative of a Lorentzian Signal"
#' plot(x, y, type = "l", ylab = ylabel, main = myt,
#' 	xlim = c(-5, 5), ylim = c(-1, 1))
#' lines(x, d2Theo/1e3, col = "red")
#' lines(x, d2Emp*1e3, col = "blue")
#' abline(h = 0.0, col = "green")
#'
#' # Create a Lorentzian curve with two overlapping peaks
#' loren2 <- data.frame(x0 = c(0, 1.0), area = c(1, 0.75), gamma = c(0.5, 0.5))
#' lorentz2 <- makeSpec(loren2, plot = FALSE, type = "lorentz",
#' dd = 100, x.range = c(-10, 10))
#' x <- lorentz2[1,]
#' y <- lorentz2[2,]
#' d2Theo <- SDL(x = x, x0 = x[1000], sigma = 0.25)
#' d2Emp <- SDL(y = y, method = "empirical")
#' plot(x, y, type = "l", ylab = ylabel, main = myt,
#' 	xlim = c(-3, 3), ylim = c(-1, 1))
#' lines(x, d2Theo/1e3, col = "red")
#' lines(x, d2Emp*1e3, col = "blue")
#' abline(h = 0.0, col = "green")
#'
#' # Create some "real" NMR data using a function from SpecHelpers that works
#' # on a typical NMR scale (plot is labeled with ppm, and input in is ppm,
#' # but internally, the units are Hz relative to TMS at 0.0)
#' # This is ethyl acetate: CH3CO2CH2CH3 @ 300 MHz (default field strength)
#' peaks <- data.frame(
#' 	delta = c(1.26, 2.04, 4.12),
#' 	mult = c(3, 1, 4),
#' 	J = c(14, 0, 14),
#' 	area = c(3, 3, 2),
#' 	pw = c(5, 5, 5))
#' nmr <- plotNMRspec(peaks = peaks, plot = TRUE, x.range = c(5, 0))
#' #
#' # Take 2nd derivative relative to the center of the
#' # triplet at 1.26 (378 Hz).
#' # The units of the data are Hz, at 300 MHz. Sigma should be
#' # 0.0005 ppm according to the reference.  plotNMRspec computes
#' # dd (passed to makeSpec) as MHz * ppHz and we must multiply by
#' # 0.0005 as above (ppHz = 1 by default; fine for 1H NMR):
#' sig <- 300 * 0.0005
#' d2Theo <- SDL(x = nmr[1,], x0 = nmr[1, 379], sigma = sig)
#' d2Emp <- SDL(y = nmr[2,], method = "empirical")
#' x <- nmr[1,]
#' y <- nmr[2,]
#'
#' # Detailed view of triplet
#' plot(x, y, type = "l", xlim = c(379-25, 379+25), ylim = c(-0.5, 0.5),
#' 	 xlab = "Hz", ylab = ylabel, main = myt)
#' lines(x, rev(d2Theo), col = "red")
#' lines(x, d2Emp*10, col = "blue")
#' abline(h = 0.0, col = "green")

SDL <- function(x = NULL, x0 = NULL, y = NULL, sigma = 0.0005, method = "theoretical"){
	
	if (method == "theoretical") {
		
		if (is.null(x)) stop("You must supply a vector of signal x values")
		if (is.null(x0)) stop("You must supply x0")
		if (!((min(x) <= x0) & (x0 <= max(x)))) stop("x0 appears out of range")
		if (!sigma > 0) stop("sigma must be greater than zero.")
		
		num <- 16 * sigma * ((12 * (x-x0)^2) - sigma^2)
		denom <- pi * ((4 * (x - x0)^2) + sigma^2)^3
		sdl <-  num/denom
		}
		
	if (method == "empirical") {
		
		if (!requireNamespace("signal", quietly = TRUE)) {
			stop("You need to install package signal to use this option")
			}	
		if (is.null(y)) stop("You must supply a vector of signal intensities")
		
		sdl <- sgolayfilt(y, m = 2)
		}
		
	return(sdl)
	}

