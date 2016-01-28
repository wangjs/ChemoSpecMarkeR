#'
#' Convolution Product of a Spectrum and its Second Derivative.
#' 
#' Computes the convolution product of a spectrum and its second derivative.
#' There are two options for the second derivative.  It may be calculated
#' from the theoretical form or
#' using the Savitzky-Golay method, which is empirical.
#'
#' @param S Numeric vector of spectral intensities.
#'
#' @param X Numeric vector of frequencies.
#'
#' @param method One of \code{c("theoretical", "empirical")}, giving the method for
#' computing the second derivative.
#'
#' @param sigma Only applies to \code{method = "theoretical"}.  Numeric.
#' Peak width at half height of the curve.  Passed to
#' SDL; See the documentation there for details.
#'
#' @param W Only applies to \code{method = "theoretical"}. Integer.
#' The number of data points in the window used to
#' compute the second derivative.  Due to the shape of a Lorentzian
#' curve, the second derivative declines dramatically beyond a
#' certain point and time can be saved by not computing values
#' outside a certain range.  The reference recommends a value of 2000
#' for a typical NMR data set.
#'
#' @return The convolution product as a numeric vector with \code{length(S)}.
#'
#' @section Note:
#' \code{method = "empirical"} currently returns the Savitzky-Golay
#' second derivative, NOT the convolution product, due to a bug.
#'
#' @author Brian K. Saulnier, Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @references Daniel Jacob, Catherine Deborde, and Annick Moing.
#' "An Efficient Spectra Processing Method for Metabolite Identification
#' from 1H-NMR Metabolomics Data."
#' Analytical and Bioanalytical Chemistry vol. 405 (2013) pgs. 5049-5061.
#' DOI: 10.1007/s00216-013-6852-y
#'
#' @keywords utilities
#' @export
#'
#' @examples
#' require("SpecHelpers")
#' require("signal")
#'
#' # Create a Lorentzian curve with a single peak
#' loren1 <- data.frame(x0 = 0, area = 1, gamma = 0.5)
#' lorentz1 <- makeSpec(loren1, plot = FALSE, type = "lorentz", dd = 100, x.range = c(-10, 10))
#'
#' # Compute convolution using each method
#' x <- lorentz1[1,] # Frequency values
#' y <- lorentz1[2,] # Intensity values
#' sig <- 100 * 0.0005 # See SDL documentation
#' cpTheo <- CP(S = y, X = x, sigma = sig)
#' cpEmp <- CP(S = y, method = "empirical")
#'
#' # Plot the original data, compare to convolution product
#' ylabel <- "data (black), Theo. Conv. Prod. (red), S-G. Conv. Prod. (blue)"
#' myt <- "Convolution Product"
#' plot(x, y, type = "l", ylab = ylabel, main = myt,
#' 	ylim = c(-1, 1), xlim = c(-2, 2))
#' lines(x, cpTheo/1000, col = "red")
#' lines(x, cpEmp*100, col = "blue")
#' abline(h = 0.0, col = "green")
#'
#' # Create a Lorentzian curve with two overlapping peaks
#' loren2 <- data.frame(x0 = c(0, 1.0), area = c(1, 0.75), gamma = c(0.5, 0.5))
#' lorentz2 <- makeSpec(loren2, plot = FALSE, type = "lorentz", dd = 100, x.range = c(-10, 10))
#' x <- lorentz2[1,]
#' y <- lorentz2[2,]
#' sig <- 100 * 0.0005 # See SDL documentation
#' cpTheo <- CP(S = y, X = x, sigma = sig)
#' cpEmp <- CP(S = y, method = "empirical")
#'
#' # Plot the original data, compare to convolution product
#' plot(x, y, type = "l", ylim = c(-0.75, 0.75), ylab = ylabel, main = myt)
#' lines(x, cpTheo/1000, col = "red")
#' lines(x, cpEmp*100, col = "blue")
#'
#' # A tougher test
#' require("ChemoSpec")
#' data(metMUD2)
#' y <- colSums(metMUD2$data)
#' cpTheo <- CP(S = y, X = metMUD2$freq, sigma = 0.01)
#' SG <- sgolayfilt(x = y, m = 2)
#'
#' # Plot a region of interest
#' xl <- c(400, 450)
#' plot(cpTheo, type = "l", xlim = xl, col = "red", ylab = ylabel, main = myt)
#' lines(y*1e6, col = "black")
#' lines(SG*1e6, col = "blue")
#' abline(h = 0, col = "green")

CP <- function(S = NULL, X = NULL, method = "theoretical", sigma = 0.0005, W = 2000) {

    if (method == "theoretical") {
    	
    		if (is.null(X)) stop("You must supply a vector of x values")

        cp <- rep(NA_real_, length(X))
		P <- floor(W/2)

        for(i in 1:length(X)) {
 			if ((i + P) > length(X)) P <- (length(X) - i + 1)
			if (i < P) P <- i
			# cat("P = ", P, "\n\n")
			# Assemble the indices corresponding to the window
			idx <- seq(i - P + 1, i + P - 1, 1)
            # Now compute the sdl
            cp[i] <- sum(SDL(x = X[idx], x0 = X[i], sigma = sigma) * S[idx])
			P <- floor(W/2) # reset at the end of each iteration
            }
        }

    if (method == "empirical") {
    	
    		if (is.null(S)) stop("You must supply a vector of signal intensities")

        sdl <- SDL(y = S, method = "empirical")
        # This next code doesn't give the right peak shape, needs work
		#cp <- convolve(S, rev(sdl), type = "open")
		#offset <- floor(length(S)/2)
		#cp <- cp[offset:(length(cp)-offset)]
		cp <- sdl # temporary
      	}

	return(cp)
	}




