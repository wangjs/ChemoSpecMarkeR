#' Sort Rows of a Covariance Matrix of a Spectra Object
#' 
#' Given a covariance matrix, the maximum absolute value of each row is
#' computed after excluding values near the diagonal.  By default, the results
#' are ordered by the absolute value of the covariance, or optionally by
#' frequency.  The results may be limited to a selected quantile of top
#' results.
#' 
#' The diagonal of a covariance matrix holds the variances, which are by
#' definition the largest value in the row.  Because of the shape of an NMR
#' peak, values near the diagonal essentially represent the variance as well.
#' We are interested in the values of the cross peaks, hence it is necessary to
#' remove a number of data points near the diagonal from the computations.  The
#' exact number to remove is subject to experimentation.
#' 
#' If the data were sorted by frequency without selecting the top values (i.e.
#' \code{Quan = NULL}), there would be no point in collapsing the values as
#' they are all equally spaced.  But when you have selected the top covariance
#' values, the frequencies are no longer consecutive and equally spaced, so
#' collapsing them replaces trivially close values with their averages.
#' 
#' @param spectra An object of S3 class \code{\link{Spectra}}.
#'
#' @param V A numeric covariance matrix, corresponding to
#' \code{cov(spectra$data)}.
#'
#' @param window Numeric.  A value in the interval (0...1). This is converted
#' to a number of data points near the diagonal which will be removed.  If
#' \code{NULL}, the default, approximately 10 percent of the points are
#' removed.
#'
#' @param Quan Numeric.  A value in the interval (0...1) giving the quantile to
#' be selected.  For instance, \code{Quan = 0.1} selects the top 10 percent of
#' the returned values.
#'
#' @param byFreq Logical.  Shall the results be sorted by frequency?  Only
#' applies if \code{Quan != NULL}.  Frequencies closer than \code{freqThres}
#' will be collapsed and replaced with their average.
#'
#' @param freqThres Numeric.  Only relevant if \code{Quan != NULL} and
#' \code{byFreq = TRUE}. Frequencies closer than \code{freqThres} will be
#' collapsed and replaced with their average.
#'
#' @param ... Arguments to be passed down stream.  In particular, you may
#' wish to set \code{verbose = TRUE} which will pass through to
#' \code{collapseRowsOrCols} and give information about which rows and columns
#' which are collapsed if \code{Quan != NULL} and \code{byFreq = TRUE}.
#'
#' @return A data frame containing the frequencies from the
#' \code{\link{Spectra}} object, the maximum covariance at that frequency, the
#' absolute value of the maximum covariance, and the relative absolute maximum
#' covariance.
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @importFrom stats cov quantile
#' @importFrom plyr arrange desc
#'
#' @export
#'
#' @keywords utilities
#'
#' @examples
#' require("ChemoSpec")
#' data(metMUD1)
#' V <- cov(metMUD1$data)
#' # Look at the top 1%
#' res <- sortCrossPeaks(metMUD1, V, Quan = 0.01)
#' res
#' 
sortCrossPeaks <- function(spectra, V = NULL, window = NULL, Quan = NULL,
	byFreq = FALSE, freqThres = 0.01, ...) {
	
	# Function to select rows with high covariance, ignoring the region near the diagonal
	# Part of the ChemoSpec package
	# Bryan Hanson, DePauw University, January 2015
	
	message("Beta version of sortCrossPeaks, check your results carefully")
	
	abs_max_cov = NULL # needed for checking / non-standard evaluation
	freq = NULL # needed for checking / non-standard evaluation
	
	if (is.null(V)) V <- cov(spectra$data)
	
	# Remove entries in V near the diagonal (diagonal is variance)
	# window as passed is a percentage of the rows. After conversion, it is the number of
	# entries on either side of the diagonal to remove.
	
	if (!is.null(window)) {
		if ((window < 0) | (window > 1)) stop("The value of window must be in (0..1)")
		window <- ceiling(window*nrow(V))
		}
	if (is.null(window)) window <- ceiling(0.1*nrow(V))
		
	if (window <= 1) {
		window <- 0
		message("Only removing the diagonal, you may wish to increase window")
		}
	
	# Next trick inspired by the code for upper.tri
	wh <- col(V) > row(V) + window # assumes square matrix! (which we have)
	V[!wh] <- NA
	
	# And this one is based upon http://stackoverflow.com/a/29857455/633251
	
	V <- pmax(V, t(V), na.rm = TRUE) # Now it is symmetric
	
	# Now locate and organize by the extreme values
	
	cv <- rep(NA_real_, nrow(V))
	for (i in 1:nrow(V)) {
		cvmin <- min(V[i,], na.rm = TRUE)  # min with sign
		cvmax <- max(V[i,], na.rm = TRUE)  # max with sign
		if (abs(cvmax) >= abs(cvmin)) cv[i] <- cvmax # save the most extreme value
		if (abs(cvmin) > abs(cvmax)) cv[i] <- cvmin
		}
	cva <- abs(cv)	
	
	df <- data.frame(freq = spectra$freq, max_cov = cv, abs_max_cov = cva)
	df <- arrange(df, desc(abs_max_cov))
	df$rel_abs_max_cov <- df$abs_max_cov/df$abs_max_cov[1]
		
	# Return only a portion if requested
	
	if (!is.null(Quan)) {
		if (!length(Quan) == 1) stop("Quan should be a single number")
		Quan <- 1 - Quan
		Q <- quantile(cva, Quan)
		keep <- which(df$abs_max_cov >= Q)
		df <- df[keep, ]

	# byFreq = TRUE sort instead by freq, then collapse rows
	
		if (byFreq) {
			df <- arrange(df, freq)
			df <- collapseRowsOrCols(as.matrix(df), 1, thres = freqThres,...)
			} # end of byFreq = TRUE
		
	} # end of !is.null(Quan)
	
	return(as.data.frame(df))
	}
