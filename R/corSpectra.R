#' @title
#' Tools for STOCSY Analysis of a Spectra Object
#'
#' @description
#' These functions provide tools for STOSCY analysis of a \code{\link{Spectra}}
#' object, using the method developed by Nicholson.  STOCSY is Statistical
#' Total Correlation Spectroscopy.  Briefly, the correlation matrix of an NMR
#' data set of \code{n} samples and \code{p} frequencies is computed (the
#' matrix dimensions are \code{p x p}).  Peaks arising from the same compound
#' are intrinsically positively correlated.  This is much like a 1D or 2D TOCSY
#' NMR spectrum. However, peaks that are correlated, positively or negatively,
#' due to metabolic processes, will also appear in the STOCSY plot.
#' \code{corSpectra} computes the correlation and covariance matrices, and can
#' display the correlation matrix in several formats.  Detailed
#' inspection/interpretation of this plot is tedious, and producing it can be
#' slow for large data sets, so it's most useful as an overview.
#' \code{covSpectra} will display a single frequency (i.e. chemical shift) from
#' the covariance matrix, but color it according to the correlation matrix
#' values.  This is point where detailed interpretation is done.  See the
#' example.
#'
#' @section Details:
#' If \code{spectra$freq} is in decreasing order, it and \code{spectra$data}
#' are silently re-ordered to be increasing before plotting.  For this reason,
#' if you are providing pre-computed covariance and correlation matrices, be
#' sure that the frequency axis is in ascending order.  See the example for
#' sample code.
#' 
#' The calculation of the correlation and covariance matrices may take quite
#' some time for large data sets.  It is possible to pre-compute these and pass
#' them into the functions to save time and avoid repetition.
#' 
#' Plotting in \code{corSpectra} can be extremely slow for large data sets.
#' The base graphics options (\code{pmode = "contour"} or \code{"image"}) are
#' much faster than the \code{lattice} options.  These plots are probably best
#' for an overall sense of the data and for publication rather than detailed
#' interpretation. If using \code{pmode = "contour"} drawing fewer contours is
#' of course faster for both drawing and computation of the contours.  Note too
#' that contour style plots have \code{n} colors for \code{n} contour levels
#' but image style plots have \code{n-1} colors for \code{n} levels.
#' 
#' The color scale for the plots is blue/low correlation to red/high
#' correlation, anchored at a shade of green for zero correlation.  The example
#' shows how to see the color scale.
#' 
#' For \code{covSpectra}, the x and y limits can be set simply by passing
#' \code{xlim} and \code{ylim} via the \ldots{}.
#' 
#' @aliases corSpectra covSpectra covSpectraJS
#'
#' @param spectra An object of S3 class \code{\link{Spectra}}.
#'
#' @param plot Logical.  Should a plot be made?  Applies to \code{corSpectra}
#' only.
#'
#' @param limX Numeric vector of length 2.  The x limits.  Applies to
#' \code{corSpectra} only.
#'
#' @param limY Numeric vector of length 2.  The y limits.  Applies to
#' \code{corSpectra} only.
#'
#' @param nticks Integer.  The number of ticks to be drawn.  Applies to
#' \code{corSpectra} only.
#'
#' @param levels Numeric.  A vector of values at which to draw the contours or
#' levels.  Applies to \code{corSpectra} only.  The default is to use
#' \code{\link{chooseLvls}} to compute 5 evenly spaced levels.  For most data
#' sets this should only be considered a starting point. A histogram of the
#' correlation matrix can be very helpful in choosing levels.
#'
#' @param pmode Character.  The plot mode.  Applies to \code{corSpectra} only.
#' One of \code{c("contour", "image", "contourplot", "levelplot", "rgl",
#' "exCon")}.  The last two are interactive.  See Details.
#'
#' @param drawGrid Logical.  Shall a faint gray grid be added to the plot?
#' Applies to \code{corSpectra} only.
#'
#' @param R Matrix.  Optional. A precomputed correlation matrix.
#'
#' @param V Matrix.  Optional. A precomputed covariance matrix.
#'
#' @param freq Numeric.  The desired frequency to be plotted.  Applies to
#' \code{covSpectra} only.  In STOCSY terminology, this is the driver peak and
#' it will be marked by a gray dotted line.
#'
#' @param yFree Logical. Applies to \code{covSpectra} only.  If \code{FALSE},
#' the y axis is scaled to the range of the covariance matrix, so that the
#' particular data row can be appreciated in terms of the overall data set.  If
#' \code{TRUE}, the y axis is scaled to the data row specified by \code{freq},
#' so it fills the vertical dimension.
#'
#' @param browser Character.  Applies to \code{covSpectra} only. See
#' \code{\link{plotSpectraJS}} for details.
#'
#' @param minify Logical.  Applies to \code{covSpectraJS} only.  Shall the
#' JavaScript be minified?  This improves performance.  However, it requires
#' package \code{js} which in turn requires package \code{V8}.  The latter is
#' not available on all platforms.  Details may be available at
#' \url{https://github.com/jeroenooms/v8}
#'
#' @param \dots Other parameters to be passed to the plotting functions.
#'
#' @return A list giving the covariance and correlation matrices.  A plot will
#' also be made if \code{plot = TRUE} for \code{corSpectra}.
#' \code{covSpectraJS} gives an interactive web page.
#'
#' @section Warning: For \code{corSpectra}, the labeling of the axis will be
#' wrong if there is a gap in the data.
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#' @seealso \code{\link{chooseLvls}} to select levels automatically.
#' 
#' \code{\link{sortCrossPeaks}} is a utility to sort the covariance matrix by
#' the maximum covariance at a given frequency.  This can be used to identify
#' the strongest interactions.
#'
#' @references O. Cloarec et. al. Analytical Chemistry vol. 77 pgs. 1282-1289
#' (2005).
#' 
#' @importFrom stats cov cov
#' @importFrom utils browseURL
#' @importFrom grDevices rainbow
#' @importFrom graphics image axis contour
#' @importFrom exCon exCon2
#' @importFrom rgl open3d surface3d
#' @importFrom jsonlite toJSON
#' @importFrom js uglify_optimize
#'
#' @export corSpectra covSpectra covSpectraJS
#'
#' @keywords hplot
#'
#' @examples
#' 
#' # The color scale used:
#' # cscale <- c(rev(rainbow(4, start = 0.45, end = 0.66)), rev(rainbow(5, start = 0.0, end = 0.25)))
#' # pie(rep(1, 9), col = cscale)
#' #
#' # This data set is a mixture of compounds, some of
#' # which have correlated concentrations, both
#' # positively and negatively.  The contour plot shows
#' # correlations w/i a spin system and among compounds.
#' # The 2nd plot shows one particular frequency.
#' # For more info about the data set see ?MUD
#' #
#' require("ChemoSpec")
#' data(metMUD2)
#' lvs <- c(-0.99, -0.95, -0.9, 0.9, 0.95, 0.99)
#' lim <- c(0.6, 4.3)
#' res <- corSpectra(metMUD2, levels = lvs, limX = lim, limY = lim, main = "metMUD2 STOCSY Analysis")
#' jnk <- covSpectra(metMUD2, freq = 1.030, C = res[[2]], V = res[[1]])
#' if (interactive()) {
#' require("jsonlite")
#' jnk <- covSpectraJS(metMUD2, freq = 1.030, C = res[[2]], V = res[[1]], minify = TRUE)
#' }
#' # Code to sort the frequency axis
#' # Run prior to pre-computing the cor/cov matrices
#' \dontrun{
#' if (is.unsorted(Spectra_object$freq)) {
#'   Spectra_object$freq <- rev(Spectra_object$freq)
#'   Spectra_object <- Spectra_object[,ncol(Spectra_object$data):1]
#'   }
#' 
#' }
#'

########## This file contains corSpectra, covSpectra & covSpectraJS

########## corSpectra

corSpectra <- function(spectra, plot = TRUE,
	limX = NULL, limY = NULL,
	nticks = 10, levels = NULL,
	pmode = "contour", drawGrid = TRUE,
	R = NULL, V = NULL, ...) {

# For large data sets, there are extreme challenges with cor()
# and in the graphical display.

# NOTE: Cannot subset before computing cor() as this gives the wrong numerical answer

	if ((pmode == "contourplot") | (pmode == "levelplot")) {
		if (!requireNamespace("lattice", quietly = TRUE)) {
			stop("You need to install package lattice to use this option")
			}
		}
	
	if (pmode == "exCon") {
		if (!requireNamespace("exCon", quietly = TRUE)) {
			stop("You need to install package exCon to use this option")
			}
		}

	if (pmode == "rgl") {
		if (!requireNamespace("rgl", quietly = TRUE)) {
			stop("You need to install package rgl to use this option")
			}
		}

	
	if (missing(spectra)) stop("No spectral data provided")
	chkSpectra(spectra)

	# Check to see if spectra$freq is increasing - if not, findInterval will fail
	# Silently reverse things
	if (is.unsorted(spectra$freq)) {
		spectra$freq <- rev(spectra$freq)
		spectra$data <- spectra$data[,ncol(spectra$data):1]
		}
	
	if (is.null(R)) { # user did not provide pre-computed correlation matrix
		X <- spectra$data
		if (ncol(X) > 10000) message("Calculating cor() may take a moment or longer")
		R <- cor(X) # same as (t(X) %*% X)/(nrow(spectra$data) - 1)
		}
	
	if (is.null(V)) { # user did not provide pre-computed covariance matrix
		X <- spectra$data
		if (ncol(X) > 10000) message("Calculating cov() may take a moment or longer")
		V <- cov(X) # same as (t(X) %*% X)/(nrow(spectra$data) - 1)
		}

	if (ncol(R) > 8000) message("Graphical output will take some time with this many data points")
	
	# Helper function to compute ticks, labels & colors
	
	decor <- function(spectra, V, R, limX, limY, levels, pmode) {
				
		# For base functions, the default levels are different than for lattice
		# For base plots, the axes range from [0...1] (rgl too, at least for surface3d)
		# whereas for lattice, range is [1...ncol/nrow]
		# Hence limits must be expressed in different units in each case
		
		# 1.  Fix axis labels
		# If user gives limits in native units, must
		# translate into appropriate units to create labels
		
		# WARNING: labeling won't work when there is a gap
		
		if ((pmode == "contour") | (pmode == "image")) { # base functions
			LX <- c(0, 1) # default values
			LY <- c(0, 1)
			tickposX <- pretty(LX, n = nticks)
			ticklabX <- (diff(range(spectra$freq)) * tickposX) + min(spectra$freq)
			tickposY <- pretty(LY, n = nticks) 
			ticklabY <- (diff(range(spectra$freq)) * tickposY) + min(spectra$freq)

			if (!is.null(limX)) { # override when limX is given
				l <- findInterval(limX[1], spectra$freq)
				r <- findInterval(limX[2], spectra$freq)
				limX<- c(l, r)
				limX <- limX/ncol(R)
				tickposX <- pretty(limX, n = nticks)
				ticklabX <- (diff(range(spectra$freq)) * tickposX) + min(spectra$freq)
				LX <- limX
				}

			if (!is.null(limY)) { # override when limY is given
				l <- findInterval(limY[1], spectra$freq)
				r <- findInterval(limY[2], spectra$freq)
				limY<- c(l, r)
				limY <- limY/ncol(R)
				tickposY <- pretty(limY, n = nticks)
				ticklabY <- (diff(range(spectra$freq)) * tickposY) + min(spectra$freq)
				LY <- limY
				}
			}


		if ((pmode == "contourplot") | (pmode == "levelplot") | (pmode == "rgl")) { # lattice functions + rgl
			LX <- c(1, ncol(R)) # See notes above
			LY <- c(1, ncol(R))
			tickposX <- seq(LX[1], LX[2], length.out = nticks) 
			tickposX <- round(tickposX)
			ticklabX <- spectra$freq[tickposX]		
			tickposY <- seq(LY[1], LY[2], length.out = nticks) 
			tickposY <- round(tickposY)
			ticklabY <- spectra$freq[tickposY]
					
			if (!is.null(limX)) {
				l <- findInterval(limX[1], spectra$freq)
				r <- findInterval(limX[2], spectra$freq)
				limX <- c(l, r)
				tickposX <- seq(limX[1], limX[2], length.out = nticks) 
				tickposX <- round(tickposX)
				ticklabX <- spectra$freq[tickposX]		
				LX <- limX
				}

			if (!is.null(limY)) {
				l <- findInterval(limY[1], spectra$freq)
				r <- findInterval(limY[2], spectra$freq)
				limY <- c(l, r)
				tickposY <- seq(limY[1], limY[2], length.out = nticks) 
				tickposY <- round(tickposY)
				ticklabY <- spectra$freq[tickposY]		
				LY <- limY
				}
			}

		if (pmode == "exCon") {
			LX <- range(spectra$freq)
			LY <- LX
			ticklabX <- NA
			ticklabY <- NA
			tickposX <- NA
			tickposY <- NA
					
			if (!is.null(limX)) {
				LX <- limX
				}

			if (!is.null(limY)) {
				LY <- limY
				}			
			
			}
			
		# 2.  Set levels (contours need a different default than images)
		#     Contours have one color for each level/break/cut
		#	  Image plots have n breaks and n-1 colors
		
		# Color scale for each level
		# blue/low -> red/high, anchored at zero (index 5, a shade of green)
		# max and min will come from the data (i.e., red will be at max of V)
		cscale <- c(rev(rainbow(4, start = 0.45, end = 0.66)), rev(rainbow(5, start = 0.0, end = 0.25)))
		# view with:
		# pie(rep(1, 9), col = cscale)
		
		refscale <- seq(-1, 1, length.out = 9)
				
		if ((pmode == "contour") | (pmode == "contourplot") | (pmode == "exCon")) {
			if (!requireNamespace("exCon", quietly = TRUE)) {
				stop("You need to install package exCon to use this plotting option")
				}
			if (is.null(levels)) {
				levels <- chooseLvls(M = R, n = 5L, mode = "even")
				msg <- paste("The levels chosen are:\n", paste(round(levels, 5), collapse = " "), sep = " ")
				message(msg)
				}	
			myc <- cscale[findInterval(levels, refscale)]
			}

		if ((pmode == "image") | (pmode == "levelplot")) {

			if (!is.null(levels)) {
				myc <- cscale[findInterval(levels, refscale)]
				# need to remove one color:
				nc <- length(myc)
				if ((nc %% 2) == 1) myc <- myc[-ceiling(nc/2)]
				if ((nc %% 2) == 0) myc <- myc[-floor(nc/2)]
				}

			if (is.null(levels)) { # must have one less color than breaks
				levels <- chooseLvls(M = R, n = 5L, mode = "even") # Gives 5 levels
				msg <- paste("The levels chosen are:\n", paste(round(levels, 5), collapse = " "), sep = " ")
				message(msg)
				myc <- cscale[c(1, 3, 7, 9)] # remove colors near zero
				}
			}
				
		if (pmode == "rgl") {
			# Levels don't apply here, simply assign color based upon value
			myc <- cscale[findInterval(R, refscale)] # the colors to be used/color selection			
			}
		
		# 3.  Labeling
		
		lab <- spectra$unit[1]
		if (lab == "ppm") {
			ticklabX <- as.character(round(ticklabX, 2))
			ticklabY <- as.character(round(ticklabY, 2))
			}
		if (lab == "wavenumber") {
			ticklabX <- as.character(round(ticklabX, 0))
			ticklabY <- as.character(round(ticklabY, 0))
			}
					
		L <- list(myc = myc, lab = lab, limX = LX, limY = LY,
			tickposX = tickposX, tickposY = tickposY,
			ticklabX = ticklabX, ticklabY = ticklabY,
			levels = levels,
			refscale = refscale, cscale = cscale)
			
		return(L)
		} # end of decor
	
	# Ready to plot
	
	if (plot) {
		d <- decor(spectra, V, R, limX, limY, levels, pmode)
		# Elements returned by decor:
		# 1. myc
		# 2. lab
		# 3. limX
		# 4. limY
		# 5. tickposX
		# 6. tickposY,
		# 7. ticklabX
		# 8. ticklabY,
		# 9. levels
		#10. reference scale
		#11. color scale
		
		# First two are lattice functions
		
		if (pmode == "levelplot") {
			p <- lattice::levelplot(R, xlab = d[[2]], ylab = d[[2]],
				col.regions = d[[1]],
				scales = list(
					x = list(at = d[[5]], labels = d[[7]]),
					y = list(at = d[[6]], labels = d[[8]])),
				xlim = d[[3]], ylim = d[[4]],
				at = d[[9]],
				colorkey = list( # fixed key, regardless of levels actually used
					at = d[[10]],
					col = d[[11]],
					labels = list(
					at = seq(-1.0, 1.0, by = 0.2), 
                         labels = as.character(seq(-1.0, 1.0, by = 0.2)))),
				...)
			print(p)
			}
		
		if (pmode == "contourplot") {
			p <- lattice::contourplot(R, xlab = d[[2]], ylab = d[[2]],
				col.regions = d[[1]],
				region = TRUE,
				scales = list(
					x = list(at = d[[5]], labels = d[[7]]),
					y = list(at = d[[6]], labels = d[[8]])),
				xlim = d[[3]], ylim = d[[4]],
				labels = FALSE,
				at = d[[9]], # passes through to panel function
				colorkey = list( # fixed key, regardless of levels actually used
					at = d[[10]],
					col = d[[11]],
					labels = list(
						at = seq(-1.0, 1.0, by = 0.2), 
                        labels = as.character(seq(-1.0, 1.0, by = 0.2)))),
				...)
			print(p)
			}

		# Next two are base functions (these are much faster)
		
		if (pmode == "image") {
			image(R, xlab = d[[2]], ylab = d[[2]],
				col = d[[1]], xlim = d[[3]], ylim = d[[4]],
				breaks = d[[9]],
				xaxt = "n", yaxt = "n", useRaster = TRUE, ...)			
			axis(1, at = d[[5]], labels = d[[7]])			
			axis(2, at = d[[6]], labels = d[[8]])			
			}

		if (pmode == "contour") {
			# Next 2 lines needed to get grid under contours
			plot(x = 0, y = 0, xlim = d[[3]], ylim = d[[4]], type = "n",
				xlab = "", ylab = "", xaxt = "n", yaxt = "n")
			if (drawGrid) abline(v = d[[6]], h = d[[5]], col = "gray95")		
			contour(R, xlab = d[[2]], ylab = d[[2]],
				xlim = d[[3]], ylim = d[[4]],
				col = d[[1]],
				levels = d[[9]],				
				drawlabels = FALSE,
				axes = FALSE, frame.plot = TRUE,
				xaxs = "i", yaxs = "i", add = TRUE, ...)
			axis(1, at = d[[5]], labels = d[[7]])			
			axis(2, at = d[[6]], labels = d[[8]])
			}
		
		# Interactive versions
		
		if (pmode == "exCon") {

			x1 <- findInterval(d[[3]][1], spectra$freq)
			x2 <- findInterval(d[[3]][2], spectra$freq)
			y1 <- findInterval(d[[4]][1], spectra$freq)
			y2 <- findInterval(d[[4]][2], spectra$freq)
			Z <- R[x1:x2, y1:y2]
					
			exCon2(M = Z, levels = d[[9]],
				x = seq(d[[3]][1], d[[3]][2], length.out = ncol(Z)),
				y = seq(d[[4]][1], d[[4]][2], length.out = nrow(Z)),
				, ...)		
			}

		if (pmode == "rgl") {
			x1 <- d[[3]][1]
			x2 <- d[[3]][2]
			y1 <- d[[4]][1]
			y2 <- d[[4]][2]
			Z <- R[x1:x2, y1:y2]

			# Crude attempt to adjust aspect ratio (max(Z) always 1 or -1)
			np <- nrow(Z)
			while(abs(max(Z)/np) < 0.40) {
				Z <- Z * 1.5	
				}

			open3d(windowRect = c(0, 0, 800, 800))
			surface3d(x = 1:np, y = 1:np, z = Z, color = d[[1]])
			}

		}
		
	L <- list(cov = V, cor = R)
	invisible(L)
	
	}
		
########## covSpectra

#' @rdname corSpectra

covSpectra <- function(spectra, freq = spectra$freq[1],
	R = NULL, V = NULL, yFree = TRUE, ...) {

# For large data sets, there are extreme challenges with cor()

# NOTE: Cannot subset before computing cor() as this gives the wrong numerical answer

	if (missing(spectra)) stop("No spectral data provided")
	chkSpectra(spectra)

	# Check to see if spectra$freq is increasing - if not, findInterval will fail
	# Silently reverse things
	if (is.unsorted(spectra$freq)) {
		spectra$freq <- rev(spectra$freq)
		spectra$data <- spectra$data[,ncol(spectra$data):1]
		}
		
	row <- findInterval(freq, spectra$freq)
	
	if (is.null(R)) { # user did not provide pre-computed correlation matrix
		X <- spectra$data
		if (ncol(X) > 10000) message("Calculating cor() may take a few minutes")
		R <- cor(X) # same as (t(X) %*% X)/(nrow(spectra$data) - 1)
		}
	
	if (is.null(V)) { # user did not provide pre-computed covariance matrix
		X <- spectra$data
		if (ncol(X) > 10000) message("Calculating cov() may take a few minutes")
		V <- cov(X) # same as (t(X) %*% X)/(nrow(spectra$data) - 1)
		}

	# Color scale for each level
	# blue/low -> red/high, anchored at zero (index 5, a shade of green)
	# max and min will come from the data (i.e., red will be at max of V)
	cscale <- c(rev(rainbow(4, start = 0.45, end = 0.66)), rev(rainbow(5, start = 0.0, end = 0.25)))
	# view with:
	# pie(rep(1, 9), col = cscale)

	refscale <- seq(-1, 1, length.out = 9)
	
	# Now average every contiguous pair of values in R[row,] so that there is one
	# less value, and use the mean value of the pair to assign colors
	# e.g. the mean of points n & n+1 determines the color used to plot that segment

	cd <- diff(R[row,])
	cr <- 0.5 * cd + R[row,][-length(R[row,])] # this will have one value less than the data
	myc <- cscale[findInterval(cr, refscale)] # color based upon cor, not cov

	# Ready to plot
	
	np <- length(spectra$freq)
	ind1 <- 1:(np-1)
	ind2 <- 2:np
	
	if (yFree) {
		plot(spectra$freq, V[row,], type = "n",
			xlab = spectra$unit[1], ylab = "covariance",
			main = paste("Frequency = ", sprintf("%5.5f", spectra$freq[row]), sep = ""),
			...)
		segments(spectra$freq[ind1], V[row, ind1], spectra$freq[ind2], V[row, ind2],  col = myc, ...)
		abline(v = spectra$freq[row], lty = 2, col = "gray")
		}
		
	if (!yFree) {
		plot(spectra$freq, V[row,], type = "n", ylim = range(V),
			xlab = spectra$unit[1], ylab = "covariance",
			main = paste("Frequency = ", sprintf("%5.5f", spectra$freq[row]), sep = ""),
			...)
		segments(spectra$freq[ind1], V[row, ind1], spectra$freq[ind2], V[row, ind2],  col = myc, ...)
		abline(v = spectra$freq[row], lty = 2, col = "gray")
		}
		
	L <- list(cov = V, cor = R)
	invisible(L)
	
	}
		
########## covSpectraJS

#' @rdname corSpectra

covSpectraJS <- function(spectra, freq = spectra$freq[1],
	R = NULL, V = NULL, browser = NULL, minify = TRUE, ...) {

# Function to carry out Nicholson's STOCSY analysis
# Part of the ChemoSpec package
# Bryan Hanson, DePauw University

# This is the interactive JS version May 2015

	if (missing(spectra)) stop("No spectral data provided")
	chkSpectra(spectra)

	# Check to see if spectra$freq is increasing - if not, findInterval will fail
	# Silently reverse things
	if (is.unsorted(spectra$freq)) {
		spectra$freq <- rev(spectra$freq)
		spectra$data <- spectra$data[,ncol(spectra$data):1]
		}

	row <- findInterval(freq, spectra$freq)

	if (is.null(R)) { # user did not provide pre-computed correlation matrix
		X <- spectra$data
		if (ncol(X) > 10000) message("Calculating cor() may take a few minutes")
		R <- cor(X)
		}

	if (is.null(V)) { # user did not provide pre-computed covariance matrix
		X <- spectra$data
		if (ncol(X) > 10000) message("Calculating cov() may take a few minutes")
		V <- cov(X)
		}

	# Color scale for each level
	# blue/low -> red/high, anchored at zero (index 5, a shade of green)
	# max and min will come from the data (i.e., red will be at max of V)
	cscale <- c(rev(rainbow(4, start = 0.45, end = 0.66)), rev(rainbow(5, start = 0.0, end = 0.25)))
	# view with:
	# pie(rep(1, 9), col = cscale)
	cscale <- substr(cscale, 1, 7) # needed as JS doesn't accept alpha channel

	refscale <- seq(-1, 1, length.out = 9)

	# Now average every contiguous pair of values in R[row,] so that there is one
	# less value, and use the mean value of the pair to assign colors
	# e.g. the mean of points n & n+1 determines the color used to plot that segment

	cd <- diff(R[row,])
	cr <- 0.5 * cd + R[row,][-length(R[row,])] # this will have one value less than the data
	myc <- cscale[findInterval(cr, refscale)] # color based upon cor, not cov

	myc <- substr(myc, 1, 7) # needed as JS doesn't accept alpha channel

	#cat("There are", length(myc), "colors.  The first 10 are:\n")
	#print(myc[1:10])
	#cat("There are", length(spectra$freq), "data points\n")

	ylab <- paste("Covariance at Frequency = ", sprintf("%5.5f", spectra$freq[row]), sep = "")

	# Ready to plot

	if (requireNamespace("jsonlite", quietly = TRUE)) {

		# Break the pieces of the Spectra object out into
		# separate JSON entities
		# These will be global variables in the JavaScript

		Freq <- toJSON(spectra$freq)
		Y <- toJSON(V[row,])
		xUnit <- toJSON(spectra$unit[1])
		yLabel <- toJSON(ylab)
		Desc <- toJSON(spectra$desc)
		Dx <- toJSON(range(spectra$freq))
		Dy <- toJSON(range(V[row,]))
		myc <- toJSON(myc)
		keyScale <- toJSON(cscale)
		driver <- toJSON(freq)

		# Prepare for writing
		# Groups commented out as it is not currently used

		data1 <- paste("var Freq = ", Freq, sep = "")
		data2 <- paste("var Y = ", Y, sep = "")
		data3 <- paste("var xUnit = ", xUnit, sep = "")
		data4 <- paste("var yLabel = ", yLabel, sep = "")
		data5 <- paste("var Desc = ", Desc, sep = "")
		data6 <- paste("var Dx = ", Dx, sep = "")
		data7 <- paste("var Dy = ", Dy, sep = "")
		data8 <- paste("var myc = ", myc, sep = "")
		data9 <- paste("var keyScale = ", keyScale, sep = "")
		data10 <- paste("var driver = ", driver, sep = "")

		# Get the JavaScript modules & related files

#		td <- getwd()
		td <- tempdir()
		fd <- system.file("extdata", package = "ChemoSpec")
		cSfiles <- c("cS.css", "cS_globals.js", "cS_controls.js",
		"cS_brushNguides.js", "cS_main.js", "covSpectraJS.html", "cS_spectra.js")
		chk2 <- file.copy(from=file.path(fd, cSfiles), to=file.path(td, cSfiles),
			 overwrite = TRUE)
		if (!all(chk2)) stop("Copying to temporary directory failed")

		js1 <- readLines(con = file.path(td,"cS_globals.js"))
		js2 <- readLines(con = file.path(td,"cS_brushNguides.js"))
		js3 <- readLines(con = file.path(td,"cS_controls.js"))
		js4 <- readLines(con = file.path(td,"cS_spectra.js"))
		js5 <- readLines(con = file.path(td,"cS_main.js"))

		# The following are used to wrap the entire JS code in a
		# scoping function so that performance is improved,

		scopeFunHeader <- "(function() {"
		scopeFunTail <- "})();"

		# Combine, then optionally minify for faster performance, write

		text = c(scopeFunHeader, data1, data2, data3, data4,
			data5, data6, data7, data8, data9, data10,
			js1, js2, js3, js4, js5, scopeFunTail)

		if (minify) {
			if (requireNamespace("js", quietly = TRUE)) {
				text <- uglify_optimize(text, unused = FALSE)
				}
			if (!requireNamespace("js", quietly = TRUE)) {
				stop("You need install package js to minify the JavaScript code")
				}
			}


		writeLines(text, sep  = "\n", con = file.path(td,"cS.js"))

		# Open the file in a browser

		pg <-  file.path(td,"covSpectraJS.html")
		if (!is.null(browser)) {
		    browseURL(pg, browser = browser)
			} else {
			# open in RStudio if viewer is not null
		    # similar to htmltools::html_print
				viewer <- getOption("viewer")
			  	if (is.null(browser) && !is.null(viewer)) {
		      		viewer(pg)
			  		} else {
			    		browseURL(pg)
			  			}
			}

		message("The covSpectraJS web page is in the following\ntemp directory which is deleted when you quit R: ")
		message(td)
		return(invisible())
	}

	if (!requireNamespace("jsonlite", quietly = TRUE)) {
		stop("You need install package jsonlite to use this function")
		}

	L <- list(cov = V, cor = R)
	invisible(L)

	}
