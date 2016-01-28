#' Process a List of Cross Peaks
#' 
#' This function assists the direct analysis of a STOCSY plot.  During manual
#' inspection of the plot, one can record a driver peak and the corresponding
#' cross peaks.  Since a STOCSCY plot is a contour plot, the exact position of
#' a peak is subject to some interpretation.  Further, since a STOCSY plot is
#' typically inspected in sectors, driver peaks may vary slightly across
#' sectors.  This function will take the list made during inspection and clean
#' it up, returning a data frame of the resulting peaks.
#' 
#' The larger the values of the thresholds, the more collapsing will occur.
#' Collapsing means rows or columns separated by less than the threshold will
#' be averaged.  For instance, two driver peaks at 4.05 and 4.07 would be
#' collapsed if \code{rowThres} were 0.02 but not if it were 0.01.  Since the
#' ultimate goal is to get unique peaks corresponding to biomarkers via look up
#' in a database, one would generally want to collapse peaks to a degree. See
#' the warnings.
#' 
#' @param spectra An object of S3 class \code{\link{Spectra}}.
#'
#' @param file The name of a csv file containing the cross peak information.
#' The first column should be the driver peaks.  Other columns should give the
#' frequencies of the cross peaks observed for a given driver (so for a given
#' row, the first entry is the driver peak, and the rest of the entries are the
#' cross peaks).  May have missing values.  Should have one value per cell.
#'
#' @param rowThres Numeric.  A value in the frequency domain.  Frequencies
#' closer than this value will be collapsed.  See details.
#'
#' @param colThres Numeric.  A value in the frequency domain.  Frequencies
#' closer than this value will be collapsed.  See details.
#'
#' @param collapseRows Logical.  Should the rows be collapsed if possible?
#'
#' @param collapseCols Logical.  Should the columns be collapsed if possible?
#'
#' @param verbose Logical.  Shall diagnostic information be printed at the
#' console?
#'
#' @return A data frame containing the driver peaks in the first column, and
#' the aligned cross peaks in the other columns.  The dimensions of this data
#' frame will depend on the options to collapse and the thresholds.
#'
#' @section Warning: Collapsing peaks can be tricky and may lead to very
#' undesirable results.  First, threshold values below the resolution of the
#' data will not collapse anything.  Second, if you have a series of peaks that
#' are closely spaced, for instance a messy multiplet that arises from the
#' overlap of several peaks, collapsing the series will replace the series with
#' one peak when there were several originally.  You should look carefully at
#' your 1D spectra to see if this could be a problem.  I recommend some
#' experimentation to see if collapsing is wise for your data set.
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @importFrom utils read.csv
#'
#' @export
#'
#' @keywords utility
#'
#' @examples
#' 
#' ##---- None at this time !! ----
#' 
cleanSTOCSYpeaks <- function(spectra, file = NULL, rowThres = 0.01, colThres = 0.01,
	collapseRows = TRUE, collapseCols = TRUE, verbose = TRUE) {
		
	# Step 0: Get the data
	
	pks <- read.csv(file, header = FALSE)
	names(pks)[1] <- "driver" # other columns will be V1, V2 ...
	driver <- NULL # NSE/CRAN checks
	pks <- plyr::arrange(pks, driver) # sort it
	
	# Step 1a: Remove any columns that are all NA
	
	throw <- NA_integer_
	
	if (verbose) message("Removing rows that are all NA")
	
	for (i in 1:nrow(pks)) {
		chk <- all(is.na(pks[i,2:ncol(pks)]))
		if (chk) throw <- c(throw, i)
		}
	
	throw <- throw[-1]
	if (verbose) cat("Removed ", length(throw), " rows\n")
	if (length(throw) > 0) pks <- pks[-throw,]
	
	# Step 1: Average & collapse any group of rows that are within the specified threshold

	if (collapseRows) {
		pks <- collapseRowsOrCols(pks, ind = 1, rows = TRUE, thres = rowThres, verbose = verbose)
		} # end of collapseRows = TRUE
		
	# Step 2: Align the peak positions with the known frequencies
	# Set up an empty matrix in which we'll align the peaks
	# But first, check the order of spectra$freq
	
	if (is.unsorted(spectra$freq)) spectra$freq <- rev(spectra$freq)

	if (verbose) message("Mapping frequencies onto spectra$freq")

	pks2 <- matrix(NA_real_, ncol = length(spectra$freq), nrow = nrow(pks))
	
	# Add a row with the frequency values, to be used later for diff calc
	
	pks2 <- rbind(spectra$freq, pks2)
	
	# dimnames(pks) <- list(
		# paste("D_", as.character(sprintf("%.3f", pks[,1])), sep = ""),
		# paste("F_", as.character(sprintf("%.3f", spectra$freq)), sep = ""))
	
	# Now fill the matrix
	
	for (i in 1:nrow(pks)) {
		for (j in 2:ncol(pks)) {
			if (is.na(pks[i,j])) next
			if (is.nan(pks[i,j])) next
			loc <- findInterval(pks[i,j], spectra$freq)
			pks2[i+1, loc] <- pks[i, j]
			}
		}
	
	# Step 3: Remove any columns that are all NA
	
	throw <- NA_integer_
	
	if (verbose) message("Removing columns that are all NA")
	
	for (i in 1:ncol(pks2)) {
		chk <- all(is.na(pks2[-1,i])) # avoiding reference row
		if (chk) throw <- c(throw, i)
		}
	
	throw <- throw[-1]

	if (length(throw) > 0) {
		pks2 <- pks2[,-throw]
		if (verbose) cat("Removed ", length(throw), " columns\n")
		}
		
	# Add back the driver peak info, which is not present in pks2
	# (needed as the source of collDiff in collapseRowOrCols)
	
	pks2 <- cbind(c(NA, pks[,1]), pks2)

	# Step 4: Repeat Step 1 except over the columns

	if (collapseCols) {
		pks2 <- collapseRowsOrCols(pks2, ind = 1, rows = FALSE, thres = colThres, verbose = verbose)
		}
	
	# Final clean up: convert NaN to NA (not sure where these come from)
	
	nan <- which(is.nan(pks2))
	pks2[nan] <- NA
	return(pks2[-1,]) # remove the reference row
	
	} # end of function

	
