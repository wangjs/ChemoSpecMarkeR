#' @title
#' Analyze the Effect of Correlation on Clustering Criteria
#' 
#' @description
#' This function carries out an analysis designed to select the optimal
#' correlation value for use in a CLV analysis.  It is based roughly on Figure 5
#' of the reference.  The correlation matrix of the intensities in the
#' \code{Spectra} object are first filtered by the values in \code{corlvls}.
#' The resulting \code{Spectra} object is smaller as a result.  It is then
#' subjected to CLV analysis via \code{\link{clvSpectra}}.  In the process, the
#' size of the largest cluster is compared to the number of clusters, a criterion
#' recommended in the reference.  Please see Details.
#'
#' @section Details:
#' There are two main ways to run this function.  (1) If \code{clusters} is \code{NULL},
#' and if \code{corlvls} has several values, each value
#' produces a delta plot with the optimal number of clusters marked (chosen by
#' \code{\link{findElbow}}, which may not be ideal for your data).  A final plot
#' shows the criterion vs. the correlation levels.  If \code{corlvls} has only one
#' value, \code{altPlot} is automatically set to \code{TRUE} and a plot of the
#' maximum cluster size vs. the number of clusters is made instead.  This approach
#' allows you to survey the results.  If you don't like the results or wish to
#' specify the number of clusters directly, then use the next approach.  (2)
#' If \code{clusters} is not \code{NULL}, then the specified number of clusters
#' is used (instead of looking for the optimal number) and a single plot is made,
#' showing the criterion vs. the correlation levels.  This is faster because only
#' one partition is computed for each value of \code{corlvls} instead of checking
#' a number of partitions (which defaults to 20 for each value of \code{corlvls}).
#'
#' For the plots showing the criterion vs. the correlation levels, two vertical
#' dashed gray lines are show.  These mark the recommended minimum and maximum
#' values of the criterion, per the reference.  If they are the same value, a
#' message is sent to the console.
#' 
#' @param spectra An object of class \code{Spectra}.
#' Should be shrinkwrapped, see the examples.
#' 
#' @param corlvls Numeric. A vector of values in (0...1).
#' The correlation matrix of the \code{Spectra} object will be filtered for each
#' of these values before further processing.
#'
#' @param clusters Integer. An optional vector of values as long as
#' \code{corlvls}.  If not \code{NULL}, these values are
#' taken as the optimal number of clusters.
#'
#' @param altPlot Logical. Should the alternate plot be made?
#' Only applies if \code{clusters == NULL}.
#' The default plot shows the criterion vs. the specified correlation values.
#' The alternative plot shows the size of the largest cluster vs. the 
#' number of clusters.
#' \code{altPlot} is forced to \code{TRUE} if only one correlation value
#' is given (since it is not possible to make the default plot in this case).
#'
#' @section Note:
#' This function does a lot of heavy lifting and it is probably
#' best to run it in batch mode, especially for large data sets.
#' 
#' @return A data frame with four columns:
#' \describe{
#' \item{cor}{the specified correlation levels}
#' \item{no.clust}{either the
#' optimal number of clusters found or the number specified by
#' \code{clusters}}
#' \item{maxclust}{the size of the largest cluster}
#' \item{crit}{the value of the criterion, which is maxclust/no.clust}
#' }
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
#' @importFrom plyr count
#' @importFrom ClustVarLV get_partition
#' @importFrom graphics text mtext
#' @importFrom grDevices heat.colors
#'
#' @examples
#' \dontrun{
#' require("plyr")
#' require("ClustVarLV")
#' require("ChemoSpec")
#' data(metMUD2)
#' 
#' sw <- shrinkwrapSpectra(spectra = metMUD2, sigma = 0.02,
#' 	thres = 4.5, xlim = c(1.25, 1.74), xlab = "ppm", plot = FALSE,
#' 	main = expression(Shrinkwrap~of~metMUD2~sigma == 0.02))
#' 
#' cl <- seq(0.9, 0.99, 0.01) # 10 values, good for testing
#' tst <- checkClusterCrit(sw$Spectra, corlvls = cl)
#' tst2 <- checkClusterCrit(sw$Spectra, corlvls = cl, clusters = rep(4, length(cl)))
#' }

checkClusterCrit <- function(spectra, corlvls = seq(0.9, 0.999, 0.002),
	clusters = NULL, altPlot = FALSE) {
	
	chkSpectra(spectra)
	if (!is.null(clusters)) {
		if (length(clusters) != length(corlvls)) stop("Length of corlvls and clusters must match")
		}
	
	nr <- length(corlvls) # no. of runs
	res <- vector("list", nr) # temp storage
	names(res) <- paste("cor", corlvls, sep = "_")
	DF <- data.frame(cor = corlvls, # this will be returned
		no.clust = rep(NA_integer_, nr),
		maxclust = rep(NA_integer_, nr))
	
	# Filter the spectra
	for (i in 1:nr) {
		res[[i]] <- corFilterSpectra(spectra, lvl = corlvls[i])
		}
	
	# Compute the criterion (two main methods)
	# First is to check all partitions, which plots
	# the results as you go (slow)
	# If nr = 1, we an only make the alternate plot
	
	if (is.null(clusters)) {
		
		if (nr == 1) altPlot <- TRUE
		
		if (!altPlot) { # Criterion plot
			for (i in 1:nr) {
				clv <- clvSpectra(res[[i]])
				tmp <- getOptimumPartition(clv)
				mtext(paste("Round", i, "cor =", corlvls[i], sep = " "))
				DF$no.clust[i] <- tmp
				part <- get_partition(clv, K = tmp)
				cnt <- count(part)
				DF$maxclust[i] <- max(cnt$freq)
				}
		
			DF$crit <- DF$maxclust/DF$no.clust
					
			plot(crit ~ cor, data = DF, type = "l",
				xlab = "correlation filter", ylab = "size of largest cluster / no. of clusters",
				main = "Effect of Correlation Filter")
			minclust <- min(DF$maxclust)
			if (minclust <= 40) mincrit <- min(which(DF$maxclust <= 40))
			if (minclust > 40) {
				mincrit <- min(which(DF$maxclust == minclust))
				message("The smallest cluster size was ", minclust, sep = " ")
				}
			maxcrit <- which.min(DF$crit)
			if (mincrit == maxcrit) message("min and max crit are the same")
			abline(v = DF$cor[c(mincrit, maxcrit)], lty = 2, col = "gray")
			}
		
		if (altPlot) { # Max clust size vs no. clust plot
			
			# Need to build a larger DF
			nm <- 20 # corresponds to nmax in CLV
			DF <- data.frame(cor = rep(corlvls, each = nm), # this will be returned
				no.clust = rep(1:nm, nr),
				maxclust = rep(NA_integer_, nr*nm),
				group = rep(letters[1:nr], each = nm))
				
			ind <- 1		
			for (i in 1:nr) {
				clv <- clvSpectra(res[[i]])
				for (j in 1:nm) {
					part <- get_partition(clv, K = j)
					cnt <- count(part)
					DF$maxclust[ind] <- max(cnt$freq)	
					ind <- ind + 1			
					}
				}
			
			myc <- rev(heat.colors(nr))
			plot(maxclust ~ no.clust, data = DF, type = "n",
				xlab = "no. of clusters (partitions)", ylab = "size of largest cluster",
				main = "Largest Cluster vs. No. of Clusters",
				ylim = c(0, max(DF$maxclust)), xlim = c(1, nm+5))
			for (i in 1:length(unique(DF$group))) {
				DF2 <- subset(DF, DF$group == levels(DF$group)[i])
				lines(maxclust ~ no.clust, data = DF2, col = myc[i])
				}
			y <- subset(DF, DF$no.clust == nm)$maxclust
			x <- rep(c(nm, nm + 2, nm + 4), length.out = length(y))
			text(x, y, labels = sprintf("%4.3f", corlvls), cex = 0.6, pos = 4)
			abline(h = 40, lty = 2, col = "gray")
			DF <- DF[,-4] # remove temporary group
			}
		
		}

	# 2nd option only checks the requested partitions (much quicker)
	
	if (!is.null(clusters)) {
		for (i in 1:nr) {
			clv <- clvSpectra(res[[i]], mode = "kmeans", clust = clusters[i])
			DF$no.clust[i] <- clusters[i]
			part <- get_partition(clv)
			cnt <- count(part)
			DF$maxclust[i] <- max(cnt$freq)
			}
	
		DF$crit <- DF$maxclust/DF$no.clust
		
		plot(crit ~ cor, data = DF, type = "l",
			xlab = "correlation filter", ylab = "size of largest cluster / no. of clusters",
			main = "Effect of Correlation Filter")
		minclust <- min(DF$maxclust)
		if (minclust <= 40) mincrit <- min(which(DF$maxclust <= 40))
		if (minclust > 40) {
			mincrit <- min(which(DF$maxclust == minclust))
			message("The smallest cluster size was ", minclust, sep = " ")
			}
		maxcrit <- which.min(DF$crit)
		if (mincrit == maxcrit) message("min and max crit are the same")
		abline(v = DF$cor[c(mincrit, maxcrit)], lty = 2, col = "gray")
		}


	return(DF)
	}
