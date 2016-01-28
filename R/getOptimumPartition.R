#' @title
#' Find the Optimum Number of Partitions in a CLV Analysis
#'
#' @description
#' This function provides an automated means of selecting the
#' optimum number of clusters in a clv analysis, which is
#' useful in certain applications.  The procedure is to
#' identify the elbow in the delta plot derived from a 
#' call to \code{\link{CLV}}.  The elbow is the delta value
#' which is farthest (orthogonally) from a line connecting
#' the delta values for the largest number of partitions to the
#' one for the smallest number of partitions.
#'
#' @param clv An object of class \code{"clv"}.
#'
#' @param plot Logical. Should a plot be made?
#' The plot is identical to that produced by
#' \code{plot(clv_object, type = "delta")}
#' except that a red diamond marks the
#' optimal partition.
#' 
#' @return Integer.  The optimal partition.
#'
#' @references E. Vigneau and E. Qannari. "Clustering of variables around latent components."
#' \emph{Communications in Statistics - Simulation and Computation}, 32(4):1131-1150, 2003. 
#' For a good introduction, see E. Vigneau, M. Chen and E. Qannari
#' "ClustVarLV: An R Package for the Clustering of Variables Around Latent 
#' Variables" \emph{The R Journal} Vol. 7/2, December 2015 
#'
#' @seealso \code{\link{findElbow}} for more about how the process works.
#'
#' @section Warning:
#' As an automated method, this function returns a consistent answer.
#' However, that answer may not be ideal for any given data set.  It
#' is advisable to inspect the delta plots when using this function
#' to ensure that the concept here works well for your data sets.
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @export
#'
#' @importFrom graphics barplot
#'
#' @examples
#'
#' require("ClustVarLV")
#' data("AUPA_psycho")
#' res1<- CLV(AUPA_psycho, method = "directional", sX = TRUE)
#' plot(res1, type = "delta")
#' part <- getOptimumPartition(res1)
#' 
#' data("apples_sh")
#' res2 <- CLV(apples_sh$pref, method = "local")
#' plot(res2, type = "delta")
#' part <- getOptimumPartition(res2)
#' 
#' require("FactoMineR")
#' data(wine)
#' X.quanti <- wine[, 3:29]
#' X.quali <- wine[, 1:2]
#' Xbig <- cbind(scale(X.quanti), stand_quali((X.quali)))
#' res3 <- CLV(Xbig, method = "directional", sX = FALSE)
#' plot(res3, type = "delta")
#' part <- getOptimumPartition(res3)
#' 
getOptimumPartition <- function(clv, plot = TRUE) {

	if (!inherits(clv, "clv")) stop("non convenient object")
	if(is.null(clv$param$nmax)) stop("clv object must have been created by CLV, not CLV_kmeans.")
	
	# First, extract the delta values from clv object
	# Code to find delta values (called tempo here)
	# modified from plot.clv.R in ClustVarLV
	
	p <- clv$param$p
	nmax <- clv$param$nmax
	sbegin <- clv$param$sbegin
	results <- clv$tabres
	if (p > nmax) gpmax <- nmax
	if (p <= nmax) gpmax <- p
	
	tempo <- (results[(p-2):(p-gpmax+1),7] - results[(p-1):(p-gpmax+2),7])
	if (results[1,7] > 0) tempo <- c(tempo, sbegin-results[1,7])
	if (results[1,7] == 0) tempo <- c(tempo, results[p-gpmax,7] - results[p-gpmax+1,7])
	tempo[which(tempo < 0)] <- 0
	tempo <- tempo[(gpmax-1):1]
	names(tempo) <- paste((gpmax):2,"->",(gpmax-1):1)
	parts <- gpmax:2
	
	# Now, find the elbow index
	
	elbow <- findElbow(tempo, plot = FALSE, returnIndex = FALSE)
	elb <- which.max(elbow$dist)
	
	if (plot) {
		mid <- barplot(tempo, col = 4, xlab = "Nb clusters", ylab = "delta", 
			main = "Variation of criterion (after consolidation)",
			axisnames = TRUE, names.arg = names(tempo),
			las = 2, cex.names = 0.8, cex.main = 0.8)
		x <- mean(c(mid[elb], mid[elb + 1]))
		points(x, elbow$y[elb], pch = 18,
			col = "red", bg = "red", cex = 2)
		}
	
	return(parts[elb+1])
 	} # end of getOptimumPartition


