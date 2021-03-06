% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findElbow.R
\name{findElbow}
\alias{findElbow}
\title{Find the Elbow in a Curve}
\usage{
findElbow(y, plot = FALSE, returnIndex = TRUE)
}
\arguments{
\item{y}{Numeric vector of y values for the curve.}

\item{plot}{Logical. Should a plot be made?}

\item{returnIndex}{Logical. Should the return value
be the index of the elbow point?}
}
\value{
If \code{returnIndex = TRUE}, the index of
the elbow point.  If \code{returnIndex = FALSE},
a data frame containing an index values (x),
the y values passed to the function, and the
the orthogonal distances of the y values from
the line connecting the first and last points.
\code{which.max(data_frame_name$dist)} will give the index of 
the elbow point.
}
\description{
This utility function finds the elbow in a curve which is concave
relative to a line drawn between the first and last points.
The elbow is defined as the point with the greatest
orthogonal distance from that line.
}
\section{Warning}{

This function makes some simple checks that the data is concave
as defined above.  Even so, it may give
answers in some cases that are not valid.  Please check
on typical data that you encounter to verify that it works
in your cases.
}
\examples{
tmp <- findElbow(c(0.9, 1.1, 1.1, 1.9, 2.5, 2.8, 4.9, 8.5),
	plot = TRUE) # wandering
tmp <- findElbow(c(0.9, 1.0, 1.2, 1.3, 1.5, 1.5, 10.0, 22.0),
	plot = TRUE) # late rise
tmp <- findElbow(c(2, 4, 6, 8, 10, 12, 14, 16)^2,
	plot = TRUE) # gradual, no obvious break

# Not the usual way to choose the number of PCs:
library("chemometrics")
data(glass)
pca <- prcomp(glass)
eigensum <- sum(pca$sdev * pca$sdev)
vv <- 100 * (pca$sdev * pca$sdev/eigensum)
cs <- cumsum(vv)
tmp <- findElbow(vv, plot = TRUE)
tmp <- findElbow(cs, plot = TRUE)

}
\author{
Bryan A. Hanson, DePauw University. \email{hanson@depauw.edu}
}
\references{
The concept of this function is based on the
clever idea in the
Stackoverflow post at stackoverflow.com/a/2022348/633251
and relies on code posted at
paulbourke.net/geometry/pointlineplane/pointline.r
(referenced in the SO post).  Minor modifications
to the code were made to that code in order to vectorize it.
}

