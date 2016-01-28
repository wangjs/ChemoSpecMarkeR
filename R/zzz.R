.onAttach <- function (libname, pkgname){

  desc <- utils::packageDescription("ChemoSpecMarkeR")

  packageStartupMessage ("Package ",  desc$Package, ", version ", desc$Version, "\n\n",
       "\tYou've been warned: everything in this package\n\tis experimental and your results should be\n\tchecked very carefully.  There may be big errors!\n\tWe appreciate your feedback via Github issues.",
       sep = "")
}

