
## ChemoSpecMarkeR WORK IN PROGRESS

`ChemoSpecMarkeR` is an `R` package for identifying (bio)markers in NMR data.  It builds on the `Spectra` class in package `ChemoSpec`.  It is a work in progress and still needs a great deal of testing.  Don't trust it!

### How to install ChemoSpecMarkeR

#### To install from Github using R:

````r
install.packages("devtools")
library("devtools")
install_github(repo = "bryanhanson/ChemoSpecMarkeR@master")
library("ChemoSpecMarkeR")
````
If you use `@devel` you can get the development branch if it is available (and there may be other branches out there too).  Development branches have new, possibly incompletely tested features.  They may may also not be ready to pass checks and thus may not install automatically using the method above.  Check the news file to see what's up.

=======
Questions?  hanson@depauw.edu
