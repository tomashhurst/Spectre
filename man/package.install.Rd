% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package.install.R
\name{package.install}
\alias{package.install}
\title{package.install - a function to install packages required for Spectre.}
\usage{
See \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
}
\arguments{
\item{type}{DEFAULT = "general". If "general", then checks for the packages required for general Spectre usage. If "spatial", then checks for additional packages required for spatial analysis. If "ML", then checks for additional packages required for machine-learing functionality.}
}
\value{
returns an error message if one of the common use packages are not installed. Proceeds in order of package importance, and only the first error message encountered will be returned.
}
\description{
This function allows you to check to see if all the common use packages dependencies for Spectre are installed. Will only install if the package has not been installed, will not update packages.
}
\examples{
package.install()
package.install(type = 'spatial')
package.install(type = 'ML')

}
\references{
\url{https://sydneycytometry.org.au/spectre}
}
\author{
Thomas M Ashhurst, \email{thomas.ashhurst@sydney.edu.au}
}
