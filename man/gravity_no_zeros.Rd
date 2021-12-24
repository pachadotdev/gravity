\name{gravity_no_zeros}
\alias{gravity_no_zeros}
\docType{data}
\title{
  Gravity dataset without zero trade flows
}
\description{
  An edited version of the full gravity dataset: The "square" gravity dataset for all possible pairs of Countries worldwide, 1948-2006, which is used in the article: Head, K., T. Mayer and J. Ries, 2010, "The erosion of colonial linkages after independence". Journal of International Economics, 81(1):1-14 (lead article).
  Dataset \code{gravity_no_zeros} corresponds to the dataset without zero trade flows, \code{gravity_zeros}, on the other hand, includes zero trade flows.
}
\usage{data("gravity_no_zeros")}
\format{
  A data frame with 17088 observations on the following 10 variables.
  \describe{
    \item{\code{iso_o}}{ISO-Code of country of origin}
    \item{\code{iso_d}}{ISO-Code of country of destination}
    \item{\code{distw}}{weighted distance}
    \item{\code{gdp_o}}{GDP of country of origin in million dollars}
    \item{\code{gdp_d}}{GDP of country of destination in million dollars}
    \item{\code{rta}}{regional trade agreement dummy}
    \item{\code{flow}}{trade flow in million dollars}
    \item{\code{contig}}{contiguity dummy}
    \item{\code{comlang_off}}{common official language dummy}
    \item{\code{comcur}}{common currency dummy}
  }
}
\source{
  \url{https://sites.google.com/site/hiegravity/data-sources}
}
\references{
  Head, K. and Mayer, T. (2014). Chapter 3 - gravity equations: 
  Workhorse,toolkit, and cookbook. 
  In Gita Gopinath, E. H. and Rogoff, K., editors, 
  Handbook of International Economics, volume 4 of Handbook of International 
  Economics, pages 131-195. Elsevier.
  (\href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook})

  Head, K., T. Mayer and J. Ries, 2010, "The erosion od colonial linkages after 
  independence". Journal of International Economics, 81(1):1-14 
  (lead article).
}
\examples{
  data(gravity_no_zeros)
  str(gravity_no_zeros)
}
\keyword{datasets}
