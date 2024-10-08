% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_summary_functions.R
\name{extract_summary_functions}
\alias{extract_summary_functions}
\title{Extract Summary Functions}
\usage{
extract_summary_functions(
  mxFDAobject,
  r_vec = seq(0, 100, by = 10),
  extract_func = c(univariate, bivariate),
  summary_func = c(Kest, Lest, Gest),
  markvar,
  mark1,
  mark2 = NULL,
  edge_correction,
  empirical_CSR = FALSE,
  permutations = 1000
)
}
\arguments{
\item{mxFDAobject}{object of class \code{mxFDA}}

\item{r_vec}{Numeric vector of radii over which to evaluate spatial summary functions. Must begin at 0.}

\item{extract_func}{Defaults to univariate, which calculates univariate spatial summary functions. Choose bivariate for bivariate spatial summary functions.}

\item{summary_func}{Spatial summary function to calculate. Options are c(Kest, Lest, Gest) which denote Ripley's K, Besag's L, and nearest neighbor G function, respectively.}

\item{markvar}{The name of the variable that denotes cell type(s) of interest. Character.}

\item{mark1}{Character string that denotes first cell type of interest.}

\item{mark2}{Character string that denotes second cell type of interest for calculating bivariate summary statistics. Not used when calculating univariate statistics.}

\item{edge_correction}{Character string that denotes the edge correction method for spatial summary function. For Kest and Lest choose one of c("border", "isotropic", "Ripley", "translate", "none"). For Gest choose one of c("rs", "km", "han")}

\item{empirical_CSR}{logical to indicate whether to use the permutations to identify the sample-specific complete spatial randomness (CSR) estimation. If there are not enough levels present in \code{markvar} column for permutations, the theoretical will be used.}

\item{permutations}{integer for the number of permtuations to use if empirical_CSR is \code{TRUE} and exact CSR not calculable}
}
\value{
an object of class \code{mxFDA} containing the corresponding spatial summary function slot filled. See \code{\link[=make_mxfda]{make_mxfda()}} for object structure details.
}
\description{
Function to extract spatial summary functions from the \code{Spatial} slot of an \code{mxFDA} object
}
\details{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

Complete spatial randomness (CSR) is the estimation or measure of a spatial summary function when the points or cells in a sample are randomly distributed,
following no clustering or dispersion pattern. Some samples do have artifacts that may influence what CSR is under the distribution of points as they
are found in the sample such as large regions of missing points or possibly in the case of tissue sections, necrotic tissue where cells are dead. Theoretical
CSR requires points have an equal chance of occurring anywhere in the sample that these artifacts violate, necessitating the need to estimate or
calculate what this CSR would be for each sample independently. Previously Wilson et al. had demonstrated cases in which sample-specific
CSR was important over the use of the theoretical in calculating how much the observed deviates from expected.
}
\examples{
#load ovarian FDA object
data('ovarian_FDA')

#run function
ovarian_FDA = extract_summary_functions(ovarian_FDA, r_vec = 0:100,
                                        extract_func = univariate,
                                        summary_func = Gest,
                                        markvar = "immune",
                                        mark1 = "immune",
                                        edge_correction = "rs")

}
\references{
Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
Fast covariance estimation for high-dimensional functional data.
\emph{Statistics and Computing}, 26, 409-421.
DOI: 10.1007/s11222-014-9485-x.

Wilson, C., Soupir, A. C., Thapa, R., Creed, J., Nguyen, J., Segura, C. M.,
Gerke, T., Schildkraut, J. M., Peres, L. C., & Fridley, B. L. (2022).
Tumor immune cell clustering and its association with survival in African
American women with ovarian cancer. PLoS computational biology, 18(3),
e1009900. https://doi.org/10.1371/journal.pcbi.1009900

Creed, J. H., Wilson, C. M., Soupir, A. C., Colin-Leitzinger, C. M., Kimmel, G. J.,
Ospina, O. E., Chakiryan, N. H., Markowitz, J., Peres, L. C., Coghill, A., & Fridley, B. L. (2021).
spatialTIME and iTIME: R package and Shiny application for visualization and analysis of
immunofluorescence data. \emph{Bioinformatics} (Oxford, England), 37(23), 4584–4586.
https://doi.org/10.1093/bioinformatics/btab757

\code{\link[spatstat.explore:Kest]{spatstat.explore::Kest()}}

\code{\link[spatstat.explore:Gest]{spatstat.explore::Gest()}}

\code{\link[spatstat.explore:Lest]{spatstat.explore::Lest()}}

\code{\link[spatstat.explore:Kcross]{spatstat.explore::Kcross()}}

\code{\link[spatstat.explore:Gcross]{spatstat.explore::Gcross()}}

\code{\link[spatstat.explore:Lcross]{spatstat.explore::Lcross()}}
}
\author{
Julia Wrobel \email{julia.wrobel@emory.edu}

Alex Soupir \email{alex.soupir@moffitt.org}
}
