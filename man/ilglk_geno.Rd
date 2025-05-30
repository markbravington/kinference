\name{ilglk_geno}
\alias{ilglk_geno}
\title{Check individual multilocus genotypes for typicality}
\description{
\code{ilglk_geno} computes the per-sample log-likelihood across its entire genotype, i.e. sum log Pr g(i,l); and compares the distribution across individuals to the theoretical distribution given allele frequencies. Some mismatch is normal (and can arise just from noise in allele-frequency estimates), but substantial mismatch is bad. You get to define "substantial". \code{ilglk_geno} can also detect outlier individuals, usually with lglks that are much too low rather than too high; I'm not sure what could generate "too typical a genome" at the individual level.

Genotype encoding (see \code{\link{diplos}} etc) must be one of 4-way, 3-way, or 6-way. The \code{useN} field is honoured.
}
\usage{ilglk_geno(snpg, hist_pars = list(), showPlot = TRUE)
}
\arguments{
\item{  snpg}{a \code{snpgeno} (6-way genotype)}
\item{  hist_pars}{\code{list()} passed to \code{hist} for controlling histogram, e.g. \code{hist_pars=list(xlim=c(-12000, -6000))}, or use \code{FALSE} to not plot.}
\item{  showPlot}{show the histogram? Defaults to TRUE, but overrideen by \code{hist_pars=FALSE}.}
}
\details{
You can use \code{locator(1)} to click the histogram to figure out where to adjust the \code{xlim/ylim} values to change the range of the data to inspect more closely- ie you then re-run the function with its \code{...hist_par} argument set accordingly.

Currently, the SPA calcs are a wee bit slow because of heavy use of \code{vecless} which in version 1.0 is slightly sluggish. The lglks themselves are computed in C and are blisteringly fast. If the SPA line (expected distro) doesn't appear, let us know; something will need fixing! There might e.g. be too many loci, so that the calculation is falling over.

We haven't added any formal uh-oh criteria yet; that could be done via the SPA, as in \code{dump_badhetz_fish}. However, reading off from the graph is probably fine. In practice, the observed and predicted \code{ilglk_geno} distributions seldom match exactly anyway (whereas they are often pretty close for \code{\link{hetzminoo_fancy}}), so that theoretical tail-probability criteria don't make sense.
}
\value{Vector of log-likelihood for each individual; also usually (but optionally), a histogram of log-likelihood values across individuals.
}
\examples{
## get rid of really bad loci
pvals <- check6and4( bluefin, thresh_pchisq_6and4 = c( 0.001, 0.0001))
bluefin_1 <- bluefin[ , pvals$pval4 > 0.01] # drastic QC!
## check for samples that are not like the others
ilglks <- ilglk_geno( bluefin_1)
## looks like anything with a lglk < -1030 is definitely abnormal
bluefin_2 <- bluefin_1[ ilglks > -1030,]
ilglks <- ilglk_geno( bluefin_2) ## much better, but not perfect -
## see other cleaning steps in the vignette
}
\keyword{misc}
