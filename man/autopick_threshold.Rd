\name{autopick_threshold}
\alias{autopick_threshold}
\title{PLOD threshold for HSPs}
\description{This function proposes a PLOD threshold for excluding almost all 3rd-order kin, and computes the associated False-Negative Probability (i.e., that a true HSP will have a PLOD below that threshold).
}
\usage{
autopick_threshold(
  x,
  kin,
  fitrange_PLOD,
  FPtol_pairs,
  use4th,
  selecto = c("ML", "paranoid"),
  NVAR = 10,
  plot_bins = NULL,
  shading_density = 10,
  want_all_results = FALSE,
  ...
)
}
\arguments{
\item{  x}{a \code{snpgeno} or its \code{locinfo} attribute. Must already have been prepared by running \code{kin_power}.}
\item{  kin}{a dataframe of "close-ish" kin-pairs and their PLODs, presumably from running \code{find_HSPs}; must have a column "PLOD".}
\item{  fitrange_PLOD}{two numbers, specifying the range of PLODs from \code{kin} to use in \emph{fitting} (though all are \emph{plotted}, by default)}
\item{  FPtol_pairs}{how many expected False-Positive 3rd-order kin should the threshold exclude?}
\item{  use4th}{whether to allow for 4th-order kin when fitting.}
\item{  selecto}{whether to choose the threshold based on the best mixture fit ("ML"), or the most conservative ("paranoid").}
\item{  NVAR}{how many variances to try, between the limits set by \code{var_PLOD_kin}}
\item{  plot_bins}{bin-width for histogram plotting. Default NULL means no plot.}
\item{  shading_density}{By default, all PLODs in \code{kin} will be included in the histogram, even though only a subset are used in fitting. The histogram bars \emph{not} used in fitting (i.e., below \code{fitrange_PLOD[1]}) will be lightened in colour, according to this parameter. Results are graphics-device-dependent, so you may need to experiment away from the default; larger numbers usually mean lighter shading. Setting \code{shading_density=NA} should result in a light transparent rectangle covering the entire LHS of the graph, which you might prefer. You can also set \code{xlim} as usual, to remove those left-hand bars altogether.}
\item{  want_all_results}{if TRUE, return dataframe(s) containing results for each variance explored. This lets you examine "sensitivity".}
\item{  ...}{other parameters passed to \code{hist}, eg \code{xlim}, \code{ylim}, \code{col}. Many others will be ignored, and some will cause problems.}
}
\details{
The rationale comes from fitting a mixture distribution to observed PLODs within some range that is expected to contain only 2nd, 3rd, and \emph{perhaps} a few 4th order kin. The threshold is then "chosen" (or proposed; it's really up to you) so that the expected number of false-positives from 3rd-order kin-pairs (i.e., with PLODs above the threshold) matches whatever you decide. A histogram with expected values is plotted (unless you tell it not to).

Means and variances of the mixture components are automatically set in advance, so the mixture-fit only has to estimate the proportion of kin-pairs of each type. The means are easily calculated from kinship coefficients, based on allele frequencies. Variances, however, also depend strongly on the degree of linkage between loci, and to some extent on the \emph{nature} of the linkage (more chromosomes, or more crossovers?). This is handled internally by the function \code{\link{var_PLOD_kin}} (qv), which uses the observed "overdispersion" of PLODs for a subset of \emph{definite} 2nd-order kin to place bounds on the variances of 3rd and 4th order kin (based on two extreme assumptions about the \emph{nature} of linkage). The code of \code{autopick_threshold} then explores different variances within those bounds and

Despite the name, \emph{you} still have to supply sensible values for a couple of parameters, based on looking at your data and understanding what you are trying to do. So it's not \emph{completely} automated- and never will be! Choosing a threshold is \emph{not} an "optimization process" with explicit bias/variance tradeoffs; rather, it's about ensuring that you have adequate "engineering tolerance" in the next stage of CKMR. The False-Negative Probability will, to a great extent, compensate for the choice of threshold (i.e. removing any bias in the fitted CKMR model) \emph{unless} you set the threshold too low, and thus end up with some 3rd-order kin-pairs in your set of "definite 2nd-orders".
\subsection{Fitrange and use4th}{If you are using an HSP-oriented PLOD, then the range of PLODs you fit to should extend from somewhere above 0 (which is verrry close to the expected PLOD for 3rd-order kin), up to the RHS of the HSP bump, but clearly not so far as to include any FSPs and POPs. If you set the lower limit high enough, then you don't need to worry about 4th-order kin intruding (their contribution would be negligible), so you can get away with fitting a 2-component mixture (simpler, less to go wrong...) by setting \code{use4th=FALSE}. But if you push the lower range closer to 0 (which does give you a larger sample size for fitting), then you might need to set \code{use4th=TRUE}. The (substantial) downside of doing that, is that there are often more PLODs close to 0 than near the HSP mean, so the mixture-fit (which has make more assumptions when also using 4th-orders- and those assumptions may not be perfect) will "concentrate its efforts" on getting a good fit near 0, rather than near the 2nd-order mean which is what we really need. It is worth experimenting.
}
}
\value{The proposed threshold, with lots of attributes. You can use those to calculate False-Neg Probabilities for \emph{other} possible thresholds, as per \emph{Examples}; you \emph{don't} have to accept the one that is proposed here. Threshold choice is \emph{up to you} (and not the fault of kinference)!
}
\seealso{\code{\link{kin_power}}, \code{\link{var_PLOD_kin}}
}
\examples{
dropbears1 <- kin_power( dropbears, k = 0.5)
hsps <- find_HSPs( dropbears1, keep_thresh = -10)
histoPLOD( hsps, log=TRUE) ## observe HSP - POP gap centred ~ PLOD = 120.
## Use that air-gap to set fitrange_PLOD.
thresh <- autopick_threshold( x= dropbears1, kin= hsps,
  fitrange_PLOD= c(0, 120), FPtol_pairs= 1, use4th= TRUE, plot_bins= 5)
thresh ## the 2nd order / 3rd order threshold value
attr(thresh, "info")["Pr_FNeg"] ## ...
## ... the estimated 2KP false-neg rate, given that threshold
}
\keyword{misc}
