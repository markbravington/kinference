\name{var_PLOD_kin}
\alias{var_PLOD_kin}
\title{Predict variance of PLOD for HCPs and HTPs}
\description{Aim is to work out how much your putative 2KPs (half-sibs, grandparent-grandchild, and full-thiatics--- here we'll say "HSP" for all of those) might be contaminated by 3KPs (eg half-thiatic pairs) and/or by 4KPs (eg half-cousins, etc) or, theoretically, by more remote kin). HSP-selection is presumably based on the pairwise PLODs for HSP:UP, taking all pairs where that PLOD exceeds some threshold. Given the allele freqs, the \emph{mean} PLOD is predictable when the truth is UP, HCP, HTP, or HSP. The variance is only predictable for UPs, though, because linkage makes loci non-independent for kin. However, an empirical variance can be estimated for HSPs based on the observed PLODs above some safe threshold (to exclude weaker kin), typically the mean PLOD when truth is HSP. Based on the empirical variance for HSPs and the analytical variance for UPs, we basically know how much linkage there might be, so we can predict the PLOD variances for the other kin-pair types. The wrinkle is that those more-remote variances also depend somewhat on the finer-scale organization of the genome, i.e. whether it's lots of chromosomes with no crossover, or a few chromosomes with lots of crossover. \code{var_PLOD_kin} therefore calculates two versions, one assuming the genome is entirely made up of equal-sized chromosome with zero crossover, and the other assuming the genome is a single chromosomes with crossover according to a memoryless random process. The output (basically, two variance estimates which ought to bound the true variance for the specified "contaminating" kin-type such as 3KP-- subject to statistical noise) can be fed into \code{\link{autopick_threshold}} (qv) to do what it says.
}
\usage{
var_PLOD_kin(
  linfo,
  emp_V_HSP = V_noX(C_equiv, 2),
  n_meio,
  debug = FALSE,
  C_equiv = NULL
)
}
\arguments{
\item{  linfo}{either a \code{snpgeno} object, or its "locinfo" attribute (or a fake one). The "locinfo" should be a dataframe with columns \code{e0}, \code{e1}, \code{v0}, \code{v1}, \code{count}. Each row is one "type" of locus, i.e., with roughly the same values of e/v 0/1, and \code{count} says how many such loci there are. e/v 0/1 are means and variances of the per-locus LOD (note no P) when the locus is or isn't co-inherited. See \emph{Details}}
\item{  emp_V_HSP}{empirical variance of PLOD for deffo HSPs. You're supposed to be running this on real data, so that \code{emp_V_HSP} is an actual number; however, for testing purposes, you can set up an artificial version via \code{C_equiv} below.}
\item{  n_meio}{Target number of meioses:2 for 2nd-order kin (e.g. HSPs; also GGPs and FTPs), 3 for 3rd-order (e.g. HTPs), etc. This is by far the main driver of variance, but technically the only one; see \emph{Subtypes of kin}.}
\item{  debug}{Logical flag. Defaults to FALSE.}
\item{  C_equiv}{for artificial test, with \code{emp_V_HSP} set to the no-crossover variance from \code{C_equiv} chromos (need not be integer). Ignored if \code{emp_V_HSP} is set.}
}
\details{
The "per-locus LOD" (whose properties are stored in the columns \code{e0}, \code{e1}, \code{v0}, \code{v1} in \code{linfo}) is created by calling \code{hsp_power} (qv). The normal use-case would be that you've done so with \code{k=0.5}, so that the (P)LOD pertains to HSP::UP comparisons. However, if you called it with \code{k=0.25} then the (P)LOD would be designed for HTP::UP comparisons, and so on. In fact, you could even hand-tweak the calculations to contain LODs for HTP::HSP comparisons, which \emph{might} in principle improve the resolution (but you'd have to fiddle manually; you could actually do it based on two calls to \code{hsp_power}, one with \code{k=0.5} and one with \code{k=0.25}, and manipulating the results). The other calculations in this function are "agnostic WRTO", ie not intrinsically dependendent on, the values of \code{e0/e1/v0/v1}, so the rest of the calcs should just work.

It's assumed that lots of loci are being used, so that the mix of loci on each "chromo", or the splatter of loci along the single "megachromo", always matches the overall population, on law-of-large-numbers grounds.

Stuff like uncertainly in allele frequencies, and in the PLOD variance for HSPs, needs to be accounted for externally, by repeatedly drawing from the posteriors and re-calculating the PLODs and re-running this function.

If the variance estimates show really good separation between the kin-pair types, then one could refine the "preliminary variance" step by reducing the super-high threshold (and assuming a truncated-Normal distribution). This might be worthwhile if the preliminary variance otherwise has to be based on a very small number of no-brainer HSPs. The "logical conclusion" of That Kind Of Thing is some kind of MLE involving estimating the population of different types of kin, and we really don't want to go there for now (since that should include the population dynamics shebang). In other words, we'd end up linking the genetic kin-finding model to the population dynamics model, which makes life statistically harder. And god knows it's hard enough. Anyway, if we were taking that approach, it might well be better to avoid PLODs altogether and instead go for inferences about the actual ppn of co-inherited loci, from which estimates-of-co-inherited-variance and inferences about kin-ppns can be made.
\subsection{Subtypes of kin}{Given the loci and the crossover rates, the PLOD variance for different kin-types is mainly determined by the number of meioses. However, at least for the with-crossover version, there is also \emph{some} effect of the \emph{type} of kin within a given order: GGPs and HSPs would have \emph{slightly} different variances. For CKMR purposes, the commonest type of kin of given order are those born closest in time, so the algorithm always uses the type with \emph{single} shared ancestor and minimax number of generations since shared ancestor. This means HSPs for \code{n_meio=2} (FTPs have _two shared ancestors; GGPs entail 2 generations of gap whereas HSPs have only 1), HTPs for \code{n_meio=3}, HC1Ps for \code{n_meio=4}, etc. If you really wanted to look at that, you could use the \code{V_allX} code inside this function, which takes two arguments \code{short} and \code{long} for the length of chains since the shared ancestor: for HSP, these are both 1, but for GGPs, one is 0 and the other is 2. But since the single-chromo equal-linkage-distance model is highly approximate anyway, do you really care?
}
}
\value{Matrix with two rows \code{V0} and \code{Vx}, and one column for each element of \code{n_meio} (which is always augmented to include 2), named "M2" etc. The two rows pertain respectively to the no-crossover multiple-chromosome scenario, and the single-chromosome multiple-crossover scenario. The matrix also has an attribute \code{info}, which is a numeric vector of elements named \code{V_UP}, \code{V_HSP}, \code{C_hat}, and \code{rho_hat} (note that \code{V_HSP} should duplicate the first column of the matrix) \code{C_hat} is estimated equivalent number of chromosomes for the no-crossover scenario, and \code{rho_hat} is per-locus crossover rate for the all-crossover scenario.
}
\examples{
# COMPLETELY MADE-UP e/v values! Nothing to do with genetics :)
var_PLOD_kin( data.frame( count=45, ev01= I( cbind( e0=-1, e1=2, v0=0.03, v1=0.02))), C_equiv=22, n_meio=3:4)
#        M2    M3    M4
#  V0 208.2 156.6  91.9
#  Vx 208.2 186.2 101.8
#  attr(,"info")
#      V_UP    V_HSP    C_hat  rho_hat
#    1.3500 208.2273  22.0000   0.2616
}
\keyword{misc}
