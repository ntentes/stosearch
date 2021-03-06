% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stosearch.R
\name{stochastic_search}
\alias{stochastic_search}
\title{Stochastic Search with resetting}
\usage{
stochastic_search(x0 = 40, t = 1000, searchers = 1,
  reset_type = "none", reset_rate = "optimal", n = 1e+05,
  output_data = TRUE, output_visual = TRUE)
}
\arguments{
\item{x0}{positive numeric initial value for all searchers (vector of length 1) (searchers will also reset to this value if resetting is chosen)}

\item{t}{time units alotted to searchers to complete the search}

\item{searchers}{number of searchers (integer from 1 to N) (warning: simulation will take a long time with a large number of searchers)}

\item{reset_type}{type of resetting, argument accepts "none", "stochastic" (each searcher is independently and randomly reset according to a poisson process with specified rate) or "deterministic" (all searchers are reset deterministically according to the specified rate).}

\item{reset_rate}{reciprocal of mean reset time. Accepts "optimal" or a positive numeric vector of length 1. "optimal" is rate = 1/(2*(x0)^2) which is the optimal resetting rate according to Bhat, De Bacco, Redner (2016).}

\item{n}{precision of brownian motion, as $n \rightarrow \inf$, brownian motion is calculated in continuous time.}

\item{output_data}{accepts TRUE or FALSE, if false, does not output list of collected data re: simulation.}

\item{output_visual}{accepts TRUE or FALSE, if false, does not output visual plot of the search}
}
\value{
if output_data==TRUE, then function returns a list of named vectors: "TargetFound" = T or F, "TimeUnits" = time elapsed until a searcher finds the target, "SearchCost" = TimeUnits * # of searchers, "SearchPaths" = a matrix of numerics where the rows (of length n) are the search path of each searcher, "NumResets" = number of times a searcher was reset (provides NA is no resetting is selected), "ResetTimes" = a vector of reset times (NA is provided if no resetting is selected). If output_visual==TRUE, a plot of the search is printed, with each searcher a different colour, a reset marked by a vertical dotted blue line, and a completed search marked by a vertical dashed green line.
}
\description{
Allows a user to specify a positive starting point and simulate (data or visually or both) 1 to N searchers moving according to standard brownian motion ($\mu = 0$ and $\sigma^2 = t$) with or without deterministic or stochastic resetting according to rates set by the user. Search will stop when any searcher finds the target (the origin).
}
