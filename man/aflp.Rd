\name{aflp.fst}
\alias{aflp}
\docType{data}
\title{
AFLP marker dataset
}
\description{
An example dataset created from simulations of 15 neutral populations under an island model. Individuals are diploid with genotypes at 100 dominant, AFLP loci.
}
\usage{data(aflp)}
\format{
  A data frame with 100 observations on the following 31 variables.
  \describe{
    \item{\code{locus}}{a numeric vector}
    \item{\code{q_hat_pop1}}{a numeric vector}
    \item{\code{q_hat_pop2}}{a numeric vector}
    \item{\code{q_hat_pop3}}{a numeric vector}
    \item{\code{q_hat_pop4}}{a numeric vector}
    \item{\code{q_hat_pop5}}{a numeric vector}
    \item{\code{q_hat_pop6}}{a numeric vector}
    \item{\code{q_hat_pop7}}{a numeric vector}
    \item{\code{q_hat_pop8}}{a numeric vector}
    \item{\code{q_hat_pop9}}{a numeric vector}
    \item{\code{q_hat_pop10}}{a numeric vector}
    \item{\code{q_hat_pop11}}{a numeric vector}
    \item{\code{q_hat_pop12}}{a numeric vector}
    \item{\code{q_hat_pop13}}{a numeric vector}
    \item{\code{q_hat_pop14}}{a numeric vector}
    \item{\code{q_hat_pop15}}{a numeric vector}
    \item{\code{q_var_pop1}}{a numeric vector}
    \item{\code{q_var_pop2}}{a numeric vector}
    \item{\code{q_var_pop3}}{a numeric vector}
    \item{\code{q_var_pop4}}{a numeric vector}
    \item{\code{q_var_pop5}}{a numeric vector}
    \item{\code{q_var_pop6}}{a numeric vector}
    \item{\code{q_var_pop7}}{a numeric vector}
    \item{\code{q_var_pop8}}{a numeric vector}
    \item{\code{q_var_pop9}}{a numeric vector}
    \item{\code{q_var_pop10}}{a numeric vector}
    \item{\code{q_var_pop11}}{a numeric vector}
    \item{\code{q_var_pop12}}{a numeric vector}
    \item{\code{q_var_pop13}}{a numeric vector}
    \item{\code{q_var_pop14}}{a numeric vector}
    \item{\code{q_var_pop15}}{a numeric vector}
  }
}
\details{
	\describe{
		\item{locus}{each row is a locus}
		\item{q_hat}{columns 2 through 16 are the q_hat values for each locus in each population}
		\item{q_var}{columns 17 through 31 are the q_hat variances for each locus in each population}
	}
}
\examples{
data(aflp)
## maybe str(aflp) ; plot(aflp) ...
}
\keyword{datasets}
