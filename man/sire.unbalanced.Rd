\name{sire.unbalanced.qst}
\alias{sire.unbalanced}
\docType{data}
\title{
Trait dataset for unbalanced half-sib sire model
}
\description{
An example dataset of a neutral trait created from simulations of 15 populations under an island model. 900 individuals total come from an average of 4 dams per sire and an average of 12 sires per population for an overall average of 60 individuals per population.
}
\usage{data(sire.unbalanced)}
\format{
  A data frame with 900 observations on the following 4 variables.
  \describe{
    \item{\code{pop}}{a numeric vector}
    \item{\code{sire}}{a numeric vector}
    \item{\code{dam}}{a numeric vector}
    \item{\code{trait}}{a numeric vector}
  }
}
\details{
  \describe{
    \item{pop}{a column of population identifiers}
    \item{sire}{a column of unique sire identifiers}
    \item{dam}{a column of unique dam identifiers}
    \item{trait}{column of trait values per individual}
  }
}
\examples{
data(sire.unbalanced)
## maybe str(sire.unbalanced) ; plot(sire.unbalanced) ...
}
\keyword{datasets}
