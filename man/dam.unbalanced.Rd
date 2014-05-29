\name{dam.unbalanced.qst}
\alias{dam.unbalanced}
\docType{data}
\title{
Trait dataset for unbalanced half-sib dam model
}
\description{
An example dataset of a neutral trait created from simulations of 15 populations under an island model. Offspring per dam are unbalanced across the breeding design with 6 dams per population, an average of 10 offspring per dam, and 900 individuals total.
}
\usage{data(dam.unbalanced)}
\format{
  A data frame with 900 observations on the following 3 variables.
  \describe{
    \item{\code{pop}}{a numeric vector}
    \item{\code{dam}}{a numeric vector}
    \item{\code{trait}}{a numeric vector}
  }
}
\details{
  \describe{
    \item{pop}{a column of population identifiers}
    \item{dam}{a column of unique dam identifiers}
    \item{trait}{column of trait values per individual}
  }
}
\examples{
data(dam.unbalanced)
## maybe str(dam.unbalanced) ; plot(dam.unbalanced) ...
}
\keyword{datasets}
