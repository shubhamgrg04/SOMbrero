\name{quality}
\alias{quality.somRes}
\alias{quality}

\title{Compute SOM algorithm quality criteria}
\description{The \code{quality} function computes several quality criteria for
the result of a SOM algorithm.}
\usage{
\method{quality}{somRes}(sommap,
        quality.type=c("all", "quantization", "topographic", 
                       "kaski.lagus"), \ldots)
}

\arguments{
  \item{sommap}{A \code{somRes} object (see \code{\link{trainSOM}} for 
  details).}
  \item{quality.type}{The quality type to compute. Three types are implemented: 
  \code{quantization}, \code{topographic} and \code{kaski.lagus}. The output of 
  the function is one of the three or all of them (using the option 
  \code{"all"}). Default value is the latter.}
  \item{\ldots}{Not used.}
}

\value{The \code{quality} function returns either a numeric value (if only one 
type is computed) or a list a numeric values (if all types are computed)

The quantization error calculates the mean (squared euclidean) distance between 
the sample vectors and their respective cluster prototypes. It is a decreasing 
function of the size of the map.

The topographic error is the simplest of the topology preservation measure: it
calculates the ratio of sample vectors for which the second best matching unit
is not in the direct neighborhood of the best matching unit.

The Kaski-Lagus error combines the quantization and topographic aspects into a
single measure. It is the sum of the mean (euclidean) distance between the 
sample vectors and their respective cluster prototypes, and the mean shortest 
path distance between the first and second best matching units.
}

\author{Laura Bendhaiba \email{laurabendhaiba@gmail.com}\cr
Madalina Olteanu \email{madalina.olteanu@univ-paris1.fr}\cr
Nathalie Villa-Vialaneix \email{nathalie.villa@toulouse.inra.fr}
}

\references{
Polzlbauer, G. (2004) Survey and comparison of quality measures for
self-organizing maps. In: \emph{Proceedings of the Fifth Workshop on Data
Analysis (WDA'04), Paralic, J., Polzlbauer, G., Rauber, A. (eds)} Sliezsky dom,
Vysoke Tatry, Slovakia: Elfa Academic Press, 67--82.

Kaski, S., Lagus, K. (1996) Comparing Self-Organizing Maps. in: 
\emph{Proceedings of International Conference on Artificial Neural Networks 
(ICANN) 96}, Springer, 809--814
}

\seealso{\code{\link{trainSOM}}, \code{\link{plot.somRes}}}

\examples{
my.som <- trainSOM(x.data=iris[,1:4])
quality(my.som, quality.type="all")
quality(my.som, quality.type="topographic")
}
\keyword{cluster}