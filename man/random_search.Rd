\name{random_search}
\alias{random_search}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
random_search(fun, maxFE = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fun}{
%%     ~~Describe \code{fun} here~~
}
  \item{maxFE}{
%%     ~~Describe \code{maxFE} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (fun, maxFE = NULL) 
{
    if (!"smoof" \%in\% rownames(installed.packages())) 
        install.packages("smoof")
    require(smoof)
    ub = 5
    lb = -5
    dimensions = length(getLowerBoxConstraints(fun))
    result = "Start"
    if (is.null(maxFE)) 
        maxFE = 10000
    bestFitness = Inf
    for (i in 1:maxFE) {
        solution = runif(dimensions, lb, ub)
        fitness = fun(solution)
        if (fitness < bestFitness) {
            bestFitness = fitness
        }
        Fopt = getGlobalOptimum(fun)$value
        if ((i\%\%10) == 0) 
            result = c(result, paste(i, i, (bestFitness - Fopt)))
    }
    result = c(result, "End")
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
