\name{omnibus}
\alias{omnibus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

Implement omnibus method
}
\description{
omnibus is an integrative multi-omics analytical framework based on p-value weight adjustment in order to incorporate observations with incomplete intermediate variables in the analysis. By splitting the data into a complete set with full information and an incomplete set with missing measurements, this mechanism can be used to derive weights and weight-adjusted p-values from the two sets.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
omnibus(Y, E, M, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A vector of phenotypic outcome
%%     ~~Describe \code{x} here~~
}
\item{E}{The matrix of intermediate variables, missing values should be filled with NA.}
\item{M}{The matrix of variables of interest.}
\item{X}{The matrix of clinical covariates. This input is set to NULL by default.}
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
A list object.
\item{general}{p-values for M generated from general weighting scheme.}
\item{reverse}{p-values for M generated from reverse weighting scheme.}
\item{omnibus}{p-values for M generated from omnibus method combining the general weighting scheme and the reverse weighting scheme.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Wenda Zhang
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
The number of observations of Y,E,M,X should be the same. Each observation should come from the same subjects. NA can be used to represent missing values.
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
omnibus(Y,E,M,X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
