% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output_interpreter.R
\name{load_results}
\alias{loadAllResults}
\alias{loadAllResultsParallel}
\alias{load_results}
\title{Load All Single BBoB Results From a Folder}
\usage{
loadAllResults(usedFunctions, usedDimensions, path, algorithmName)

loadAllResultsParallel(usedFunctions, usedDimensions, path, algorithmName)
}
\arguments{
\item{usedFunctions}{the function identifiers for which bbob data exists and interpretation is desired.}

\item{usedDimensions}{the problem dimensions that have been used for optimization.}

\item{path}{the path of the folder containing the output data of a bbob experiment}

\item{algorithmName}{The identifier of the optimizer for detecting single output files in the specified folder}
}
\value{
\code{loadAllResults} returns an object \code{allResults}, which is a list of single bbob results as generated by \code{readOutput}
(see \code{\link{readOutput}} for details on the lists' elements).
}
\description{
\code{loadAllResults} loads the bbob results of a a number of single bbob output files generated by \code{\link{bbob_custom}} or 
\code{\link{bbob_custom_parallel}}.
}
\details{
\code{loadAllResults} detects and loads the files to be loaded/interpreted, 
i.e. \code{\link{readOutput}} is applied to every single file specified by \code{usedFunctions, usedDimensions, path, algorithmName}.
For parallelisation of the loading process, the user might invoke \code{loadAllResultsParallel} with the same parameter setup.
}

