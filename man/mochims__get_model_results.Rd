% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mochims__get_model_results.R
\name{mochims__get_model_results}
\alias{mochims__get_model_results}
\title{mochims__get_model_results}
\usage{
mochims__get_model_results(
  input_folder,
  normalisation_file,
  RT = 0.001987 * (273 + 24)
)
}
\arguments{
\item{input_folder}{path to input folder with model results (required)}

\item{normalisation_file}{path to normalisation file (required)}

\item{RT}{constant (default:0.001987*(273+24))}
}
\value{
data.table with model results
}
\description{
Get model results and coefficients.
}
