% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliary.R
\name{download_PWM}
\alias{download_PWM}
\title{Download PWM from Jaspar or Hocomoco}
\usage{
download_PWM(id, pseudocount = NULL)
}
\arguments{
\item{id}{string; The ID of the motif to download}

\item{pseudocount}{numeric; value of pseudocount to add to every entry of the PWM (to avoid 0 or 1 counts)}
}
\value{
A list of two items, a PWM matrix and a name
}
\description{
Download PWM from Jaspar or Hocomoco
}
\examples{
download_PWM("ALX1_MOUSE.H11MO.0.B")

}
