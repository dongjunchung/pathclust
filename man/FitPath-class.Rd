% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R, R/methods_FitPath.R
\docType{class}
\name{FitPath-class}
\alias{FitPath-class}
\alias{predict,FitPath-method}
\alias{plot,FitPath,missing-method}
\alias{plot,FitPath-method}
\title{An S4 class for pathway-level analysis results}
\usage{
\S4method{predict}{FitPath}(object, newx = NULL, cuts = NULL)

\S4method{plot}{FitPath,missing}(x, y, type = "KM")
}
\arguments{
\item{object}{output of selectPath function}

\item{newx}{test gene expression data}

\item{cuts}{cut points to determine risk groups}

\item{x}{output of selectPath function}

\item{y}{(unused argument)}

\item{type}{type of plot. one of Kaplan-Meier curves (\code{type="KM"}), ROC (\code{type="ROC"}), or hazard ratio (\code{type="HR"}).}
}
\description{
An S4 class for pathway-level analysis results

Risk group prediction.

Plot Kaplan-Meier curves (\code{type="KM"}), ROC (\code{type="ROC"}), or hazard ratio (\code{type="HR"}).
}
\section{Slots}{

\describe{
\item{\code{PathAll}}{all pathways}

\item{\code{PathSelected}}{pathways selected by LASSO}

\item{\code{coef}}{LASSO coefficient estimates for pathway latent components}

\item{\code{fitGene}}{gene-level analysis results}
}}

