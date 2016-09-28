\name{popdemo-package}
\alias{popdemo}
\docType{package}
\title{
    Demographic Modelling Using Projection Matrices
}
\description{
    This package provides a set of tools designed for modelling population demography
    using population projection matrix (PPM) models.  The focus is mainly on indices of transient dynamics,
    and the use of control theory (specifically transfer functions).  However,
    \code{popdemo} includes a number of tools focussing on other areas of demographic
    modelling.  Some tools may prove useful to anyone working with
    projection matrix models of any type, or indeed matrices in general.\cr
    \cr
    The package is designed to build on those tools already available in \code{popbio}.
}
\details{
  \tabular{ll}{
    Package: \tab popdemo\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2011-05-12\cr
    License: \tab GPL-2\cr
    LazyLoad: \tab yes\cr
  }
    Selected functions (see bottom of page for a link to the full index):\cr
    \cr
    Functions useful for working with matrices:\cr
    \code{popdemo} contains a number of tools to ease working with matrices
    (and specifically PPMs) in \R.
    \code{\link{Matlab2R}} allows coding of matrices in a Matlab style, which
    also facilitates import of multiple matrices simultaneously if comma-seperated
    files are used to import dataframes.  Its analogue, \code{\link{R2Matlab}},
    converts \R matrices to Matlab-style strings, for easier export.\cr
    \code{\link{is.matrix_primitive}}, \code{\link{is.matrix_irreducible}} and
    \code{\link{is.matrix_ergodic}} facilitate diagnosis of matrix properties
    pertaining to ergodicity.\cr
    \cr
    Population projection:\cr
    \code{popdemo} provides a simple means of projecting and plotting PPM models.
    \code{\link{project}} will project population dynamics and a plotting method
    is available via \code{\link{plot.projection}}.\cr
    \cr
    Indices of transient density:\cr
    Transient dynamics are important to study in PPM models, and \code{popdemo}
    provides the means to work with transient indices of PPM models.
    \code{\link{reactivity}}, \code{\link{firststepatt}} measure immediate transient
    density of a population (within the first time step).  \code{\link{maxamp}} and
    \code{\link{maxatt}} are near-term indices that measure the largest and smallest 
    transient dynamics a population may exhibit overall.  
    \code{\link{inertia}} measures asymptotic population density relative to stable
    state, and has many perturbation methods in the package (see below).  All indices
    can be calculated using specific population structures, as well as bounds on
    population size.\cr
    \cr
    Perturbation analysis:\cr
    Perturbation analysis is a key part of population studies. \code{popdemo} provides
    methods for nonlinear perturbation analysis of both asymptotic dynamics, using \code{\link{tfa}} and
    \code{\link{tfamatrix}}, and transient dynamics, using \code{\link{inertia.tfa}} and
    \code{\link{inertia.tfamatrix}}.  These all have associated plotting methods linked to them: see
    \code{\link{plot.tfa}} and \code{\link{plot.tfamatrix}}).  Sensitivity analyses are also available
    using transfer function methods: see \code{\link{tfsens}}, \code{\link{tfsensmatrix}}, 
    \code{\link{inertia.tfsens}} and \code{\link{inertia.tfsensmatrix}}.  (Traditional perturbation
    analyses are also available: \code{\link{sens}}, \code{\link{elas}}).\cr
    \cr
    Indices of convergence:\cr
    Information on the convergence of populations to stable state can be useful, and
    \code{popdemo} provides several means of analysing convergence.
    \code{dr} measures the damping ratio, and there are several distance measures available
    (see \code{\link{Keyfitz.delta}}, \code{\link{projection.distance}} and 
    \code{\link{Cohen.cumulative}}).  There is also a means of calculating convergence time
    through simulation: \code{\link{convergence.time}}.\cr
}
\author{
  Iain Stott <iainmstott@gmail.com>, 
  Dave Hodgson <D.J.Hodgson@exeter.ac.uk>, 
  Stuart Townley <S.B.Townley@exeter.ac.uk>
  Maintainer: Iain Stott <iainmstott@gmail.com>
}
