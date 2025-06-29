% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synpar_filter.R
\name{synpar_filter}
\alias{synpar_filter}
\title{The SyNPar procedure}
\usage{
synpar_filter(
  X,
  y,
  fdr_value,
  best_lambda = NULL,
  B_reps = NULL,
  dist_type = "normal",
  model_type
)
}
\arguments{
\item{X}{n-by-p matrix or data frame of predictors.}

\item{y}{n-vector of response variables. For cox models, there should be two columns:
the first column is the survival time and the second column is the censoring indicator. For graphical models,
there is no y.}

\item{fdr_value}{Numeric value FDR levels, must be between 0 and 1.}

\item{best_lambda}{Regularization parameter (lambda) used in the model.
If not specified by the user (default: \code{NULL}), it will be automatically
determined based on the chosen \code{model_type} using cross-validation.
For linear models, the default is \code{0.5*fit$lambda.min}.
When \eqn{n} is small and \eqn{p > n}, it is recommended to use \code{fit$lambda.min}.
For generalized linear models, the default is \code{0.5*fit$lambda.1se}. For cox models, the default is
\code{fit$lambda.min}. For graphical models, the default is the one corresponding to the maximum AIC from cross-validation}

\item{B_reps}{Number of repetitions for computing the correction factor. The correction factor is calculated B_reps times using different random seeds,
and the maximum value across all repetitions is used. If not specified by the user (default: \code{NULL}),
it will be automatically determined based on the chosen \code{model_type}. For linear models,
the default is 5. For generalized linear models, the default is 2. For cox models, the default is 2.
For graphical models, the default is 2.}

\item{dist_type}{Type of distribution for the error term. Applicable only to linear models.
The default is \code{"normal"}.}

\item{model_type}{Type of model, must be one of "linear", "glm", "cox", "graphical".}
}
\value{
Returns a list containing the following elements:
\describe{
\item{threshold}{computed selection threshold}
\item{selected}{named vector of selected variables}
\item{statistic}{computed test statistics}
}
}
\description{
This function runs the full SyNPar procedure to identify significant variables under a user-specified FDR threshold.
}
\examples{
n <- 2000
p <- 1000
s <- 30
Amp <- 0.25
rho <- 0.8
Theta.8 <- toeplitz(rho^(0:(p - 1)))
X <- mvrnorm(n, rep(0, p), Sigma = Theta.8)
X <- scale(X)
beta <- rep(0, p)
beta[1:s] <- sample(c(-Amp, Amp), s, replace = TRUE)
true_labels <- beta != 0
y <- X \%*\% beta + rnorm(n)
y <- y - mean(y)
result_normal <- synpar_filter(
  X, y, fdr_value = 0.1, best_lambda = NULL,
  dist_type = "normal", model_type = "linear"
)
y <- X \%*\% beta + rt(n, df = 3)
y <- y - mean(y)
result_t <- synpar_filter(
  X, y, fdr_value = 0.1, best_lambda = NULL, B_reps = NULL,
  dist_type = "t", model_type = "linear"
)
result_normal$selected
result_t$selected

}
