# ---------- Khattree-Bahuguna's Univariate Skewness ---------------
#' Khattree-Bahuguna's Univariate Skewness
#' @name kbSkew
#' @description Compute Khattree-Bahuguna's Univariate Skewness.
#'
#' @param x a vector of original observations.
#'
#' @details
#' Given a univariate random sample of size \eqn{n} consist of observations \eqn{x_1, x_2, \ldots, x_n}, let \eqn{x_{(1)} \le x_{(2)} \le \cdots \le x_{(n)}} be the order statistics of \eqn{x_1, x_2, \ldots, x_n} after being centered by their mean. Define
#' \deqn{y_ i = \frac{x_{(i)} + x_{(n - i + 1)}}{2}}
#' and
#' \deqn{w_ i = \frac{x_{(i)} - x_{(n - i + 1)}}{2}}
#' The sample Khattree-Bahuguna's univariate skewness is defined as
#' \deqn{\hat{\delta} = \frac{\sum y_i^2}{\sum y_i^2 + \sum w_i^2}.}
#' Clearly, \eqn{0 \le \hat{\delta} \le 1}. Values close to zero indicate, low skewness while those close to one indicate the presence of high degree of skewness.
#'
#' @references
#' Khattree, R. and Bahuguna, M. (2019). An alternative data analytic approach to measure the univariate and multivariate skewness. \emph{International Journal of Data Science and Analytics}, Vol. 7, No. 1, 1-16.
#'
#' @return \code{kbSkew} gives the Khattree-Bahuguna's univariate skewness of the data.
#'
#' @examples
#' # Compute Khattree-Bahuguna's univariate skewness
#'
#' set.seed(2019)
#' x <- rnorm(1000) # Normal Distribution
#' kbSkew(x)
#'
#' set.seed(2019)
#' y <- rlnorm(1000, meanlog = 1, sdlog = 0.25) # Log-normal Distribution
#' kbSkew(y)
#'
#' @export
kbSkew <- function(x) {
  x <- as.vector(x)
  orderX <- x[order(x)] - mean(x)
  revOrdrX <- x[order(x, decreasing = TRUE)] - mean(x)
  even <- (orderX + revOrdrX)/2
  odd <- (orderX - revOrdrX)/2
  combine <- c(even, odd)
  est_sigma <- sum(even^2)/sum(combine^2)
  return(est_sigma)
}

# ---------- Khattree-Bahuguna's Multivariate Skewness ---------------
#' Khattree-Bahuguna's Multivariate Skewness
#' @name kbMvtSkew
#' @description Compute Khattree-Bahuguna's Multivariate Skewness.
#'
#' @param x a matrix of original observations.
#'
#' @details Let \eqn{\mathbf{X}=(X_1,\ldots,X_p)'} be the multivariate random vector and \eqn{(X_{i_1}, X_{i_2}, \ldots, X_{i_p})'} be one of the \eqn{p!} permutations of \eqn{(X_1,\ldots,X_p)'}. We predict \eqn{X_{i_j}} conditionally on subvector \eqn{(X_{i_1}, \ldots,X_{i_{j-1}})} and compute the corresponding residual \eqn{V_{i_j}} through a linear regression model for \eqn{j = 2, \cdots, p}. For \eqn{j=1}, we define \eqn{V_{i_1} = X_{i_1} - \bar{X}_{i_1}}, where \eqn{\bar{X}_{i_1}} is the mean of \eqn{X_{i_1}}. For \eqn{j \ge 2}, we have
#' \deqn{\hat{X}_{i_2}  = \beta_0 + \beta_1 X_{i_1}, \quad V_{i_2} = X_{i_2} - \hat{X}_{i_2}}
#' \deqn{\hat{X}_{i_3} = \beta_0 + \beta_1 X_{i_1} + \beta_2 X_{i_2}, \quad V_{i_3} = X_{i_3} - \hat{X}_{i_3}}
#' \deqn{\vdots}
#' \deqn{\hat{X}_{i_p} = \beta_0 + \beta_1 X_{i_1} + \beta_2 X_{i_2} + \cdots + \beta_{p-1} X_{i_{p-1}}, \quad V_{i_p} = X_{i_p} - \hat{X}_{i_p}.}
#'
#' We calculate the sample skewness \eqn{\hat{\delta}_{i_j}} of \eqn{V_{i_j}} by the sample Khattree-Bahuguna's univariate skewness formula (see details of \code{\link{kbSkew}} that follows) respectively for \eqn{j=1,\cdots,p} and define \eqn{\hat{\Delta}_{i} = \sum_{j=1}^{p} \hat{\delta}_{i_j}, i = 1, 2, \ldots, P} for all \eqn{P = p!} permutations of \eqn{(X_1,\ldots,X_p)'}. The sample Khattree-Bahuguna's multivariate skewness is defined as
#' \deqn{\hat{\Delta} = \frac{1}{P} \sum_{i=1}^{P} \hat{\Delta}_{i}.}
#'
#' @references
#' Khattree, R. and Bahuguna, M. (2019). An alternative data analytic approach to measure the univariate and multivariate skewness. \emph{International Journal of Data Science and Analytics}, Vol. 7, No. 1, 1-16.
#'
#' @return \code{kbMvtSkew} computes the Khattree-Bahuguna's multivairate skewness for a \eqn{p}-dimensional data.
#'
#' @seealso \code{\link{kbSkew}} for Khattree-Bahuguna's univariate skewness.
#'
#' @examples
#' # Compute Khattree-Bahuguna's multivairate skewness
#'
#' data(OlymWomen)
#' kbMvtSkew(OlymWomen[, c("m800","m1500","m3000","marathon")])
#'
#' @export
kbMvtSkew <- function(x) {
  x <- as.matrix(x)
  p <- dim(x)[2]

  perm <- function(v) {
    n <- length(v)
    if (n == 1) v
    else {
      X <- NULL
      for (i in 1:n) X <- rbind(X, cbind(v[i], perm(v[-i])))
      X
    }
  }

  cardP <- factorial(p)
  permX <- perm(seq(p))
  Delta_Hat <- rep(0, cardP)

  for (k in 1:cardP) {
    delta_hat <- rep(0, p)
    X = x[, permX[k, ]]
    meanX = apply(X, 2, mean)
    deMeanX = sweep(X, 2, meanX, FUN = "-")
    delta_hat[1] = kbSkew(X[, 1])

    for (j in 1:(p-1)) {
      V = meanX[j+1] + cov(X[,(j+1)], X[,1:j]) %*% solve(cov(X[,1:j], X[,1:j])) %*% t(deMeanX[,1:j])
      Z = X[,j+1] - as.vector(V)
      delta_hat[j+1] = kbSkew(Z)
    }
    Delta_Hat[k] = sum(delta_hat)
  }

  return(mean(Delta_Hat))
}

# ---------- Principal-component-based Khattree-Bahuguna's Multivariate Skewness ---------------
#' Principal-component-based Khattree-Bahuguna's Multivariate Skewness
#' @name pcKbSkew
#' @description Compute Principal-component-based Khattree-Bahuguna's Multivariate Skewness.
#'
#' @param x a matrix of original scale observations.
#' @param cor a logical value indicating whether the calculation should use the correlation matrix (\code{cor = TRUE}) or the covariance matrix (\code{cor = FALSE}). The default value is \code{cor = FALSE}.
#'
#' @details
#' Let \eqn{\mathbf{X} = X_1, \ldots, X_p} be a \eqn{p}-dimensional multivariate random vector. We compute the sample skewness for \eqn{p} principal components of \eqn{\mathbf{X}} respectively by the sample Khattree-Bahuguna's univariate skewness formula (see details of \code{\link{kbSkew}} that follows). Let \eqn{\eta_1, \eta_2, \ldots, \eta_p} be the \eqn{p} univariate skewnesses for \eqn{p} principal components. Principal-component-based Khattree-Bahuguna's multivariate skewness for a sample is then defined as
#' \deqn{\eta = \sum_{i=1}^{p} \eta_i.}
#'
#' @seealso \code{\link{kbSkew}} for Khattree-Bahuguna's univariate skewness.
#'
#' @references
#' Khattree, R. and Bahuguna, M. (2019). An alternative data analytic approach to measure the univariate and multivariate skewness. \emph{International Journal of Data Science and Analytics}, Vol. 7, No. 1, 1-16.
#'
#'
#' @return \code{pcKbSkew} gives the sample principal-component-based Khattree-Bahuguna's multivairate skewness.
#' @importFrom stats princomp
#'
#' @examples
#' # Compute principal-component-based Khattree-Bahuguna's multivairate skewness
#'
#' data(OlymWomen)
#' pcKbSkew(OlymWomen[, c("m800","m1500","m3000","marathon")])
#'
#' @export
pcKbSkew <- function(x, cor = FALSE) {
  x_pca <- princomp(x, cor = cor)
  pcSkew <- apply(x_pca$scores, 2, kbSkew)
  # print(pcSkew)
  return(sum(pcSkew))
}


# ---------- Mardia's Multivariate Skewness ---------------
#' Mardia's Multivariate Skewness
#' @name MardiaMvtSkew
#' @description Compute Mardia's Multivariate Skewness.
#'
#' @param x a matrix of original observations.
#'
#' @details
#' Given a \eqn{p}-dimensional multivariate random vector with mean vector \eqn{\boldsymbol{\mu}} and positive definite variance-covariance matrix \eqn{\boldsymbol{\Sigma}}, Mardia's multivariate skewness is defined as
#' \deqn{\beta_{1,p} = E[(\boldsymbol{X}_1 - \boldsymbol{\mu})' \boldsymbol{\Sigma}^{-1} (\boldsymbol{X}_2 - \boldsymbol{\mu})]^3,}
#' where \eqn{\boldsymbol{X}_1} and \eqn{\boldsymbol{X}_2} are independently and identically distributed copies of \eqn{\boldsymbol{X}}. For a multivariate random sample of size \eqn{n}, \eqn{\boldsymbol{x}_1, \boldsymbol{x}_1, \ldots, \boldsymbol{x}_n}, its sample version is defined as
#' \deqn{\hat{\beta}_{1,p} = \frac{1}{n^2} \sum_{i=1}^{n} \sum_{j=1}^{n} [(\boldsymbol{x}_i - \bar{\boldsymbol{x}})'\boldsymbol{S}^{-1} (\boldsymbol{x}_j - \bar{\boldsymbol{x}})]^3,}
#' where the sample mean \eqn{\bar{\boldsymbol{x}} = \frac{1}{n}\sum_{i=1}^{n} \boldsymbol{x}_i} and the sample variance-covariance matrix \eqn{\boldsymbol{S} = \frac{1}{n} \sum_{i=1}^{n} (\boldsymbol{x}_i - \bar{\boldsymbol{x}}) (\boldsymbol{x}_i - \bar{\boldsymbol{x}})'}. It is assumed that \eqn{n \ge p}.
#'
#' @references
#'
#' Mardia, K.V. (1970). Measures of multivariate skewness and kurtosis with applications. \emph{Biometrika}, 57(3), 519–530.
#'
#' @return \code{MardiaMvtSkew} gives the sample Mardia's multivairate skewness.
#' @importFrom stats cov
#' @examples
#' # Compute Mardia's multivairate skewness
#'
#' data(OlymWomen)
#' MardiaMvtSkew(OlymWomen[, c("m800","m1500","m3000","marathon")])
#'
#' @export
MardiaMvtSkew <- function(x) {
  x <- as.matrix(x)
  n <- dim(x)[1]
  sample.mean <- apply(x, 2, mean)
  cx <- sweep(x, 2, sample.mean, FUN = "-")
  S_Cov <- cov(x)*(n - 1)/n
  inv_S <- solve(S_Cov)
  quad_mat <- cx %*% inv_S %*% t(cx)
  quad_cube <- quad_mat^3
  beta_hat <- 1/(n^2) * sum(quad_cube)
  return(beta_hat)
}

# ---------- Pearson's Univariate Skewness ---------------
#' Pearson's coefficient of skewness
#' @name PearsonSkew
#' @description Compute Pearson's coefficient of skewness.
#'
#' @param x a vector of original observations.
#'
#' @details
#' Pearson's coefficient of skewness is defined as
#' \deqn{\gamma_1 = \frac{E[(X - \mu)^3]}{(\sigma^3)}}
#' where \eqn{\mu = E(X)} and \eqn{\sigma^2 = E[(X - \mu)^2]}. The sample version based on a random sample \eqn{x_1,x_2,\ldots,x_n} is defined as
#' \deqn{\hat{\gamma_1} = \frac{\sum_{i=1}^n (x_i - \bar{x})^3}{n s^3}}
#' where \eqn{\bar{x}} is the sample mean and \eqn{s} is the sample standard deviation of the data, respectively.
#'
#' @references
#'
#' Pearson, K. (1894). Contributions to the mathematical theory of evolution. \emph{Philos. Trans. R. Soc. Lond.} A 185, 71-110.
#'
#' Pearson, K. (1895). Contributions to the mathematical theory of evolution II: skew variation in homogeneous material. \emph{Philos. Trans. R. Soc. Lond.} A 86, 343-414.
#'
#' @return \code{PearsonSkew} gives the sample Pearson's univariate skewness.
#' @importFrom stats sd
#' @examples
#' # Compute Pearson's univariate skewness
#'
#' set.seed(2019)
#' x <- rnorm(1000) # Normal Distribution
#' PearsonSkew(x)
#'
#' set.seed(2019)
#' y <- rlnorm(1000, meanlog = 1, sdlog = 0.25) # Log-normal Distribution
#' PearsonSkew(y)
#'
#' @export
PearsonSkew <- function(x) {
  if (!is.vector(x, mode = "numeric")) {
    stop(sQuote("x"), " must be a vector of numeric values")
  }
  gamma <- mean((x - mean(x))^3)/(sd(x)^3)
  return(gamma)
}

# ---------- Bowley's Univariate Skewness ---------------
#' Bowley's Univariate Skewness
#' @name BowleySkew
#' @description Compute Bowley's Univariate Skewness.
#'
#' @param x a vector of original observations.
#'
#' @details
#' Bowley's skewness is defined in terms of quantiles as
#' \deqn{\hat{\gamma} = \frac{Q_3 + Q_1 - 2 Q_2}{Q_3 - Q_1}}
#' where \eqn{Q_i} is the \eqn{i}th quartile \eqn{i=1,2,3} of the data.
#'
#' @references
#' Bowley, A. L. (1920). \emph{Elements of Statistics}. London : P.S. King & Son, Ltd.
#'
#' @return \code{BowleySkew} gives the Bowley's univariate skewness of the data.
#' @importFrom stats sd
#' @importFrom stats quantile
#' @examples
#' # Compute Bowley's univariate skewness
#'
#' set.seed(2019)
#' x <- rnorm(1000) # Normal Distribution
#' BowleySkew(x)
#'
#' set.seed(2019)
#' y <- rlnorm(1000, meanlog = 1, sdlog = 0.25) # Log-normal Distribution
#' BowleySkew(y)
#'
#' @export
BowleySkew <- function(x) {
  if (!is.vector(x, mode = "numeric")) {
    stop(sQuote("x"), " must be a vector of numeric values")
  }
  Q <- as.vector(quantile(x, prob = c(0.25, 0.50, 0.75)))
  gamma = (Q[3] + Q[1] - 2 * Q[2]) / (Q[3] - Q[1])
  return(gamma)
}
