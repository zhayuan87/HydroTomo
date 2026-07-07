# =============================================================================
# Theis & Cooper-Jacob 解析解
#
# 承压含水层非稳定流抽水试验的经典解析解。
#
# 核心公式：
#   Theis:        s = (Q / 4πT) * W(u)
#   Cooper-Jacob: s = (Q / 4πT) * (-0.5772 - ln(u))  (u < 0.03)
#   其中 u = r²S / (4Tt)
#
# 符号：
#   s = 降深 drawdown [L]
#   Q = 抽水流量 pumping rate [L³/T]（正值 = 抽水）
#   T = 导水系数 transmissivity [L²/T]
#   S = 储水系数 storativity [-]
#   r = 观测孔距抽水井距离 [L]
#   t = 时间 [T]
# =============================================================================


# ---------------------------------------------------------------------------
# 1. Theis 井函数 W(u) —— 分段多项式近似
#    参考：Abramowitz & Stegun (1964), Handbook of Mathematical Functions
# ---------------------------------------------------------------------------

#' Theis well function (exponential integral)
#'
#' Computes the Theis well function \eqn{W(u) = \int_u^\infty e^{-x}/x \,dx}
#' using piecewise polynomial approximations from Abramowitz & Stegun (1964).
#'
#' @param u numeric vector, dimensionless parameter \eqn{u = r^2 S / (4 T t)}.
#' @return numeric vector of \eqn{W(u)} values.
#' @keywords internal
theis_well_function_small <- function(u) {
  a <- c(-0.57721566, 0.99999193, -0.24991055, 0.05519968,
         -0.00976004, 0.00107857)
  W <- -log(u) + a[1] + a[2]*u + a[3]*u^2 + a[4]*u^3 + a[5]*u^4 + a[6]*u^5
  return(W)
}

theis_well_function_large <- function(u) {
  a <- c(2.334733, 0.250621)
  b <- c(3.330657, 1.681534)
  num <- exp(-u) * (u^2 + a[1]*u + a[2])
  den <- u * (u^2 + b[1]*u + b[2])
  W <- num / den
  return(W)
}

theis_well_function <- function(u) {
  if (length(u) == 1) {
    if (u < 1) return(theis_well_function_small(u))
    else return(theis_well_function_large(u))
  }
  W <- numeric(length(u))
  idx_small <- u < 1
  idx_large <- u >= 1
  if (any(idx_small)) W[idx_small] <- theis_well_function_small(u[idx_small])
  if (any(idx_large)) W[idx_large] <- theis_well_function_large(u[idx_large])
  return(W)
}


# ---------------------------------------------------------------------------
# 2. Cooper-Jacob 井函数 —— 对数近似
# ---------------------------------------------------------------------------

cooper_jacob_well_function <- function(u) {
  return(-0.5772157 - log(u))
}


# ---------------------------------------------------------------------------
# 3. 降深计算
# ---------------------------------------------------------------------------

#' Theis drawdown
#'
#' Compute transient drawdown in a confined aquifer using the Theis (1935)
#' solution.
#'
#' @param Q pumping rate \eqn{[L^3/T]} (positive = extraction).
#' @param r radial distance from pumping well \eqn{[L]}.
#' @param T aquifer transmissivity \eqn{[L^2/T]}.
#' @param S aquifer storativity \eqn{[-]}.
#' @param t numeric vector of elapsed times \eqn{[T]}.
#' @return numeric vector of drawdown values \eqn{[L]}.
#' @export
#' @examples
#' t <- 10^seq(1, 5, length.out = 100)
#' s <- theis_drawdown(Q = 1.3e-3, r = 200, T = 1.5e-3, S = 2e-5, t = t)
#' plot(t, s, type = "l", log = "x", xlab = "t", ylab = "s")
theis_drawdown <- function(Q, r, T, S, t) {
  u <- (S * r^2) / (4 * T * t)
  W <- theis_well_function(u)
  s <- (Q / (4 * pi * T)) * W
  return(s)
}


#' Cooper-Jacob drawdown
#'
#' Compute transient drawdown using the Cooper-Jacob (1946) logarithmic
#' approximation to the Theis solution.  Valid only when \eqn{u < 0.03}.
#'
#' @param Q pumping rate \eqn{[L^3/T]} (positive = extraction).
#' @param r radial distance from pumping well \eqn{[L]}.
#' @param T aquifer transmissivity \eqn{[L^2/T]}.
#' @param S aquifer storativity \eqn{[-]}.
#' @param t numeric vector of elapsed times \eqn{[T]}.
#' @return numeric vector of drawdown values \eqn{[L]}.
#' @export
#' @examples
#' t <- 10^seq(2, 5, length.out = 100)
#' s <- cooper_jacob_drawdown(Q = 1.3e-3, r = 200, T = 1.5e-3, S = 2e-5, t = t)
#' plot(t, s, type = "l", log = "x")
cooper_jacob_drawdown <- function(Q, r, T, S, t) {
  u <- (S * r^2) / (4 * T * t)
  if (any(u > 0.03)) {
    warning("部分 u 值 > 0.03，Cooper-Jacob 近似可能不准确。建议使用 theis_drawdown()")
  }
  W <- cooper_jacob_well_function(u)
  s <- (Q / (4 * pi * T)) * W
  return(s)
}


# ---------------------------------------------------------------------------
# 4. 对数导数 (Bourdet 方法)
# ---------------------------------------------------------------------------

#' Log-derivative via central differences
#'
#' Compute \eqn{ds/d(\ln t)} using three-point central differences.
#'
#' @param t numeric vector of times.
#' @param s numeric vector of drawdown values (same length as \code{t}).
#' @param d half-window size in points (default 2).
#' @return list with \code{x} (times) and \code{y} (derivative values).
#' @keywords internal
log_derivative_central <- function(t, s, d = 2) {
  n <- length(t)
  if (n < 2*d + 1) stop("数据点不足以计算导数")

  x <- log(t)
  y <- s

  dt <- x[(d+1):(n-d)]
  dy <- numeric(length(dt))

  for (i in seq_along(dt)) {
    idx <- i + d
    dy[i] <- (y[idx+d] - y[idx-d]) / (x[idx+d] - x[idx-d])
  }

  return(list(x = exp(dt), y = dy))
}


#' Bourdet log-derivative with L-smoothing
#'
#' Compute the Bourdet diagnostic derivative \eqn{ds/d(\ln t)} with Gaussian
#' \eqn{L}-smoothing, commonly used in well-test analysis.
#'
#' @param t numeric vector of times.
#' @param s numeric vector of drawdown values.
#' @param L smoothing parameter in \eqn{\log_{10}} units (default 0.2).
#' @return list with \code{x} (times) and \code{y} (smoothed derivatives).
#' @export
#' @examples
#' t <- 10^seq(1, 5, length.out = 100)
#' s <- theis_drawdown(Q = 1.3e-3, r = 200, T = 1.5e-3, S = 2e-5, t = t)
#' d <- log_derivative_bourdet(t, s)
#' plot(t, s, type = "l", log = "xy"); lines(d$x, d$y, col = "red")
log_derivative_bourdet <- function(t, s, L = 0.2) {
  n <- length(t)
  x <- log10(t)

  dy <- numeric(n)
  dy[1] <- (s[2] - s[1]) / (x[2] - x[1])
  dy[n] <- (s[n] - s[n-1]) / (x[n] - x[n-1])

  for (i in 2:(n-1)) {
    dXi  <- x[i] - x[i-1]
    dXi1 <- x[i+1] - x[i]
    dy[i] <- ((s[i]-s[i-1])*dXi1/dXi + (s[i+1]-s[i])*dXi/dXi1) / (dXi + dXi1)
  }

  # Gaussian L-smoothing
  ws <- numeric(n)
  for (i in 1:n) {
    w <- exp(-((x[i] - x) / L)^2)
    ws[i] <- sum(w * dy) / sum(w)
  }

  return(list(x = t, y = ws))
}


# ---------------------------------------------------------------------------
# 5. Cooper-Jacob 直线法参数估算
# ---------------------------------------------------------------------------

#' Estimate aquifer parameters via Cooper-Jacob straight-line method
#'
#' Fits a straight line to \eqn{s} vs \eqn{\log_{10}(t)} and computes
#' \eqn{T} and \eqn{S} from the slope and intercept.
#'
#' @param t numeric vector of observed times \eqn{[T]}.
#' @param s numeric vector of observed drawdowns \eqn{[L]}.
#' @param Q pumping rate \eqn{[L^3/T]}.
#' @param r distance from pumping well \eqn{[L]}.
#' @return list with components \code{T}, \code{S}, \code{a} (slope),
#'   \code{t0} (intercept on log-time axis), \code{radius_influence},
#'   and \code{lm_fit} (the \code{\link{lm}} object).
#' @export
#' @examples
#' set.seed(123)
#' t_obs <- 10^seq(1.5, 4.5, length.out = 30)
#' s_true <- theis_drawdown(1.3e-3, 200, 1.5e-3, 2e-5, t_obs)
#' s_obs <- s_true + rnorm(length(s_true), 0, 0.02)
#' est <- estimate_parameters(t_obs, s_obs, Q = 1.3e-3, r = 200)
#' print(est$T)
#' print(est$S)
estimate_parameters <- function(t, s, Q, r) {
  fit <- lm(s ~ log10(t))
  a <- unname(coef(fit)[2])
  intercept <- unname(coef(fit)[1])
  t0 <- 10^(-intercept / a)

  T <- 0.1832339 * Q / a
  S <- 2.2458394 * T * t0 / r^2

  n <- length(t)
  Ri <- 2 * sqrt(T * t[n] / S)

  return(list(
    T = T, S = S,
    a = a, t0 = t0,
    radius_influence = Ri,
    lm_fit = fit
  ))
}
