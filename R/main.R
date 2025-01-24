#' Conduct MR-GMM method
#'
#' Main function for MR-GMM
#'
#' @import stats
#'
#' @param beta_exp An m-dimensional vector of effect sizes from the exposure GWAS data.
#' @param beta_out An m-dimensional vector of effect sizes from the outcome GWAS data.
#' @param se_exp An m-dimensional vector of standard errors for `beta_exp` from the exposure GWAS data.
#' @param se_out An m-dimensional vector of standard errors for `beta_out` from the outcome GWAS data.
#' @param sel.p The p-value threshold for selecting SNPs. Default is 1.
#' @param inits A list containing the initial parameter estimates for the Gaussian Mixture Model (GMM).
#' The list must include four elements: `beta`, `tau`, `eta`, and `pi`.
#' - `beta`: The causal effect.
#' - `tau`: The inverse variance of latent pleiotropic effects on the outcome.
#' - `eta`: The inverse variance of latent SNP effects on the exposure.
#' - `pi`: The probabilities for the four clusters: invalid, valid, invalid & null, and null IVs.
#' The default is a random start.
#' @param lambda The penalty on `pi`. Default is the square root of m.
#' @param DA A logical value indicating whether the deterministic annealing procedure should be applied. Default is TRUE.
#' @param xi The initial temperature for the deterministic annealing (DA) procedure.
#' If `DA = TRUE`, `xi` defaults to 0.1; otherwise, it is set to 1.
#' @param rate Temperature scheduling parameter. In each iteration, `xi` is increased by multiplying `rate`. Default is 1.1.
#' @param maxit Maximum number of iterations of the (DA)EM algorithm. Default is 10,000.
#' @param tol The convergence tolerance. Default is 1e-7.
#'
#' @return A list containing the following components:
#' * `data` An m × 4 data frame of the input GWAS data.
#' * `beta.hat` The causal effect estimate.
#' * `beta.se` The standard error of `beta_hat`.
#' * `beta.p.value` The p-value from the two-tailed Wald test for \eqn{H_0: \beta=0}.
#' * `tau.hat` The estimated inverse variance of latent pleiotropic effects on the outcome.
#' * `eta.hat` The estimated inverse variance of latent SNP effects on the exposure.
#' * `rho` An m × 4 matrix indicating the estimated latent probabilities of each SNP belonging to the four clusters.
#' * `pi` The estimated probabilities for the four clusters.
#' * `ELBO` A vector recording the Evidence Lower BOund (ELBO) at each iteration, useful for assessing convergence.
#' * `cluster` An m × 2 data frame. The first column contains the most probable cluster assignments for each SNP, and the second column contains the corresponding estimated latent probabilities.
#' * `iter` The total number of iterations performed by the DAEM algorithm.
#'
#' @export
#'
#' @references Yuquan Wang, Yunlong Cao, Dong Chen, Dapeng Shi, Liwan Fu, et al. Pseudo-p-Value-Based Clumping Enhanced Proteome-wide Mendelian Randomization with Application in Identifying Coronary Heart Disease-Associated Plasma Proteins. \url{https://doi.org/10.1101/2025.01.13.25320450}.
#'
#' @examples
#' data(bmi.bmi)
#' attach(bmi.bmi)
#' bmi.bmi.res <- mr.gmm(beta.exposure, beta.outcome, se.exposure, se.outcome)
#' plot_wald(bmi.bmi.res)
#' plot_scatter(bmi.bmi.res)
#'
mr.gmm <- function(beta_exp,
                   beta_out,
                   se_exp,
                   se_out,
                   sel.p = 1,
                   inits = NULL,
                   lambda = NULL,
                   DA = T,
                   xi = ifelse(DA, 0.1, 1),
                   rate = ifelse(DA, 1.1, 1),
                   maxit = 10000,
                   tol = 1e-7){
  if (missing(beta_exp) || missing(beta_out) || missing(se_exp) || missing(se_out)) {
    stop("Missing input effect sizes or standard errors")
  }
  if (!is.numeric(beta_exp) || !is.vector(beta_exp) || length(beta_exp) < 10) {
    stop("Argument 'beta_exp' must be a numeric vector of length greater than 10.")
  }
  if (!is.numeric(beta_out) || !is.vector(beta_out) || length(beta_out) != length(beta_exp)) {
    stop("Argument 'beta_out' must be a numeric vector of the same length as 'beta_exp'.")
  }
  if (!is.numeric(se_exp) || !is.vector(se_exp) || length(se_exp) != length(beta_exp)) {
    stop("Argument 'se_exp' must be a numeric vector of the same length as 'beta_exp'.")
  }
  if (!is.numeric(se_out) || !is.vector(se_out) || length(se_out) != length(beta_exp)) {
    stop("Argument 'se_out' must be a numeric vector of the same length as 'beta_exp'.")
  }
  m <- length(beta_exp)
  t <- -qnorm(sel.p/2)
  if (is.null(lambda)) {
    lambda <- sqrt(m)
  }

  ## Random start
  if (is.null(inits)) {
    if (!requireNamespace("gtools", quietly = TRUE)) {
      stop("gtools package must be installed when initial estimate is missing.")
    }
    inits$beta <- 0
    inits$tau <- rgamma(1, 2, scale = 10000)
    inits$eta <- rgamma(1, 2, scale = 100)
    inits$pi <- gtools::rdirichlet(1, rep(1, 4))
  }
  required_names <- c("beta", "tau", "eta", "pi")

  if (!is.list(inits)) {
    stop("'inits' must be a list.")
  }

  if (!identical(sort(names(inits)), sort(required_names))) {
    stop("The list must contain exactly 'beta', 'tau', 'eta', and 'pi' as elements.")
  }

  if (!is.numeric(inits$beta) || length(inits$beta) != 1) {
    stop("Initial estimate for 'beta' must be a numeric value.")
  }

  if (!is.numeric(inits$tau) || length(inits$tau) != 1 || inits$tau < 0) {
    stop("Initial estimate for 'tau' must be a positive numeric value.")
  }

  if (!is.numeric(inits$eta) || length(inits$eta) != 1 || inits$eta < 0) {
    stop("Initial estimate for 'eta' must be a positive numeric value.")
  }

  if (any(inits$pi <= 0 | inits$pi >= 1)) {
    stop("Each element of probability 'pi' must be in the interval (0, 1).")
  }

  if (sum(inits$pi) != 1) {
    stop("Summation of probability 'pi' must be equal to 1.")
  }

  if (!is.numeric(inits$pi) || length(inits$pi) != 4) {
    stop("Initial estimate for probability 'pi' must be a numeric vector of length 4.")
  }

  beta_hat <- inits$beta
  tau_hat <- inits$tau
  eta_hat <- inits$eta
  pi_hat <- inits$pi

  ELBO <- 0

  for (i in 1:maxit) {
    beta_hat_old <- beta_hat
    tau_hat_old <- tau_hat
    eta_hat_old <- eta_hat
    pi_hat_old <- pi_hat

    ## V-E

    ve.res <- .Estep(beta_exp, beta_out, se_exp, se_out, t,
                    beta_hat_old, tau_hat_old, eta_hat_old, pi_hat_old,
                    xi)
    rho_jk <- ve.res$rho_jk
    Sigma_j1_12 <- ve.res$Sigma_j1_12
    mu_jk_alpha <- ve.res$mu_jk_alpha
    mu_jk_gamma <- ve.res$mu_jk_gamma
    Sigma_jk_alpha <- ve.res$Sigma_jk_alpha
    Sigma_jk_gamma <- ve.res$Sigma_jk_gamma
    det_Sigma <- ve.res$det_Sigma
    quad <- ve.res$quad
    Phi <- ve.res$Phi

    ## VB-M

    ### beta
    sigma_hat <- 1/(sum(rho_jk[, 1]/se_out^2*(mu_jk_gamma[, 1]^2+Sigma_jk_gamma[, 1])) +
                      sum(rho_jk[, 2]/se_out^2*(mu_jk_gamma[, 2]^2+Sigma_jk_gamma[, 2])))

    beta_hat <- sigma_hat*(sum(rho_jk[, 1]/se_out^2*(beta_out*mu_jk_gamma[, 1]-
                                                       mu_jk_alpha[, 1]*mu_jk_gamma[, 1]-
                                                       Sigma_j1_12))+
                             sum(rho_jk[, 2]/se_out^2*beta_out*mu_jk_gamma[, 2]))

    ### tau
    tau_hat <- m/sum(rho_jk*(mu_jk_alpha^2+Sigma_jk_alpha))

    ### eta
    eta_hat <- ((m/sqrt(eta_hat_old)+sum((rho_jk[, 1]+rho_jk[, 2])*t*se_exp*(1+se_exp^2*eta_hat_old)^(-3/2)*exp(dnorm(-t*se_exp/sqrt(1/eta_hat_old+se_exp^2), log = T) - pnorm(-t*se_exp/sqrt(1/eta_hat_old+se_exp^2), log.p = T))))/
                  sum(rho_jk*(mu_jk_gamma^2+Sigma_jk_gamma)))^2

    ### pi
    pi_hat = (colSums(rho_jk) + lambda)/(m + 4*lambda)

    ### ELBO
    ELBO_j <- rowSums(rho_jk*(0.5*(log(det_Sigma)+quad)+
                                matrix(log(pi_hat_old), m, 4, byrow = T)-
                                log(rho_jk) -
                                Phi)) + 0.5*log(tau_hat_old*eta_hat_old)

    ELBO <- c(ELBO, sum(ELBO_j) + lambda*sum(log(pi_hat_old)))
    xi <- min(c(1, xi*rate))

    if(abs(ELBO[i+1]-ELBO[i])/m < tol) {
      break
    }
  }

  ## Reparametrize
  theta0 <- c(beta_hat, log(tau_hat), log(eta_hat), log(pi_hat[1:3])-log(pi_hat[4]))

  ### Jacobian for VFE
  jac0 <- .jac(theta0, lambda, beta_exp, beta_out, se_exp, se_out, t)
  A0 <- jac0$fi
  A0_inv <- solve(A0)
  B0 <- crossprod(jac0$score)/m
  S0 <- colMeans(jac0$score)
  cov_mat0 <- A0_inv%*%B0%*%A0_inv/m
  cat("####################MR-GMM##################\n")
  cat("--------------------------------------------\n")
  cat(sprintf("Causal effect estimate: %.3f", theta0[1]),"\n")
  cat(sprintf("Standard error: %.3f", sqrt(cov_mat0[1, 1])),"\n")
  cat(sprintf("Wald test statistic: %.3f", theta0[1]/sqrt(cov_mat0[1, 1])),"\n")
  cat(sprintf("Two-tailed P-value: %.3f", 2*pnorm(-abs(theta0[1]/sqrt(cov_mat0[1, 1])))),"\n")
  cat(sprintf("Average IV strength = %.3f", mean(beta_exp^2/se_exp^2)-1), "\n")
  cat("--------------------------------------------")

  invisible(list(data = data.frame(beta_exp, beta_out, se_exp, se_out),
                 beta.hat = theta0[1],
                 beta.se = sqrt(cov_mat0[1, 1]),
                 beta.p.value = 2*pnorm(-abs(theta0[1]/sqrt(cov_mat0[1, 1]))),
                 tau.hat = exp(theta0[2]),
                 eta.hat = exp(theta0[3]),
                 rho = rho_jk,
                 pi = exp(c(theta0[4:6], 0))/(sum(exp(theta0[4:6]))+1),
                 ELBO = ELBO,
                 cluster = data.frame(cluster = apply(rho_jk, 1, which.max), prob = apply(rho_jk, 1, max)),
                 iter = i))
}

