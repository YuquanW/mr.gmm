.Estep <- function(beta_exp,
                   beta_out,
                   se_exp,
                   se_out,
                   t,
                   beta_hat_old,
                   tau_hat_old,
                   eta_hat_old,
                   pi_hat_old,
                   xi) {

  m <- length(beta_exp)
  nomi1 <- beta_out/se_out^2
  nomi2 <- beta_hat_old*beta_out/se_out^2 + beta_exp/se_exp^2

  ### alpha, gamma | L = 1

  Sigma_j1_inv11 <- 1/se_out^2 + tau_hat_old
  Sigma_j1_inv12 <- beta_hat_old/se_out^2
  Sigma_j1_inv22 <- beta_hat_old^2/se_out^2 + 1/se_exp^2 + eta_hat_old

  det_Sigma_j1_inv <- Sigma_j1_inv11*Sigma_j1_inv22 - Sigma_j1_inv12^2
  det_Sigma_j1 <- 1/det_Sigma_j1_inv


  Sigma_j1_11 <- Sigma_j1_inv22/det_Sigma_j1_inv
  Sigma_j1_12 <- -Sigma_j1_inv12/det_Sigma_j1_inv
  Sigma_j1_22 <- Sigma_j1_inv11/det_Sigma_j1_inv

  mu_j1_alpha <- Sigma_j1_11*nomi1 + Sigma_j1_12*nomi2
  mu_j1_gamma <- Sigma_j1_12*nomi1 + Sigma_j1_22*nomi2

  quad_j1 <- Sigma_j1_inv11*mu_j1_alpha^2 +
    2*Sigma_j1_inv12*mu_j1_alpha*mu_j1_gamma +
    Sigma_j1_inv22*mu_j1_gamma^2
  Phi_j1 <- pnorm(-t*se_exp/sqrt(1/eta_hat_old+se_exp^2), log.p = T)

  ### alpha, gamma | L = 2

  Sigma_j2_11 <- 1/tau_hat_old
  Sigma_j2_12 <- 0
  Sigma_j2_22 <- 1/Sigma_j1_inv22

  det_Sigma_j2 <- Sigma_j2_11*Sigma_j2_22

  mu_j2_alpha <- 0
  mu_j2_gamma <- Sigma_j2_22*nomi2

  quad_j2 <- 1/Sigma_j2_22*mu_j2_gamma^2
  Phi_j2 <- Phi_j1

  ### alpha, gamma | L = 3

  Sigma_j3_11 <- 1/Sigma_j1_inv11
  Sigma_j3_12 <- 0
  Sigma_j3_22 <- 1/eta_hat_old

  det_Sigma_j3 <- Sigma_j3_11*Sigma_j3_22

  mu_j3_alpha <- Sigma_j3_11*nomi1
  mu_j3_gamma <- 0

  quad_j3 <- 1/Sigma_j3_11*mu_j3_alpha^2
  Phi_j3 <- pnorm(-t, log.p = T)

  ### alpha, gamma | L = 4

  Sigma_j4_11 <- 1/tau_hat_old
  Sigma_j4_12 <- 0
  Sigma_j4_22 <- 1/eta_hat_old

  det_Sigma_j4 <- Sigma_j4_11*Sigma_j4_22

  mu_j4_alpha <- 0
  mu_j4_gamma <- 0

  quad_j4 <- 0
  Phi_j4 <- Phi_j3

  ### Update rho

  det_Sigma <- cbind(det_Sigma_j1,
                     det_Sigma_j2,
                     det_Sigma_j3,
                     det_Sigma_j4)

  quad <- cbind(quad_j1,
                quad_j2,
                quad_j3,
                quad_j4)

  Phi <- cbind(Phi_j1,
               Phi_j2,
               Phi_j3,
               Phi_j4)

  dist <- 0.5*(log(det_Sigma)+quad)+
    matrix(log(pi_hat_old), m, 4, byrow = T) -
    Phi

  dist_scale <- dist - apply(dist, 1, max)

  rho_jk_prop <- exp(dist_scale)
  rho_jk <- (rho_jk_prop+1e-300)^xi/rowSums((rho_jk_prop+1e-300)^xi)

  mu_jk_alpha <- cbind(mu_j1_alpha, mu_j2_alpha, mu_j3_alpha, mu_j4_alpha)
  mu_jk_gamma <- cbind(mu_j1_gamma, mu_j2_gamma, mu_j3_gamma, mu_j4_gamma)
  Sigma_jk_alpha <- cbind(Sigma_j1_11, Sigma_j2_11, Sigma_j3_11, Sigma_j4_11)
  Sigma_jk_gamma <- cbind(Sigma_j1_22, Sigma_j2_22, Sigma_j3_22, Sigma_j4_22)

  return(list(rho_jk = rho_jk,
              Sigma_j1_12 = Sigma_j1_12,
              mu_jk_alpha = mu_jk_alpha,
              mu_jk_gamma = mu_jk_gamma,
              Sigma_jk_alpha = Sigma_jk_alpha,
              Sigma_jk_gamma = Sigma_jk_gamma,
              det_Sigma = det_Sigma,
              quad = quad,
              Phi = Phi))
}

.jac <- function(theta0, lambda, beta_exp, beta_out, se_exp, se_out, t){
  m <- length(beta_exp)
  beta <- theta0[1]
  tau <- exp(theta0[2])
  eta <- exp(theta0[3])
  p <- exp(c(theta0[4:6], 0))/(sum(exp(theta0[4:6]))+1)

  logf1 <- rep(0, m)
  logf2 <- rep(0, m)
  dlogf1_dbeta <- rep(0, m)
  dlogf1_dtau <- rep(0, m)
  dlogf1_deta <- rep(0, m)
  d2logf1_d2beta <- rep(0, m)
  d2logf1_dbeta_dtau <- rep(0, m)
  d2logf1_dbeta_deta <- rep(0, m)
  d2logf1_d2tau <- rep(0, m)
  d2logf1_dtau_deta <- rep(0, m)
  d2logf1_d2eta <- rep(0, m)

  dlogf2_dbeta <- rep(0, m)
  dlogf2_deta <- rep(0, m)
  d2logf2_d2beta <- rep(0, m)
  d2logf2_dbeta_deta <- rep(0, m)
  d2logf2_d2eta <- rep(0, m)
  for (i in 1:m) {
    y <- c(beta_out[i], beta_exp[i])

    ## f1
    Sigma1 <- matrix(c(1/tau+beta^2/eta+se_out[i]^2,
                       beta/eta,
                       beta/eta,
                       1/eta+se_exp[i]^2), 2, 2)
    Sigma1_inv <- solve(Sigma1)
    logf1[i] <- -log(2*pi) - 0.5*log(det(Sigma1)) - 0.5*t(y)%*%Sigma1_inv%*%y

    dSigma1_dbeta <- matrix(c(2*beta/eta,
                              1/eta,
                              1/eta,
                              0), 2, 2)
    dSigma1_dtau <- matrix(c(-1/tau^2, 0,
                             0, 0), 2, 2)
    dSigma1_deta <- matrix(c(-beta^2/eta^2, -beta/eta^2,
                             -beta/eta^2, -1/eta^2), 2, 2)
    d2Sigma1_d2beta <- matrix(c(2/eta, 0,
                                0, 0), 2, 2)
    d2Sigma1_dbeta_deta <- matrix(c(-2*beta/eta^2, -1/eta^2,
                                    -1/eta^2, 0), 2, 2)
    d2Sigma1_d2tau <- matrix(c(2/tau^3, 0,
                               0, 0), 2, 2)
    d2Sigma1_d2eta <- matrix(c(2*beta^2/eta^3, 2*beta/eta^3,
                               2*beta/eta^3, 2/eta^3), 2, 2)
    dlogf1_dbeta[i] <- -0.5*sum(diag(Sigma1_inv%*%dSigma1_dbeta)) +
      0.5*t(y)%*%Sigma1_inv%*%dSigma1_dbeta%*%Sigma1_inv%*%y
    dlogf1_dtau[i] <- -0.5*sum(diag(Sigma1_inv%*%dSigma1_dtau)) +
      0.5*t(y)%*%Sigma1_inv%*%dSigma1_dtau%*%Sigma1_inv%*%y
    dlogf1_deta[i] <- -0.5*sum(diag(Sigma1_inv%*%dSigma1_deta)) +
      0.5*t(y)%*%Sigma1_inv%*%dSigma1_deta%*%Sigma1_inv%*%y

    d2logf1_d2beta[i] <- 0.5*sum(diag(Sigma1_inv%*%dSigma1_dbeta%*%Sigma1_inv%*%dSigma1_dbeta - Sigma1_inv%*%d2Sigma1_d2beta))-
      t(y)%*%Sigma1_inv%*%dSigma1_dbeta%*%Sigma1_inv%*%dSigma1_dbeta%*%Sigma1_inv%*%y +
      0.5*t(y)%*%Sigma1_inv%*%d2Sigma1_d2beta%*%Sigma1_inv%*%y
    d2logf1_dbeta_dtau[i] <- 0.5*sum(diag(Sigma1_inv%*%dSigma1_dtau%*%Sigma1_inv%*%dSigma1_dbeta))-
      t(y)%*%Sigma1_inv%*%dSigma1_dtau%*%Sigma1_inv%*%dSigma1_dbeta%*%Sigma1_inv%*%y
    d2logf1_dbeta_deta[i] <- 0.5*sum(diag(Sigma1_inv%*%dSigma1_deta%*%Sigma1_inv%*%dSigma1_dbeta - Sigma1_inv%*%d2Sigma1_dbeta_deta))-
      t(y)%*%Sigma1_inv%*%dSigma1_deta%*%Sigma1_inv%*%dSigma1_dbeta%*%Sigma1_inv%*%y +
      0.5*t(y)%*%Sigma1_inv%*%d2Sigma1_dbeta_deta%*%Sigma1_inv%*%y
    d2logf1_d2tau[i] <- 0.5*sum(diag(Sigma1_inv%*%dSigma1_dtau%*%Sigma1_inv%*%dSigma1_dtau - Sigma1_inv%*%d2Sigma1_d2tau))-
      t(y)%*%Sigma1_inv%*%dSigma1_dtau%*%Sigma1_inv%*%dSigma1_dtau%*%Sigma1_inv%*%y +
      0.5*t(y)%*%Sigma1_inv%*%d2Sigma1_d2tau%*%Sigma1_inv%*%y
    d2logf1_dtau_deta[i] <- 0.5*sum(diag(Sigma1_inv%*%dSigma1_deta%*%Sigma1_inv%*%dSigma1_dtau))-
      t(y)%*%Sigma1_inv%*%dSigma1_deta%*%Sigma1_inv%*%dSigma1_dtau%*%Sigma1_inv%*%y
    d2logf1_d2eta[i] <- 0.5*sum(diag(Sigma1_inv%*%dSigma1_deta%*%Sigma1_inv%*%dSigma1_deta - Sigma1_inv%*%d2Sigma1_d2eta))-
      t(y)%*%Sigma1_inv%*%dSigma1_deta%*%Sigma1_inv%*%dSigma1_deta%*%Sigma1_inv%*%y +
      0.5*t(y)%*%Sigma1_inv%*%d2Sigma1_d2eta%*%Sigma1_inv%*%y

    ## f2

    Sigma2 <- Sigma1 - matrix(c(1/tau,0,
                                0,0), 2, 2)
    Sigma2_inv <- solve(Sigma2)
    logf2[i] <- -log(2*pi) - 0.5*log(det(Sigma2)) - 0.5*t(y)%*%Sigma2_inv%*%y

    dSigma2_dbeta <- dSigma1_dbeta
    dSigma2_deta <- dSigma1_deta
    d2Sigma2_d2beta <- d2Sigma1_d2beta
    d2Sigma2_dbeta_deta <- d2Sigma1_dbeta_deta
    d2Sigma2_d2eta <- d2Sigma1_d2eta
    dlogf2_dbeta[i] <- -0.5*sum(diag(Sigma2_inv%*%dSigma2_dbeta)) +
      0.5*t(y)%*%Sigma2_inv%*%dSigma2_dbeta%*%Sigma2_inv%*%y
    dlogf2_deta[i] <- -0.5*sum(diag(Sigma2_inv%*%dSigma2_deta)) +
      0.5*t(y)%*%Sigma2_inv%*%dSigma2_deta%*%Sigma2_inv%*%y

    d2logf2_d2beta[i] <- 0.5*sum(diag(Sigma2_inv%*%dSigma2_dbeta%*%Sigma2_inv%*%dSigma2_dbeta - Sigma2_inv%*%d2Sigma2_d2beta))-
      t(y)%*%Sigma2_inv%*%dSigma2_dbeta%*%Sigma2_inv%*%dSigma2_dbeta%*%Sigma2_inv%*%y +
      0.5*t(y)%*%Sigma2_inv%*%d2Sigma2_d2beta%*%Sigma2_inv%*%y
    d2logf2_dbeta_deta[i] <- 0.5*sum(diag(Sigma2_inv%*%dSigma2_deta%*%Sigma2_inv%*%dSigma2_dbeta - Sigma2_inv%*%d2Sigma2_dbeta_deta))-
      t(y)%*%Sigma2_inv%*%dSigma2_deta%*%Sigma2_inv%*%dSigma2_dbeta%*%Sigma2_inv%*%y +
      0.5*t(y)%*%Sigma2_inv%*%d2Sigma2_dbeta_deta%*%Sigma2_inv%*%y
    d2logf2_d2eta[i] <- 0.5*sum(diag(Sigma2_inv%*%dSigma2_deta%*%Sigma2_inv%*%dSigma2_deta - Sigma2_inv%*%d2Sigma2_d2eta))-
      t(y)%*%Sigma2_inv%*%dSigma2_deta%*%Sigma2_inv%*%dSigma2_deta%*%Sigma2_inv%*%y +
      0.5*t(y)%*%Sigma2_inv%*%d2Sigma2_d2eta%*%Sigma2_inv%*%y


  }

  ## f3

  logf3 <- dnorm(beta_out, 0, sd = sqrt(1/tau+se_out^2), log = T) +
    dnorm(beta_exp, 0, sd = se_exp, log = T)
  dlogf3_dtau <- 0.5/(tau+se_out^2*tau^2)-0.5*beta_out^2/(1+se_out^2*tau)^2

  d2logf3_d2tau <- -0.5/(tau+se_out^2*tau^2)^2*(1+2*se_out^2*tau)+beta_out^2/(1+se_out^2*tau)^3*se_out^2

  ## f4

  logf4 <- dnorm(beta_out, 0, sd = se_out, log = T) +
    dnorm(beta_exp, 0, sd = se_exp, log = T)

  ## Phi
  q <- t*se_exp/sqrt(1/eta+se_exp^2)
  logphi1 <- pnorm(-q, log.p = T)
  logphi2 <- logphi1
  logphi3 <- pnorm(-t, log.p = T)
  logphi4 <- logphi3

  logimr <- dnorm(-q, log = T) - logphi1
  dphi1_deta <- 0.5*q*1/(1/eta + se_exp^2)*eta^(-2)*exp(logimr)
  dphi2_deta <- dphi1_deta
  d2phi1_d2eta <- 2*dphi1_deta^2 - dphi1_deta^2*q/exp(logimr) +
    0.75*exp(logimr)*q*1/(1/eta + se_exp^2)^2*eta^(-4) -
    2*dphi2_deta/eta
  d2phi2_d2eta <- d2phi1_d2eta

  ## 1st order derivate

  f_all <- exp(cbind(logf1 - logphi1,
                     logf2 - logphi2,
                     logf3 - logphi3,
                     logf4 - logphi4))

  pi_all <- matrix(p, m, 4, byrow = T)

  lik <- rowSums(pi_all*f_all)
  jac_beta <- 1/lik*(p[1]*f_all[, 1]*dlogf1_dbeta +
                       p[2]*f_all[, 2]*dlogf2_dbeta)
  jac_logtau <- tau/lik*(p[1]*f_all[, 1]*dlogf1_dtau +
                           p[3]*f_all[, 3]*dlogf3_dtau)
  jac_logeta <- eta/lik*(p[1]*f_all[, 1]*(dphi1_deta + dlogf1_deta) +
                           p[2]*f_all[, 2]*(dphi2_deta + dlogf2_deta))
  jac_logpi1 <- (f_all[, 1]/lik - 1)*p[1] + lambda/m*(1-4*p[1])
  jac_logpi2 <- (f_all[, 2]/lik - 1)*p[2] + lambda/m*(1-4*p[2])
  jac_logpi3 <- (f_all[, 3]/lik - 1)*p[3] + lambda/m*(1-4*p[3])

  ## 2nd order derivate
  jac_beta2 <- 1/lik*(p[1]*f_all[, 1]*(dlogf1_dbeta^2 + d2logf1_d2beta) +
                        p[2]*f_all[, 2]*(dlogf2_dbeta^2 + d2logf2_d2beta)) -
    jac_beta^2
  jac_beta_logtau <- tau/lik*p[1]*f_all[, 1]*(dlogf1_dtau*dlogf1_dbeta + d2logf1_dbeta_dtau) -
    jac_beta*jac_logtau
  jac_beta_logeta <- eta/lik*(p[1]*f_all[, 1]*((dlogf1_deta + dphi1_deta)*dlogf1_dbeta + d2logf1_dbeta_deta) +
                                p[2]*f_all[, 2]*((dlogf2_deta + dphi2_deta)*dlogf2_dbeta + d2logf2_dbeta_deta)) -
    jac_beta*jac_logeta
  jac_beta_logpi1 <- 1/lik*(p[1]*(1-p[1])*f_all[, 1]*dlogf1_dbeta - p[1]*p[2]*f_all[, 2]*dlogf2_dbeta) -
    jac_beta*jac_logpi1
  jac_beta_logpi2 <- 1/lik*(-p[1]*p[2]*f_all[, 1]*dlogf1_dbeta + p[2]*(1-p[2])*f_all[, 2]*dlogf2_dbeta) -
    jac_beta*jac_logpi2
  jac_beta_logpi3 <- 1/lik*(-p[1]*p[3]*f_all[, 1]*dlogf1_dbeta - p[2]*p[3]*f_all[, 2]*dlogf2_dbeta) -
    jac_beta*jac_logpi3

  jac_logtau2 <- tau^2/lik*(p[1]*f_all[, 1]*(dlogf1_dtau^2 + d2logf1_d2tau) +
                              p[3]*f_all[, 3]*(dlogf3_dtau^2 + d2logf3_d2tau)) -
    jac_logtau^2 + jac_logtau
  jac_logtau_logeta <- tau*eta/lik*p[1]*f_all[, 1]*((dlogf1_deta + dphi1_deta)*dlogf1_dtau+d2logf1_dtau_deta) -
    jac_logtau*jac_logeta
  jac_logtau_logpi1 <- tau/lik*(p[1]*(1-p[1])*f_all[, 1]*dlogf1_dtau - p[3]*p[1]*f_all[, 3]*dlogf3_dtau) -
    jac_logtau*jac_logpi1
  jac_logtau_logpi2 <- tau/lik*(-p[1]*p[2]*f_all[, 1]*dlogf1_dtau - p[3]*p[2]*f_all[, 3]*dlogf3_dtau) -
    jac_logtau*jac_logpi2
  jac_logtau_logpi3 <- tau/lik*(-p[1]*p[3]*f_all[, 1]*dlogf1_dtau + p[3]*(1-p[3])*f_all[, 3]*dlogf3_dtau) -
    jac_logtau*jac_logpi3

  jac_logeta2 <- eta^2/lik*(p[1]*f_all[, 1]*((dlogf1_deta + dphi1_deta)^2 - dphi1_deta^2 + d2logf1_d2eta + d2phi1_d2eta) +
                              p[2]*f_all[, 2]*((dlogf2_deta + dphi2_deta)^2 - dphi2_deta^2 + d2logf2_d2eta + d2phi2_d2eta)) -
    jac_logeta^2 + jac_logeta
  jac_logeta_logpi1 <- eta/lik*(p[1]*(1-p[1])*f_all[, 1]*(dphi1_deta + dlogf1_deta) - p[2]*p[1]*f_all[, 2]*(dphi2_deta + dlogf2_deta)) -
    jac_logeta*jac_logpi1
  jac_logeta_logpi2 <- eta/lik*(-p[1]*p[2]*f_all[, 1]*(dphi1_deta + dlogf1_deta) + p[2]*(1-p[2])*f_all[, 2]*(dphi2_deta + dlogf2_deta)) -
    jac_logeta*jac_logpi2
  jac_logeta_logpi3 <- eta/lik*(-p[1]*p[3]*f_all[, 1]*(dphi1_deta + dlogf1_deta) - p[2]*p[3]*f_all[, 2]*(dphi2_deta + dlogf2_deta)) -
    jac_logeta*jac_logpi3

  jac_logpi12 <- p[1]*(f_all[, 1]/lik - 1) - p[1]^2*(f_all[, 1]/lik - 1)*(f_all[, 1]/lik + 1) -4*lambda/m*p[1]*(1-p[1])
  jac_logpi1_logpi2 <- -p[1]*p[2]*(f_all[, 1]*f_all[, 2]/lik^2 - 1 - 4*lambda/m)
  jac_logpi1_logpi3 <- -p[1]*p[3]*(f_all[, 1]*f_all[, 3]/lik^2 - 1 - 4*lambda/m)


  jac_logpi22 <- p[2]*(f_all[, 2]/lik - 1) - p[2]^2*(f_all[, 2]/lik - 1)*(f_all[, 2]/lik + 1) -4*lambda/m*p[2]*(1-p[2])
  jac_logpi2_logpi3 <- -p[2]*p[3]*(f_all[, 2]*f_all[, 3]/lik^2 - 1 - 4*lambda/m)

  jac_logpi32 <- p[3]*(f_all[, 3]/lik - 1) - p[3]^2*(f_all[, 3]/lik - 1)*(f_all[, 3]/lik + 1) -4*lambda/m*p[3]*(1-p[3])

  jac2 <- c(mean(jac_beta2),
            mean(jac_beta_logtau),
            mean(jac_beta_logeta),
            mean(jac_beta_logpi1),
            mean(jac_beta_logpi2),
            mean(jac_beta_logpi3),
            mean(jac_logtau2),
            mean(jac_logtau_logeta),
            mean(jac_logtau_logpi1),
            mean(jac_logtau_logpi2),
            mean(jac_logtau_logpi3),
            mean(jac_logeta2),
            mean(jac_logeta_logpi1),
            mean(jac_logeta_logpi2),
            mean(jac_logeta_logpi3),
            mean(jac_logpi12),
            mean(jac_logpi1_logpi2),
            mean(jac_logpi1_logpi3),
            mean(jac_logpi22),
            mean(jac_logpi2_logpi3),
            mean(jac_logpi32))
  fi <- matrix(0, 6, 6)
  fi[lower.tri(fi, diag = T)] <- jac2
  fi <- fi + t(fi*lower.tri(fi))
  return(list(score = cbind(jac_beta, jac_logtau, jac_logeta, jac_logpi1, jac_logpi2, jac_logpi3),
              fi = -fi))

}


