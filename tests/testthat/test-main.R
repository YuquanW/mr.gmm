test_that("Input validation", {
  attach(bmi.bmi)
  expect_error(mr.gmm(), "Missing input effect sizes or standard errors")
  expect_error(mr.gmm("a", beta.outcome, se.exposure, se.outcome), "Argument 'beta_exp' must be a numeric vector of length greater than 10.")
  expect_error(mr.gmm(beta.exposure[1:5], beta.outcome[1:5], se.exposure[1:5], se.outcome[1:5]), "Argument 'beta_exp' must be a numeric vector of length greater than 10.")
  expect_error(mr.gmm(beta.exposure, "b", se.exposure, se.outcome), "Argument 'beta_out' must be a numeric vector of the same length as 'beta_exp'")
  expect_error(mr.gmm(beta.exposure, beta.outcome[-1], se.exposure, se.outcome), "Argument 'beta_out' must be a numeric vector of the same length as 'beta_exp'")
  expect_error(mr.gmm(beta.exposure, beta.outcome, "c", se.outcome), "Argument 'se_exp' must be a numeric vector of the same length as 'beta_exp'")
  expect_error(mr.gmm(beta.exposure, beta.outcome, se.exposure[-1], se.outcome), "Argument 'se_exp' must be a numeric vector of the same length as 'beta_exp'")
  expect_error(mr.gmm(beta.exposure, beta.outcome, se.exposure, "d"), "Argument 'se_out' must be a numeric vector of the same length as 'beta_exp'")
  expect_error(mr.gmm(beta.exposure, beta.outcome, se.exposure, se.outcome[-1]), "Argument 'se_out' must be a numeric vector of the same length as 'beta_exp'")
})

test_that("Check initial estimate: gtools not installed",{
  skip_if(requireNamespace("gtools", quietly = TRUE))
  attach(bmi.bmi)
  expect_error(mr.gmm(beta.exposure[1:50], beta.outcome[1:50], se.exposure[1:50], se.outcome[1:50]), "gtools package must be installed when initial estimate is missing.")
})

test_that("Check initial estimate: gtools installed",{
  skip_if_not_installed("gtools")
  valid_inits <- list(beta = 1, tau = 200, eta = 200, pi = c(0.25, 0.25, 0.25, 0.25))
  invalid_inits_1 <- list(beta = 1, tau = 2, eta = 3)
  invalid_inits_2 <- list(beta = 1, tau = 2, eta = 3, pi = c(0.25, 0.25, 0.25, 0.25), extra = 5)
  invalid_inits_3 <- list(beta = "a", tau = 2, eta = 3, pi = c(0.25, 0.25, 0.25, 0.25))
  invalid_inits_4 <- list(beta = 1, tau = -2, eta = 3, pi = c(0.25, 0.25, 0.25, 0.25))
  invalid_inits_5 <- list(beta = 1, tau = 2, eta = -3, pi = c(0.25, 0.25, 0.25, 0.25))
  invalid_inits_6 <- list(beta = 1, tau = 2, eta = 3, pi = c(-1, 1, 0.5, 0.5))
  invalid_inits_7 <- list(beta = 1, tau = 2, eta = 3, pi = c(0.5, 0.5, 0.5, 0.5))
  invalid_inits_8 <- c(beta = 1, tau = 2, eta = 3, pi = c(0.25, 0.25, 0.25, 0.25))
  invalid_inits_9 <- list(beta = c(1, 2), tau = 2, eta = 3, pi = c(0.25, 0.25, 0.25, 0.25))
  invalid_inits_10 <- list(beta = 1, tau = c(1, 2), eta = 3, pi = c(0.25, 0.25, 0.25, 0.25))
  invalid_inits_11 <- list(beta = 1, tau = 2, eta = c(1, 2), pi = c(0.25, 0.25, 0.25, 0.25))
  invalid_inits_12 <- list(beta = 1, tau = 2, eta = 3, pi = c(0.5, 0.25, 0.25))
  attach(bmi.bmi)
  expect_invisible(mr.gmm(beta.exposure[1:50],
                          beta.outcome[1:50],
                          se.exposure[1:50],
                          se.outcome[1:50],
                          inits = valid_inits))
  expect_error(mr.gmm(beta.exposure[1:50],
                      beta.outcome[1:50],
                      se.exposure[1:50],
                      se.outcome[1:50],
                      inits = invalid_inits_1), "The list must contain exactly 'beta', 'tau', 'eta', and 'pi' as elements.")
  expect_error(mr.gmm(beta.exposure[1:50],
                      beta.outcome[1:50],
                      se.exposure[1:50],
                      se.outcome[1:50],
                      inits = invalid_inits_2), "The list must contain exactly 'beta', 'tau', 'eta', and 'pi' as elements.")
  expect_error(mr.gmm(beta.exposure[1:50],
                      beta.outcome[1:50],
                      se.exposure[1:50],
                      se.outcome[1:50],
                      inits = invalid_inits_3), "Initial estimate for 'beta' must be a numeric value.")
  expect_error(mr.gmm(beta.exposure[1:50],
                      beta.outcome[1:50],
                      se.exposure[1:50],
                      se.outcome[1:50],
                      inits = invalid_inits_4), "Initial estimate for 'tau' must be a positive numeric value.")
  expect_error(mr.gmm(beta.exposure[1:50],
                      beta.outcome[1:50],
                      se.exposure[1:50],
                      se.outcome[1:50],
                      inits = invalid_inits_5), "Initial estimate for 'eta' must be a positive numeric value.")
  expect_error(mr.gmm(beta.exposure[1:50],
                      beta.outcome[1:50],
                      se.exposure[1:50],
                      se.outcome[1:50],
                      inits = invalid_inits_6), "Each element of probability 'pi' must be in the interval \\(0, 1\\).")
  expect_error(mr.gmm(beta.exposure[1:50],
                      beta.outcome[1:50],
                      se.exposure[1:50],
                      se.outcome[1:50],
                      inits = invalid_inits_7), "Summation of probability 'pi' must be equal to 1.")
  expect_error(mr.gmm(beta.exposure[1:50],
                      beta.outcome[1:50],
                      se.exposure[1:50],
                      se.outcome[1:50],
                      inits = invalid_inits_8), "'inits' must be a list.")
  expect_error(mr.gmm(beta.exposure[1:50],
                      beta.outcome[1:50],
                      se.exposure[1:50],
                      se.outcome[1:50],
                      inits = invalid_inits_9), "Initial estimate for 'beta' must be a numeric value.")
  expect_error(mr.gmm(beta.exposure[1:50],
                      beta.outcome[1:50],
                      se.exposure[1:50],
                      se.outcome[1:50],
                      inits = invalid_inits_10), "Initial estimate for 'tau' must be a positive numeric value.")
  expect_error(mr.gmm(beta.exposure[1:50],
                      beta.outcome[1:50],
                      se.exposure[1:50],
                      se.outcome[1:50],
                      inits = invalid_inits_11), "Initial estimate for 'eta' must be a positive numeric value.")
  expect_error(mr.gmm(beta.exposure[1:50],
                      beta.outcome[1:50],
                      se.exposure[1:50],
                      se.outcome[1:50],
                      inits = invalid_inits_12), "Initial estimate for probability 'pi' must be a numeric vector of length 4.")
})
