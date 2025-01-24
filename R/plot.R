#' Visualize the results of MR-GMM
#'
#' @description
#' Visualization functions for MR-GMM analysis. Note that the clusters shown in the plots represent the most probable clusters of SNPs, not the true clusters.
#'
#' * `plot_wald`: Plots the Wald ratio estimates against IV strength.
#' * `plot_scatter`: Plots the IV-outcome associations against IV-exposure associations.
#'
#' @import stats grDevices graphics
#'
#' @param res A list obtained from `mr.gmm`.
#'
#' @return R plots.
#' @export
#'
#'
plot_wald <- function(res) {
  x <- res$data$beta_exp^2/res$data$se_exp^2-1
  y <- res$data$beta_out/res$data$beta_exp
  xlim <- c(min(x), max(x))
  ylim <- quantile(y, c(0.025, 0.975), na.rm = T)
  cluster1 <- res$cluster$cluster == 1
  cluster2 <- res$cluster$cluster == 2
  cluster3 <- res$cluster$cluster == 3
  cluster4 <- res$cluster$cluster == 4
  par(mar = c(4.5, 4.5, 1, 1))
  plot(x[cluster4], y[cluster4],
       xlab = "IV strength", ylab = "Wald ratio estimate",
       col = rgb(0, res$cluster$prob[cluster4], 0), pch = 19,
       xlim = xlim,
       ylim = ylim,
       cex = 1.5,
       cex.lab = 1.5,
       cex.axis = 1.5)
  abline(h = res$beta.hat, lwd = 3, lty = 2, col = "blue")
  points(x[cluster2], y[cluster2],cex = 1.5, col = rgb(0, 0, res$cluster$prob[cluster2]), pch = 19)
  points(x[cluster1], y[cluster1],cex = 1.5, col = rgb(res$cluster$prob[cluster1], 0, 0), pch = 19)
  points(x[cluster3], y[cluster3],cex = 1.5, col = rgb(res$cluster$prob[cluster3], 165/255*res$cluster$prob[cluster3], 0), pch = 19)
  legend("topright",
         legend = c("Invalid", "Valid", "Invalid&Null", "Null"),
         col = c(rgb(1, 0, 0),
                 rgb(0, 0, 1),
                 rgb(1, 165/255, 0),
                 rgb(0, 1, 0)),
         pch = 19,
         bty = "n",
         text.width = rep((xlim[2]-xlim[1])/4.3, 4),
         cex = 1.5,
         pt.cex = 1.5)
}

#' @rdname plot_wald
#' @export

plot_scatter <- function(res) {
  x <- res$data$beta_exp
  y <- res$data$beta_out
  xlim <- c(min(x), max(x))
  ylim <- c(min(y), max(y))
  cluster1 <- res$cluster$cluster == 1
  cluster2 <- res$cluster$cluster == 2
  cluster3 <- res$cluster$cluster == 3
  cluster4 <- res$cluster$cluster == 4
  par(mar = c(4.5, 4.5, 1, 1))
  plot(x[cluster4], y[cluster4],
       xlab = "IV-exposure association", ylab = "IV-outcome association",
       col = rgb(0, res$cluster$prob[cluster4], 0), pch = 19,
       xlim = xlim,
       ylim = ylim,
       cex = 1.5,
       cex.lab = 1.5,
       cex.axis = 1.5)
  abline(a = 0, b = res$beta.hat, col = "blue", lwd = 3, lty = 2)
  points(x[cluster2], y[cluster2], col = rgb(0, 0, res$cluster$prob[cluster2]), pch = 19, cex = 1.5)
  points(x[cluster1], y[cluster1], col = rgb(res$cluster$prob[cluster1], 0, 0), pch = 19, cex = 1.5)
  points(x[cluster3], y[cluster3], col = rgb(res$cluster$prob[cluster3], 165/255*res$cluster$prob[cluster3], 0), pch = 19, cex = 1.5)
  legend("topright",
         legend = c("Invalid", "Valid", "Invalid&Null", "Null"),
         col = c(rgb(1, 0, 0),
                 rgb(0, 0, 1),
                 rgb(1, 165/255, 0),
                 rgb(0, 1, 0)),
         pch = 19,
         bty = "n",
         text.width = rep((xlim[2]-xlim[1])/4.3, 4),
         cex = 1.5,
         pt.cex = 1.5)
}
