# mr.gmm: Mendelian randomization using a four-cluster Gaussian mixture model

## Installation

`mr.gmm` performs two-sample Mendelian randomization analysis using a four-cluster Gaussian mixture model and generate visualizations of the results. The deterministic annealing expectation-maximization (DAEM) algorithm is applied to fit the mixture model. To install the package, run:

```
if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("YuquanW/mr.gmm")
```

## Usage

We recommend preprocessing your data using `TwoSampleMR::harmonise_data()`. You can test the example datasets included in the package by running the following:

```
## Use p<0.05 as IV selection threshold
data(abo.chd)
abo.chd.selected <- subset(abo.chd, pval.exposure > -log10(0.05))
attach(abo.chd.selected)
res <- mr.gmm(beta.exposure,
              beta.outcome,
              se.exposure,
              se.outcome,
              sel.p = 0.05)
```

To visualize the results, simply run:

```
plot_wald(res)
plot_scatter(res)
```

To reproduce the simulation results or perform pseudo-$p$-value-based LD clumping, please refer to the [Clumping](https://github.com/YuquanW/Clumping) repository, specifically the `Batch_simu.py` and `LDclumping.sh` scripts.

## Reference

Yuquan Wang, Yunlong Cao, Dong Chen, Dapeng Shi, Liwan Fu, Anqi Chen, Siyuan Shen, Yue-Qing Hu. Pseudo-p-Value-Based Clumping Enhanced Proteome-wide Mendelian Randomization with Application in Identifying Coronary Heart Disease-Associated Plasma Proteins. doi:
[https://doi.org/10.1101/2025.01.13.25320450](https://doi.org/10.1101/2025.01.13.25320450)

## Contact

If you encounter any issues, please report them on the [Issues](https://github.com/YuquanW/mr.gmm/issues) page or contact Dr. Yuquan Wang (21110700112\@m.fudan.edu.cn) and Prof. Yue-Qing Hu (yqhu\@fudan.edu.cn).
