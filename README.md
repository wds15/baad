# Source for "Bayesian aggregation of average data: An application in drug development"

Date: 13th May 2020

This repository includes updated R and Stan source files used to
generate the main outputs in our paper on Bayesian aggregation of
average data. The updated codes ensure that the source runs with the
RStan 2.19 release.

The code is distributed under the BSD-3 license, see LICENSE.

Please cite our work when you use these codes or apply our approach:

_S. Weber, A. Gelman, D. Lee, M. Betancourt, A. Vehtari, A. Racine
(2018), "Bayesian aggregation of average data: An application in drug
development", Ann. Appl. Stat., Volume 12, Number 3 (2018), 1583-1604._

URL: https://projecteuclid.org/euclid.aoas/1536652966  
doi: 10.1214/17-AOAS1122  
arXiv: TBD

These sources are also available in their original form as
supplemental material to the original article under the DOI

doi: 10.1214/17-AOAS1122SUPP

To use the code, we recommend to start with one of the examples below
as template and modify it to your needs.

The repository includes a Docker image defintion which is configured
to use R 4.0.0 and R packages as dated from the 20th April of 2020. To
run the scripts with this environment please first build the
respective Docker image with the `build-rimage.sh` script. Once the
Docker image is created you can run the analysis in the source
directory via

- Analytical linear example: `./run-docker-cmd.sh ./linear_baad.R`
- Non-linear example: `./run-docker-cmd.sh ./pkpd_baad.R`

Outputs from these scripts will be saved in the current source
directory directly. Note that the scripts assume a bash shell as
commonly used on Linux or Mac.

Updates will be posted to https://github.com/wds15/baad

# Contents

- `build-rimage.sh` Build Docker image
- `run-docker-cmd.sh` Run analysis with Docker image as container
- `linear_baad.R` Linear example
- `pkpd_badd.R` Non-linear example
- `utils_baad.R` Various utility functions

# Software used

R session info:

```
R version 4.0.0 (2020-04-24)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux bullseye/sid

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] plyr_1.8.6         mvtnorm_1.1-0      abind_1.4-5        rstan_2.19.3
[5] ggplot2_3.3.0      StanHeaders_2.19.2

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.4.6        RColorBrewer_1.1-2  pillar_1.4.3
 [4] compiler_4.0.0      prettyunits_1.1.1   tools_4.0.0
 [7] digest_0.6.25       pkgbuild_1.0.6      lifecycle_0.2.0
[10] tibble_3.0.0        gtable_0.3.0        lattice_0.20-41
[13] pkgconfig_2.0.3     rlang_0.4.5         Matrix_1.2-18
[16] cli_2.0.2           parallel_4.0.0      loo_2.2.0
[19] gridExtra_2.3       withr_2.1.2         vctrs_0.2.4
[22] stats4_4.0.0        grid_4.0.0          glue_1.4.0
[25] inline_0.3.15       R6_2.4.1            BH_1.72.0-3
[28] processx_3.4.2      fansi_0.4.1         farver_2.0.3
[31] callr_3.4.3         magrittr_1.5        codetools_0.2-16
[34] scales_1.1.0        ps_1.3.2            ellipsis_0.3.0
[37] matrixStats_0.56.0  assertthat_0.2.1    colorspace_1.4-1
[40] labeling_0.3        munsell_0.5.0       crayon_1.3.4
[43] RcppEigen_0.3.3.7.0
```
