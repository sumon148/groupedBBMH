# groupedBBMH

**Grouped Beta-Binomial Model with Metropolis-Hastings**

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

`groupedBBMH` is an R package for analyzing group-testing biosecurity surveillance data using a nested beta-binomial model. 
It implements maximum likelihood estimation (MLE) and a Metropolis-Hastings (MH) algorithm to handle censored data from group tests, where only batch-level or group-level outcomes are observed rather than individual item contamination.



> **Das, S., Clark, R., Parsa, M., &  Barnes, B. (2025). Inferring Leakage in Imports of Frozen Seafood allowing for Censoring, Testing Accuracy and a Minimum Positive Prevalence.**
> *Journal of the American Statistical Association*, 

The methodology is described in detail in the accompanying paper and supplementary material.

> 📘 **Repository:**
> [https://github.com/sumon148/groupedBBMH](https://github.com/sumon148/groupedBBMH)

## Overview

The **groupedBBMH** R package provides statistical tools for modeling 
**group-testing data in biosecurity surveillance** using a 
**nested beta-binomial framework**. 
It explicitly accounts for **imperfect diagnostic accuracy** and the 
**biologically meaningful minimum prevalence threshold**.

## Key Highlights

- **Hierarchical Modeling under Imperfect Testing**  
  Implements a nested beta-binomial model that captures both within- and between-group variability, 
  accounting for imperfect test **sensitivity** and **specificity**.  
  Supports exact inference via a *Metropolis–Hastings (MH)* algorithm and 
  approximate estimation via *Maximum Likelihood Estimation (MLE)*.

- **Flexible Integration with Standard R Packages**  
  Supports Bayesian modeling and diagnostic visualization using familiar tools such as  
  **`coda`** for MCMC diagnostics and **`bayesplot`** for posterior summaries and convergence assessment.
  
- **Risk and Decision Support for Biosecurity Surveillance**  
  Provides methods to estimate contamination prevalence, quantify 
  **leakage risk** (undetected contamination), and evaluate sampling strategies.  
  Users can also specify a **minimum biologically plausible prevalence threshold** `k`, 
  ensuring realistic contamination estimates and robust, risk-informed decision-making.


## Installation

The package is **not** available on CRAN. You can install it directly from GitHub.

### Install without vignettes

This is faster if you do not need the tutorial vignette:

```r
# install.packages("devtools") if needed
devtools::install_github("sumon148/groupedBBMH",dependencies=TRUE)
```

### Install with vignette (recommended)

To install **and build the vignette**, use:

```r
devtools::install_github("sumon148/groupedBBMH", build_vignettes = TRUE)
```

> **Note:** Building vignettes requires `rmarkdown` and `knitr`:
>
> ```r
> install.packages(c("rmarkdown", "knitr"))
> ```

---

## Getting Started

After installation, load the package:

```r
library(groupedBBMH)
```

To see the vignette tutorial:

1. **List available vignettes:**

```r
vignette(package = "groupedBBMH")
```

2. **Open the vignette:**

```r
vignette("groupedBBMH")
```

### The vignette includes:

- **Mathematical background** of the grouped beta-binomial (BB) model, including extensions to handle test measurement errors.  
- **Generalized and truncated grouped BB models**, incorporating a biologically meaningful minimum positive prevalence threshold.  
- **Implementation of Maximum Likelihood Estimation (MLE)** and the **Metropolis–Hastings (MH)** algorithm in a practical case study.  
- **Leakage risk estimation** under different testing and prevalence scenarios.  
- **Model-based simulation** to evaluate biosecurity or community-level risk under varying sampling and diagnostic test conditions.


---

## Example (Quick Preview)

Below is a minimal illustration. For detailed use cases, refer to the vignette.

```r
# Load library
library(groupedBBMH)

# Develop Beta-binomial model using MLE method
ty=c(0:13)
freq=c(2815, 9, 10, 6, 1, 3, 2, 0, 1, 2, 1, 0, 0, 0)
trbb <- fit_trGroupedBB(ty, freq, b = 13, m = 5, cutoff = 0.01,
                       sensitivity = 0.8, specificity = 1)
trbb
```


## Methods

The package implements the analytical results described in:

* Das, S., Clark, R., Parsa, M., & Barnes, B. (2025). Inferring Leakage in Imports of Frozen Seafood allowing for Censoring, Testing Accuracy and a Minimum Positive Prevalence. *Journal of the American Statistical Association* (under review).



## Citation

If you use this package in publications, please cite:

* Das, S., Clark, R., Parsa, M., & Barnes, B. (2025). *groupedBBMH*: R package implementing Beta-binomial model to group-testing biosecurity data with imperfect tests. Version 0.1.0. https://github.com/sumon148/groupedBBMH


## License

This package is released under the [Creative Commons CC0 1.0 Universal Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/).

[![License: CC0-1.0](https://img.shields.io/badge/License-CC0%201.0-lightgrey.svg)](https://creativecommons.org/publicdomain/zero/1.0/)

