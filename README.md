proteoQ
================
true
2026-09-30

- [Introduction to proteoQ](#introduction-to-proteoq)
- [Installation](#installation)
- [](#section)

The family currently contains members of

- [Mzion](https://github.com/qzhang503/mzion) for database searches.
- [proteoQ](https://github.com/qzhang503/proteoQ) for quality metrics
  and informatics.

The document is mainly about proteoQ.

## Introduction to proteoQ

Tandem mass tag ([TMT](https://en.wikipedia.org/wiki/Tandem_mass_tag))
and label-free quantitation (LFQ) of proteomics data. It can utilize
results from various search engines, including:

- TMT workflows: [Mascot](https://http://www.matrixscience.com/),
  [MaxQuant](https://www.maxquant.org/),
  [MSFragger](http://msfragger.nesvilab.org/) and
  [Mzion](https://github.com/qzhang503/mzion)

- LFQ workflows: DIA-NN, Mzion

## Installation

To install this package, start R (latest version) and enter:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("qzhang503/proteoQ")
```

## 
