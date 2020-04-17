# ShiftConvolve
<!-- badges: start -->
[![Travis build status](https://travis-ci.org/andrew12678/ShiftConvolve.svg?branch=master)](https://travis-ci.org/andrew12678/ShiftConvolve)
<!-- badges: end -->
`ShiftConvolve` is a R package which uses exponential shifting and the Fast Fourier Transformations (FFT) to compute the (right) tail of distribution of the Poisson Binomial Distribution. 
This package makes use of the minimalist Fast Fourier Transform library known as [minFFT](https://github.com/aimukhin/minfft) to perform the necessary DFT and Inverse DFT computations.

For the `ShiftConvolve` implementation which uses [FFTW3](http://www.fftw.org/) to perform the Fourier Transformations please go to [https://github.com/andrew12678/ShiftConvolveFFTW](https://github.com/andrew12678/ShiftConvolveFFTW).

In the interest of speed we all significant computational aspects of our procedure are executed by code written in C and the R code mainly acts as a wrapper around that compiled C code.

We have successfully installed `ShiftConvolveFFTW` on `Windows 10: Professional Version 1909`, `macOS: Mojave 10.14.6/Catalina 10.15.4` and `Linux: Ubuntu 18.04/Manjaro 19.02 KDE Plasma`. A full installation guide will be provided below for all 3 Operating Systems.

`ShiftConvolveFFTW` has also been installed on both `AMD` and `Intel` CPU systems. Specifically, a  `Intel Core i7 MacMini with 32GB of RAM`, `AMD Ryzen 5 2600 with 16GB of RAM` and a `Macbook Pro i7 with 16GB of RAM` 

## Dependencies

There are no external dependencies required to be installed for this version of `ShiftConvolve`. The required `minFFT` files have been compiled and included into the `ShiftConvolve` project.    

`Windows` users may be required (instructed by compilation errors later anyways) to install `Rtools` from the following website [https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/). Make sure to tick the box saying 'Add to system path variables' when installing. 

## Installation

The most simple installation involves simply cloning this repository and installing `ShiftConvolve` from source. 

```R
install.packages('devtools') # If you have not installed devtools
devtools::install_github('https://github.com/andrew12678/ShiftConvolve.git')
```

An alternative installation procedure involves cloning the repository, creating a `RStudio` project and then building. 

## Examples

An simple example with the uniform distribution

```R
library(ShiftConvolvePoibin)
set.seed(18)
n=1000
probs <- runif(n)
x <- c(200, 500, 800)
p <- seq(0, 1, 0.01)
dpoisbin(x,probs,method="ShiftConvolve",log.p=FALSE)
ppoisbin(x,probs,method="ShiftConvolve",lower.tail=FALSE,log.p=TRUE)
qpoisbin(p,probs,method="ShiftConvolve",lower.tail=TRUE,log.p=FALSE)
rpoisbin(n,probs)
```

## References
<a id="1">[1]</a> 
Peres, N., Lee, A., and Keich, U. (2020). Exactly computing the tail of the Poisson-Binomial Distribution. 
[arXiv:2004.07429](https://arxiv.org/abs/2004.07429)

