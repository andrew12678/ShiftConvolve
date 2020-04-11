# ShiftConvolve
`ShiftConvolve` is a R package which uses exponential shifting and the Fast Fourier Transformations (FFT) to compute the (right) tail of distribution of the Poisson Binomial Distribution. 
This package makes use of the minimalist Fast Fourier Transform library known as [minFFT](https://github.com/aimukhin/minfft) to perform the necessary DFT and Inverse DFT computations.

For the `ShiftConvolve` implementation which uses [FFTW3](http://www.fftw.org/) to perform the Fourier Transformations please go to [https://github.com/andrew12678/ShiftConvolveFFTW](https://github.com/andrew12678/ShiftConvolveFFTW).

In the interest of speed we all significant computational aspects of our procedure are executed by code written in C and the R code mainly acts as a wrapper around that compiled C code.

## Dependencies

There are no external dependencies required to be installed for this version of `ShiftConvolve`. The required `minFFT` files have been compiled and included into the `ShiftConvolve` project.    

## Installation

The most simple installation involves simply cloning this repository and installing `ShiftConvolve` from source. 

```console
git clone https://github.com/andrew12678/ShiftConvolve.git
cd ShiftConvolve
install.packages('ShiftConvolvePoibin_2.5.0.tar.gz', repos = NULL, type="source")
```

An alternative installation procedure involves cloning the repository, creating a `RStudio` project in the `ShiftConvolvePoibin` folder and then building. 

## Examples

An simple example with the uniform distribution

```R
library(ShiftConvolvePoibin)
set.seed(18)
n = 10000
p = runif(n)
s0 = 5200
shiftpval(p, s0)	# compute the p-value, or right tail at s0
shiftpval(1-p, n-s0)	# compute the p-value, or left tail at s0
```
