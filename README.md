# ShiftConvolve
`ShiftConvolve` is a R package which uses exponential shifting and Fast Fourier Transformations to compute the distribution of the Poisson Binomial Distribution. 
This package makes use of the minimalist Fast Fourier Transform library known as [minFFT](https://github.com/aimukhin/minfft) to perform the necessary DFT and Inverse DFT computations in the ShiftConvolve approach.

For the `ShiftConvolve` implementation which uses [FFTW3](http://www.fftw.org/) to perform the Fourier Transformations please go to [https://github.com/andrew12678/ShiftConvolveFFTW](https://github.com/andrew12678/ShiftConvolveFFTW).

All heavy computational features and convolution functions are written in C, compiled and called by R.

## Dependencies

There are no external dependencies required to be installed for this version of `ShiftConvolve`. The required `minFFT` files have been compiled and included into the `ShiftConvolve` project.    

## Installation

The most simple installation involves simply cloning this repository and installing `ShiftConvolve` from source. 

```console
git clone https://github.com/andrew12678/ShiftConvolve.git
cd ShiftConvolve
install.packages('ShiftConvolvePoibin_2.2.0.tar.gz', repos = NULL, type="source")
```

An alternative installation procedure involves cloning the repository, creating a `RStudio` project in the `ShiftConvolvePoibin` folder and then building. 

## Examples

An simple example with the uniform distribution

```R
library(ShiftConvolvePoibin)
set.seed(18)
a = runif(10000)
shiftpval(a,50000)
```