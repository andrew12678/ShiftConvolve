#' Computing the Poisson Binomial Distribution using Shifted FFT
#'
#' A package which uses exponential shifting and Fast Fourier Transformations to
#' compute the distribution of the Poisson Binomial Distribution
#'
#' @rdname ShiftedConvolvePB
#' @author Andrew Ray Lee, Noah Peres and Uri Keich
#' @title ShiftConvolve
#'
process = function(p){
  n0 = sum(p == 0)
  n1 = sum(p == 1)
  if((n0 != 0) | (n1 != 0)){
    p = p[(p != 0) & (p != 1)]
  }
  return(c(n0,n1,p))
}
#'
#' fftconv
#' Performs the 'streamlined' approach of convoluting a given probability vector in the Fourier domain,
#' this function prepares neccessary vectors, including the final result vector, which will be passed
#' to the a C program to attain the probability mass function.
#' @param p Accepts a list of probabilities
#' @return The final probability mass function from the convolution
fftconv = function(p){
  #Calculating the sizes of the vectors so they can be declared and passed to C
  input = as.vector(p)
  n = length(input)/2
  r <- 4
  c <-n
  while(c != 2){
    c <- ceiling(c/2)
    r <- r*2
  }
  s = r*2

  #Declaring vectors to be passed to C
  result = vector("numeric", s/2) #The vector for the final probability mass function
  A = vector("complex", s) #The vector for the probabilities transformed to the Fourier domain
  B = vector("complex", s) #The vector to store pointwise multiplications
  U = vector("complex", s) #The vector to store retrieve pmf before further pointwise multiplication
  in_vec = vector("numeric", s)
  w_vec = vector("numeric", s)
  out_vec = vector("numeric", s)

  pmf = .C("fftconvPairs",
           as.numeric(input),
           as.integer(n),
           as.complex(A),
           as.complex(B),
           as.complex(U),
           as.numeric(result),
           as.numeric(in_vec),
           as.numeric(w_vec),
           as.numeric(out_vec))[[6]]
  return(pmf)
}
#'
#'shiftMean
#'
#'The shiftMean function evaluates the mean of a PB distribution shifted by t (minus S0)
shiftedmean <- function(t, p, S0 = 0){
  v <- exp(t)*(p)/((1-p)+exp(t)*(p))
  return(sum(v)-S0)
}
#'
#' shiftConvolveLog
#' @description
#' Convolutes a vector of probabilities with exponential shift and returns the log of the distribution
#' @param p Accepts a list of probabilities
#' @param S0 An observed value, default set to -1, which will set the shifting paramater t0 to 0
#' @param process Whether we want to process for 0s and 1s, always set to TRUE unless called from shiftpval
#' @examples
#' set.seed(18)
#' a = runif(100)
#' shiftConvolveLog(a, 50)
#' @export
shiftConvolveLog<-function(p, S0 = -1, process = TRUE){
  if(process){
    p = process(p)
    n0 = p[1]
    n1 = p[2]
    p = p[3:length(p)]
  }

  if(S0 == -1){
    t0 <-0
  }else{
    U <- uniroot(shiftedmean, interval = c(-50,50), p, S0, tol = 10^-16) #finds a suitable shifting parameter
    t0 <- U$root
  }


  lp <- log(c(1-p,p))
  n = length(p)
  M <- vector("double",length = n)
  S <- vector("double", length = n*2)

  xmax <-vector("double", length = 1)
  ret = .C("computeMGF",
           as.numeric(lp),
           as.integer(n),
           as.numeric(t0),
           as.numeric(M),
           as.numeric(S),
           as.numeric(xmax))
  M<-ret[[4]]
  S <- ret[[5]]
  PBshift <- fftconv(exp(S)) #convolve using FFT
  lPBshift <- log(PBshift)
  lMPB <- -sum(M) #(log)MGF of shifted PB at -t0
  l<- length(lPBshift)
  k <- c(0:(l-1))
  PB <- lPBshift + -t0*k + sum(M)
  return(PB)
}
#'
#' logsum
#' This function calculates the a sum of a vector (in log space, returns a log)
#' @param x A vector
logsum <- function(x) {
  xmax <- max(x)
  v <- x-xmax
  v <- exp(v)
  ans <- xmax + log(sum(v))
  return(ans)
}
#'
#'
dcConvolvePaired <- function(p){
  n = length(p);
  result <- vector("double", length = n+1)
  ret = .C("fullconvolvePaired",
           as.numeric(p),
           as.integer(n),
           as.numeric(result))
  #pmf = ret[[1]]
  return(ret[[3]])
}
#'
dcConvolvePairedLog <- function(p){
  n = length(p);
  result <- vector("double", length = n+1)
  ret = .C("fullconvolvePairedLog",
           as.numeric(p),
           as.integer(n),
           as.numeric(result))
  #pmf = ret[[1]]
  return(ret[[3]])
}
#' shiftpval
#'
#' @description
#' The function takes in a vector of (Bernoulli) probabilities, an observed value s0, and calculates (log(p-value), p-value)
#' of the PB distribution using shiftConvolveLog
#'
#' @param pb A vector of Bernoulli probabilities
#' @param S0 An observed value
#' @param right.tail Logical determining whether to sum the right tail (default) or left tail
#'
#' @return A vector containing the log of the p-value and the p-value
#'
#' @examples
#' set.seed(18)
#' a = runif(100)
#' shiftpval(a, 50)
#' shiftpval(a, 50, right.tail=FALSE)
#'
#' @export
shiftpval <- function(pb,S0,right.tail=TRUE){
  pb = process(pb)
  n0 = pb[1]
  n1 = pb[2]
  pb = pb[3:length(pb)]
  n <- length(pb)+1

  S1 = S0 - n1
  if (S1 <= 0){
    if (right.tail){
      return(c(0,1))
    } else {
      return(c(-Inf,0))
    }
  }

  if (S1 >= n){
    if (right.tail){
      return(c(-Inf,0))
    } else {
      return(c(0,1))
    }
  }

  if (S1 == n-1){
    ldist <- shiftConvolveLog(pb, (S1-1), process = FALSE)
    U <- uniroot(shiftedmean, interval = c(-50,50), pb, S1-1, tol = 10^-16) #find a suitable shifting parameter
    t0 <- U$root
  }else{
    ldist <- shiftConvolveLog(pb, S1, process = FALSE)
    U <- uniroot(shiftedmean, interval = c(-50,50), pb, S1, tol = 10^-16) #find a suitable shifting parameter
    t0 <- U$root
  }

  ldist = c(rep(-Inf, n1), ldist)
  ldist = c(ldist, rep(-Inf,n0))

  if (right.tail){
    if (t0 > 0){
      lp <- logsum(ldist[c((S0+1):length(ldist))])
      p <- exp(lp)
    }
    else{
      llefttail <- logsum(ldist[c(0:S0)])
      p <- 1-exp(llefttail)
      lp <- log(p)
    }
  } else {
    if (t0 > 0){
      lrighttail <- logsum(ldist[c((S0+2):length(ldist))])
      p <- 1-exp(lrighttail)
      lp <- log(p)
    }
    else{
      lp <- logsum(ldist[c(0:(S0+1))])
      p <- exp(lp)
    }
  }
  return(c(lp,p))
}
#' dcpval
#'
#' @description
#' The function takes in a vector of (Bernoulli) probabilities, an observed value s0, and calculates (log(p-value), p-value)
#' of the PB distribution using a direct convolution
#'
#' @param pb A vector of Bernoulli probabilities
#' @param S0 An observed value
#' @param right.tail Logical determining whether to sum the right tail (default) or left tail
#' @param log.p Logical determining whether to return the log of the p-value or the p-value
#' @return A vector containing the log of the p-value and the p-value
#'
#' @examples
#' set.seed(18)
#' a = runif(100)
#' dcpval(a, 50)
#' dcpval(a, 50, right.tail=FALSE)
#' dcpval(a, 50, right.tail=TRUE, log.p=TRUE)
#' @export
dcpval <- function(pb,S0,right.tail=TRUE,log.p=FALSE){
  pb = process(pb)
  n0 = pb[1]
  n1 = pb[2]
  pb = pb[3:length(pb)]
  n <- length(pb)+1

  S1 = S0 - n1
  if (S1 <= 0){
    if (right.tail){
      return(c(0,1))
    } else {
      return(c(-Inf,0))
    }
  }

  if (S1 >= n){
    if (right.tail){
      return(c(-Inf,0))
    } else {
      return(c(0,1))
    }
  }

  if(log.p){
    ldist <- dcConvolvePairedLog(pb)
    ldist = c(rep(-Inf, n1), ldist)
    ldist = c(ldist, rep(-Inf,n0))
    if(right.tail){
      lp <- logsum(ldist[c((S0+1):length(ldist))])
    }else{
      lp <- logsum(ldist[c(0:(S0+1))])
    }
    return(lp)
  }else{
    dist<- dcConvolvePaired(pb)
    dist = c(rep(0, n1), dist)
    dist = c(dist, rep(0,n0))
    if(right.tail){
      p <- sum(dist[c((S0+1):length(dist))])
    }else{
      p <- sum(ldist[c(0:(S0+1))])
    }
    return(p)

  }
}
#dyn.load("../src/fftconv_min.so")
#dyn.load("../src/ShiftConvolvePoibin.so")
#dyn.load("../src/fftconv.so")
#set.seed(18)
#a = runif(5)
#shiftpval(a, 4)
