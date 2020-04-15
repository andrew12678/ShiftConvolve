#' Computing the Poisson Binomial Distribution using ShiftConvolve
#'
#' A package which uses exponential shifting and Fast Fourier Transformations with the minFFT library to
#' compute the distribution of the Poisson Binomial Distribution
#'
#'@rdname ShiftedConvolvePB
#'@useDynLib ShiftConvolvePoibin
#'@title ShiftConvolve Poisson Binomial
#'
#'@name ShiftedConvolvePB
#'@author Andrew Ray Lee, Noah Peres and Uri Keich
#'@description
#'Density, distribution function, quantile function and random generation for
#'the Poisson binomial distribution with the option of using the ShiftConvolve method.
#'
#'
#'@param x           Either a vector of observed numbers of successes
#'                   (or vector of quantiles as dbinom/pbinom refers to) or NULL.
#'                   If NULL, probabilities of all possible observations are returned.
#'@param p           Vector of probabilities for computation of quantiles.
#'@param n           Number of observations. If \code{length(n) > 1}, the
#'                   length is taken to be the number required.
#'@param probs       Vector of probabilities of success of each Bernoulli
#'                   trial.
#'@param method      Character string that specifies the method of computation
#'                   and must be either \code{"ShiftConvolve"} or \code{"DC"}
#'@param log.p       logical; if TRUE, probabilities p are given as log(p).
#'@param lower.tail  Logical value indicating if results are \eqn{P[X \le x]}
#'                   (if \code{TRUE}; default) or \eqn{P[X > x]} (if
#'                   \code{FALSE}).
#'
#'@examples
#'set.seed(18)
#'n=1000
#'probs <- runif(n)
#'x <- c(200, 500, 800)
#'p <- seq(0, 1, 0.01)
#'dpoisbin(x,probs,method="ShiftConvolve",log.p=FALSE)
#'ppoisbin(x,probs,method="ShiftConvolve",lower.tail=FALSE,log.p=TRUE)
#'qpoisbin(p,probs,method="ShiftConvolve",lower.tail=TRUE,log.p=FALSE)
#'rpoisbin(n,probs,method="ShiftConvolve")
#'@export
dpoisbin <- function(x, probs, method = "ShiftConvolve", log.p = FALSE){
  #Check if input x is NULL,
  #If so we return all possible probabilities
  #Otherwise we ensure X is all integers
  if(is.null(x)){
    n<-length(probs)
    x <- 0:n
  }else{
    x <- as.vector(x)
    if(!all(x == floor(x))){
      stop("'x' must be either NULL or all integers")
    }
  }

  #Check if we are given valid probabilities
  if(is.null(probs) || any(is.na(probs) | probs < 0 | probs > 1)){
    stop("'probs' must contain valid probabilities")
  }

  #Ensure we are using a valid method
  if(!(method %in% c("ShiftConvolve", "DC"))){
    stop("'method' provided can only be 'ShiftConvolve' or 'DC'")
  }

  #Ensure we receive a legimate boolean for log.p
  if(!(log.p %in% c(TRUE, FALSE))){
    stop("'log.p' provided can only be TRUE or FALSE")
  }

  if(method == "ShiftConvolve"){
    return(shift_pb_density(x = x, probs = probs, log.p = log.p))
  }else if(method == "DC"){
    return(dc_pb_density(x = x, probs = probs, log.p = log.p))
  }


}

#'@rdname ShiftedConvolvePB
#'@export
ppoisbin <- function(x, probs, method = "ShiftConvolve", lower.tail = TRUE, log.p = FALSE){
  #Check if input x is NULL,
  #If so we return all possible probabilities
  #Otherwise we ensure X is all integers
  if(is.null(x)){
    n<-length(probs)
    x <- 0:n
  }else{
    x <- as.vector(x)
    if(!all(x == floor(x))){
      stop("'x' must be either NULL or all integers")
    }
  }

  #Check if we are given valid probabilities
  if(is.null(probs) || any(is.na(probs) | probs < 0 | probs > 1)){
    stop("'probs' must contain valid probabilities")
  }

  #Ensure we are using a valid method
  if(!(method %in% c("ShiftConvolve", "DC"))){
    stop("'method' provided can only be 'ShiftConvolve' or 'DC'")
  }

  #Ensure we receive a legimate boolean for lower.tail
  if(!(lower.tail %in% c(TRUE, FALSE))){
    stop("'lower.tail' provided can only be TRUE or FALSE")
  }

  #Ensure we receive a legimate boolean for log.p
  if(!(log.p %in% c(TRUE, FALSE))){
    stop("'log.p' provided can only be TRUE or FALSE")
  }
  if(method == "ShiftConvolve"){
    return(shift_pb_dist(x = x, probs = probs, lower.tail = lower.tail, log.p = log.p))
  }else if(method == "DC"){
    return(dc_pb_dist(x = x, probs = probs, lower.tail = lower.tail, log.p = log.p))
  }

}

#'@rdname ShiftedConvolvePB
#'@export
qpoisbin <- function(p, probs, method = "ShiftConvolve", lower.tail = TRUE, log.p = FALSE){
  #Ensure we receive a legimate boolean for log.p
  if(!(log.p %in% c(TRUE, FALSE))){
    stop("'log.p' provided can only be TRUE or FALSE")
  }

  # Check that 'p' contains only probabilities
  if(!log.p){
    if(is.null(p) || any(is.na(p) | p < 0 | p > 1))
      stop("'p' must contain real numbers between 0 and 1 if log.p is FALSE")
  }else{
    if(is.null(p) || any(is.na(p) | p > 0))
      stop("'p' must contain real numbers between -Inf and 0 if log.p is TRUE")
  }

  ## Setting x=NULL gives us the CDF (also checkse other variables), log.p also set to FALSE
  cdf <- ppoisbin(NULL, probs = probs, method = method, lower.tail = lower.tail, log.p = FALSE)

    #Matching with CDF, log.p is TRUE
  if(log.p){
    p <- exp(p)
  }

  return_quantiles = vector(mode = "integer", length(p))

  curr_cdf_pos <- 1
  cdf_length <- length(cdf)
  if(lower.tail){
    for(prob in sort(unique(p))){
      indexes <- which(p == prob)
      while(curr_cdf_pos < cdf_length && cdf[curr_cdf_pos] <= prob){
        curr_cdf_pos <- curr_cdf_pos + 1
      }
      return_quantiles[indexes] <- curr_cdf_pos - 1
    }
  }else{
    for(prob in sort(unique(p), decreasing = TRUE)){
      indexes <- which(p == prob)
      while(curr_cdf_pos < cdf_length && cdf[curr_cdf_pos] >= prob){
        curr_cdf_pos <- curr_cdf_pos + 1
      }
      return_quantiles[indexes] <- curr_cdf_pos - 1
    }
  }

  return(return_quantiles)

}

#'@rdname ShiftedConvolvePB
#'@export
rpoisbin <- function(n, probs){
  # check if 'n' is NULL
  if(is.null(n)){
    stop("'n' must not be NULL!")
  }

  #Check if we are given valid probabilities
  if(is.null(probs) || any(is.na(probs) | probs < 0 | probs > 1)){
    stop("'probs' must contain valid probabilities")
  }

  len <- length(n)
  if(len > 1){
    n <- len
  }

  return_probs <- rep(0, n)

  for(i in probs){
    return_probs <- return_probs + rbinom(n,1,i)
  }

  return(return_probs)

}

fftconv = function(p){
  input = as.vector(p)   #Calculating the sizes of the vectors so they can be declared and passed to C
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

process = function(p){
  n0 = sum(p == 0)
  n1 = sum(p == 1)
  if((n0 != 0) | (n1 != 0)){
    p = p[(p != 0) & (p != 1)]
  }
  return(c(n0,n1,p))
}

shiftedmean <- function(t, p, S0 = 0){
  v <- exp(t)*(p)/((1-p)+exp(t)*(p))
  return(sum(v)-S0)
}

shiftConvolveLog<-function(p, S0 = -1){
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
  l<- length(lPBshift)
  k <- c(0:(l-1))
  PB <- lPBshift + -t0*k + sum(M)
  return(PB)
}

logsum <- function(x) {
  xmax <- max(x)
  v <- x-xmax
  v <- exp(v)
  ans <- xmax + log(sum(v))
  return(ans)
}

dcConvolvePaired <- function(p){
  n = length(p);
  result <- vector("double", length = n+1)
  ret = .C("fullconvolvePaired",
           as.numeric(p),
           as.integer(n),
           as.numeric(result))
  return(ret[[3]])
}

dcConvolvePairedLog <- function(p){
  n = length(p);
  result <- vector("double", length = n+1)
  ret = .C("fullconvolvePairedLog",
           as.numeric(p),
           as.integer(n),
           as.numeric(result))
  return(ret[[3]])
}

shift_pb_dist <- function(x, probs, lower.tail = TRUE, log.p=FALSE){
  probs = process(probs)
  n0 = probs[1]
  n1 = probs[2]
  probs = probs[3:length(probs)]
  n <- length(probs)+1

  x_ordinary = c()
  x_return_index = c()
  pmf_computed = NULL
  return_dist = vector(mode = "numeric", length = length(x))
  for(i in 1:length(x)){
    edge <- edge_prob(x[i], n,  n1, lower.tail, log.p)
    if(!is.null(edge)){
      return_dist[i] <- edge
    }else{
      x_ordinary = c(x_ordinary, x[i])
      x_return_index = c(x_return_index, i)
    }
  }

  n_x_ordinary = length(x_ordinary)
  if(n_x_ordinary == 0){
    return(return_dist)
  }else if(n_x_ordinary == 1){
    y = x_ordinary[1] - n1
    if (y == n-1){
      if(log.p){
        pmf_computed <- shiftConvolveLog(probs, (y-1))
      }else{
        pmf_computed <- exp(shiftConvolveLog(probs, (y-1)))
      }

      U <- uniroot(shiftedmean, interval = c(-50,50), probs, y-1, tol = 10^-16) #find a suitable shifting parameter
      t0 <- U$root
    }else{
      if(log.p){
        pmf_computed <- shiftConvolveLog(probs, y)
      }else{
        pmf_computed <- exp(shiftConvolveLog(probs, y))
      }
      U <- uniroot(shiftedmean, interval = c(-50,50), probs, y, tol = 10^-16) #find a suitable shifting parameter
      t0 <- U$root
    }
  }else{
    pbd_mean <- sum(probs)
    diff_x_pbd_mean <- x_ordinary - pbd_mean
    count_pos = sum(diff_x_pbd_mean > 0)
    count_neg = sum(diff_x_pbd_mean < 0)
    if((count_pos > 0) & (count_neg > 0)){
      if(log.p){
        pmf_computed <- shiftConvolveLog(probs, S0 = -1)
      }else{
        pmf_computed <- exp(shiftConvolveLog(probs, S0 = -1))
      }
      t0 <- 0
    }else{
      closest <- x_ordinary[which.min(diff_x_pbd_mean)]
      if(log.p){
        pmf_computed <- shiftConvolveLog(probs, closest)
      }else{
        pmf_computed <- exp(shiftConvolveLog(probs, closest))
      }
      U <- uniroot(shiftedmean, interval = c(-50,50), probs, closest, tol = 10^-16) #find a suitable shifting parameter
      t0 <- U$root
    }
  }

  pmf_computed = c(rep(-Inf, n1), pmf_computed)
  pmf_computed = c(pmf_computed, rep(-Inf,n0))

  for(i in 1:length(x_ordinary)){
    if(lower.tail){
      if(t0 > 0){
        if(log.p){
          log_right_tail <- logsum(pmf_computed[c((x_ordinary[i]+2):length(pmf_computed))])
          p <- 1-exp(log_right_tail)
          lp <- log(p)
          return_dist[x_return_index[i]]<-lp
        }else{
          right_tail <- sum(pmf_computed[c((x_ordinary[i]+2):length(pmf_computed))])
          return_dist[x_return_index[i]]<-1-right_tail
        }
      }else{
        if(log.p){
          lp <- logsum(pmf_computed[c(0:(x_ordinary[i]+1))])
          return_dist[x_return_index[i]]<-lp
        }else{
          p <- sum(pmf_computed[c(0:(x_ordinary[i]+1))])
          return_dist[x_return_index[i]]<-p
        }
      }
    }else{
      if (t0 > 0){
        if(log.p){
          lp <- logsum(pmf_computed[c((x_ordinary[i]+2):length(pmf_computed))])
          return_dist[x_return_index[i]]<-lp
        }else{
          p <- sum(pmf_computed[c((x_ordinary[i]+2):length(pmf_computed))])
          return_dist[x_return_index[i]]<-p
        }
      }
      else{
        if(log.p){
          log_left_tail <- logsum(pmf_computed[c(0:x_ordinary[i]+1)])
          p <- 1-exp(log_left_tail)
          lp <- log(p)
          return_dist[x_return_index[i]]<-lp
        }else{
          left_tail <- sum(pmf_computed[c(0:x_ordinary[i]+1)])
          return_dist[x_return_index[i]]<-1-left_tail
        }
      }
    }
  }
  return(return_dist)
}

shift_pb_density <- function(x, probs, log.p=FALSE){
  probs = process(probs)
  n0 = probs[1]
  n1 = probs[2]
  probs = probs[3:length(probs)]
  n <- length(probs)+1

  x_ordinary = c()
  x_return_index = c()
  pmf_computed = NULL
  return_density = vector(mode = "numeric", length = length(x))
  for(i in 1:length(x)){
    edge <- edge_density(x[i], n, n1, log.p)
    if(!is.null(edge)){
      return_density[i] <- edge
    }else{
      x_ordinary = c(x_ordinary, x[i])
      x_return_index = c(x_return_index, i)
    }
  }

  n_x_ordinary = length(x_ordinary)
  if(n_x_ordinary == 0){
    return(return_density)
  }else if(n_x_ordinary == 1){
    y = x_ordinary[1] - n1
    if (y == n-1){
      pmf_computed <- shiftConvolveLog(probs, (y-1))
      U <- uniroot(shiftedmean, interval = c(-50,50), probs, y-1, tol = 10^-16) #find a suitable shifting parameter
      t0 <- U$root
    }else{
      pmf_computed <- shiftConvolveLog(probs, y)
      U <- uniroot(shiftedmean, interval = c(-50,50), probs, y, tol = 10^-16) #find a suitable shifting parameter
      t0 <- U$root
    }
  }else{
    pbd_mean <- sum(probs)
    diff_x_pbd_mean <- x_ordinary - pbd_mean
    count_pos = sum(diff_x_pbd_mean > 0)
    count_neg = sum(diff_x_pbd_mean < 0)
    if((count_pos > 0) & (count_neg > 0)){
      pmf_computed <- shiftConvolveLog(probs, S0 = -1)
      t0 <- 0
    }else{
      closest <- x_ordinary[which.min(diff_x_pbd_mean)]
      pmf_computed <- shiftConvolveLog(probs, closest)
      U <- uniroot(shiftedmean, interval = c(-50,50), probs, closest, tol = 10^-16) #find a suitable shifting parameter
      t0 <- U$root
    }
  }

  pmf_computed = c(rep(-Inf, n1), pmf_computed)
  pmf_computed = c(pmf_computed, rep(-Inf,n0))

  for(i in 1:length(x_ordinary)){

    log_density = pmf_computed[x_ordinary[i]+1]
    density = exp(log_density)
    if(log.p){
      return_density[x_return_index[i]]<-log_density
    }else{
      return_density[x_return_index[i]]<-density
    }
  }
  return(return_density)
}

edge_density <- function(x, n, n1, log.p){
  y <- x - n1
  if(y <= 0){
    if(log.p){
      return(-Inf)
    }else{
      return(0)
    }
  }

  if(y >= n){
    if(log.p){
      return(-Inf)
    }else{
      return(0)
    }
  }
  return(NULL)
}

edge_prob <- function(x, n, n1, lower.tail, log.p){
  y <- x - n1
  if(y <= 0){
    if(lower.tail){
      if(log.p){
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      if(log.p){
        return(0)
      }else{
        return(1)
      }
    }

  }

  if(y >= n){
    if(lower.tail){
      if(log.p){
        return(0)
      }else{
        return(1)
      }
    }else{
      if(log.p){
        return(-Inf)
      }else{
        return(0)
      }
    }

  }
  return(NULL)
}

dc_pb_dist <- function(x, probs, lower.tail = TRUE,log.p=FALSE){
  probs = process(probs)
  n0 = probs[1]
  n1 = probs[2]
  probs = probs[3:length(probs)]
  n <- length(probs)+1

  pmf_computed = NULL
  return_dist = vector(mode = "numeric", length = length(x))
  for(i in 1:length(x)){
    edge <- edge_prob(x[i], n, n1, lower.tail, log.p)
    if(!is.null(edge)){
      return_dist[i] <- edge

    }else{
      if(is.null(pmf_computed)){
        #Computing distribution first time if not computed already
        if(log.p){
          pmf_computed <- dcConvolvePairedLog(probs)
          pmf_computed <- c(rep(-Inf, n1), pmf_computed)
          pmf_computed <- c(pmf_computed, rep(-Inf,n0))
        }else{
          pmf_computed <- dcConvolvePaired(probs)
          pmf_computed <- c(rep(0, n1), pmf_computed)
          pmf_computed <- c(pmf_computed, rep(0,n0))
        }
      }

      if(log.p){
        if(lower.tail){
          return_dist[i] <- logsum(pmf_computed[c(0:(x[i]+1))])
        }else{
          return_dist[i] <- logsum(pmf_computed[c((x[i]+2):length(pmf_computed))])
        }
      }else{
        if(lower.tail){
          return_dist[i] <- sum(pmf_computed[c(0:(x[i]+1))])
        }else{
          return_dist[i] <- sum(pmf_computed[c((x[i]+2):length(pmf_computed))])
        }
      }

    }
  }
  return(return_dist)
}

dc_pb_density <- function(x, probs, log.p=FALSE){
  probs = process(probs)
  n0 = probs[1]
  n1 = probs[2]
  probs = probs[3:length(probs)]
  n <- length(probs)+1

  pmf_computed = NULL
  return_density = vector(mode = "numeric", length = length(x))
  for(i in 1:length(x)){
    edge <- edge_density(x[i], n, n1, log.p)
    if(!is.null(edge)){
      return_density[i] <- edge
    }else{
      if(is.null(pmf_computed)){
        #Computing distribution first time if not computed already
        if(log.p){
          pmf_computed <- dcConvolvePairedLog(probs)
          pmf_computed <- c(rep(-Inf, n1), pmf_computed)
          pmf_computed <- c(pmf_computed, rep(-Inf,n0))
        }else{
          pmf_computed <- dcConvolvePaired(probs)
          pmf_computed <- c(rep(0, n1), pmf_computed)
          pmf_computed <- c(pmf_computed, rep(0,n0))
        }
      }
      return_density[i] <- pmf_computed[x[i]+1] # Previous if block would have already decided log.p or not
    }
  }
  return(return_density)
}
