#' Distributional/Tail Structural Change Test
#'
#' @description
#' The function supisd() performs the distributional/tail structural change test
#' proposed in Lu and Ker (2024).
#' Note: This version is for test only. Please contact the authors if have any
#' questions. Thank you.
#'
#'
#' @param data data in chronological order to be tested.
#' @param breakstart a numeric value giving the starting period of all possible
#'        breaks. Default leaves 20 observations from the start of the data.
#' @param breakend a numeric value giving the ending period of all possible
#'        breaks. Default leaves 20 observations from the end of the data.
#' @param evalrange.l lower bound of the evaluation support. Either a numeric
#'        value from 0 to 1 indicating the percentage of data range from the
#'        minimum of the data, or a character sting "mean" indicating the mean
#'        of the data. Default is 10 percent of data range to the left of the
#'        minimum of the data.
#' @param evalrange.r upper bound of the evaluation support. Either a numeric
#'        value from 0 to 1 indicating the percentage of data range from the
#'        minimum of the data, or a character sting "mean" indicating the mean
#'        of the data. Default is 10% of the data range to the right of the
#'        maximum of the data.
#' @param width space between the points in the grid for kernel estimation.
#'        Default results in 500 points at which the density is estimated.
#' @param null method to recover null distribution of the test. Either "ran"
#'        indicating randomization method or "asy" indicating asymptotic
#'        distribution. Default is randomization method.
#' @param S number of simulations for generating the null distribution of the
#'        test statistic. Default is 1000.
#'
#' @details
#' A Sup-type nonparametric test for structural change with an unknown break
#' in time-series data as proposed in Lu and Ker (2024). The test has power
#' against structural change in any moment and can be tailored to a specific
#' range of the underlying distribution. The function uses Gaussian kernel
#' function and leave-one-out cross validation for finding the optimal bandwidth
#' for kernel estimation. Test p value is found via randomization method by
#' default setting as it performs better empirically than the asymptotic null
#' distribution. However, the latter is also made available in cases of large
#' samples. It is suggested by simulations that the possible breaks to be
#' tested should not be too close to the start/end of data (leaving at least 20
#' observations from the start/end is recommended).
#'
#' @return test statistic, estimated break and p value
#'
#' @import stats
#'
#' @author Hanjun Lu and Alan Ker
#'
#' @references Hanjun Lu, Alan P Ker, Testing for distributional structural
#' change with unknown breaks: application to pricing crop insurance contracts,
#' Journal of the Royal Statistical Society Series C: Applied Statistics, 2024;,
#' qlae020, https://doi.org/10.1093/jrsssc/qlae020
#'
#' @examples
#' #width and S in examples below are set to reduce running time.
#' #Still, it may take a few seconds.
#' data <- rnorm(100)
#' #Tests structural change over the entire domain of data
#' supisd(data, width=0.1, S=200)
#' #Test structural change within the lower tail (below mean) of data
#' supisd(data, evalrange.r="mean", width=0.1, S=200)
#' #Test structural change within the upper tail (above mean) of data
#' supisd(data, evalrange.l="mean", width=0.1, S=200)
#'
#' @export

supisd <- function(data, breakstart = NULL, breakend = NULL, evalrange.l = NULL,
                   evalrange.r = NULL, width = NULL, null = "ran", S = NULL){

  T <- length(data)
  if (T <= 50) {warning("Small sample size")}

  if (!is.null(breakstart) == FALSE) {breakstart <- 20}
  if (!is.null(breakend) == FALSE) {breakend <- T-20}
  if (breakstart > breakend) {stop("Wrong possible breaks")}
  if (breakstart < 20 ) {warning("Possible break is close to the start of data")}
  if (breakstart > T-20 ) {warning("Possible break is close to the end of data")}

  breakset<-seq(breakstart,breakend,1)

  range.data<-max(data)-min(data)


  if (!is.null(evalrange.l) == TRUE & !is.null(evalrange.r) == TRUE){ #both have values
    if (is.numeric(evalrange.l) == TRUE & is.numeric(evalrange.r) == TRUE){
      if (evalrange.l < 0 | evalrange.l >= 1) {stop("Inaccurate lower evaluation range")}
      if (evalrange.r <= 0 | evalrange.r > 1) {stop("Inaccurate lower evaluation range")}
      if (evalrange.l >= evalrange.r) {stop("Inaccurate lower/upper evaluation range")}
      grid.l <- min(data)+evalrange.l*range.data
      grid.r <- min(data)+evalrange.r*range.data
    } else if (is.numeric(evalrange.l) == TRUE & is.character(evalrange.r) == TRUE){
      if (evalrange.l < 0 | evalrange.l >= 1) {stop("Inaccurate lower evaluation range")}
      if (evalrange.r != "mean") {stop("Inaccurate upper evaluation range")}
      if (mean(data) <= (min(data)+evalrange.l*range.data)) {stop("Inaccurate lower/upper evaluation range")}
      grid.l <- min(data)+evalrange.l*range.data
      grid.r <- mean(data)
    } else if (is.character(evalrange.l) == TRUE & is.numeric(evalrange.r) == TRUE){
      if (evalrange.r <= 0 | evalrange.r > 1) {stop("Inaccurate upper evaluation range")}
      if (evalrange.l != "mean") {stop("Inaccurate lower evaluation range")}
      if (mean(data) >= (min(data)+evalrange.r*range.data)) {stop("Inaccurate lower/upper evaluation range")}
      grid.l <- mean(data)
      grid.r <- min(data)+evalrange.r*range.data
    } else if (is.character(evalrange.l) == TRUE & is.character(evalrange.r) == TRUE){
      stop("Inaccurate lower/upper evaluation range")
    }
  } else if (!is.null(evalrange.l) == TRUE & !is.null(evalrange.r) == FALSE){ #lower has value
    if (is.numeric(evalrange.l) == TRUE){
      if (evalrange.l < 0 | evalrange.l >= 1) {stop("Inaccurate lower evaluation range")}
      grid.l <- min(data)+evalrange.l*range.data
      grid.r <- max(data)+0.1*range.data
    } else if (is.character(evalrange.l) == TRUE){
      if (evalrange.l != "mean") {stop("Inaccurate lower evaluation range")}
      grid.l <- mean(data)
      grid.r <- max(data)+0.1*range.data
    }
  } else if (!is.null(evalrange.l) == FALSE & !is.null(evalrange.r) == TRUE){ #upper has value
    if (is.numeric(evalrange.r) == TRUE){
      if (evalrange.r <= 0 | evalrange.r > 1) {stop("Inaccurate upper evaluation range")}
      grid.l <- min(data)-0.1*range.data
      grid.r <- min(data)+evalrange.r*range.data
    } else if (is.character(evalrange.r) == TRUE){
      if (evalrange.r != "mean") {stop("Inaccurate upper evaluation range")}
      grid.l <- min(data)-0.1*range.data
      grid.r <- mean(data)
    }
  } else if (!is.null(evalrange.l) == FALSE & !is.null(evalrange.r) == FALSE){ #neither has value
    grid.l <- min(data)-0.1*range.data
    grid.r <- max(data)+0.1*range.data
  }




  if (!is.null(width) == FALSE) {width <- (grid.r-grid.l)/499}
  grid <- seq(grid.l,grid.r,width)
  est_f <- matrix(0, length(grid))
  est_g <- matrix(0, length(grid))

  ISDresult <- matrix(0, length(breakset))
  rownames(ISDresult) <- seq(breakstart,breakend,1)

  if (!is.null(S) == FALSE) {S <- 1000}

  for (t in breakset){
    f <- data[1:t]
    g <- data[(t+1):T]

    n_f <- length(f)
    n_g <- length(g)

    #optimal h for f
    mlcv.f <- matrix(0, length(f))
    MLCV.f <- function(h_f){
      for(i in 1:length(f)){
        mlcv.f[i] <- (sum(dnorm((f[i] - f)/h_f))-dnorm(0))/((n_f-1)*h_f)
      }
      likelihood.f <- sum(log(mlcv.f))/n_f
      return(likelihood.f)
    }
    H_f <- optimize(MLCV.f,c(0.01, 100),maximum = TRUE)
    h_f <- H_f[["maximum"]]

    #optimal h for g
    mlcv.g <- matrix(0, length(g))
    MLCV.g <- function(h_g){
      for(i in 1:length(g)){
        mlcv.g[i] <- (sum(dnorm((g[i] - g)/h_g))-dnorm(0))/((n_g-1)*h_g)
      }
      likelihood.g <- sum(log(mlcv.g))/n_g
      return(likelihood.g)
    }
    H_g <- optimize(MLCV.g,c(0.01, 100),maximum = TRUE)
    h_g <- H_g[["maximum"]]

    #kernel estimates for f
    for(i in 1:length(grid)){
      est_f[i] <- sum(dnorm((grid[i]-f)/h_f))/(n_f*h_f)
    }

    #kernel estimates for g
    for(i in 1:length(grid)){
      est_g[i] <- sum(dnorm((grid[i]-g)/h_g))/(n_g*h_g)
    }

    #ISD estimates
    ISD_hat <- sum((est_f-est_g)^2)*width

    t <-toString(t)
    ISDresult[t,] <- ISD_hat
  }

  #estimate results
  ISD_sup <- max(ISDresult)
  t_star <- as.numeric(rownames(ISDresult)[which.max(apply(ISDresult,MARGIN=1,max))])


  ################################### Bootstrap Null
  if (null == "ran"){
    start.time<-Sys.time()
    ISD_sup.null <- matrix(0, S)
    est_f.B <- matrix(0, length(grid))
    est_g.B <- matrix(0, length(grid))

    ISDresult.null <- matrix(0, length(breakset))
    rownames(ISDresult.null) <- seq(breakstart,breakend,1)

    f <- data[1:t_star]
    g <- data[(t_star+1):T]
    n_f <- length(f)
    n_g <- length(g)
    #optimal h for f
    mlcv.f <- matrix(0, length(f))
    MLCV.f <- function(h_f){
      for(i in 1:length(f)){
        mlcv.f[i] <- (sum(dnorm((f[i] - f)/h_f))-dnorm(0))/((n_f-1)*h_f)
      }
      likelihood.f <- sum(log(mlcv.f))/n_f
      return(likelihood.f)
    }
    H_f <- optimize(MLCV.f,c(0.001, 100),maximum = TRUE)
    h_f <- H_f[["maximum"]]
    #optimal h for g
    mlcv.g <- matrix(0, length(g))
    MLCV.g <- function(h_g){
      for(i in 1:length(g)){
        mlcv.g[i] <- (sum(dnorm((g[i] - g)/h_g))-dnorm(0))/((n_g-1)*h_g)
      }
      likelihood.g <- sum(log(mlcv.g))/n_g
      return(likelihood.g)
    }
    H_g <- optimize(MLCV.g,c(0.001, 100),maximum = TRUE)
    h_g <- H_g[["maximum"]]

    for (s in 1:S){
      temp.var <- sample(data,T,replace = FALSE)
      for (t in breakset){
        f.B <- temp.var[1:t]
        g.B <- temp.var[(t+1):T]

        n_f <- length(f.B)
        n_g <- length(g.B)

        #optimal h for f
        h_f.B <- h_f

        #optimal h for g
        h_g.B <- h_g

        for(i in 1:length(grid)){
          est_f.B[i] <- sum(dnorm((grid[i]-f.B)/h_f.B))/(n_f*h_f.B)
        }

        for(i in 1:length(grid)){
          est_g.B[i] <- sum(dnorm((grid[i]-g.B)/h_g.B))/(n_g*h_g.B)
        }

        #ISD estimates
        ISD_hat <- sum((est_f.B-est_g.B)^2)*width

        t <-toString(t)
        ISDresult.null[t,] <- ISD_hat
      }
      ISD_sup.null[s] <- max(ISDresult.null)

      end.time<-Sys.time()
      time.taken<-as.numeric(round(end.time - start.time))
      if (time.taken!=0 & time.taken%%60==0){
        print("Recovering the null distribution using randomization method could take some time. Thank you for your patience.")
      }
    }

    #p value
    PV <- rep(0,S)
    PV[ISD_sup.null>ISD_sup] <- 1

    pvalue <- sum(PV)/S

  } else if (null == "asy"){
    if (T < 1000){warning("Randomization method is recommended for recovering
                          the null distribution")}
    ISD_sup.std <- ISD_sup*sqrt(T)
    pvalue <- 1-asy.null[which.min(abs(ISD_sup.std-asy.grid))]
  }

  SupISD_results <- list("test statistic" = ISD_sup, "estimated break" = t_star, "p value" = pvalue)

  return(SupISD_results)
}


