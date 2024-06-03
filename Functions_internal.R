

ifrm <- function(obj, env = globalenv()) {
  obj <- deparse(substitute(obj))
  if(exists(obj, envir = env)) {
    rm(list = obj, envir = env)
  }
}

prepare_data <- function(Mat, #matrix of each observation with left/right censoring and left/right truncation. 
                         #if only 2 columns assume to truncation
                         expr = as.formula(Mat[,1]~1), # formula for the covariates 
                         data = NULL, # dataset for covariates (only allows NULL if intercept)
                         exact = FALSE, # whether observations are allowed to be exact observations
                         drop = "soft", # 'no' stops function if observation is invalid, 'soft' drops invalid and 'hard' drops                                                #meaningless, defaults to no for consistency
                         message = TRUE #return message about data
){
  if(ncol(Mat) ==2){
    trunc <- FALSE 
    Mat <- cbind(Mat,Mat)
    Mat[,3] <- 0
    Mat[,4] <- Inf
  }else if(ncol(Mat) ==3){stop("observation data has only 3 columns")
  }else{
    Mat <- Mat[,1:4]
    trunc <- FALSE
    if( any(c((Mat[,3] > 0), (Mat[,4] < Inf)))){
      trunc <- TRUE
    }
  }
  text <- ""
  if(!(drop %in% c("no","soft","hard")) ){
    text <- cat(text, "\n ", ("\"drop\" variable not declared properly, defaults to \"no\". "))
    drop <- "no"
  }
  
  MM <- model.matrix(expr, data = data)
  
  dat <- as.data.frame(cbind(Mat,MM))
  colnames(dat) <- c("CL","CR","TL","TR",colnames(MM))
  
  if(trunc){
    sub_test <- TRUE
    if(min(Mat[,3]) > 0 | max(Mat[,4]) < Inf){
      sub_test <- FALSE
    } else{
      right <- 0
      TRUNC <- cbind(Mat[,3],Mat[,4])
      TRUNC <- TRUNC[which(!duplicated(TRUNC)),]
      while(TRUE){
        right <- max(TRUNC[which(TRUNC[,1] <= right),2])
        if(right == Inf){
          break
        } else if(max(TRUNC[,1] <= right) == 0){
          sub_test <- FALSE
          break
        }
      }
    }
    if(!sub_test){text <- paste(text, "\n ","Warning: The Union of the support of the truncated random variables
          does not cover the positive reals. The estimate of the survival function will not
          correspond to an untruncated distribution, but the distribution conditioned on being
          on the joint support.")
    }
  }
  test_good <- TRUE
  if(any(Mat[,1] > Mat[,2])){
    text <- paste(text, "\n ","Warning: The left endpoint of the censoring interval
          is larger than the right endpoint in ", sum(Mat[,1] > Mat[,2]), "datapoints.", sep = "")
    test_good <- FALSE}
  if(any(Mat[,3] > Mat[,1])){
    text <- paste(text, "\n ","Warning: The left endpoint of the truncation interval
          is larger than the left endpoint of censoring in ", sum(Mat[,3] > Mat[,1]), "datapoints.", sep = "")
    test_good <- FALSE}
  if(any(Mat[,4] < Mat[,2])){
    text <- paste(text, "\n ","Warning: The right endpoint of the censoring interval
          is larger than the right endpoint of truncation in ", sum(Mat[,4] < Mat[,2]), "datapoints.", sep = "")
    test_good <- FALSE}
  if((!exact) & any(Mat[,1] == Mat[,2]) ){
    text <- paste(text, "\n ","Warning: ", sum( Mat[,1] == Mat[,2]),
                  " datapoints are exact observation, when this is not wanted.",sep = "")
    test_good <- FALSE
  }
  if(drop == "no" & !test_good){
    if(message){print(text)}
    stop()}
  dat <- dat[which(!(dat$CL > dat$CR)),]
  dat <- dat[which(!(dat$CL < dat$TL)),]
  dat <- dat[which(!(dat$CR > dat$TR)),]
  if(!exact){
    dat <- dat[which(!(dat$CL == dat$CR)),]
  }
  if(drop == "soft"){
    if(message){print(text)}
    return( list( data = dat, trunc = trunc))}
  dat <- dat[which( !((dat$CL == dat$TL) & (dat$CR == dat$TR))), ]
  return( list( data = dat, trunc = trunc, message = text))
}

print.turnbull_est <- function(obj){
  print(obj$s)
}

turnbull.internal.IMat <- function(vec_l,vec_r, q , p){
  if(any(vec_l>vec_r)){
    stop(paste("turnbull.internal.IMat recieved incompatible values ", sum(vec_l>vec_r), " times!", sep = ""))
  }
  m <- length(q)
  MAT <- cbind(vec_l,vec_r)
  MAT <- as.matrix(MAT[which(!duplicated(MAT)),])
  if(dim(MAT)[2] == 1){MAT <- t(MAT)}
  MAT <- as.matrix(MAT[order(MAT[,1],MAT[,2]),])
  if(dim(MAT)[2] == 1){MAT <- t(MAT)}
  IMat_n <- nrow(MAT)
  IMat_w <- numeric(IMat_n)
  for(i in 1:IMat_n){
    IMat_w[i] <- sum(MAT[i,1] == vec_l & MAT[i,2] ==vec_r)
  }
  
  IMat <- matrix(0,nrow = IMat_n,ncol  = m)
  
  for( i in 1:m){
    IMat[,i] <-( MAT[,1] <=q[i] ) & (p[i] <= MAT[,2])
  }
  colnames(IMat) <- paste(q,p,sep = ",")
  L <- list( IMat = IMat,IMat_w=IMat_w,IMat_n=IMat_n)
  return(L)
}

turnbull.Obs.FI <- function(obj, cov = TRUE){
  CL <- obj$left_censor
  CR <- obj$right_censor
  TL <- obj$left_trunc
  TR <- obj$right_trunc
  q <- obj$q
  p <- obj$p
  s <- obj$s
  m <- length(q)
  
  if(all(TL == 0 & TR == Inf)){trunc <- FALSE}else{trunc <- TRUE}
  foo <- turnbull.internal.IMat(CL,CR,q,p)
  Alpha <- foo[[1]]
  w_c <- foo[[2]]
  n_c <- foo[[3]]
  if(trunc){
    foo <- turnbull.internal.IMat(TL,TR,q,p)
    Beta <- foo[[1]]
    w_t <- foo[[2]]
    n_t <- foo[[3]]
  }else{
    Beta <- 0
  }
  
  OFI <-  NULL
  
  theta <- s[-m]
  index <- rowSums(Alpha) != m
  Alpha_m <- Alpha[index,m]
  Alpha <- Alpha[index,-m]
  
  Alpha_res_nm <- rowSums((rep(1,sum(index))%*%t(theta)*Alpha))
  Alpha_res_nm <- (1-Alpha_m)*w_c[index]*ifelse(Alpha_res_nm == 0,0,Alpha_res_nm^-2)
  
  Alpha_res_m <- 1-rowSums((rep(1,sum(index))%*%t(theta)*(1-Alpha)))
  Alpha_res_m <- Alpha_m*w_c[index]*Alpha_res_m^-2
  if(trunc){
    index <- rowSums(Beta) != m
    Beta_m <- Beta[index,m]
    Beta <- Beta[index,-m]
    Beta_res_nm <- rowSums((rep(1,sum(index))%*%t(theta)*Beta))
    Beta_res_nm <- (1-Beta_m)*w_t[index]*ifelse(Beta_res_nm == 0,0,Beta_res_nm^-2)
    
    Beta_res_m <- 1-rowSums((rep(1,sum(index))%*%t(theta)*(1-Beta)))
    Beta_res_m <- Beta_m*w_t[index]*Beta_res_m^-2
  }
  AFI <- BFI <- matrix(0,m-1,m-1)
  for(i in 1:(m-1)){
    for(j in i:(m-1)){
      AFI[i,j] <- AFI[j,i] <-  sum(Alpha[,i]*Alpha[,j]*Alpha_res_nm)+sum((1-Alpha[,i])*(1-Alpha[,j])*Alpha_res_m)
    }
  }
  if(trunc){
    for(i in 1:(m-1)){
      for(j in i:(m-1)){
        BFI[i,j] <- BFI[j,i] <-  sum(Beta[,i]*Beta[,j]*Beta_res_nm)+sum((1-Beta[,i])*(1-Beta[,j])*Beta_res_m)
      }
    }
  }
  result <- AFI-BFI
  colnames(result) <- rownames(result) <-names(s)[-m]
  return(result)
}

turnbull.CI <- function(obj,
                        type = "transSurv",
                        alpha = 0.05
){
  s <- obj$s
  try(COV <- solve(turnbull.Obs.FI(obj)))
  if(!exists("COV")){
    stop("Observed Fisher Information not invertible.")
  }else if(min(diag(COV)) <= 0){
    stop("Diagonal entry in Observed Fisher Information non-positive.")
  }else{
    if(type == "NormSurv"){
      J <- matrix(0,ncol = nrow(COV),nrow = nrow(COV))
      J[upper.tri(J,diag = T)] <- -1
      se <- sqrt(diag(t(J)%*%COV %*%J))
      lower <- s - qnorm(1-alpha/2)*se
      upper <- s + qnorm(1-alpha/2)*se
      
    }
    else{
      if(type != "transSurv"){
        message("Did not recognise CI type, default is \"transSurv\" ")
      }
      n <- nrow(COV)
      logit <- function(x){exp(x)/(1+exp(x))}
      invlogit <- function(x){-log((1-x)/x)}
      dinvlogit <- function(x){1/(x*(1-x))}
      # J <- matrix(0,ncol = nrow(COV),nrow = nrow(COV))
      # J[upper.tri(J,diag = T)] <- -1
      # COV <- sqrt(diag(t(J)%*%COV %*%J))
      # COV <- c(0,0,rep(COV,each = 2),1,1)
      
      Nablag <- matrix(0,ncol = n,nrow = n)
      for(i in 1:n){
        Nablag[1:i,i] <- -dinvlogit(1-sum(s[1:i]))
      }
      phidiag <- sqrt(diag(t(Nablag)%*%COV %*%Nablag))
      lower <- logit(invlogit(1-cumsum(s[1:n]))-qnorm(1-alpha/2)*phidiag)
      upper <- logit(invlogit(1-cumsum(s[1:n]))+qnorm(1-alpha/2)*phidiag)
    }
    
    names(lower) <- names(upper) <- names(s[1:n])
    L <- list(CI.lower = lower, CI.upper = upper)
    return(L)
  }
}

plot.turnbull_est <- function(obj,type = "S",add = FALSE,...){
  L <- list(...)
  if( !("type" %in% names(L))){
    L$type = "l"
  }
  
  
  s <- obj$s
  p <- obj$p
  q <- obj$q
  y <- c(1,1,rep(1-cumsum(s), each = 2))
  
  L$y <- c(y,y[length(y):1])
  x1 <- c(0,sort(rep(p,each = 2)),p[length(p)])
  x2 <- c(0,sort(rep(q,each = 2)),q[length(q)])
  L$x <- c(x1,x2[length(x2):1])
  if(tolower(type) %in% c("cumulativehazard","cumhaz","h")){
    L$y <- -log(L$y)
    if(!("xlab" %in% names(L))){L$xlab = "time"}
    if(!("ylab" %in% names(L))){L$ylab = "H(t)"}
  }else if(tolower(type) %in% c("logcumulativehazard","logcumhaz","logh","weibullcheck","weibull","weib")){
    L$y <- ifelse(L$y > 0,L$y,NA)
    L$y <- log(-log(L$y))
    if(!("xlab" %in% names(L))){L$xlab = "time"}
    if(!("ylab" %in% names(L))){L$ylab = "log H(t)"}
  }else if(tolower(type) %in% c("loglogcumulativehazard","loglogcumhaz","loglogh","gompertzcheck","gomp")){
    L$y <- ifelse(L$y > 0,L$y,NA)
    L$y <- log(-log(L$y))
    L$x <- log(L$x)
    if(!("xlab" %in% names(L))){L$xlab = "log time"}
    if(!("ylab" %in% names(L))){L$ylab = "log H(t)"}
  }else if(tolower(type) %in% c("exppowercheck","exppower","ep")){
    L$y <- ifelse(L$y > 0,L$y,NA)
    L$y <- log(log(1-log(L$y)))
    L$x <- log(L$x)
    if(!("xlab" %in% names(L))){L$xlab = "log time"}
    if(!("ylab" %in% names(L))){L$ylab = "log log (1 - H(t))"}
  }
  else{
    if(!("ylim" %in% names(L))){L$ylim = c(0,1)}
    if(!("xlab" %in% names(L))){L$xlab = "time"}
    if(!("ylab" %in% names(L))){L$ylab = "S(t)"}
  }
  
  if(add){
    do.call(lines,c(L))
  }else{
    do.call(plot,c(L))
  }
}

turnbull <- function(CL, # Left endpoint of censoring interval (0 for left censored)
                     CR, # Right endpoint of censoring interval (Inf for Right censored)
                     TL =rep(0,length(CL)) , # Left endpoint of truncation interval
                     TR = rep(Inf,length(CL)), # right endpoint of truncation interval
                     tol = 10^-4, # relative tolerance in iteration. (number of turnbull intervals * max abs coordinate difference)
                     verbose = FALSE, # print every n'th statement for n = verbose
                     message = TRUE, # print ending message
                     start = "equal",# initialise probability vector
                     include0 = FALSE #whether S(0+) = 1 or S(0+) < 1 
){
  start_time <- Sys.time()
  
  foo <- prepare_data(cbind(CL,CR,TL,TR), drop = "hard",message = message, exact = FALSE)
  
  if(message){print(foo$message)}
  trunc <- foo[[2]]
  
  foo <- foo[[1]]
  
  CL <- foo[,1]
  CR <- foo[,2]
  TL <- foo[,3]
  TR <- foo[,4]
  
  L <- sort(unique(c(CL,TR)))
  L <- as.data.frame(L[L<Inf])
  R <- sort(unique(c(CR,TL)))
  if(include0){
    R <- as.data.frame(R)
  }else{
    R <- as.data.frame(R[R>0])
  }
  colnames(L) <-colnames(R) <- "numb"
  L$name <- "L"
  R$name <- "R"
  
  I <- rbind(L,R)
  I <- I[order(I$numb),]
  
  Match <- stringr::str_locate_all(paste(I$name, collapse = ""), "LR")[[1]]
  
  #create turnbull (inner) intervals 
  q <- p <- numeric(nrow(Match))
  for(i in 1:nrow(Match)){
    q[i] <- I$numb[Match[i,1]]
    p[i] <- I$numb[Match[i,2]]
  }
  
  n <- length(CL)
  
  m <- length(q)
  if(m > 200){
    message(paste("warning: estimator uses ",m," inner (turnbull) intervals. It is recomended to discretise the input by interval censoring."))
  }
  if(m > 1000){
    errorCondition("warning: Too many inner (turnbull) intervals created. discretise the input by interval censoring.")
  }
  
  if(!(start %in% c("random","both","tail","begin","equal"))){
    message("variable \"start\"  did not match any options, reverted to default: \"equal\"")
  }
  if(start =="random"){
    s <- runif(m)
  }else if(start =="both"){
    s <- c(m/2,rep(1,m-2),m/2)
  }else if(start =="tail"){
    s <- c(rep(1,m-1),m)
  }else if(start == "begin"){
    s <- c(m,rep(1,m-1))
  }else{
    s <- rep(1,m) 
  }
  s <- s/sum(s)
  
  # C_MAT <- cbind(CL,CR)
  # C_MAT <- C_MAT[which(!duplicated(C_MAT)),]
  # if(is.null(dim(C_MAT))){
  #   stop("Warning: All censoring intervals are identical.")
  # }
  # C_MAT <- C_MAT[order(C_MAT[,1],C_MAT[,2]),]
  # n_c <- nrow(C_MAT)
  # w_c <- numeric(n_c)
  # for(i in 1:n_c){
  #   w_c[i] <- sum(C_MAT[i,1] == CL & C_MAT[i,2] ==CR)
  # }
  # 
  # Alpha <- matrix(0,nrow = n_c,ncol  = m)
  # 
  # for( i in 1:m){
  #   Alpha[,i] <-( C_MAT[,1] <=q[i] ) & (p[i] <= C_MAT[,2])
  # }
  # colnames(Alpha) <- paste(q,p,sep = ",")
  foo <- turnbull.internal.IMat(CL,CR,q,p)
  Alpha <- foo[[1]]
  w_c <- foo[[2]]
  n_c <- foo[[3]]
  if(n_c == 1){
    stop("All censoring interval are identical!")
  }
  if(trunc){
    foo <- turnbull.internal.IMat(TL,TR,q,p)
    Beta <- foo[[1]]
    w_t <- foo[[2]]
    n_t <- foo[[3]]
  }else{
    Beta <- 0
  }
  
  itt <-1
  while(TRUE){
    mu_. <- as.numeric(w_c%*%((Alpha*(rep(1,n_c)%*%t(s)))/((Alpha%*%s)%*%t(rep(1,m)))))
    
    if(trunc){
      nu_. <- as.numeric(w_t%*%(((1-Beta)*rep(1,n_t)%*%t(s))/((Beta%*%s)%*%t(rep(1,m)))))
      M <- sum(nu_.) + sum(mu_.)
    }else{
      nu_. <- rep(0,m)
      M <- n
    }
    
    pi <- (mu_. + nu_.)/M
    if((verbose>0) &(itt%%verbose == 0) ){print(max(abs(s-pi)*m))}
    if(max(abs(s-pi)*m) < tol){break}
    s <- pi
    itt <- itt +1
  }
  names(s) <- paste(q,p,sep = ",")
  
  if(message){print(Sys.time() - start_time)}
  
  obj <- list( s = s, turnbull_int = names(s),q = q ,p = p,left_censor = CL,
               right_censor = CR, left_truncation = TL , right_truncation = TR,
               n_obs = n, n_int = m, start = start, include0 = include0)
  class(obj) <- "turnbull_est"
  return(obj)
}

turnbull.CI.boot <- function(obj,boot.sample = 1000, alpha = 0.05,tol = 10^-4){
  CL <- obj$left_censor
  CR <- obj$right_censor
  TL <- obj$left_trunc
  TR <- obj$right_trunc
  
  Index <- sort(unique(c(CL,CR,TL,TR)))
  if(!obj$include0){
    Index <- Index[Index>0]
  }
  I_n <- length(Index) 
  
  rev_sort <- function(x){
    x <- sort(x,decreasing = T)
    return(x)
  }
  
  quant <- ceiling(boot.sample*alpha)
  
  Upper.quant <-  matrix(0,nrow = quant +1, ncol = I_n)
  Lower.quant <-  matrix(Inf,nrow = quant +1, ncol = I_n)
  
  new_upper <- new_lower <- rep(NA,I_n)
  names(new_upper) <- names(new_lower) <- Index
  
  pb = txtProgressBar(min = 0, max = boot.sample, initial = 1, style = 3) 
  
  for(i in 1:boot.sample){
    setTxtProgressBar(pb,i)
    new_upper[] <- new_lower[] <- rep(NA,I_n)
    res <- sample(length(CL),length(CL),replace = T)
    
    fit_temp <- turnbull(CL[res],CR[res],TL[res],TR[res],tol = tol,message = F, start = obj$start, include0 = obj$include0)
    
    s <- fit_temp$s
    
    if(length(s) == I_n){
      Upper.quant[1,] <- Lower.quant[1,] <- s
      
    }else{
      name_ind <- rbind(matrix(as.numeric(unlist(strsplit(names(s),","))),ncol = length(s)),1:length(s))
      name_ind_eq <- as.matrix(name_ind[,which(name_ind[1,]==name_ind[2,])])
      name_ind_uneq <- as.matrix(name_ind[,which(name_ind[1,]!=name_ind[2,])])
      for(j in 1:I_n){
        value <- match(Index[j],name_ind_eq[1,]) 
        if(!is.na(value)){
          foo <- which(names(new_lower) == Index[j])
          new_lower[foo] <- new_upper[foo] <- s[value]
        }
      }
      if(dim(name_ind_uneq)[2] > 0){
        for(j in 1:ncol(name_ind_uneq)){
          value <- match(name_ind_uneq[1,j],Index)
          new_upper[which(names(new_upper) == name_ind_uneq[2,j])] <- s[name_ind_uneq[3,j]]
          value <- match(name_ind_uneq[2,j],Index)
          new_lower[which(names(new_lower) == name_ind_uneq[1,j])] <- s[name_ind_uneq[3,j]]
        }
      } 
      if(any(is.na(new_lower))){
        miss <- is.na(new_lower)
        new_lower[miss] <- 0
      }
      if(any(is.na(new_upper))){
        miss <- is.na(new_upper)
        new_upper[miss] <- 0
      }
    }
    Upper.quant[1,] <- 1-cumsum(new_upper[1:(I_n)])
    Upper.quant <- apply(Upper.quant,2,sort)
    
    Lower.quant[1,] <- 1-cumsum(new_lower[1:(I_n)])
    Lower.quant <- apply(Lower.quant,2,rev_sort)
  } 
  close(pb)
  CI.upper <- Upper.quant[2,]
  names(CI.upper) <- names(new_upper)
  CI.lower <- Lower.quant[2,]
  names(CI.lower) <- names(new_lower)
  return(list(CI.upper = CI.upper, CI.lower = CI.lower))
}

S <- function(obj,time){
  if(time > obj$p[length(obj$p)]){
    return(0)
  }else if(time < obj$p[1]){
    return(1)
  }else{
    x <- which.min(time-obj$p>0)-1
    return(1-unname(sum(obj$s[1:x])))
  }
}

S <- Vectorize(S,"time")


LL_AFM_std <- function(par,dist){
  #likelihood for models where log(X)=beta*Z+sigma*W, for W ~ f(.).
  #just need survival function for W to add new distribution to model
  
  if(tolower(dist) %in% c("lognorm","ln","lognormal")){ 
    # standard normal for lognormal
    S <- function(x){1-pnorm(x)}
    
  }else if(tolower(dist) %in% c("loglogistic","loglog","ll")){
    S <- function(x){1/(1+exp(x))}
  }else if(tolower(dist)%in% c("ep","exponentialpower","exp-pow")){
    S <- function(x){exp(1-exp(exp(x)))}
  }else{
    
    # standard extreme for Weibull
    S <- function(x){exp(-exp(x))}
  }
  
  #unbound <- function(x){ ifelse(x>1,x,exp(x-1))}
  
  n_par <- length(par)
  
  inner <- function(index,entry){
    X <- fit_data[index,entry]
    Y <- fit_data[index,5:ncol(fit_data)]
    return(list( y = (log(X) - as.numeric(par[1:(n_par-1)]%*%t(Y)))/par[n_par], X = Y))
  }
  
  S_mat <- matrix(0,nrow = nrow(fit_data),ncol = 4)
  for(i in 1:4){
    foo <- inner(1:nrow(fit_data),i)[[1]]
    S_mat[,i] <- ifelse(foo == -Inf,1,ifelse(foo == Inf,0,S(foo)))
    
  }
  
  L <- (S_mat[,1]-S_mat[,2])/(S_mat[,3]-S_mat[,4])
  
  LogL <- 0-sum(log(L))
  
  LogL <- ifelse(is.na(LogL) |LogL == Inf, 10^10,LogL )
  
  return(LogL)
}

LL_AFM_std_grad <- function(par,dist){
  #likelihood for models where log(X)=beta*Z+sigma*W, for W ~ f(.).
  #just need survival function for W to add new distribution to model
  
  n_par <- length(par)
  
  # if(tolower(dist) %in% c("lognorm","ln","lognormal")){ 
  # # standard normal for lognormal
  # S <- function(x){1-pnorm(x)}
  # dS <- function(x){-dnorm(x)}
  # }else
  
  if(tolower(dist) %in% c("loglogistic","loglog","ll")){
    S <- function(x){1/(1+exp(x))}
    dH <- function(x){1/(1+exp(-x))}
  }else if(tolower(dist)%in% c("ln","logn","lognormal")){
    S <- function(x){1-pnorm(x)}
    dH <- function(x){dnorm(x)/(S(x))}
  }else if(tolower(dist)%in% c("ep","exponentialpower","exp-pow")){
    S <- function(x){exp(1-exp(exp(x)))}
    dH <- function(x){exp(x)*exp(exp(x))}
  }else{
    # standard extreme for Weibull
    S <- function(x){exp(-exp(x))}
    dH <- function(x){exp(x)}
  }
  
  Y <- fit_data[,5:ncol(fit_data)]
  eta <- exp(as.numeric(par[1:(n_par-1)]%*%t(Y)))
  sigma <- par[n_par]
  
  S_mat <- matrix(0,nrow = nrow(fit_data),ncol = 4)
  for(i in 1:4){
    foo <- fit_data[,i]
    S_mat[,i] <- ifelse(foo == 0,1,ifelse(foo == Inf,0,S(log(foo/eta)/sigma)))
  }
  
  d_H_d_eta_mat <- matrix(0,nrow(fit_data), ncol = 4)
  for(i in 1:4){
    foo <- fit_data[,i]
    d_H_d_eta_mat[,i] <- ifelse(foo == 0,0,ifelse(foo == Inf,0,dH(log(foo/eta)/sigma)))
  }
  
  # vector for {(d/d(eta)H1)S1 - (d/d(eta)H2)S2}/{S1-S2} både censurering og trunkering
  vec <-  ifelse(S_mat[,1] ==0, (-d_H_d_eta_mat[,1]),
                 ((-d_H_d_eta_mat[,1])*S_mat[,1]-(-d_H_d_eta_mat[,2])*S_mat[,2])/(S_mat[,1]-S_mat[,2]))
  vec <- vec - ifelse(S_mat[,3]==0,(-d_H_d_eta_mat[,3]),
                      ((-d_H_d_eta_mat[,3])*S_mat[,3]-(-d_H_d_eta_mat[,4])*S_mat[,4])/(S_mat[,3]-S_mat[,4]))
  
  grad <- numeric(n_par)
  
  grad[1:(n_par-1)] <- colSums(outer(vec, rep(1,n_par-1))*(Y)/sigma)
  
  if(any(!is.finite(grad))){
    index <- which(!is.finite(grad))
    for(i in index){
      grad[i] <- sign(par[i])*10^10
    }
  }
  
  d_H_d_sigma_mat <- matrix(0,nrow(fit_data), ncol = 4)
  for(i in 1:4){
    foo <- fit_data[,i]
    d_H_d_sigma_mat[,i] <- ifelse(foo == 0,0,ifelse(foo == Inf,0,dH(log(foo/eta)/sigma)*log(foo/eta)))
  }
  
  vec <-  ifelse(S_mat[,1]-S_mat[,2] ==0, 0,
                 ((-d_H_d_sigma_mat[,1])*S_mat[,1]-(-d_H_d_sigma_mat[,2])*S_mat[,2])/(S_mat[,1]-S_mat[,2]))
  vec <- vec - ifelse(S_mat[,3]-S_mat[,4]==0,0,
                      ((-d_H_d_sigma_mat[,3])*S_mat[,3]-(-d_H_d_sigma_mat[,4])*S_mat[,4])/(S_mat[,3]-S_mat[,4]))
  
  grad[n_par] <- sum(vec/sigma^2)
  
  if(!is.finite(grad[n_par])){
    grad[n_par] <- ifelse(par[n_par] < 1,-10^5,10^5)
  }
  
  return(grad)
}

LL_AFM_gen <- function(par){
  #likelihood for models where X=exp(beta*Z)W, for W ~ f(.; sigma).
  #need surival function for W as a function of x and sigma
  
  n_par <- length(par)
  
  #if(tolower(dist) %in% c("g","gom")){ 
  S <- function(x){exp((1-exp(x))/par[n_par])}
  
  
  
  inner <- function(index,entry){
    X <- fit_data[index,entry]
    Y <- fit_data[index,5:ncol(fit_data)]
    return(list( y = exp(-as.numeric(par[1:(n_par-1)]%*%t(Y)))*X, X = Y))
  }
  
  S_mat <- matrix(0,nrow = nrow(fit_data),ncol = 4)
  for(i in 1:4){
    foo <- inner(1:nrow(fit_data),i)[[1]]
    S_mat[,i] <- ifelse(foo == 0,1,ifelse(foo == Inf,0,S(foo)))
  }
  
  L <- (S_mat[,1]-S_mat[,2])/(S_mat[,3]-S_mat[,4])
  
  LogL <- 0-sum(log(L))
  
  LogL <- ifelse(is.na(LogL) |LogL == Inf, 10^10,LogL )
  
  return(LogL)
}

LL_AFM_gen_grad <- function(par){
  #likelihood for proportional hazard models
  #fitted as = exp(beta*Z)*H(u|sigma,Z)), where H is cumulative hazard function
  n_par <- length(par)
  #if(tolower(dist) %in% c("g","gom")){ 
  
  Y <- fit_data[,5:ncol(fit_data)]
  eta <- exp(as.numeric(par[1:(n_par-1)]%*%t(Y)))
  sigma <- par[n_par]
  
  S <- function(x){exp((1-exp(x))/sigma)}
  
  S_mat <- matrix(0,nrow = nrow(fit_data),ncol = 4)
  for(i in 1:4){
    foo <- fit_data[,i]
    S_mat[,i] <- ifelse(foo == 0,1,ifelse(foo == Inf,0,S(foo/eta)))
    
  }
  
  d_H_d_eta <- function(x){x*exp(x/eta)/(eta*sigma)}
  
  d_H_d_sigma <- function(x){(exp(x/eta)-1)/sigma^2}
  
  d_H_d_eta_mat <- matrix(0,nrow(fit_data), ncol = 4)
  for(i in 1:4){
    foo <- fit_data[,i]
    d_H_d_eta_mat[,i] <- ifelse(foo == 0,0,ifelse(foo == Inf,0,d_H_d_eta(foo)))
  }
  
  # vector for {(d/d(eta)H1)S1 - (d/d(eta)H2)S2}/{S1-S2} både censurering og trunkering
  vec <-  ifelse(S_mat[,1] ==0, (-d_H_d_eta_mat[,1]),
                 ((-d_H_d_eta_mat[,1])*S_mat[,1]-(-d_H_d_eta_mat[,2])*S_mat[,2])/(S_mat[,1]-S_mat[,2]))
  vec <- vec - ifelse(S_mat[,3]==0,(-d_H_d_eta_mat[,3]),
                      ((-d_H_d_eta_mat[,3])*S_mat[,3]-(-d_H_d_eta_mat[,4])*S_mat[,4])/(S_mat[,3]-S_mat[,4]))
  
  grad <- numeric(n_par)
  
  grad[1:(n_par-1)] <- colSums(outer(vec, rep(1,n_par-1))*(Y))
  
  if(any(!is.finite(grad))){
    index <- which(!is.finite(grad))
    for(i in index){
      grad[i] <- sign(par[i])*10^10
    }
  }
  
  d_H_d_sigma_mat <- matrix(0,nrow(fit_data), ncol = 4)
  for(i in 1:4){
    foo <- fit_data[,i]
    d_H_d_sigma_mat[,i] <- ifelse(foo == 0,0,ifelse(foo == Inf,0,d_H_d_sigma(foo)))
  }
  
  vec <-  ifelse(S_mat[,1]-S_mat[,2] ==0, 0,
                 ((-d_H_d_sigma_mat[,1])*S_mat[,1]-(-d_H_d_sigma_mat[,2])*S_mat[,2])/(S_mat[,1]-S_mat[,2]))
  vec <- vec - ifelse(S_mat[,3]-S_mat[,4]==0,0,
                      ((-d_H_d_sigma_mat[,3])*S_mat[,3]-(-d_H_d_sigma_mat[,4])*S_mat[,4])/(S_mat[,3]-S_mat[,4]))
  
  grad[n_par] <- sum(vec)
  
  if(!is.finite(grad[n_par])){
    grad[n_par] <- ifelse(par[n_par] < 1,-10^5,10^5)
  }
  
  return(grad)
}

LL_PH_gen <- function(par){
  #likelihood for proportional hazard models
  #fitted as = exp(beta*Z)*H(u|sigma,Z)), where H is cumulative hazard function
  n_par <- length(par)
  #if(tolower(dist) %in% c("g","gom")){ 
  H <- function(x){(exp(x/par[n_par])-1)}
  
  Y <- fit_data[,5:ncol(fit_data)]
  eta <- exp(-as.numeric(par[1:(n_par-1)]%*%t(Y)))
  
  S_mat <- matrix(0,nrow = nrow(fit_data),ncol = 4)
  for(i in 1:4){
    foo <- fit_data[,i]
    S_mat[,i] <- ifelse(foo == 0,1,ifelse(foo == Inf,0,exp(-eta*H(foo))))
    
  }
  
  L <- (S_mat[,1]-S_mat[,2])/(S_mat[,3]-S_mat[,4])
  
  LogL <- 0-sum(log(L))
  
  LogL <- ifelse(is.na(LogL) |LogL == Inf, 10^10,LogL )
  
  return(LogL)
}

LL_PH_gen_grad <- function(par){
  #likelihood for proportional hazard models
  #fitted as = exp(beta*Z)*H(u|sigma,Z)), where H is cumulative hazard function
  n_par <- length(par)
  #if(tolower(dist) %in% c("g","gom")){ 
  
  Y <- fit_data[,5:ncol(fit_data)]
  eta <- exp(as.numeric(par[1:(n_par-1)]%*%t(Y)))
  sigma <- par[n_par]
  
  H <- function(x){(exp(x/sigma)-1)}
  
  S_mat <- matrix(0,nrow = nrow(fit_data),ncol = 4)
  for(i in 1:4){
    foo <- fit_data[,i]
    S_mat[,i] <- ifelse(foo == 0,1,ifelse(foo == Inf,0,exp(-1/eta*H(foo))))
    
  }
  
  d_H_d_eta <- function(x){(1-exp(x/sigma))/eta}
  
  d_H_d_sigma <- function(x){-x*exp(x/sigma)/(sigma^2*eta)}
  
  d_H_d_eta_mat <- matrix(0,nrow(fit_data), ncol = 4)
  for(i in 1:4){
    foo <- fit_data[,i]
    d_H_d_eta_mat[,i] <- ifelse(foo == 0,0,ifelse(foo == Inf,0,d_H_d_eta(foo)))
  }
  
  # vector for {(d/d(eta)H1)S1 - (d/d(eta)H2)S2}/{S1-S2} både censurering og trunkering
  vec <-  ifelse(S_mat[,1] ==0, (-d_H_d_eta_mat[,1]),
                 ((-d_H_d_eta_mat[,1])*S_mat[,1]-(-d_H_d_eta_mat[,2])*S_mat[,2])/(S_mat[,1]-S_mat[,2]))
  vec <- vec - ifelse(S_mat[,3]==0,(-d_H_d_eta_mat[,3]),
                      ((-d_H_d_eta_mat[,3])*S_mat[,3]-(-d_H_d_eta_mat[,4])*S_mat[,4])/(S_mat[,3]-S_mat[,4]))
  
  grad <- numeric(n_par)
  
  grad[1:(n_par-1)] <- colSums(outer(vec, rep(1,n_par-1))*(Y))
  
  if(any(!is.finite(grad))){
    index <- which(!is.finite(grad))
    for(i in index){
      grad[i] <- sign(par[i])*10^10
    }
  }
  
  d_H_d_sigma_mat <- matrix(0,nrow(fit_data), ncol = 4)
  for(i in 1:4){
    foo <- fit_data[,i]
    d_H_d_sigma_mat[,i] <- ifelse(foo == 0,0,ifelse(foo == Inf,0,d_H_d_sigma(foo)))
  }
  
  vec <-  ifelse(S_mat[,1]-S_mat[,2] ==0, 0,
                 ((-d_H_d_sigma_mat[,1])*S_mat[,1]-(-d_H_d_sigma_mat[,2])*S_mat[,2])/(S_mat[,1]-S_mat[,2]))
  vec <- vec - ifelse(S_mat[,3]-S_mat[,4]==0,0,
                      ((-d_H_d_sigma_mat[,3])*S_mat[,3]-(-d_H_d_sigma_mat[,4])*S_mat[,4])/(S_mat[,3]-S_mat[,4]))
  
  grad[n_par] <- sum(vec)
  
  if(!is.finite(grad[n_par])){
    grad[n_par] <- ifelse(par[n_par] < 1,-10^5,10^5)
  }
  
  return(0-grad)
}


AFM_fit <- function(data,dist,hessian = FALSE, start = NULL){
  n_par <- ncol(data)-3
  
  if(dist %in% c("gompertzph","gph","gomph")){
    LL_foo <- function(x){LL_PH_gen(x)}
    gr <- function(x){LL_PH_gen_grad(x)}
    scale <- 100
  }else if(dist %in% c("g","gom","gompertzafm","gafm","gomafm")){
    LL_foo <- function(x){LL_AFM_gen(x)}
    gr <- function(x){LL_AFM_gen_grad(x)}
  }else{
    LL_foo <- function(x){LL_AFM_std(x,dist)}
    gr <- function(x){LL_AFM_std_grad(x,dist)}
  }
  
  choose_start <-function(){
    best_LL<- 10^10
    best_stard <- numeric(n_par)
    for(i in c(1,(1:10)*10)){
      for(j in 1:10){
        start <- c(j,rep(0,n_par-2),i)
        new_LL <-LL_foo(start)
        if(new_LL < best_LL){
          best_LL <- new_LL
          best_start <- start
        }
        start <- c(rep(j,n_par-1),i)
        new_LL <-LL_foo(start)
        if(new_LL < best_LL){
          best_LL <- new_LL
          best_start <- start
        }
      }
    }
    return(best_start)
  }
  if(is.null(start)){
  start <- choose_start()
  }else if(length(start) == n_par & LL_foo(start) < 10^10){
    print(LL_foo(start))
  }else(stop())
  
  obj <- optim(start, LL_foo,method = "L-BFGS-B", lower = c(rep(-Inf,n_par-1),0.001),
               gr = gr, hessian = hessian,control = list(maxit = 1e4, pgtol = 0, ndeps = rep(1e-4, n_par)))
  
  names(obj$par) <- c(colnames(data)[5:(3+n_par)],"scale")
  
  if(hessian){
    colnames(obj$hessian) <- rownames(obj$hessian) <- names(obj$par)
  }
  
  return(obj)
}


ssurv <- function(x,eta,sigma,dist = "weib"){
  if(dist %in% c("gompertzph","gph","gomph")){
    return(exp((1-exp(x/sigma))/eta))
  }else if(dist %in% c("g","gom","gompertzafm","gafm","gomafm")){
    return(exp((1-exp(x/eta))/sigma))
  }else if(tolower(dist) %in% c("lognorm","ln","lognormal")){ 
    # standard normal for lognormal
    S <- function(x){1-pnorm(x)}
  }else if(tolower(dist) %in% c("loglogistic","loglog","ll")){
    S <- function(x){1/(1+exp(x))}
  }else if(tolower(dist)%in% c("ep","exponentialpower","exp-pow")){
    S <- function(x){exp(1-exp(exp(x)))}
  }else{
    # standard extreme for Weibull
    S <- function(x){exp(-exp(x))}
  }
  
  return(S((log(x)-log(eta))/sigma))
}

dsurv <- function(x,eta,sigma,dist = "weib"){
  if(tolower(dist) %in% c("lognorm","ln","lognormal")){ 
    # standard normal for lognormal
    dS <- function(x){-dnorm(x)}
  }else if(tolower(dist) %in% c("gafm","gomafm")){
    return(exp(x/eta)*exp((1-exp(x/eta))/sigma)/(eta*sigma))
  }else if(tolower(dist) %in% c("loglogistic","loglog","ll")){
      dS <- function(x){-exp(x)/(1+exp(x))^2}
  }else if(tolower(dist)%in% c("ep","exponentialpower","exp-pow")){
    dS <- function(x){-exp(x)*exp(exp(x))*exp(1-exp(exp(x)))}
  }else{
    # standard extreme for Weibull
    dS <- function(x){-exp(x-exp(x))}
  }
  
  return(-dS((log(x)-log(eta))/sigma)/(x*sigma))
}

dsurv <- Vectorize(dsurv,c("x","eta","sigma"))

hsurv <- function(x,eta,sigma,dist = "weib"){
  if(dist %in% c("gompertzph","gph","gomph")){
    return(exp(x/sigma)/eta)
  }else if(dist %in% c("g","gom","gompertzafm","gafm","gomafm")){
    return(exp(x/eta)/sigma/eta)
  }else if(tolower(dist) %in% c("lognorm","ln","lognormal")){ 
    # standard normal for lognormal
    S <- function(x){1-pnorm(x)}
    dS <- function(x){-dnorm(x)}
  }else if(tolower(dist) %in% c("loglogistic","loglog","ll")){
    S <- function(x){1/(1+exp(x))}
    dS <- function(x){-exp(x)/(1+exp(x))^2}
    
  }else if(tolower(dist)%in% c("ep","exponentialpower","exp-pow")){
    S <- function(x){exp(1-exp(exp(x)))}
    dS <- function(x){-exp(x)*exp(exp(x))*exp(1-exp(exp(x)))}  
  }else{
    # standard extreme for Weibull
    S <- function(x){exp(-exp(x))} 
    dS <- function(x){-exp(x-exp(x))}
  }
  
  return(-dS((log(x)-log(eta))/sigma)/(x*sigma*S((log(x)-log(eta))/sigma)))
}

qsurv <- function(x,eta,sigma,dist = "weib"){
  if(!(0<x & x < 1)){stop()}
  if(dist %in% c("gompertzph","gph","gomph")){
    return(log(-log(1-x)*eta+1)*sigma)
  }else if(dist %in% c("g","gom","gompertzafm","gafm","gomafm")){
    return(log(-log(1-x)*sigma+1)*eta)
  }else if(tolower(dist)%in% c("ep","exponentialpower","exp-pow")){
    return(exp(log(log(1-log(1-x))*eta)*sigma))
  }else{
    return(exp(log(-log(1-x)*eta)*sigma))
  }
}

qsurv <- Vectorize(qsurv, vectorize.args = c("x","eta","sigma"))

rsurv <- function(n, eta,sigma,dist = "weib"){
  return(qsurv(runif(n),eta,sigma,dist))
}

rsurv <- Vectorize(rsurv, vectorize.args = c("n","eta","sigma"))


AFM_fit_misc <- function(data,dist = "weib",hessian = FALSE){
  LL_AFM_std_misc <- function(par,dist){
    #likelihood for models where log(X)=beta*Z+sigma*W, for W ~ f(.).
    #just need survival function for W to add new distribution to model
    
    if(tolower(dist) %in% c("lognorm","ln","lognormal")){ 
      # standard normal for lognormal
      S <- function(x){1-pnorm(x)}
      
    }else if(tolower(dist) %in% c("loglogistic","loglog","ll")){
      S <- function(x){1/(1+exp(x))}
    }else if(tolower(dist)%in% c("ep","exponentialpower","exp-pow")){
      S <- function(x){exp(1-exp(exp(x)))}
    }else{
      
      # standard extreme for Weibull
      S <- function(x){exp(-exp(x))}
    }
    
    #unbound <- function(x){ ifelse(x>1,x,exp(x-1))}
    
    n_par <- length(par)
    
    inner <- function(index,entry){
      X <- fit_data[index,entry]
      Y <- fit_data[index,9:ncol(fit_data)]
      return(list( y = (log(X) - as.numeric(par[1:(n_par-1)]%*%t(Y)))/par[n_par], X = Y))
    }
    
    S_mat <- matrix(0,nrow = nrow(fit_data),ncol = 8)
    for(i in 1:8){
      foo <- inner(1:nrow(fit_data),i)[[1]]
      S_mat[,i] <- ifelse(foo == -Inf,1,ifelse(foo == Inf,0,S(foo)))
      
    }
    
    L <- ifelse(fit_data$CR == fit_data$CR_B, (S_mat[,1]-S_mat[,2]), alpha*S_mat[,5]+(beta-alpha)*S_mat[,6])
    
    LogL <- 0-sum(log(L))
    
    LogL <- ifelse(is.na(LogL) |LogL == Inf, 10^10,LogL )
    
    return(LogL)
  }
  
  LL_AFM_gen_misc <- function(par){
    if(length(par) != (length(fit_data) - 7)){stop()}
    #likelihood for models where X=exp(beta*Z)W, for W ~ f(.; sigma).
    #need surival function for W as a function of x and sigma
    
    n_par <- length(par)
    
    
    #if(tolower(dist) %in% c("g","gom")){ 
    S <- function(x){exp((1-exp(x))/par[n_par])}
    
    inner <- function(index,entry){
      X <- fit_data[index,entry]
      Y <- fit_data[index,9:ncol(fit_data)]
      return(list( y = exp(-as.numeric(par[1:(n_par-1)]%*%t(Y)))*X, X = Y))
    }
    
    S_mat <- matrix(0,nrow(fit_data),ncol = 8)
    for(i in 1:8){
      foo <- inner(1:nrow(fit_data),i)[[1]]
      S_mat[,i] <- ifelse(foo == -Inf,1,ifelse(foo == Inf,0,S(foo)))
      
    }
    
    L <- ifelse(fit_data$CR == fit_data$CR_B, (S_mat[,1]-S_mat[,2]), alpha*S_mat[,5]+(beta-alpha)*S_mat[,6])
    
    LogL <- 0-sum(log(L))
    
    LogL <- ifelse(is.na(LogL) |LogL == Inf, 10^10,LogL )
    
    return(LogL)
  }
  
  LL_PH_gen_misc <- function(par){
    #likelihood for proportional hazard models
    #fitted as = exp(beta*Z)*H(u|sigma,Z)), where H is cumulative hazard function
    n_par <- length(par)
    #if(tolower(dist) %in% c("g","gom")){ 
    H <- function(x){par[n_par]*(exp(x/par[n_par])-1)}
    
    Y <- fit_data[,9:ncol(fit_data)]
    eta <- exp(-as.numeric(par[1:(n_par-1)]%*%t(Y)))
    
    S_mat <- matrix(0,nrow = nrow(fit_data),ncol = 8)
    for(i in 1:8){
      foo <- fit_data[,i]
      S_mat[,i] <- ifelse(foo == -Inf,1,ifelse(foo == Inf,0,exp(-eta*H(foo))))
      
    }
    
    L <- ifelse(fit_data$CR == fit_data$CR_B, (S_mat[,1]-S_mat[,2]), alpha*S_mat[,5]+(beta-alpha)*S_mat[,6])
    
    LogL <- 0-sum(log(L))
    
    LogL <- ifelse(is.na(LogL) |LogL == Inf, 10^10,LogL )
    
    return(LogL)
  }
  if(dist %in% c("gompertzph","gph","gomph")){
    LL_foo <- function(x){LL_PH_gen_misc(x)}
    scale <- 100
  }else if(dist %in% c("g","gom","gompertzafm","gafm","gomafm")){
    LL_foo <- function(x){LL_AFM_gen_misc(x)}
  }else{
    LL_foo <- function(x){LL_AFM_std_misc(x,dist)}
  }
  for(i in c(1,(1:10)*10)){
    start <- c(5,rep(0,n_par-2),i)
    if(LL_foo(start) < 5*10^9){
      break
    }
    start <- c(rep(5,n_par-1),i)
    if(LL_foo(start) < 5*10^9){
      break
    }
  }
  
  obj <- optim(start, LL_foo,method = "L-BFGS-B", lower = c(rep(-Inf,n_par-1),0),hessian = hessian)
  
  names(obj$par) <- c(colnames(data)[9:(7+n_par)],"scale")
  
  if(hessian){
    colnames(obj$hessian) <- rownames(obj$hessian) <- names(obj$par)
  }
  
  return(obj)
}

BrierRes <- function(data,par, dist = "gafm", t= 70){
  n_par <- ncol(data)-3
  if(length(par) != n_par ){stop()}
  
  beta <- par[1:(length(par)-1)]
  sigma <- par[length(par)]
  if(n_par == 2){
    eta <- exp(beta*data[,5])
  }else{
    eta <- exp(rowSums(outer(rep(1,nrow(data)),beta)*data[5:(ncol(data)-1)]))
  }
  S_CL <- ssurv(data[,1],eta,sigma, dist)
  S_CR <- ssurv(data[,2],eta,sigma, dist)
  S_t <-  ssurv(t,eta,sigma, dist)
  
  P_t <- (S_t- S_CR)/(S_CL-S_CR)
  
  P_t <- ifelse(is.nan(P_t),1,P_t)
  
  forecast <- ifelse(t >= data[,2], 0, ifelse(t < data[,1], 1, P_t))
  
  total <- (forecast-S_t)^2
  
  return(total)
}

BrierCurve <- function(data,par, dist = "gafm", t= 70){
  return(mean(BrierRes(data,par,dist,t)))
}

BrierScore <- function(data,par, dist = "gafm"){
  total <- numeric(100)
  for(i in 1:100){
    total[i] <- BrierCurve(data,par,dist,i)
  }
  return(sum(total))
}
