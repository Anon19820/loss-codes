## Functions for MNV Samples
if(FALSE){

gen_symb_mats <- function(dim_1 = NULL, sd_vec = NULL, cov = F, return_Ryacas = T ){
  
  if(!is.null(sd_vec)){
    dim_1 = length(sd_vec)
  }
  n_corr_param <- dim_1*(dim_1-1)/2
  R_symb <- base::diag(dim_1)
  R_symb[upper.tri(R_symb)] <- letters[1:n_corr_param]
  R_symb[lower.tri(R_symb)] <-   t(R_symb)[lower.tri(R_symb)] 
  
  if(return_Ryacas){
  R_symb <- ysym(R_symb)
  if(cov){
    sim_mat <- ysym(base::diag(sd_vec))
    sim_mat%**%R_symb%**%sim_mat
    }
  }
  R_symb
  
}



gen_R_mat <- function(corr_param, corr_row_index, corr_col_index, mu_length){ # type declarations
    R_mat <- base::diag(mu_length)
    #R_mat <- matrix(1, nrow = mu_length, ncol = mu_length) 
    num_corr <- mu_length*(mu_length-1)/2
    for(i in 1:num_corr){
      R_mat[corr_row_index[i],corr_col_index[i]] <- corr_param[i]
      R_mat[corr_col_index[i],corr_row_index[i]] <- corr_param[i]
    }
    return(R_mat)
}

gen_pd_corr_par<- function(p_concord,corr_row_index,corr_col_index, mu_length, corr_index_pd, alpha_concord, beta_concord, plot_roots = F){
  corr_param <- sin(2*pi*(p_concord/2 -0.25))
  
  R_mat<- gen_R_mat(corr_param,corr_row_index,corr_col_index, mu_length)
  R_mat <- ysym(R_mat)
  R_mat[corr_row_index[corr_index_pd],corr_col_index[corr_index_pd]] <- R_mat[corr_col_index[corr_index_pd],corr_row_index[corr_index_pd]] <- "r" 
  
  #can be any permutation
  dependent.ind <- 1
  given.ind <- c(2:mu_length)
  
  B <- R_mat[dependent.ind, dependent.ind, drop = F] #By definition - will be 1 
  C <- R_mat[dependent.ind, given.ind, drop = FALSE]
  D <- R_mat[given.ind, given.ind]
  
  CDinv <- C %**% solve(D)
  cVar <- B - CDinv %**% t(C)
  corr_eval <-   "r"
  
  fin_var <- simplify(cVar)
  func_eval <- as.function(eval(parse(text = paste0("alist(",corr_eval," = , eval(parse(text = fin_var$yacas_cmd)))"))))
  
  if(plot_roots){
  x_seq <- seq(-1,1, length.out =100)  
  plot(x_seq, func_eval(x_seq), ylim = c(0,1))
  }

  roots<- uniroot.all(func_eval, c(-.99,.99)) # sometimes there is a third root very near 1 won't worry about that
  #sometimes if roots are always positive then the length is zero.
  if(length(roots) ==0){
    roots <- c(-1,1)
  }
  if(length(roots) ==1){ #Need to see which direction PD is
    #roots <- roots[1]
    if(func_eval(roots*1.01) >= 0){ #If it is the lower interval then the other root is 1
      roots <- c(roots,.99)
    }else{
      roots <- c(-.99,roots)
    }
  }
  
  #Convert back to concordance prob; 
  prob_bounds <- (asin(roots)/pi) +0.5
  print(paste0("Index - ",corr_index_pd," Bounds =" ,round(prob_bounds, digits = 2)))
  # #prob_bounds <-(asin(c(-1,1))/pi) +0.5
  # F.ab <- pbeta(prob_bounds, alpha_concord[corr_index_pd], beta_concord[corr_index_pd])
  # u <- runif(1, min = min(F.ab), max = max(F.ab))
  # prob_sim <-  qbeta(u, alpha_concord[corr_index_pd], beta_concord[corr_index_pd])
  # r_corr <- sin(2*pi*(prob_sim/2 -0.25))
  # 
  # #R_mat <- with_value(R_mat, "r", r_corr)
  # #Also need the correlation vals
  # corr_param[corr_index_pd] <- r_corr
  # corr_param
  # 
}

#Yacas will coerce to a scalar so need to add the variable drop to stop this if required
`[.yac_symbol` <- function (x, i, j, drop = T){
  stopifnot(methods::is(x, "yac_symbol"))
  y_res <- yac_str(x$yacas_cmd)
  y <- ysym(y_res)
  stopifnot(y$is_mat | y$is_vec)
  if (y$is_vec && !missing(j)) {
    stop("Cannot specify second dimension for a vector")
  }
  w <- as_r(y$yacas_cmd)
  z <- NULL
  if (y$is_vec) {
    z <- base::`[`(x = w, i = i)
  }
  else if (y$is_mat) {
    if (missing(j)) {
      n_args <- nargs()
      if (n_args == 2L) {
        z <- base::`[`(x = w, i = i)
      }
      else if (n_args == 3L) {
        z <- base::`[`(x = w, i = i, )
      }
    }
    else {
      z <- base::`[`(x = w, i = i, j = j, drop = drop)
    }
  }
  stopifnot(!is.null(z))
  v <- ysym(z)
  return(v)
}



tmpfun <- get("[.yac_symbol", envir = asNamespace("Ryacas"))
environment(`[.yac_symbol`) <- environment(tmpfun)
attributes(`[.yac_symbol`) <- attributes(tmpfun)
assignInNamespace("[.yac_symbol", `[.yac_symbol`, ns="Ryacas") 

`%**%` <- Ryacas:::`%*%`
}

dinvgamma <-function (x, shape, rate = 1, scale = 1/rate, log = FALSE){
  if (missing(rate) && !missing(scale)) 
    rate <- 1/scale
  log_f <- dgamma(1/x, shape, rate, log = TRUE) - 2 * log(x)
  if (log) 
    return(log_f)
  exp(log_f)
}

qinvgamma <- function (p, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, 
                       log.p = FALSE){
  if (missing(rate) && !missing(scale)) 
    rate <- 1/scale
  qgamma(1 - p, shape, rate, lower.tail = lower.tail, log.p = log.p)^(-1)
}


pt.scaled <- function (q, df, mean = 0, sd = 1, ncp, lower.tail = TRUE, log.p = FALSE){
  stats::pt((q - mean)/sd, df, ncp = ncp, log.p = log.p)
}


rinvgamma <- function (n, shape, rate = 1, scale = 1/rate) 
{
  if (missing(rate) && !missing(scale)) 
    rate <- 1/scale
  1/rgamma(n, shape, rate)
}


qt.scaled <- function (p, df, mean = 0, sd = 1, ncp, lower.tail = TRUE, log.p = FALSE){
  mean + sd * stats::qt(p, df, ncp = ncp, log.p = log.p)
}

dt.scaled <-function (x, df, mean = 0, sd = 1, ncp, log = FALSE){
  if (!log) 
    stats::dt((x - mean)/sd, df, ncp = ncp, log = FALSE)/sd
  else stats::dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - 
    log(sd)
}
student_t_error <- function(par,vals , probs,sample_size){
  scale <- sample_size
  a <- sample_size/2
  mu <- par[1]
  b <- exp(par[2])
  sum((pt.scaled(vals,df = 2*a,mean = mu,sd = sqrt(b*(scale+1)/(scale*a))) - probs)^2)
}
student_t_error_1d <- function(par,vals , probs,sample_size, mu){
  scale <- sample_size
  a <- sample_size/2
  #mu <- par[1]
  b <- exp(par)
  sum((pt.scaled(vals,df = 2*a,mean = mu,sd = sqrt(b*(scale+1)/(scale*a))) - probs)^2)
}
norm_error <- function(par,vals , probs){
  mu <- par[1]
  sd <- exp(par[2])
  sum((pnorm(vals,mean = mu,sd = sd) - probs)^2)
}

hyper_norm_gamma <- function(vals, probs, sample_size, plot_marg = FALSE){
  
  #Need to solve the two quantile equations: i.e. 
  #val_i = mu + sd*qnorm(prob_i, mu = 0, sd = 1)
  #Write out matrix
  #conflict with Ryacas
  init_res <- base:::`%*%`(solve(as.matrix(cbind(1, qnorm(probs)))),  t(t(vals)))
  init_sd <- init_res[2]
  init_mu <- init_res[1]
  
  #Find the log(b) conditional on the mean 
  Brent_opt<- optim(c(0),student_t_error_1d, vals = vals,probs= probs,sample_size=sample_size, mu =init_mu, method= "Brent", lower = 0, upper= 100)
  
  #student_t_error(c(init_mu,Brent_opt$par),vals = vals,probs= probs,sample_size=sample_size)
  t_dist_opt <- optimx::optimx(c(init_mu,Brent_opt$par),student_t_error, vals = vals,probs= probs,sample_size=sample_size)
  if(min(t_dist_opt$value) > 0.001){
    warning("Optimization may not have been sucessfull")
  }
  final_pars <- t_dist_opt[which.min(t_dist_opt$value),1:2]
  mu_final <- final_pars[[1]]
  b_final <- exp(final_pars[[2]])
  
  #Plot
  mu <- mu_final
  a <- sample_size/2
  b <- b_final
  scale <- sample_size
  
  if(plot_marg){
    range_eval<- qt.scaled(c(0.001, .99999),df = 2*a,mu,sd = sqrt(b*(scale+1)/(scale*a) ))
    seq_x_eval <- seq(from = range_eval[1], to = range_eval[2], length.out = 1000)
    
    plot(seq_x_eval, dt.scaled(seq_x_eval,df = 2*a,mu,sd = sqrt(b/(scale*a))),
         col = "blue",type = "l", xlab = "x", ylab= "Density") # Marginal distribution of the mean
    lines(seq_x_eval, dt.scaled(seq_x_eval,df = 2*a,mu,sd = sqrt(b*(scale+1)/(scale*a)) ),col = "red") 
    legend("topright", legend=c("Marginal distibution of data", expression(paste("Marginal distibution of ", mu))),
           col=c("red", "blue"), lty=1, cex=0.8)
    
    range_var_eval<- qinvgamma(c(0.001, .999),a, b)
    range_sd_eval<- sqrt(range_var_eval)
    seq_var_x_eval <- seq(from = range_sd_eval[1], to = range_sd_eval[2], length.out = 1000)
    plot(seq_var_x_eval,dinvgamma((seq_var_x_eval)^2,a,b)*2*seq_var_x_eval,
         main = expression(paste("Marginal distribution of ",sigma)), xlab = "x", ylab= "Density", type ="l")
    
  }
  return(c(mu=mu,scale=scale,alpha =a,beta=b))
  
}

`%!in%` = Negate(`%in%`)


dens_func <- function(j,alpha, beta, t) { #Function which returns the density of survival at a particular timepoint
  dinvgamma(-t/log(j), shape= alpha, rate = beta)*(t/(j*(log(j)^2)))
}
get_y_val <- function(log_y, tau, gamma, alpha, time){
  y <- exp(log_y)
  prob_vec <- pgamma(q = -alpha*y*log(gamma)/time, alpha,1)
  return(abs(prob_vec-tau))
}
dlomax <- function(x, a, b, log = F){
  dens_eval_log <-  (log(a)+ a*log(b)) - (a +1)*log(x+b)
  
  if(log){
    res <- dens_eval_log
  }else{
    res <- exp(dens_eval_log)
  }
}
pow<- function (x1, x2){
  x1^x2
}

