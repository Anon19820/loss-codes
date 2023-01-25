
#--- Make Sure Relevant Packages are loaded ---#

list.of.packages <- need<-c("ggplot2",  "nimble", "Ryacas" ) #needed libraries

res <- lapply(list.of.packages, require, character.only = TRUE)

not.loaded <-   list.of.packages[which(sapply(res, unlist) ==F)]
not.installed <-   not.loaded
#load the packages
if(length(not.installed)) install.packages(not.installed)
if(length(not.installed)) lapply(not.installed, require, character.only = TRUE)

#path_files <- #Define working directory

source(paste0(path_files,"helper functions.R"))


#--- Example 1A: Exponential Likelihood Prior Predictive Approach ---#

#Figure 1 Plot of the tertiles for different sample sizes

samp_size <- c(5,10,25,100)
a_final<- b_final<- L_quant <- U_quant <- rep(NA,length(samp_size))
x_vec <- seq(0,30, by = .0001) 

quant_test <- 0.75
alpha_eval <-samp_size
beta_eval <- 10/(pow(0.5,-1/alpha_eval)-1)

#Each distribution has the same median survival of 10
for(i in 1:length(samp_size)){
  a_test <- alpha_eval[i]
  b_test <- beta_eval[i]
  # lomax_cdf <- 1- pow(b_test/(x_vec+b_test),a_test)
  # plot(x = x_vec, y = lomax_cdf , type = "l", xlim = c(0,30),ylim = c(0,1))
  
  U_quant[i] <- b_test*((1-0.75)^(-1/a_test)-1)
  L_quant[i] <- b_test*((1-0.25)^(-1/a_test)-1)
  
  FU_sim <- 0.75#runif(1,0.5,1)
  FL_sim <- 0.25#runif(1,0,0.5)
  
  # Use the brute force approach
  a_eval <- function(log_b, FU_sim,U_quant,FL_sim, L_quant){
    b <- exp(log_b)
    a_propose_1 <- log(1-FU_sim)/(log_b-log(U_quant+b))
    a_propose_2 <-  log(1-FL_sim)/(log_b-log(L_quant+b))
    return(abs(a_propose_1 - a_propose_2))
  }
  
  b_eval <- seq(0.01, 1000, by = 0.01)

  #Seem to have to brute force grid search it
  a_eval_res <- a_eval(log(b_eval), FU_sim, U_quant[i], FL_sim, L_quant[i])
  b_final[i] <- b_eval[which.min(a_eval_res)]
  
  a_eval2 <- function(b, FU_sim,U_quant,FL_sim, L_quant){
    log_b <- log(b)
    log(1-FU_sim)/(log_b-log(U_quant+b))
  }
  a_final[i] <- a_eval2(b_final[i], FU_sim,U_quant[i],FL_sim, L_quant[i])
}

png(paste0(path_files,"plots/Quantiles-Lomax.png"),width = 7, height = 7, units = 'in',res = 480)
plot(x_vec, dlomax(x_vec, a= alpha_eval[1], b = beta_eval[1]), type= "l", col = 1,
     ylab = "Density", xlab = "Time", 
     main = "Tertiles of Lomax distribution with equal median survival", ylim = c(0,0.1))
abline(v = c(L_quant[1],U_quant[1]))
for(i in 2:length(alpha_eval)){
  lines(x_vec, dlomax(x_vec, a= alpha_eval[i], b = beta_eval[i]), col = i)
  abline(v = c(L_quant[i],U_quant[i]),col = i)
}
legend("topright", 
       legend=paste0("n = ",samp_size ),
       col=1:length(alpha_eval),
       lty=1,
       cex=0.8)
dev.off()


#Verify by simulation
index <- 3
lambda <- rgamma(10000, alpha_eval[index], beta_eval[index])
sims_final <- sapply(lambda, rexp, n = 1000)
plot(density(sims_final))
lines(x_vec, dlomax(x_vec, a= alpha_eval[index], b = beta_eval[index]),col = "red")

#--- Example 1B: Exponential Likelihood with Data Augmentation Priors ---#


#Specify two timepoints and find the hyperparameters through optimization

alpha <- 10
beta <- alpha*1

t1 <- 1
gamma1 <- 0.5

t2 <- 2
gamma2 <- .1

j_vec_surv <- seq(0, .99, by = 0.01)


#Find tau_1 and tau_2 through integration
integral1 <- stats::integrate(dens_func, lower = gamma1, upper = 1,alpha = alpha, beta = beta, t = t1)
tau_1 <- integral1$value

integral2 <- stats::integrate(dens_func, lower = gamma2, upper = 1,alpha = alpha, beta = beta, t = t2)
tau_2 <- integral2$value


get_y_val2 <- function(par, tau, gamma, time){
  
  y <- exp(par[1])
  alpha <- exp(par[2])
  prob_vec <- rep(NA, length(tau))
  for(i in 1:length(prob_vec)){
    prob_vec[i] <- pgamma(q = -alpha*y*log(gamma[i])/time[i], alpha,1)
  }
  
  return(sum(abs(prob_vec-tau)))
}

tau_vec_2 <- c(tau_1, tau_2)
gamma_vec_2 <- c(gamma1,gamma2)
time_2 <- c(t1,t2)

optim_val <- optim(c(log(1),log(1)),
                   fn = get_y_val2,
                   tau = tau_vec_2,
                   gamma = gamma_vec_2,
                   time = time_2)

optim_par <- exp(optim_val$par)

#Figure 2: Probability of surviving past time equal 1
png(paste0(path_files,"plots/Example Bedrick.png"),width = 7, height = 7, units = 'in',res = 480)
plot(j_vec_surv, dens_func(j_vec_surv, alpha =optim_par[2] , beta = prod(optim_par), t = time_2[1]),
     type = "l", ylab = "Density", 
     #main =  paste0("Probability of Surviving past Time = ",time_2[1]),
     main = expression(paste("S(t=1) with ",
                             tau, " = 0.16 ,",gamma, " = 0.5")),
     xlab = "Probability of Surviving")
abline(v = gamma_vec_2[1], lty = 2)
dev.off()


#Define hyperparamters based on median

j_vec_median <- seq(0, 5, by = 0.01)

dens_func_median <- function(j,alpha, beta) {
  dinvgamma(j/log(2), shape= alpha, rate = beta)/log(2)
}

if(F){
  #Confirm change of variables calculation by simulation
  median_St <- log(2)*rinvgamma(1000000, shape= alpha, rate = beta)
  
  plot(density(median_St))
  lines(j_vec_median, 
        dens_func_median(j_vec_median, alpha =optim_par[2],
                         beta = prod(optim_par)),
        col = "red")
}

#Figure 3
png(paste0(path_files,"plots/Example Bedrick - Median.png"),width = 7, height = 7, units = 'in',res = 480)
plot(j_vec_median, 
     dens_func_median(j_vec_median, alpha =optim_par[2],
                      beta = prod(optim_par)),
     col = "black", type = "l", ylab = "Density", xlab = "Time", 
     main = "Density for Median Survival", xlim = c(0, 3))
dev.off()




#--- Example 1C: Exponential Likelihood with Loss Function BHW ---#

# We use NIMBLE to illustrate how this can be done in the BUGS language;
# It would be computationally more efficient to define "distribtuions" for each loss 
# functions in NIMBLE, however, the "zeros" trick works in both NIMBLE and JAGS.
# Other examples are in Stan.

#Because we have an explicit expression for the density that an IG (inverse gamma)
# distribution induces on median survival we can include this as a loss function.
# We can then see that they produce the same results.


# mod_median1: We can remove the impact of the uniform prior that lambda has on 
# median survival by specifying the corresponding loss function.
mod_median1<- nimbleCode({
  a <- 0.01
  b <- 10
  cons <- 1/(b-a)
  lambda ~ dunif(a, b)
  median_St <- log(2)/lambda
  LL_tot <- alpha*log(beta)-lgamma(alpha)+(-alpha-1)*log(median_St/log(2))+((-beta*log(2))/median_St)+log(1/log(2)) - log(cons*(lambda^2)/log(2))
  zero ~ dpois(-LL_tot+C)

})

# mod_median2: We can remove the impact of the uniform prior that lambda has on 
# median survival by specifying median_St as a parameter and lambda as a 
# deterministic function of median_St.
mod_median<- nimbleCode({
  median_St ~ dunif(0, 100)
  lambda <- log(2)/median_St
  mean_St <- 1/lambda
  LL_tot <- alpha*log(beta)-lgamma(alpha)+(-alpha-1)*log(median_St/log(2))+((-beta*log(2))/median_St)+log(1/log(2))
  zero ~ dpois(-LL_tot+C)
  
})

#mod_median3: If we parameterized the model in terms of mean survival (psi)  
# we don't need to include any adjustment as a uniform prior on psi produces a
#uniform density on median survival
mod_median3<- nimbleCode({

  a <- 0.01
  b <- 10
  
  psi ~ dunif(a, b) #Mean Survival
  median_St <- log(2)*psi
  LL_tot <- alpha*log(beta)-lgamma(alpha)+(-alpha-1)*log(median_St/log(2))+((-beta*log(2))/median_St)+log(1/log(2)) 
  zero ~ dpois(-LL_tot+C)
  
})

#mod_median4: In this example we assume that the expert's opinion has a 
#log normal distribution.

mod_median4<- nimbleCode({

  median_St ~ dunif(0, 100)
  lambda <- log(2)/median_St
  LL_tot <- log(dlnorm(median_St,mean_expert, sd_expert))
  zero ~ dpois(-LL_tot+C)
  
})


model_median_eval1 <- nimbleMCMC(mod_median1, 
                                data = list(zero = 0),
                                constants = list(alpha =optim_par[2] , beta = prod(optim_par), C = 1000),
                                monitors = c("median_St", "lambda", "LL_tot"),
                                #inits =inits_nimble, 
                                niter = 100000, nchains =2, nburnin =1000, thin  = 5)


mod_res <- rbind(model_median_eval1[[1]],model_median_eval1[[2]]) %>% data.frame()

png(paste0(path_files,"plots/Comparison of Priors.png"),width = 7, height = 7, units = 'in',res = 480)
plot(density(mod_res$median_St), xlim = c(0, 5),main = "Median Survival via two specifications of the prior", xlab = "Time")
lines(j_vec_median, dens_func_median(j_vec_median, alpha =optim_par[2] , beta = prod(optim_par)),col = "red")
legend("topright", 
       legend= c("Bissiri method", "DAP"),
       col=c("black", "red"),
       lty=1,
       cex=0.8)
dev.off()


# model_median_eval2 <- nimbleMCMC(mod_median2, 
#                                 data = list(zero = 0),
#                                 constants = list(alpha =optim_par[2] , beta = prod(optim_par), C = 1000),
#                                 monitors = c("median_St", "lambda", "LL_tot"),
#                                 #inits =inits_nimble, 
#                                 niter = 100000, nchains =2, nburnin =1000, thin  = 5)
# mod_res <- rbind(model_median_eval2[[1]],model_median_eval2[[2]]) %>% data.frame()


# model_median_eval3 <- nimbleMCMC(mod_median3,
#                                 data = list(zero = 0),
#                                 constants = list(alpha =optim_par[2] , beta = prod(optim_par), C = 1000),
#                                 monitors = c("median_St", "psi", "LL_tot"),
#                                 #inits =inits_nimble,
#                                 niter = 100000, nchains =2, nburnin =1000, thin  = 5)
# mod_res <- rbind(model_median_eval3[[1]],model_median_eval3[[2]]) %>% data.frame()



median_St <- log(2)*rinvgamma(1000000, shape= alpha, rate = beta)
mu_St <- mean(median_St)
var_St <- var(median_St)
#Find the lognormal parameters - Method of Moments
mu_ln <- log(mu_St/(sqrt(var_St/(mu_St^2) +1)))
sd_ln <- sqrt(log(var_St/(mu_St)^2 +1))

model_median_eval4 <- nimbleMCMC(mod_median4, 
                                data = list(zero = 0),
                                constants = list(mean_expert = mu_ln,
                                                 sd_expert = sd_ln, 
                                                 C = 1000),
                                monitors = c("median_St", "lambda", "LL_tot"),
                                niter = 100000, nchains =2, nburnin =5000, thin  = 5)


mod_res <- rbind(model_median_eval4[[1]],model_median_eval4[[2]]) %>% data.frame()

png(paste0(path_files,"plots/Comparison_of_Priors_LogNormal.png"), ,width = 7, height = 7, units = 'in',res = 480)
plot(density(mod_res$median_St), xlim = c(0, 5),main = "Median Survival using BHW Approach", xlab = "Time")
lines(j_vec_median, dlnorm(j_vec_median, meanlog  =mu_ln, sdlog  = sd_ln),col = "red")
legend("topright", 
       legend= c("BHW method", "True Density (Log-Normal)"),
       col=c("black", "red"),
       lty=1,
       cex=0.8)
dev.off()



##--Not used in the manuscript - Plot of the different sample sizes and 
# posterior uncertainty they imply


alpha_vec <- c(100,30,10,1)
y_vec <- rep(NA, length(alpha_vec))
tau_1 = 0.16
gamma_1 = 0.5
time_1 = 1#52

#Optimization can Fail! - Brute force grid search required
for(i in 1:length(alpha_vec)){
  seq_y <- seq(0,1000, by = .01)
  deviation_y <- get_y_val(log(seq_y), tau = tau_1, gamma = gamma_1,alpha = alpha_vec[i], time = time_1)
  print(deviation_y[which.min(deviation_y)])
  y_vec[i] <- seq_y[which.min(deviation_y)]
}


#For this plot we know the values of alpha and y so we use them rather than the optimized values

png(paste0(path_files,"plots/Bedrick-DAP-multiple_n.png"),width = 7, height = 7, units = 'in',res = 480)
plot(x =j_vec_surv,
     dens_func(j_vec_surv,alpha = alpha_vec[1], beta = y_vec[1]*alpha_vec[1], t = time_1 ), col = 1, type = "l",
     xlab = "Survival S(t)", ylab = "Density", main = expression(paste("S(t=1) for different samples sizes ",
                                                                       tau, " = 0.16 ,",gamma, " = 0.5")))
for(i in 2:length(alpha_vec)){
  lines(x =j_vec_surv,
        dens_func(j_vec_surv,alpha = alpha_vec[i], beta = y_vec[i]*alpha_vec[i], t = time_1 ), col = i)
}
abline(v = gamma_1, lty = "dashed")
legend("topright", 
       legend=paste0("n = ",alpha_vec ),
       col=1:length(alpha_vec),
       lty=1,
       cex=0.8)
dev.off()




