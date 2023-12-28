
#--- Make Sure Relevant Packages are loaded ---#

list.of.packages <- need<-c("ggplot2",  "nimble", "runjags", "posterior", "ggmcmc","coda", "latex2exp", "SHELF", "optimx") #needed libraries

res <- lapply(list.of.packages, require, character.only = TRUE)

not.loaded <-   list.of.packages[which(sapply(res, unlist) ==F)]
not.installed <-   not.loaded
#load the packages
if(length(not.installed)) install.packages(not.installed)
if(length(not.installed)) lapply(not.installed, require, character.only = TRUE)

#path_files <- #Define working directory
path_files <- "~/Expert Opinion - General/"
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


#Verify density of lomax by simulation
index <- 3
lambda <- rgamma(10000, alpha_eval[index], beta_eval[index])
sims_final <- sapply(lambda, rexp, n = 1000)
plot(density(sims_final))
lines(x_vec, dlomax(x_vec, a= alpha_eval[index], b = beta_eval[index]),col = "red")

#--- Example 1B: Exponential Likelihood with Data Augmentation Priors ---#


#Define hyperparamters based on beliefs about median survival
#Optimizing both alpha and beta is difficult so we just set alpha to 5
# We then pick a quantile for the 50% percentile (i.e. median)

alpha <- 5
j_vec_median <- seq(0, 50, by = 0.01)

dens_func_median <- function(j,alpha, beta) {
  dinvgamma(j/log(2), shape= alpha, rate = beta)/log(2)
}


eval_cdf <- function(par, alpha, prob =0.5, quantile){
  par2 <- exp(par)
  beta <- par2[1]
  cdf_1<- integrate(dens_func_median, 0, quantile, alpha = alpha, beta = beta)$value
 
  return((cdf_1-prob)^2)
}

res_final <- optim(par =c(0), eval_cdf, alpha = alpha, quantile = 10)
beta <- exp(res_final$par)

if(F){
  #Confirm change of variables calculation by simulation
  median_St <- log(2)*rinvgamma(1000000, shape= alpha, rate = beta)
  
  plot(density(median_St))
  lines(j_vec_median, 
        dens_func_median(j_vec_median, alpha =alpha,
                         beta = beta),
        col = "red")
}

#Figure 2
png(paste0(path_files,"plots/Example Bedrick - Median.png"),width = 7, height = 7, units = 'in',res = 480)
plot(j_vec_median, 
     dens_func_median(j_vec_median, alpha =alpha,
                      beta = beta),
     col = "black", type = "l", ylab = "Density", xlab = "Time", 
     main = "Density for Median Survival", xlim = c(0, 60))
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
# median survival by specifying the corresponding loss function (which is based on an inverse gamma distribution).
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

  lambda ~ dunif(0, 10)
  median_St <- log(2)/lambda
  LL_tot <- log(dlnorm(median_St,mean_expert, sd_expert)) -log(lambda^2)*equals(1,jacob_adj) 
  
  zero ~ dpois(-LL_tot+C)
  
})


model_median_eval1 <- nimbleMCMC(mod_median1, 
                                 data = list(zero = 0),
                                 constants = list(alpha =alpha , beta = beta, C = 1000),
                                 monitors = c("median_St", "lambda", "LL_tot"),
                                 #inits =inits_nimble, 
                                 niter = 100000*5+10000, nchains =2, nburnin =10000, thin  = 5)

mod_res <- rbind(model_median_eval1[[1]],model_median_eval1[[2]]) %>% data.frame()

# model_median_eval2 <- nimbleMCMC(mod_median2, 
#                                 data = list(zero = 0),
#                                 constants = list(alpha = alpha , beta = beta, C = 1000),
#                                 monitors = c("median_St", "lambda", "LL_tot"),
#                                 #inits =inits_nimble, 
#                                 niter = 100000, nchains =2, nburnin =1000, thin  = 5)
# mod_res <- rbind(model_median_eval2[[1]],model_median_eval2[[2]]) %>% data.frame()


# model_median_eval3 <- nimbleMCMC(mod_median3,
#                                 data = list(zero = 0),
#                                 constants = list(alpha = alpha, beta = beta, C = 1000),
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
                                data = list(zero = 0,jacob_adj = 1),
                                constants = list(mean_expert = mu_ln,
                                                 sd_expert = sd_ln, 
                                                 C = 1000),
                                monitors = c("median_St", "lambda", "LL_tot"),
                                niter = 100000*5+10000, nchains =2, nburnin =10000, thin  = 5)

mcmc.list_eval<- lapply(model_median_eval4, as.mcmc)
ggmcmc(ggs(as.mcmc.list(mcmc.list_eval)), "output.pdf")


mod_res <- rbind(model_median_eval4[[1]],model_median_eval4[[2]]) %>% data.frame()

rstan::Rhat(cbind(model_median_eval4[[1]][,"lambda"],model_median_eval4[[2]][,"lambda"]))


model_median_eval4_wo_adj <- nimbleMCMC(mod_median4, 
                                        data = list(zero = 0,jacob_adj = 0),
                                        constants = list(mean_expert = mu_ln,
                                                         sd_expert = sd_ln, 
                                                         C = 1000),
                                        monitors = c("median_St", "lambda", "LL_tot"),
                                        niter = 100000*5+10000, nchains =2, nburnin =10000, thin  = 5)


mod_res_wo_adj <- rbind(model_median_eval4_wo_adj[[1]],model_median_eval4_wo_adj[[2]]) %>% data.frame()
dens_expo <- density(mod_res_wo_adj[,"median_St"])
x_seq <- seq(0,50, by = 0.01)
rstan::Rhat(cbind(model_median_eval4_wo_adj[[1]][,"median_St"],model_median_eval4_wo_adj[[2]][,"median_St"]))


#Figure 3
png(paste0(path_files,"plots/Comparison of Priors.png"),width = 7, height = 7, units = 'in',res = 480)
plot(density(mod_res[,"median_St"]), xlab = "Time", ylim = c(0, 0.2),xlim = c(0,50), main = "Median Survival using LAPs")
lines(x_seq,dlnorm(x_seq, meanlog = mu_ln, sdlog = sd_ln), col = "red", lty = 2)
lines(dens_expo$x,dens_expo$y, col = "black", lty = 2)
legend("topright", 
       legend= c("LAP (Loss adjusted Posterior)", "True Density (Log-Normal)",
                 "LAP without adjustment for prior" ),
       col=c("black", "red", "black"),
       lty=c(1,2,2),
       cex=0.8)
dev.off()


# Figure 4: Plot of Posterior distributions of lambda and median survival
#For completness we specify the model in JAGS
#We assume different values of w2 the weight function for the loss associated with the expert opinion.

mod.jags_expo <- " #Exponential


model{
lambda ~ dunif(0,10)

n_events <- 5
cum_time <- 100

prec_expert <- pow(sd_expert,-2)
t_med <- log(2)/lambda

Lik <- dlnorm(t_med,mu_expert, prec_expert)
LL <- log(Lik)*w2*equals(jacob_adj,1)
LL_data <- (log(lambda)*n_events-cum_time*lambda)*equals(data_adj,1)

jacobian <- log(2)/(lambda^2)
log_jacobian <- log(jacobian)*equals(jacob_adj,1)
zero.mean <- -(LL + log_jacobian + LL_data) + C

zero ~ dpois(zero.mean)

C <- 10000

}"

n_sim <- 100000
mu_expert <- mu_ln#0.5
sd_expert <- sd_ln #0.1
data.jags <- list(mu_expert = mu_expert, sd_expert = sd_expert,
                  zero= 0,jacob_adj = 1, data_adj = 0, w2 = 1)


results_exponential_expert<- runjags::run.jags(model=mod.jags_expo,
                                               data = data.jags,
                                               monitor=c("lambda", "t_med", "jacobian"),
                                               n.chains=2, method="parallel", sample  = n_sim,summarise = FALSE,
                                               inits = function(chain){list(lambda = 1)})

data.jags <- list(mu_expert = mu_expert, sd_expert = sd_expert,
                  zero= 0,jacob_adj = 1, data_adj = 1, w2 = 0)

results_exponential_data<- runjags::run.jags(model=mod.jags_expo,
                                             data = data.jags,
                                             monitor=c("lambda", "t_med", "jacobian"),
                                             n.chains=2, method="parallel", sample  = n_sim,summarise = FALSE,
                                             inits = function(chain){list(lambda = 1)})
data.jags <- list(mu_expert = mu_expert, sd_expert = sd_expert,
                  zero= 0,jacob_adj = 1, data_adj = 1, w2 = 1)

results_exponential_w1 <- runjags::run.jags(model=mod.jags_expo,
                                            data = data.jags,
                                            monitor=c("lambda", "t_med", "jacobian"),
                                            n.chains=2, method="parallel", sample  = n_sim,summarise = FALSE,
                                            inits = function(chain){list(lambda = 1)})

data.jags <- list(mu_expert = mu_expert, sd_expert = sd_expert,
                  zero= 0,jacob_adj = 1, data_adj = 1, w2 = 0.5)

results_exponential_w_half <- runjags::run.jags(model=mod.jags_expo,
                                                data = data.jags,
                                                monitor=c("lambda", "t_med", "jacobian"),
                                                n.chains=2, method="parallel", sample  = n_sim,summarise = FALSE,
                                                inits = function(chain){list(lambda = 1)})


ggmcmc(ggs(results_exponential_expert$mcmc, family  = "lambda|t_med"), 
       plot = c("histogram", "density", "traceplot", "running", "compare_partial", 
                "Rhat", "ggs_effective"),
       file = paste0(path_files,"plots/Exponential model Diagnostics.pdf"))

results_exponential_expert_mat <- as.matrix(results_exponential_expert$mcmc)
results_exponential_data_mat <- as.matrix(results_exponential_data$mcmc)
results_exponential_w1_mat <- as.matrix(results_exponential_w1$mcmc)
results_exponential_w_half_mat <- as.matrix(results_exponential_w_half$mcmc)


dens_data_lambda <- density(results_exponential_data_mat[,"lambda"])
dens_data_tmed <- density(results_exponential_data_mat[,"t_med"])

dens_expert_lambda <- density(results_exponential_expert_mat[,"lambda"])
dens_expert_tmed <- density(results_exponential_expert_mat[,"t_med"])

dens_w1_lambda <- density(results_exponential_w1_mat[,"lambda"])
dens_w1_tmed <- density(results_exponential_w1_mat[,"t_med"])

dens_w_half_lambda <- density(results_exponential_w_half_mat[,"lambda"])
dens_w_half_tmed <- density(results_exponential_w_half_mat[,"t_med"])

.pardefault <- par()


png(paste0(path_files,"plots/Comparison_of_Log_normal_opinions.png"), width = 10, height = 7, units = 'in',res = 480)
par(mfrow = c(1,2), oma = c(5, 0, 0, 0))
plot(dens_data_lambda$x, dens_data_lambda$y , col = "red",type = "l",  ylab = "Density",ylim = c(0,25),
     xlab = expression(lambda),lty = 2, main = expression(paste("Posterior for  ", lambda)))
lines(dens_expert_lambda$x,dens_expert_lambda$y,lty = 2)
lines(dens_w1_lambda$x,dens_w1_lambda$y, col = "blue")
lines(dens_w_half_lambda$x,dens_w_half_lambda$y, col = "purple")

plot(dens_data_tmed$x, dens_data_tmed$y , col = "red",type = "l",  ylab = "",ylim = c(0,0.15),
     xlim = c(0, 60),
     xlab = "Time",lty = 2, main = expression("Posterior for Median Survival"))
lines(dens_expert_tmed$x,dens_expert_tmed$y,lty = 2)
lines(dens_w1_tmed$x,dens_w1_tmed$y, col = "blue")
lines(dens_w_half_tmed$x,dens_w_half_tmed$y, col = "purple")
par(xpd=NA)
legend(-50, -0.05, 
       legend= c("Posterior with Data only", "Posterior with Expert Opinion",
                 TeX(r'(Posterior with data and expert opinion given equal weight [$w_2$ = 1])'),
                 TeX(r'(Posterior with data and expert opinion given 50% weight [$w_2$ = 0.5])')),
       col=c("red", "black", "blue", "purple"),
       lty=c(2,2,1,1),
       cex=0.8)

dev.off()

par(.pardefault)

# p <- c(0.01,0.25, 0.5, 0.75, 0.99)
# quant_vec <- quantile(median_St, p)
# myfit <- fitdist(vals = quant_vec, probs = p, lower = 0, upper = 100)

p <- c(0.01,0.1,0.2,0.3, 0.4, 0.5,0.6,0.7,0.8,0.9,0.99)
quant_vec <- quantile(median_St, p)
myfit <- fitdist(vals = quant_vec, probs = p, lower = 0, upper = 100)


# Compare against a Student-t distribution

mod.jags_expo_t <- " #Exponential


model{
lambda ~ dunif(0,10)

n_events <- 5
cum_time <- 100

prec_expert <- pow(sd_expert,-2)
t_med <- log(2)/lambda

Lik <- dt(t_med,mu_expert, prec_expert, df )
LL <- log(Lik)*w2*equals(jacob_adj,1)
LL_data <- (log(lambda)*n_events-cum_time*lambda)*equals(data_adj,1)

jacobian <- log(2)/(lambda^2)
log_jacobian <- log(jacobian)*equals(jacob_adj,1)
zero.mean <- -(LL + log_jacobian + LL_data) + C

zero ~ dpois(zero.mean)

C <- 10000

}"
data.jags <- list(mu_expert = myfit$Student.t$location, sd_expert = myfit$Student.t$scale,
                  df = myfit$Student.t$df,
                  zero= 0,jacob_adj = 1, data_adj = 1, w2 =1)


results_exponential_t <- runjags::run.jags(model=mod.jags_expo_t,
                                                data = data.jags,
                                                monitor=c("lambda", "t_med", "jacobian", "Lik"),
                                                n.chains=2, method="parallel", sample  = n_sim,summarise = FALSE,
                                                inits = function(chain){list(lambda = 1)})


results_exponential_t_mat <- as.matrix(results_exponential_t$mcmc)

results_exponential_t_mat_median <- density(results_exponential_t_mat[,"t_med"])
results_exponential_t_mat_lambda <- density(results_exponential_t_mat[,"lambda"])
# index <- 10
# results_exponential_t_mat[index,]
# dt.scaled(results_exponential_t_mat[index,"t_med"], mean = data.jags$mu_expert, 
#           sd = data.jags$sd_expert, df = data.jags$df)




# Histogram Belief

#Histogram grid

mod.jags_hist <- "

model{

lambda ~ dunif(0,10)

n_events <- 5
cum_time <- 100

t_med <- log(2)/lambda


#Continous so no need for any interval adjusment
for(i in 1:dim_vec[1]){
int_eval[i]<- step(t_med - Prob_x[i,2])*step(Prob_x[i,3]-t_med)
prob_eval[i] <- Prob_x[i,1]*int_eval[i]
}

Lik <- sum(prob_eval)

LL <- log(Lik)*w2*equals(jacob_adj,1)
LL_data <- (log(lambda)*n_events-cum_time*lambda)*equals(data_adj,1)

jacobian <- log(2)/(lambda^2)
log_jacobian <- log(jacobian)*equals(jacob_adj,1)
zero.mean <- -(LL + log_jacobian + LL_data) + C

zero ~ dpois(zero.mean)

C <- 10000

}"


Probs_x <- matrix(nrow = 12, ncol = 3)
Probs_x[1,1]   <-p[1]
Probs_x[1,2:3] <- c(min(median_St), quant_vec[1])

for(i in 2:nrow(Probs_x)){
  Probs_x[i,1]   <-p[i]-p[i-1]
  Probs_x[i,2:3] <- c(quant_vec[i-1], quant_vec[i])
  
  if(i ==nrow(Probs_x)){
    Probs_x[i,1] <-  1-p[i-1]
    Probs_x[i,2:3] <- c(quant_vec[i-1], 100)
  }
}

# 
# Probs_x <- matrix(nrow = 6, ncol = 3)
# Probs_x[1,1]   <-p[1]
# Probs_x[1,2:3] <- c(min(median_St), quant_vec[1])
# Probs_x[2,1]   <-p[2]-p[1]
# Probs_x[2,2:3] <- c(quant_vec[1], quant_vec[2])
# Probs_x[3,1]   <-p[3]-p[2]
# Probs_x[3,2:3] <- c(quant_vec[2], quant_vec[3])
# Probs_x[4,1]   <- p[4]-p[3]
# Probs_x[4,2:3] <- c(quant_vec[3], quant_vec[4])
# Probs_x[5,1] <-  p[5]-p[4]
# Probs_x[5,2:3] <- c(quant_vec[4], quant_vec[5])
# Probs_x[6,1] <-  1-p[5]
# Probs_x[6,2:3] <- c(quant_vec[5], 100)

Probs_x[,1] <- Probs_x[,1]/sum(Probs_x[,1])#normalize

Probs_x[,1]  <- Probs_x[,1]/apply(Probs_x,1, function(x){x[3]-x[2]})  #turn to pdf

n_sim <- 500000




data.jags <- list(Prob_x  = Probs_x, dim_vec = nrow(Probs_x), 
                  zero= 0,jacob_adj = 1, data_adj = 1, w2 =1)


results_exponential_hist <- runjags::run.jags(model=mod.jags_hist,
                                           data = data.jags,
                                           monitor=c("lambda", "t_med", "jacobian", "Lik"),
                                           n.chains=2, method="parallel", sample  = n_sim,summarise = TRUE,
                                           inits = function(chain){list(lambda = 0.1)})

results_exponential_hist <- as.matrix(results_exponential_hist$mcmc)

results_exponential_hist_median <- density(results_exponential_hist[,"t_med"])
results_exponential_hist_lambda <- density(results_exponential_hist[,"lambda"])

plot(results_exponential_hist_median)
lines(Probs_x[,2], Probs_x[,1], type = "s", col = "red")

png(paste0(path_files,"plots/Comparison_of_Loss_types.png"), width = 10, height = 7, units = 'in',res = 480)
par(mfrow = c(1,2), oma = c(5, 0, 0, 0))
plot(dens_w1_lambda$x, dens_w1_lambda$y , col = "blue",type = "l",  ylab = "Density",ylim = c(0,35),
     xlab = expression(lambda),lty = 1, main = expression(paste("Posterior for  ", lambda)))
lines(results_exponential_t_mat_lambda$x,results_exponential_t_mat_lambda$y, col = "purple")
lines(results_exponential_hist_lambda$x,results_exponential_hist_lambda$y, col = "red")

plot(dens_w1_tmed$x, dens_w1_tmed$y , col = "blue",type = "l",  ylab = "",ylim = c(0,0.15),
     xlim = c(0, 60),
     xlab = "Time",lty = 1, main = expression("Posterior for Median Survival"))
lines(results_exponential_t_mat_median$x,results_exponential_t_mat_median$y, col = "purple")
lines(results_exponential_hist_median$x,results_exponential_hist_median$y, col = "red")

par(xpd=NA)
legend(-50, -0.05, 
       legend= c("Posterior - Log-Normal belief", "Posterior - t distribution belief", "Posterior - Histogram belief"),
       col=c("blue", "purple", "red"),
       lty=c(1,1,1),
       cex=0.8)

dev.off()

par(.pardefault)


# lines(x_seq, dt.scaled(x_seq, mean = data.jags$mu_expert, 
#                        sd = data.jags$sd_expert, df = data.jags$df), col = "red")

# Pathological example


mod.jags_expo_mix <- " #Exponential


model{
lambda ~ dunif(0,10)

n_events <- 5
cum_time <- 100

prec_expert <- pow(sd_expert,-2)
t_med <- log(2)/lambda

Lik <- dnorm(t_med,mu_expert[1],prec_expert[1])*0.5 +dnorm(t_med,mu_expert[2],prec_expert[2])*0.5
LL <- log(Lik)*w2*equals(jacob_adj,1)
LL_data <- (log(lambda)*n_events-cum_time*lambda)*equals(data_adj,1)

jacobian <- log(2)/(lambda^2)
log_jacobian <- log(jacobian)*equals(jacob_adj,1)
zero.mean <- -(LL + log_jacobian + LL_data) + C

zero ~ dpois(zero.mean)

C <- 10000

}"

n_sim <- 100000
mu_expert <- c(5, 20)#0.5
sd_expert <- c(2,2) #0.1
pi <- c(0.5,0.5)
data.jags <- list(mu_expert = mu_expert, sd_expert = sd_expert,
                  zero= 0,jacob_adj = 1, data_adj = 0, w2 = 1, pi = pi)


results_exponential_exper_mixt<- runjags::run.jags(model=mod.jags_expo_mix,
                                               data = data.jags,
                                               monitor=c("lambda", "t_med", "jacobian"),
                                               n.chains=2, method="parallel", sample  = n_sim,summarise = TRUE,
                                               inits = function(chain){list(lambda = 1)})


ggmcmc(ggs(results_exponential_exper_mixt$mcmc, family  = "lambda|t_med"), 
       plot = c("histogram", "density", "traceplot", "running", "compare_partial", 
                "Rhat", "ggs_effective"),
       file = paste0(path_files,"plots/Mixture Diagnostics.pdf"))


mix_mat <- as.matrix(results_exponential_exper_mixt$mcmc)
plot(density(mix_mat[,"lambda"]))

seq_x <- seq(0, 30, by = 0.01)
dens_mix <- dnorm(seq_x, mu_expert[1], sd = sd_expert[1])*pi[1]+dnorm(seq_x, mu_expert[2], sd = sd_expert[2])*pi[2]
plot(density(mix_mat[,"t_med"]))
lines(seq_x, dens_mix, col = "red")


t_expert <- 10
plot(seq_x, log(2)/(t_expert^seq_x), ylab = "Optimal scale", xlab = "shape")


mod.jags <- "
#data{
#zero <- 0
#}

model{
shape ~ dunif(0,100)
scale ~ dunif(0,100)

prec_expert <- pow(sd_expert,-2)
t_med <- pow(log(2)/scale, 1/shape)
Lik <- dnorm(t_med,mu_expert[1],prec_expert[1])*0.5 +dnorm(t_med,mu_expert[2],prec_expert[2])*0.5


jacobian1 <- pow(log(2), 1/shape)*pow(scale,-1/shape-1)/shape 
jacobian2 <- pow(log(2), 1/shape)*(log(scale)-log(log(2)))/(pow(scale, 1/shape)*pow(shape,2)) 
jacobian <- abs(jacobian2) +abs(jacobian1) #+jacobian2 #abs(jacobian1*jacobian2)*equals(jacob_adj,1)+ equals(jacob_adj,0)

zero.mean <- -log(Lik*jacobian) + C

zero ~ dpois(zero.mean)

C <- 10000

}"



n_sim <- 100000
mu_expert <- c(5, 20)#0.5
sd_expert <- c(2,2) #0.1
pi <- c(0.5,0.5)
data.jags <- list(mu_expert = mu_expert, sd_expert = sd_expert,
                  zero= 0,jacob_adj = 1, data_adj = 0, w2 = 1, pi = pi)


results <- runjags::run.jags(model=mod.jags,
                             data = data.jags,
                             monitor=c("shape", "scale","t_med", "zero.mean"),
                             n.chains=2, method="parallel", sample  = 300000, summarise = TRUE,
                             inits = function(chain){list(shape = 1,#runif(1,0,2),
                                                          scale = 10#runif(1,0,2)
                             )})
mat_res <- as.matrix(results$mcmc)
#plot(x = mat_res[,"shape"], y = mat_res[,"scale"])

png(paste0(path_files,"plots/Mixture Model comparison.png"), width = 10, height = 7, units = 'in',res = 480)
plot(density(mat_res[,"t_med"]), main = "Median Survival")
lines(seq_x, dens_mix, col = "red")
legend("topright", col = c("black", "red"), legend = c("Posterior sample", "True mixture density"))
dev.off()
plot(mat_res[,"shape"],mat_res[,"scale"])




ggmcmc(ggs(results$mcmc, family  = "lambda|t_med"), 
       plot = c("histogram", "density", "traceplot", "running", "compare_partial", 
                "Rhat", "ggs_effective"),
       file = paste0(path_files,"plots/Mixture Diagnostics2.pdf"))




mod.jags <- "#Weibull PH
#data{
#zero <- 0
#}

model{
shape ~ dunif(0,50)
scale ~ dunif(0,5)

prec_expert <- pow(sd_expert,-2)
St <- exp(-scale*pow(t_expert,shape))
Lik <- dnorm(St,mu_expert[1],prec_expert[1])*0.5 +dnorm(St,mu_expert[2],prec_expert[2])*0.5

#Lik <- dbeta(St, 1,1)
jacobian1 <- pow(t_expert,shape)*St*log(t_expert)*scale
jacobian2 <- pow(t_expert,shape)*St

jacobian <- abs(jacobian2)+abs(jacobian1)
zero.mean <- -log(Lik*abs(jacobian)) + C

zero ~ dpois(zero.mean)

C <- 10000

}"


mod.jags <- " #Weibull ATF
#data{
#zero <- 0
#}

model{
shape ~ dunif(0,10)
scale ~ dunif(0,50)

prec_expert <- pow(sd_expert,-2)
St <- exp(-pow(t_expert/scale,shape))
Lik <- dnorm(St,mu_expert[1],prec_expert[1])*0.5 +dnorm(St,mu_expert[2],prec_expert[2])*0.5

#Lik <- dbeta(St, 1,1)
jacobian1 <- -log(t_expert/scale)*pow(t_expert/scale,shape)*exp(-(t_expert/scale)^shape)
jacobian2 <- (shape*(t_expert/scale)^shape*exp(-(t_expert/scale)^shape))/scale

jacobian <- abs(jacobian2)+abs(jacobian1)
zero.mean <- -log(Lik*abs(jacobian)) + C

zero ~ dpois(zero.mean)

C <- 10000

}"



n_sim <- 500000
mu_expert <- c(0.2, 0.8)#0.5
sd_expert <- c(0.05,0.05) #0.1
pi <- c(0.5,0.5)
data.jags <- list(mu_expert = mu_expert, sd_expert = sd_expert,
                  zero= 0,jacob_adj = 1, data_adj = 0, w2 = 1, pi = pi,
                  t_expert = 3)



weibull_sims <- flexsurv::rweibullPH(100, shape = 0.6 , scale = 0.6)
# library(flexsurv)
# mle.fit <- flexsurvreg(Surv(weibull_sims,rep(1, length(weibull_sims)))~1, dist = "weibullPH")

results <- runjags::run.jags(model=mod.jags,
                             data = data.jags,
                             monitor=c("shape", "scale","St"),
                             n.chains=2, method="parallel", sample  = n_sim,
                             inits = function(chain){list(shape = runif(1,0,2), scale = runif(1,0,2))},
                             summarise = FALSE)
mat_res <- as.matrix(results$mcmc)
#plot(x = mat_res[,"shape"], y = mat_res[,"scale"])
dens_mix <- dnorm(seq_x, mu_expert[1], sd = sd_expert[1])*pi[1]+dnorm(seq_x, mu_expert[2], sd = sd_expert[2])*pi[2]

png(paste0(path_files,"plots/St_mix.png"), width = 10, height = 7, units = 'in',res = 480)
plot(density(mat_res[,"St"]), xlab = "St", col = "red")
lines(seq_x, dens_mix)
dev.off()

png(paste0(path_files,"plots/shape-scale plot.png"), width = 10, height = 7, units = 'in',res = 480)

plot(mat_res[,"shape"],mat_res[,"scale"], xlab = "scale", ylab = "shape", xlim= c(0,5))
dev.off()



ggmcmc(ggs(results$mcmc, family  = "shape|scale|St"), 
       plot = c("histogram", "density", "traceplot", "running", "compare_partial", 
                "Rhat", "ggs_effective"),
       file = paste0(path_files,"plots/Mixture Diagnostics2.pdf"))




