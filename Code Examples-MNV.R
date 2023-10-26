
list.of.packages <- need<-c("ggplot2", "conflicted","rootSolve", "rstan","optimx","mvtnorm","Ryacas") #needed libraries
res <- lapply(list.of.packages, require, character.only = TRUE)

  not.loaded <-   list.of.packages[which(sapply(res, unlist) ==F)]
  not.installed <-   not.loaded
  #load the packages
  if(length(not.installed)) install.packages(not.installed)
  if(length(not.installed)) lapply(not.installed, require, character.only = TRUE)
  
  #path_files <- #Define working directory
  
source(paste0(path_files,"helper functions.R"))




conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("footnote", "kableExtra")
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("summarise", "plyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("ar", "brms")
conflict_prefer("diag", "base")
conflict_prefer("extract", "rstan")
conflict_prefer("%*%", "base")
source(paste0(path_files,"helper functions.R"))



#"True" data
mu_vec <- c(5,2, 1, 3)
sd_vec <- c(2,1, 0.5, 3)

#Expert elicated quantiles and sample size
quantile_probs <- c( 0.5,0.75)
sample_size <-10

# Elicited values generated from a normal distribution
elicited_vals <- t(apply(cbind(mu_vec,sd_vec),1, function(x){qnorm(quantile_probs,x[1], x[2])}))
hyperprior_mat<- t(apply(elicited_vals,1,  function(x){hyper_norm_gamma(x, quantile_probs,sample_size)}))


n_corr_param <- length(mu_vec)*(length(mu_vec)-1)/2

#Median is invariant to transformation (see with n = odd number so that median is unique)
median_concord <- c(0.60, 0.25, 0.4, rep(0.5,3)) 
# Calculate median correlation (from median concordance) then to Fisher's z transformation
z_concord_transform <- atanh(sin(2*pi*(median_concord/2 -0.25))) 
se_z_concord <- 1/sqrt(sample_size-3) #assume to be equal for all correlations


quantile((asin(tanh(rnorm(100000,z_concord_transform[2], se_z_concord)))/pi)+0.5, probs = c(0.025,0.975))

symb_mat <- gen_symb_mats(dim_1 = length(mu_vec), return_Ryacas = F)
symb_mat[lower.tri(symb_mat)] <-   0

corr_row_index <- rep(NA, n_corr_param)
corr_col_index <-  rep(NA, n_corr_param)

for(i in 1:n_corr_param){
  
 index_temps <-  which(symb_mat == letters[i],arr.ind = TRUE)
 corr_row_index[i] <- index_temps[1,1]
 corr_col_index[i] <- index_temps[1,2]
 symb_mat[index_temps[1,1],index_temps[1,2]] <- median_concord[i]
}
concord_mat<- cbind(corr_row_index,corr_col_index, median_concord)

seq_prob <- seq(0,1, by = 0.01)

for(i in 1:length(median_concord)){
  gen_pd_corr_par(median_concord,corr_row_index,corr_col_index, mu_length = length(mu_vec), i,  plot_roots = F)
}

model_stan <- "

functions {
  vector to_corr_vec(matrix R, int K , int num_corr) {
    vector[num_corr] y;    
    int pos = 1;
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        y[pos] = R[i, j];
        pos += 1;
      }
    }
    return y;
  }
}

data {
  int<lower=1> K; // dimension of observations
  int<lower=1> num_quantile; // Number of quantiles
  real eta;
  int num_corr;
  vector[num_corr] z_concord_transform;
  real se_z_concord;
  real pi;
  matrix[K,4] hyperprior_mat;
}
parameters {
  cholesky_factor_corr[K] Lcorr; // cholesky factor (L_u matrix for R)
  vector[K] mu_vec;
  vector[K] prec_vec;
}
transformed parameters {
  corr_matrix[K] R; // correlation matrix
  vector[num_corr] corr_vector;
  vector[num_corr] z_vector;
  vector[K] sd_vec_mu;
  vector[K] sd_vec;
  
  R = multiply_lower_tri_self_transpose(Lcorr); // R = Lcorr * Lcorr'
  corr_vector = to_corr_vec(R, K,num_corr);
  for(i in 1:num_corr){
  z_vector[i] = atanh(corr_vector[i]);
  }
  for(i in 1:K){
  sd_vec_mu[i] = (hyperprior_mat[i,2]*prec_vec[i])^(-0.5);
  sd_vec[i] = (prec_vec[i])^(-0.5);
  }

}
model {
  Lcorr ~ lkj_corr_cholesky(eta); // prior for cholesky factor of a correlation matrix
  target += -lkj_corr_cholesky_lpdf(Lcorr |  eta);
  for(i in 1:num_corr){
  target += normal_lpdf(z_vector[i]| z_concord_transform[i],se_z_concord );
  }
  
  for(i in 1:K){
  prec_vec[i] ~ gamma(hyperprior_mat[i,3],hyperprior_mat[i,4]);
  mu_vec[i] ~ normal(hyperprior_mat[i,1],sd_vec_mu[i]);
  }
  
}  

generated quantities {
vector[K] x_data_sim;
vector[num_corr]concord_vector;
cov_matrix[K] Sigma;
vector[K] x_data_sim2;

// prior predictive
  for(i in 1:K) {
    x_data_sim[i] = normal_rng(mu_vec[i],sd_vec[i]);
  }
 for(i in 1:num_corr){
  concord_vector[i] = (asin(corr_vector[i])/pi)+0.5;
  }
 Sigma = quad_form_diag(R, sd_vec); // quad_form_diag: diag_matrix(sig) * R * diag_matrix(sig)

x_data_sim2 = multi_normal_rng(mu_vec,Sigma);
}
  
"

num_corr <- length(mu_vec)*(length(mu_vec)-1)/2
 
#Generate some Latex tables
# xtable::xtable(apply(elicited_vals,c(1,2), round,digits = 2))
# xtable::xtable(concord_mat)


data.stan <- list()
data.stan$K <- length(mu_vec)
data.stan$eta <- 1
data.stan$z_concord_transform <- z_concord_transform
data.stan$se_z_concord <- se_z_concord
data.stan$num_corr <- num_corr
data.stan$pi <- pi
data.stan$hyperprior_mat <- hyperprior_mat

data.stan$quantile_probs <- quantile_probs
data.stan$num_quantile <- length(quantile_probs)
iter <- 10000
chains <- 2
model_stan_MVN <- rstan::stan_model(model_code = model_stan) # compilation



model_MVN <- rstan::sampling(model_stan_MVN, data.stan, chains = chains, 
                             iter = iter, warmup = iter/10, thin =1 )



x_data_sim2 <- rstan::extract(model_MVN, pars = c("x_data_sim2")) 
stan_res  <- rstan::extract(model_MVN, pars = c('concord_vector',"corr_vector", "x_data_sim")) 
stan_res_concord  <- rstan::extract(model_MVN, pars = c('concord_vector')) 


for(i in 1:num_corr){
  assign(paste0("dens_",i),density(stan_res$concord_vector[,i]) )
  assign(paste0("dens_corr_",i),density(stan_res$corr_vector[,i]) )
}
seq_prob <- seq(0,1, by = 0.01)

#Figure 6a and 6b

png(paste0(path_files,"plots/Concordance_Probability.png"),  width = 7, height = 7, units = 'in',res = 250)
plot(dens_1,xlim = c(0,1), ylim = c(0,6), type = "l", main  = "Concordance Probability Distributions",xlab = "Probability")
for(i in 2:num_corr){
  temp_dens <- get(paste0("dens_",i))
  lines(temp_dens$x,temp_dens$y,col = i)
}
legend("topright", 
       legend=expression('p'[12],'p'[13],'p'[23],'p'[14],'p'[24],'p'[34]),
       lty = 1,
       col=1:num_corr,
       title=NULL, text.font=4, bg='white')
dev.off()

png(paste0(path_files,"plots/Correlation_Probability.png"),  width = 7, height = 7, units = 'in',res = 250)
plot(dens_corr_1,xlim = c(-1,1), ylim = c(0,6), type = "l", main  = "Partial Correlation Distributions",xlab = "Correlation")
for(i in 2:num_corr){
  temp_dens <- get(paste0("dens_corr_",i))
  lines(temp_dens$x,temp_dens$y,col = i)
}
legend("topright", 
       legend=expression(rho[12],rho[13],rho[23],rho[14],rho[24],rho[34]),
       lty = 1,
       col=1:num_corr,
       title=NULL, text.font=4, bg='white')
dev.off()


#Confirm Prior Predicitive Dist is a Student's t distribution

index <-4
a_eval <- hyperprior_mat[index,3]
b_eval <- hyperprior_mat[index,4]
scale_eval <- hyperprior_mat[index,2]
mu_eval <- hyperprior_mat[index,1]


dens_prior_pred <- density(x_data_sim2$x_data_sim2[,index])
#dens_prior_pred <- density(stan_res$x_data_sim[,index])
plot(dens_prior_pred$x,dens_prior_pred$y,type = "l", col = "blue" )
lines(dens_prior_pred$x, dt.scaled(dens_prior_pred$x,df = 2*a_eval,mu_eval,sd = sqrt(b_eval*(scale_eval+1)/(scale_eval*a_eval)) ),col = "red") 






#Generate Plot for LJK distribution 

model_stan_LKJ <- "

data {
  int<lower=1> K; // dimension of observations
  real eta;
}
parameters {
  cholesky_factor_corr[K] Lcorr; // cholesky factor (L_u matrix for R)
}
transformed parameters {
  corr_matrix[K] R; // correlation matrix
  R = multiply_lower_tri_self_transpose(Lcorr); // R = Lcorr * Lcorr'
 }
model {
  Lcorr ~ lkj_corr_cholesky(eta); // prior for cholesky factor of a correlation matrix
  target += -lkj_corr_cholesky_lpdf(Lcorr |  eta);
 }  

"

model_stan_LKJ2 <- "

data {
  int<lower=1> K; // dimension of observations
  real eta;
}
parameters {
  cholesky_factor_corr[K] Lcorr; // cholesky factor (L_u matrix for R)
}
transformed parameters {
  corr_matrix[K] R; // correlation matrix
  R = multiply_lower_tri_self_transpose(Lcorr); // R = Lcorr * Lcorr'
 }
model {
  Lcorr ~ lkj_corr_cholesky(eta); // prior for cholesky factor of a correlation matrix
  //target += -lkj_corr_cholesky_lpdf(Lcorr |  eta);
 }  

"


stan_LKJ1 <- rstan::stan_model(model_code = model_stan_LKJ) # compilation
stan_LKJ2 <- rstan::stan_model(model_code = model_stan_LKJ2) # compilation

input = list(eta =1, K = 4)
iter <- 200000

test_1 <- rstan::sampling(stan_LKJ1, input, chains = 1, 
                          iter = iter, warmup = 10, thin =1 )
stan_res1  <- rstan::extract(test_1, pars = c('R[1,2]')) 

#plot(test_1, show_density = TRUE, pars = c('R[1,2]'))

test_2 <- rstan::sampling(stan_LKJ2, input, chains = 1, 
                          iter = iter, warmup = 10, thin =1 )
stan_res2  <- rstan::extract(test_2, pars = c('R[1,2]')) 

stan_res_df1 <- data.frame(R_vals = stan_res1$`R[1,2]`, model = "Including Loss Function")
stan_res_df2 <- data.frame(R_vals = stan_res2$`R[1,2]`, model = "Standard")


#Seems to be half as effiecnt in terms of effective samples size)
stan_final_df <- rbind(stan_res_df1,stan_res_df2)

ggplot(stan_final_df, aes(R_vals, colour = model)) +
  geom_density()+
  xlab(expression(rho))+
  ylab("Density")+
  ggtitle(expression("Marginal correlation for LKJ("~eta~"=1) and k = 4"))+
  theme_bw()

ggsave(paste0(path_files,"plots/LKJ_prior.png"))



#Verify that concordance probability works for any Normal distribution; need not be mu = zero; sd = 1
#Additionally verify that the same calculation works for a (non-standard) Student's-t distribution.
corr <- -0.5
n_sim <- 20000000
R = matrix(c(1, corr,
             corr, 1), 2)
sd_vec <- c(1,2)
mu_vec <- c(1,2)
sigma <- diag(sd_vec)%*%R%*%diag(sd_vec)

if(F){
  x <- rmvnorm(n=n_sim, mean=mu_vec, sigma=sigma)
  
  sum(x[,1]-colMeans(x)[1] < 0 & x[,2]-colMeans(x)[2] <0)/n_sim
}
0.25+ asin(corr)/(2*pi)
dof <- 10
x2 <-rmvt(n=n_sim,delta  =mu_vec,  sigma=sigma*(dof-2)/dof, df=dof)
sum(x2[,1]-colMeans(x2)[1] < 0 & x2[,2]-colMeans(x2)[2] <0)/n_sim

#Show that the marginal of a multivariate-t is a univariate-t
if(F){
  library("mvtnorm")
  R = matrix(c(1, 0.7, 0.2,
               0.7, 1, -0.5,
               0.2, -0.5, 1), 3)
  
  sd_vec <- c(.2,0.5,0.1)
  Sigma <-base::diag(sd_vec)%*%R%*%base::diag(sd_vec)
  mu_vec <- c(5, 2,1)
  dof <- 10
  x <- rmvt(1000000,mu_vec, sigma = Sigma, df = dof) # t_3(0, diag(2)) sample
  #plot(x[1:1000,])
  dens_sim <- density(x[,1])
  x_dens <- rep(NA, length(dens_sim$x))
  for(i in 1:length(dens_sim$x)){
    x_dens[i] <-  dmvt(dens_sim$x[i], delta  = mu_vec[1,drop =F],sigma = Sigma[1,1, drop = F], df = dof,log = F) 
  }
  plot(dens_sim)
  lines(dens_sim$x, x_dens, col = "blue")
  
}


