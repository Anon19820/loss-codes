
list.of.packages <- need<-c("ggplot2","dplyr","tidyr", "brms", "stringr", "posterior", "ggmcmc" ) #needed libraries

res <- lapply(list.of.packages, require, character.only = TRUE)
  not.loaded <-   list.of.packages[which(sapply(res, unlist) ==F)]
  not.installed <-   not.loaded
  #load the packages
  if(length(not.installed)) install.packages(not.installed)
  if(length(not.installed)) lapply(not.installed, require, character.only = TRUE)
  
  #path_files <- #Define working directory
  path_files <- "~/Expert Opinion - General/"
  source(paste0(path_files,"helper functions.R"))

#Load Data;



data_exercise <- read.delim(paste0(path_files,"exercise.txt"),sep = " ")
data_exercise2 <- data_exercise %>% dplyr::select(-X, )
data_exercise2$Subject <-  paste0(data_exercise2$PROGRAM, "-", data_exercise2$ID)
data_exercise2 <- dplyr::select(data_exercise2,-ID)
data_exercise3 <- gather(data_exercise2, "Outcome", "Score", -PROGRAM ,-Subject)

data_exercise3$time <- as.numeric(stringr:::str_remove_all(data_exercise3$Outcome, "S"))-1
data_exercise3$PROGRAM <- as.factor(data_exercise3$PROGRAM)


data_exercise3_orth <- data_exercise3
data_exercise3_orth$time.orth.1 <- poly(data_exercise3_orth$time, 2, raw=FALSE)[,1] # orthogonal version of time
data_exercise3_orth$time.orth.2 <- poly(data_exercise3_orth$time, 2, raw=FALSE)[,2] # orthogonal version of time

chains <- 2
iter <- 10000

multiModel_stan_random_slope<-brms::brm(Score~PROGRAM+ time.orth.1 +time.orth.2+ PROGRAM*time.orth.1 +PROGRAM*time.orth.2 +
                                        (1+time.orth.1+time.orth.2|Subject),
                                        data=data_exercise3_orth,
                                        iter =iter,
                                        warmup = iter/10,
                                        chains = chains)

plot(multiModel_stan_random_slope, variable = "b_Intercept", regex = TRUE, )

#Extract Stan Data
data.stan <- brms::standata(multiModel_stan_random_slope)


#Extract the index for the design Matrix
index_vec <-c(min(which(data_exercise3$PROGRAM == "WI" & data_exercise3$time == 0)),
              min(which(data_exercise3$PROGRAM == "WI" & data_exercise3$time == 6)))

#We Verify that the uniform prior on Beta coeffeicents have no impact on the observable outcome
n_sims <- 1000000
sims_1 <- runif(n_sims, min = -1000, max = 1000)
sims_2 <- runif(n_sims, min = -1000, max = 1000)

plot(density(6*sims_1 +36*sims_2))


# Define data for modified brms model

unique_index <-as.numeric(rownames(unique(data_exercise3[ ,-c(2,4)])))
data.stan2 <- data.stan
data.stan2$pred_mat <- data.stan$X[unique_index,-1]
data.stan2$n_pred <- nrow(data.stan2$pred_mat)

X_expert <- data.stan$X[index_vec[2],2:ncol(data.stan$X),drop = F]-data.stan$X[index_vec[1],2:ncol(data.stan$X),drop = F]
X_expert<- apply(X_expert,c(1,2),function(x){if(x <1e-15){x=0}else{x=x} })
data.stan2$X_expert <- X_expert
data.stan2$par_expert <- matrix(c(2.5,1.5), nrow= 1, ncol = 2)
data.stan2$par_expert <- matrix(c(0.5,0.5), nrow= 1, ncol = 2)

#data.stan2$prior_only <- 1 #0 if including data

#The default prior for population-level effects (including monotonic and category specific effects) is an improper flat prior over the reals.
model_test<- "// generated with brms 2.18.0
functions {
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  vector[N] Z_1_3;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
  
  matrix[1, K-1] X_expert; // Design matrix corresponding to the expert opinion (difference between final observation and baseline)
  matrix[1, 2] par_expert; // Parameters of the expert's distribution
  int<lower=1> n_pred;
  matrix[n_pred,K-1] pred_mat;
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // dispersion parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  vector[N_1] r_1_3;
  real lprior = 0;  // prior contributions to the log posterior
  vector[1] mu_diff;

  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_1 = r_1[, 1];
  r_1_2 = r_1[, 2];
  r_1_3 = r_1[, 3];
  lprior += student_t_lpdf(Intercept | 3, 81, 3);
  lprior += student_t_lpdf(sigma | 3, 0, 3)
    - 1 * student_t_lccdf(0 | 3, 0, 3);
  lprior += student_t_lpdf(sd_1 | 3, 0, 3)
    - 3 * student_t_lccdf(0 | 3, 0, 3);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1);
  
   mu_diff = X_expert*b;
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n] + r_1_3[J_1[n]] * Z_1_3[n];
    }
    target += normal_id_glm_lpdf(Y | Xc, mu, b, sigma);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(to_vector(z_1));
  
  // Expert Opinion
  target += normal_lpdf(mu_diff[1]|par_expert[1,1],par_expert[1,2]);

}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  vector[n_pred] pred_fixed;

  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
  
  pred_fixed =  pred_mat*b +b_Intercept;
}"
model_stan_custom2 <- rstan::stan_model(model_code = model_test) # compilation


model_rep_meas2 <- rstan::sampling(model_stan_custom2, data.stan2, chains = chains, 
                                   iter = iter, warmup = iter/10, thin =1, cores = 1 )



pars_eval <- rstan::extract(model_rep_meas2, pars = c("b", "Intercept","sigma", "mu_diff"), permuted = FALSE)

posterior_array <- posterior::as_draws_array(pars_eval)

posterior_summary<- summarize_draws(posterior_array)

ggmcmc(ggs(model_rep_meas2, family  = "b|sigma|mu_diff"), 
       plot = c("histogram", "density", "traceplot", "running", "compare_partial", 
                "Rhat", "ggs_effective"),
           file = paste0(path_files,"plots/Regression Model Diagnostics.pdf"))



#Some plots to assess convergence
#traceplot(model_rep_meas2, pars = "b", inc_warmup = TRUE)
#traceplot(model_rep_meas2, pars = "b_Intercept", inc_warmup = TRUE)

#Verify that I calculated the population level response values correctly:
predicted_stan <- as.numeric(colMeans(rstan::extract(model_rep_meas2, pars = "pred_fixed")[[1]]))
predicted_manual <- as.numeric(data.stan2$pred_mat%*%colMeans(rstan::extract(model_rep_meas2, pars = "b")[[1]]) +mean(rstan::extract(model_rep_meas2, pars = "b_Intercept")[[1]]))
predicted_stan-predicted_manual

# Extract the fixed effects
data_exercise4 <- data_exercise3_orth %>%
  group_by(PROGRAM, time)%>% summarise(mean_data = mean(Score))

data_exercise5 <- data_exercise4 %>% ungroup()%>%arrange(time ,PROGRAM)
data_exercise5$mean_expert <- predicted_stan

samples1 <- posterior_samples(multiModel_stan_random_slope, "^b")
data_exercise5$mean_fixed <- cbind(1,data.stan2$pred_mat)%*%colMeans(samples1)


fixed_effect_data_posterior <- apply(samples1, 1, function(x){cbind(1,data.stan2$pred_mat)%*%x})
fixed_effect_data_posterior_quantile <- t(apply(fixed_effect_data_posterior, 1, quantile, probs = c(0.025,0.975)))

pred_expert <- rstan::extract(model_rep_meas2, pars = "pred_fixed")[[1]]
pred_expert <- t(apply(pred_expert, 2, quantile, probs = c(0.025,0.975)))
colnames(pred_expert) <- paste0(colnames(pred_expert), "-Expert")

WI_compare <- cbind(data_exercise5,fixed_effect_data_posterior_quantile,pred_expert) %>% dplyr::filter(PROGRAM == "WI")

#Figure 7
ggplot(data_exercise5,
       aes(x=time, y=mean_fixed, colour =as.factor(PROGRAM)))+
  geom_line()+
  geom_point(aes(x=time, y=mean_data, colour =as.factor(PROGRAM)))+
  geom_line(aes(x=time, y=mean_expert, colour =as.factor(PROGRAM)),linetype = "dashed")+
  theme_bw()+
  xlab("Time since baseline")+
  ylab("Functional Score")+
  guides(col=guide_legend("Program"))


ggsave(paste0(path_files,"plots/Expert-Repeated-Measures.png"))

#Heuristic Estimation of ESS
mu_diff_data <- apply(samples1, 1, function(x){cbind(0,data.stan2$X_expert)%*%x})
sd(mu_diff_data) # 0.6722
dens_eval_data <- density(mu_diff_data)

n_WI <- length(which(data_exercise3$PROGRAM == "WI" & data_exercise3$time == 0))

n_WI*(sd(mu_diff_data)/data.stan2$par_expert[1,2])^2

#Figure 8

mu_diff_expert <- rstan::extract(model_rep_meas2, pars = "mu_diff[1]")[[1]]
sd(mu_diff_expert)
dens_eval_expert <- density(mu_diff_expert)
png(paste0(path_files,"plots/Change-Baseline.png"),  width = 7, height = 7, units = 'in',res = 250)
plot(dens_eval_expert$x, dens_eval_expert$y,col = "red", type = "l", xlab = "Change from baseline", ylab = "Density", 
     main = "Expected Change from Baseline - WI Group")
lines(dens_eval_data$x, dens_eval_data$y,col = "blue")
legend("topright", legend=c("Data only", "Expert + Data"),
       col=c("blue", "red"), lty=1, cex=0.8)
dev.off()



