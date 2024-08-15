data {
  int<lower=1> N; //length of data
  int<lower=1> J; //n pops
  int<lower=1> K; //n regions
  int<lower=1> pop[N]; //vector indexing populations
  int<lower=1> region[N]; //vector indexing regions
  int<lower=1> reg_pop[J]; //region-pop index
  matrix[N, J] X; //predictor matrix of spawners
  vector[N] logRS; //the response variable
  vector[J] SmaxPR; //Smax estimates for each pop
  vector[J] SmaxPR_SD; //Smax error for each pop
}

transformed data{
  vector[J] logbeta_pr;
  vector[J] logbeta_pr_sig;

  for(i in 1:J){
    logbeta_pr_sig[i] = sqrt(log(1+((1/SmaxPR_SD[i])*(1/SmaxPR_SD[i]))/((1/SmaxPR[i])*(1/SmaxPR[i])))); //this converts sigma on the untransformed scale to a log scale
    logbeta_pr[i] = log(1/SmaxPR[i])-0.5*logbeta_pr_sig[i]*logbeta_pr_sig[i]; //convert smax prior to per capita slope - transform to log scale with bias correction
  }
}

parameters {
  real log_a; //global alpha
  vector[J] z_a_dev_pop; //pop specific deviation from global
  vector[K] z_a_dev_region; //regional alpha deviations
  vector[J] log_b_pop; //slopes for each pop 
  vector<lower=0>[J] sigma; //main variance term
  vector<lower=0>[J] sd_a_pop; //variance of means for around pops
  real<lower=0> sd_a_region; //variance of means for region
}

transformed parameters {
  vector<lower=0>[K] log_a_reg; //regional alpha deviations
  vector<lower=0>[J] log_a_pop; //population alpha deviations
  vector[J] b_pop;
  
  log_a_reg = log_a + z_a_dev_region[reg_pop]*sd_a_region; //regional alpha estimated
  log_a_pop = log_a_reg[reg_pop]+ z_a_dev_pop[reg_pop].*sd_a_pop[reg_pop];//realized log_a_pop
  
  b_pop = exp(log_b_pop);
}

model {
  //priors
  log_a ~ normal(1.5,2);
  z_a_dev_region ~ normal(0,1);
  z_a_dev_pop ~ normal(0,1);
  log_b_pop ~ normal(logbeta_pr, logbeta_pr_sig); 
  
  //variance priors
  sigma ~ normal(1,1); 
  sd_a_pop ~ gamma(2,3); 
  sd_a_region ~ gamma(2,3);

  //likelihood model
  logRS ~ normal(log_a_pop[pop] - X*b_pop, sigma[pop]);
}

//helper to get log_lik from extract_log_lik help file
//generated quantities{
//  vector[N]  log_lik; 
//  for (i in 1:N){log_lik[i] = normal_lpdf(logRS[i] | log_a_pop[pop[i]] - X[i]*b_pop, sigma); //HOW REF SIGMA?
//  }
//}
