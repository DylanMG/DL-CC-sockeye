library(tidyverse)
library(here)
library(rstan) 
set.seed(123)

# wrangle and load data ------------------------------------------------------------------
SmaxPR <- read.table(here("output/data/SmaxPRs.txt"), 
                     header = TRUE) |>
  select(population, region, water_clarity, SmaxPR, SmaxPR_SD) |>
  dplyr::rename(pop = population) #to match below

SR_data <- read.table(here("data/FinalBroodSK.txt"), 
                      header = TRUE) |>
  select(-n.data, -n.stock) |> #dump useless vars
  mutate(watershed = gsub("_lake|_lakes", "", watershed)) |>
  #fix some names so it joins
  mutate(watershed = case_when(watershed == "owekino" ~ "owikeno",
                               watershed == "portjohn" ~ "port_john",
                               watershed == "freeda_brodie" ~ "freeda_lakes",
                               TRUE ~ watershed)) |> #leave the rest alone
  mutate(logRS = log(total_rec/exp_spawn)) |> #shorten names for plotting, etc.
  dplyr::rename(pop = watershed) |>
  relocate(pop, 1) |>
  filter(pop !="end_hill") |>
  arrange(pop, year) |>
  left_join(SmaxPR, by = "pop") |> #join it to have 1 object to reference
  mutate(pop = gsub("_", "-", pop)) |> #change delimiter for use later
  na.omit() |>
  group_by(pop) |>
  mutate(exp_spawn_scaled = scale(exp_spawn))

#loop pops, fit models -------------------------------------------------------------------
iter <- 4000 #how many iterations in a stan model?
m1_LL <- NULL #empty objects to cbind LLs to
m2a_LL <- NULL
m2b_LL <- NULL
m3_LL <- NULL
m4a_LL <- NULL
m4b_LL <- NULL

#SR_data <- filter(SR_data, pop %in% unique(SR_data$pop)[1:10]) #subset for testing
for(i in unique(SR_data$pop)[1:3]){
  begin <- Sys.time()
  #could add "for j in 1:n_reps" or something to repeat sampling
  sub_pop <- filter(SR_data, pop == i)
  
  thinned_pop <- sub_pop |> #reduce 1 pop's data by 50%
    sample_frac(0.5)
  
  sub_data <- filter(SR_data, pop !=i) |> #write reduced pop back in with rest of data
    bind_rows(thinned_pop) |>
    arrange(pop, year)
  
  OOS_pop <- anti_join(sub_pop, thinned_pop) #the data that was removed; to predict OOS
  
  #write new stan data ---
  logRS <- sub_data$logRS #response variable
  pop <- as.numeric(as.factor(sub_data$pop)) #pop index vector
  N <- as.numeric(length(logRS)) #n observations
  J <- max(pop) #n populations
  K <- max(as.numeric(as.factor(sub_data$region))) #n regions
  region <- as.numeric(as.factor(sub_data$region)) #region index vector
  reg_pop <- sub_data |> #region-pop key for indexing
    select(pop, region) |>
    distinct() |>
    pull(region) |>
    as.factor() |>
    as.numeric()
  
  X <- matrix(rep(0, (N*J)), 
              nrow=N, 
              ncol=J) #empty predictor matrix to populate
  
  for(j in 1:N){
    X[j, pop[j]] <- sub_data$exp_spawn[j] #populated!
  }
  int_X <- ifelse(X!=0, 1,0) #binary matrix to index when to "turn on" intercept 
  SmaxPR <- as.numeric(unique(sub_data$SmaxPR)) #SmaxPR prior mus
  SmaxPR_SD <- sub_data |> #SmaxPR SDs
    select(pop, SmaxPR_SD) |>
    distinct() |>
    pull(SmaxPR_SD) |> #got creative because a SD repeats
    as.numeric() 

  #getting the out of sample data to predict --- 
  N_OOS <- as.numeric(nrow(OOS_pop))
  y_OOS <- as.numeric(OOS_pop$logRS) 
  X_OOS <- matrix(rep(0, (N_OOS*J)), 
              nrow=N_OOS, 
              ncol=J) #empty predictor matrix to populate
  pop_OOS <- which(unique(SR_data$pop)==i)
  
  X_OOS[, pop_OOS] <- OOS_pop$exp_spawn #populated!
  
  int_X_OOS <- ifelse(X_OOS!=0, 1,0) #binary matrix to index when to "turn on" intercept 
  
  stan_data <- list(logRS=logRS,
                    N=N,
                    pop=pop,
                    J=J,
                    K=K,
                    region=region, 
                    reg_pop=reg_pop,
                    X=X,
                    int_X=int_X,
                    SmaxPR=SmaxPR,
                    SmaxPR_SD=SmaxPR_SD, 
                    N_OOS = N_OOS, 
                    y_OOS = y_OOS,
                    X_OOS = X_OOS, 
                    int_X_OOS = int_X_OOS, 
                    pop_OOS = pop_OOS)
  
  #fit new models ---
  m1 <- stan(here("Stan/CV/m1_CV.stan"), data=stan_data, model_name = "m1", iter = iter)
  m2a <- stan(here("Stan/CV/m2a_CV.stan"), data=stan_data, model_name = "m2a", iter = iter)
  m2b <- stan(here("Stan/CV/m2b_CV.stan"), data=stan_data, model_name = "m2b", iter = iter)
  m3 <- stan(here("Stan/CV/m3_CV.stan"), data=stan_data, model_name = "m3", iter = iter)
  m4a <- stan(here("Stan/CV/m4a_CV.stan"), data=stan_data, model_name = "m4a", iter = iter)
  m4b <- stan(here("Stan/CV/m4b_CV.stan"), data=stan_data, model_name = "m4b", iter = iter)

  #pull the log_lik for the OOS preds
  m1_LL <- cbind(m1_LL, extract(m1)$log_lik)
  m2a_LL <- cbind(m2a_LL, extract(m2a)$log_lik)
  m2b_LL <- cbind(m2b_LL, extract(m2b)$log_lik)
  m3_LL <- cbind(m3_LL, extract(m3)$log_lik)
  m4a_LL <- cbind(m4a_LL, extract(m4a)$log_lik)
  m4b_LL <- cbind(m4b_LL, extract(m4b)$log_lik)
  
  write.csv(m1_LL, paste0(here("output/model_fits/CV/m1_LL_"),i))
  write.csv(m2a_LL, paste0(here("output/model_fits/CV/m2a_LL_"),i))
  write.csv(m2b_LL, paste0(here("output/model_fits/CV/m2b_LL_"),i))
  write.csv(m3_LL, paste0(here("output/model_fits/CV/m3_LL_"),i))
  write.csv(m4a_LL, paste0(here("output/model_fits/CV/m4a_LL_"),i))
  write.csv(m4b_LL, paste0(here("output/model_fits/CV/m4b_LL_"),i))
  
  end <- Sys.time()
  end - begin
}

#compute the loos
loo(m1_LL)$estimates[3,1]
loo(m2a_LL)$estimates[3,1]
loo(m2b_LL)$estimates[3,1]
loo(m3_LL)$estimates[3,1]
loo(m4a_LL)$estimates[3,1]
loo(m4b_LL)$estimates[3,1]
