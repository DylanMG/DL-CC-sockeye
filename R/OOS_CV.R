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
for(i in unique(SR_data$pop)){
  #could add "for j in 1:n_reps" or something to repeat sampling
  thinned_pop <- filter(SR_data, pop == i) |> #reduce 1 pop's data by 50%
    sample_frac(0.5)
  
  sub_data <- filter(SR_data, pop !=i) |> #write reduced pop back in with rest of data
    bind_rows(thinned_pop) |>
    arrange(pop, year)
  
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
  
  for(i in 1:N){
    X[i, pop[i]] <- sub_data$exp_spawn[i] #populated!
  }
  int_X <- ifelse(X!=0, 1,0) #binary matrix to index when to "turn on" intercept 
  SmaxPR <- as.numeric(unique(sub_data$SmaxPR)) #SmaxPR prior mus
  SmaxPR_SD <- sub_data |> #SmaxPR SDs
    select(pop, SmaxPR_SD) |>
    distinct() |>
    pull(SmaxPR_SD) |> #got creative because a SD repeats
    as.numeric() 
  
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
                    SmaxPR_SD=SmaxPR_SD)
  
  #fit new models ---
  m1 <- stan(here("Stan/model1.stan"), data=stan_data, model_name = "m1", iter = iter)
  m2a <- stan(here("Stan/model2a.stan"), data=stan_data, model_name = "m2a", iter = iter)
  m2b <- stan(here("Stan/model2b.stan"), data=stan_data, model_name = "m2b", iter = iter)
  m3 <- stan(here("Stan/model3.stan"), data=stan_data, model_name = "m3", iter = iter)
  m4a <- stan(here("Stan/model4a.stan"), data=stan_data, model_name = "m4a", iter = iter)
  m4b <- stan(here("Stan/model4b.stan"), data=stan_data, model_name = "m4b", iter = iter)

  #get log_lik(?) from models we fit ---
  m1.s <- as.data.frame(rstan::summary(m1)$summary) #get log lik fo the OOS or something? 
}