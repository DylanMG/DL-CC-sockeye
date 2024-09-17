library(tidyverse)
library(here)
library(rstan)
library(loo)
#set.seed(123) #setting seed affects which pops are randomly selected, used 1, 123, n_pops(69) in analysis
set.seed(2)

# wrangle and load data ------------------------------------------------------------------
SmaxPR <- read.table(here("output/data/SmaxPRs.txt"), 
                     header = TRUE) |>
  select(population, region, water_clarity, SmaxPR, SmaxPR_SD) |>
  dplyr::rename(pop = population) #to match below

SR_data <- read.table(here("data/FinalBroodSK.txt"), 
                      header = TRUE) |>
  select(-n.data, -n.stock) |>
  mutate(watershed = gsub("_lake|_lakes", "", watershed)) |>
  mutate(watershed = case_when(watershed == "owekino" ~ "owikeno",   #fix some names so it joins
                               watershed == "portjohn" ~ "port_john",
                               watershed == "freeda_brodie" ~ "freeda_lakes",
                               TRUE ~ watershed)) |>
  mutate(logRS = log(total_rec/exp_spawn)) |>
  dplyr::rename(pop = watershed) |>
  relocate(pop, 1) |>
  filter(pop !="end_hill") |>
  arrange(pop, year) |>
  left_join(SmaxPR, by = "pop") |> #join it to have 1 object to reference
  mutate(pop = gsub("_", "-", pop)) |> #change delimiter for use later
  na.omit() |>
  group_by(pop) |>
  mutate(exp_spawn_scaled = scale(exp_spawn))

#loop pops, fit models for a single run --------------------------------------------------
new_run <- FALSE #toggle to do run. make sure you set.seed() different if so

if(new_run == FALSE){pop.index <- read.csv(here("output/CV/pop_index.csv")) |>
  pull(x)} #need this for later if not doing new run
if(new_run == TRUE){
  iter <- 4000 #how many iterations in a single stan model
  m1_LL <- NULL #empty objects to cbind LLs to
  m2a_LL <- NULL
  m2b_LL <- NULL
  m3_LL <- NULL
  m4a_LL <- NULL
  m4b_LL <- NULL
  
  track_OOS <- NULL #to track which data was pulled out for OOS
  for(i in unique(SR_data$pop)){
    begin <- Sys.time()
    sub_pop <- filter(SR_data, pop == i)
    
    thinned_pop <- sub_pop |> #reduce 1 pop's data by 50%
      sample_frac(0.5)
    
    sub_data <- filter(SR_data, pop !=i) |> #write reduced pop back in with rest of data
      bind_rows(thinned_pop) |>
      arrange(pop, year)
    
    OOS_pop <- anti_join(sub_pop, thinned_pop) #the data that was removed; to predict OOS
    track_OOS <- rbind(track_OOS, OOS_pop) #keep track of OOS obs by run

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
    
    #write likelihoods by pop*model
    write.csv(m1_LL, paste0(here("output/model_fits/CV/run 3/m1_LL_"),i,".csv"))
    write.csv(m2a_LL, paste0(here("output/model_fits/CV/run 3/m2a_LL_"),i,".csv"))
    write.csv(m2b_LL, paste0(here("output/model_fits/CV/run 3/m2b_LL_"),i,".csv"))
    write.csv(m3_LL, paste0(here("output/model_fits/CV/run 3/m3_LL_"),i,".csv"))
    write.csv(m4a_LL, paste0(here("output/model_fits/CV/run 3/m4a_LL_"),i,".csv"))
    write.csv(m4b_LL, paste0(here("output/model_fits/CV/run 3/m4b_LL_"),i,".csv"))
    
    end <- Sys.time()
    print(end - begin)
  
  #compute the loos
  loos <- as.data.frame(rbind(loo(m1_LL)$estimates,
                              loo(m2a_LL)$estimates,
                              loo(m2b_LL)$estimates,
                              loo(m3_LL)$estimates,
                              loo(m4a_LL)$estimates,
                              loo(m4b_LL)$estimates))
  loos_table <- loos |>
    mutate(Metric = rep(c("elpd_loo", "p_loo", "looic"), 6), 
           Model = c(rep("m1", 3), 
                     rep("m2a", 3), 
                     rep("m2b", 3), 
                     rep("m3", 3), 
                     rep("m4a", 3), 
                     rep("m4b", 3))) |>
    select(Model, Metric, Estimate, SE) |>
    arrange(Metric, Estimate)
  rownames(loos_table) <- NULL

  #need to get looic by pop, so make a single, long df with c(pop, model, metric, estimate, CV)
  pop.index <- as.data.frame(track_OOS) |> #how many OOS data for each pop? - stays same each run
    group_by(pop) |>
    summarise(n()) |>
    pull()

  pop_loos <- NULL #long object to populate with pop specific loos
  pop_weights <- NULL #and weights
  j <- 1 #column tracker
  for(i in unique(SR_data$pop)){
    n.cols <- pop.index[which(unique(SR_data$pop)== i)] #cols of OOS
    cols <- j:(j+n.cols-1)
    
    #calc pop's loos
    loos <- as.data.frame(rbind(loo(m1_LL[,cols])$estimates,
                                loo(m2a_LL[,cols])$estimates,
                                loo(m2b_LL[,cols])$estimates,
                                loo(m3_LL[,cols])$estimates,
                                loo(m4a_LL[,cols])$estimates,
                                loo(m4b_LL[,cols])$estimates)) |>
      mutate(Metric = rep(c("elpd_loo", "p_loo", "looic"), 6), 
             Model = c(rep("m1", 3),
                       rep("m2a", 3), 
                       rep("m2b", 3), 
                       rep("m3", 3), 
                       rep("m4a", 3), 
                       rep("m4b", 3)), 
             Pop = rep(i, 18)) |>
      select(Pop, Model, Metric, Estimate, SE)
    pop_loos <- rbind(pop_loos, loos)
    
    pop_weight <- as.data.frame(cbind(Pop = rep(i, 6),
                                      Model = c("m1", "m2a", "m2b", "m3", "m4a", "m4b"), 
                                      Weight = loo_model_weights(list(m1 = m1_LL[,cols], 
                                                                      m2a = m2a_LL[,cols], 
                                                                      m2b = m2b_LL[,cols],
                                                                      m3 = m3_LL[,cols],
                                                                      m4a = m4a_LL[,cols], 
                                                                      m4b = m4b_LL[,cols]))))
    pop_weights <- rbind(pop_weights, pop_weight)
  }
    j <- j+n.cols #advance index to start of next pop 
  }
}
#write run specific output summaries to folders. rename these depending on which is run
  #yes, looping could work better, but it takes a week-ish to do 3 runs... 
  # if worried that runs are the same compare them with all_equal() or something.
#write.csv(pop.index, here("output/CV/pop_index.csv"), row.names = FALSE) #stays same as long as new data not passed
#write.csv(loos_table, here("output/CV/loo_table_run3.csv"))
#write.csv(pop_loos, here("output/CV/pop_loos_run3.csv"))
#write.csv(pop_weights, here("output/CV/pop_weights_run3.csv"))

### aggregate the loos -------------------------------------------------------------------
# code prior to this has worked with a single run and written LLs to here("ouput/model_fits/CV") 
  # by run, now we want to summarise multiple runs of OOS CV fits to capture the variability 
  # in which data was left out of sample 
# read in stored LLs from fits by model, where yeo is the last pop with all LLs bound to it
m1_LLs<- do.call(rbind, 
                lapply(here("output/model_fits/CV", list.files(here("output/model_fits/CV"), 
                                                               pattern = "m1_LL_yeo.csv", 
                                                               recursive = T)), read.csv)) |>
  select(-X) |>
  as.matrix()

m2a_LLs<- do.call(rbind, 
                 lapply(here("output/model_fits/CV", list.files(here("output/model_fits/CV"), 
                                                                pattern = "m2a_LL_yeo.csv", 
                                                                recursive = T)), read.csv))|>
  select(-X) |>
  as.matrix()

m2b_LLs<- do.call(rbind, 
                 lapply(here("output/model_fits/CV", list.files(here("output/model_fits/CV"), 
                                                                pattern = "m2b_LL_yeo.csv", 
                                                                recursive = T)), read.csv))|>
  select(-X) |>
  as.matrix()

m3_LLs<- do.call(rbind, 
                 lapply(here("output/model_fits/CV", list.files(here("output/model_fits/CV"), 
                                                                pattern = "m3_LL_yeo.csv", 
                                                                recursive = T)), read.csv))|>
  select(-X) |>
  as.matrix()

m4a_LLs<- do.call(rbind, 
                 lapply(here("output/model_fits/CV", list.files(here("output/model_fits/CV"), 
                                                                pattern = "m4a_LL_yeo.csv", 
                                                                recursive = T)), read.csv))|>
  select(-X) |>
  as.matrix()

m4b_LLs<- do.call(rbind, 
                 lapply(here("output/model_fits/CV", list.files(here("output/model_fits/CV"), 
                                                                pattern = "m4b_LL_yeo.csv", 
                                                                recursive = T)), read.csv))|>
  select(-X) |>
  as.matrix()

# what was the best model overall?
#compute the loos
loos_all <- as.data.frame(rbind(loo(m1_LLs)$estimates,
                            loo(m2a_LLs)$estimates,
                            loo(m2b_LLs)$estimates,
                            loo(m3_LLs)$estimates,
                            loo(m4a_LLs)$estimates,
                            loo(m4b_LLs)$estimates))
#get looic weights
weights <- as.data.frame(cbind(Model= c("m1", "m2a", "m2b", "m3", "m4a", "m4b"),
                      Metric = rep("looic", 6),
                      Weight = loo_model_weights(list(m1 = m1_LLs, 
                                                       m2a = m2a_LLs,
                                                       m2b = m2b_LLs,
                                                       m3 = m3_LLs,
                                                       m4a = m4a_LLs, 
                                                       m4b = m4b_LLs))))
                       

loos_table <- loos_all |>
  mutate(Metric = rep(c("elpd_loo", "p_loo", "looic"), 6), 
         Model = c(rep("m1", 3), 
                   rep("m2a", 3), 
                   rep("m2b", 3), 
                   rep("m3", 3), 
                   rep("m4a", 3), 
                   rep("m4b", 3))) |>
  left_join(weights, by = c("Model", "Metric")) |>
  select(Model, Metric, Estimate, SE, Weight) |>
  mutate(Weight = round(as.numeric(Weight), 4)) |>
  arrange(Metric, Estimate)
rownames(loos_table) <- NULL

#write.csv(loos_table, here("output/CV/loos_table.csv"))

# by population which performed best overall? (add error among models?)
pop_loos <- NULL #long object to populate with pop specific loos
pop_weights <- NULL #and weights
pop_weights_hier <- NULL #weights between m2a and m4a to see w. and w/o prior
j <- 1 #column tracker
for(i in unique(SR_data$pop)){
  n.cols <- pop.index[which(unique(SR_data$pop)== i)] #cols of OOS; same among runs
  cols <- j:(j+n.cols-1)
  
  #calc pop's loos
  pop_loo <- as.data.frame(rbind(loo(m1_LLs[,cols])$estimates,
                              loo(m2a_LLs[,cols])$estimates,
                              loo(m2b_LLs[,cols])$estimates,
                              loo(m3_LLs[,cols])$estimates,
                              loo(m4a_LLs[,cols])$estimates,
                              loo(m4b_LLs[,cols])$estimates)) |>
    mutate(Metric = rep(c("elpd_loo", "p_loo", "looic"), 6), 
           Model = c(rep("m1", 3),
                     rep("m2a", 3), 
                     rep("m2b", 3), 
                     rep("m3", 3), 
                     rep("m4a", 3), 
                     rep("m4b", 3)), 
           Pop = rep(i, 18)) |>
    select(Pop, Model, Metric, Estimate, SE)
  pop_loos <- rbind(pop_loos, pop_loo)
  
  pop_weight <- as.data.frame(cbind(Pop = rep(i, 6),
                                    Model = c("m1", "m2a", "m2b", "m3", "m4a", "m4b"), 
                                    Weight = loo_model_weights(list(m1 = m1_LLs[,cols], 
                                                                    m2a = m2a_LLs[,cols], 
                                                                    m2b = m2b_LLs[,cols],
                                                                    m3 = m3_LLs[,cols],
                                                                    m4a = m4a_LLs[,cols], 
                                                                    m4b = m4b_LLs[,cols]))))
  pop_weights <- rbind(pop_weights, pop_weight)
  
  #and weights between m2a and m4a
  pop_weight_hier <- as.data.frame(cbind(Pop = rep(i, 2),
                                    Model = c("m2a", "m4a"), 
                                    Weight = loo_model_weights(list(m2a = m2a_LLs[,cols], 
                                                                    m4a = m4a_LLs[,cols]))))
  pop_weights_hier <- rbind(pop_weights_hier, pop_weight_hier)
  
  j <- j+n.cols #advance index to start of next pop 
}

write.csv(pop_loos, here("output/CV/pop_loos.csv"), row.names = FALSE)
write.csv(pop_weights, here("output/CV/pop_weights.csv"), row.names = FALSE)
write.csv(pop_weights_hier, here("output/CV/pop_weights_hier.csv"), row.names = FALSE)