#thinning code from pervious exercise
#And also look at $S_{MSY}$ varied between the top 2 models. This is some data we'll use for figure 1.  
fig1_data <- all_draws |>
  filter(model %in% c("m2b", "m4a")) |> #THIS is what you change to compare different models 
  group_by(model, pop) |>
  summarise(CV = round(sd(Smsy)/mean(Smsy), 2)) |>
  pivot_wider(names_from = model, 
              values_from = CV) |>
  mutate(perc_diff_CV = perc_diff(m2b, m4a), 
         diff_CV = m2b - m4a, 
         diff_CV_percentage = diff_CV*100) |>
  left_join(select(pop_data, pop, region, est_origin, water_clarity, years, sr_pairs),
            by = "pop")

## Thinning exercise  
#We'll also take different amounts of data from Meziadin, which is a heavily monitored 
#wild population that has high #quality age and escapement data.
if(FALSE){
  thin <- c(.2, .4, .6, .8, 1)
  
  thinning_pops <- "meziadin"  #set up robustly in case one wants multiple pops
  
  thinning_data <- vector(mode="list", length(thin)*length(thinning_pops))
  n <- NULL #empty obj to populate w/ n_data by percent for plotting later
  k <- 1
  for(i in thinning_pops){
    for(j in thin){
      
      other_data <- SR_data |>
        filter(pop != i)
      
      thinned_data <- SR_data |>
        filter(pop==i) |>
        slice_sample(prop = j)
      
      thin_df <-  bind_rows(other_data, thinned_data) |>
        arrange(pop, year) 
      
      n_data <- thinned_data |>
        summarise(n_data = as.character(n())) |>
        mutate(pop = str_to_title(i), 
               thin = gsub(".", "-", j, fixed = TRUE))
      
      n <- rbind(n, n_data)
      
      t_logRS <- thin_df$logRS #response variable
      t_pop <- as.numeric(as.factor(thin_df$pop)) #pop index vector
      t_N <- length(t_logRS) #n observations
      t_J <- max(t_pop) #n populations
      t_K <- max(as.numeric(as.factor(thin_df$region))) #n regions
      t_region <- as.numeric(as.factor(thin_df$region)) #region index vector
      t_reg_pop <- thin_df |> #region-pop key for indexing
        select(pop, region) |>
        distinct() |>
        pull(region) |>
        as.factor() |>
        as.numeric()
      
      t_X <- matrix(rep(0, (t_N*t_J)), 
                    nrow=t_N, 
                    ncol=t_J) #empty predictor matrix to populate
      
      for(l in 1:t_N){
        t_X[l, t_pop[l]] <- thin_df$exp_spawn[l] #populated!
      }
      t_int_X <- ifelse(t_X!=0, 1,0) #binary matrix to index when to "turn on" intercept 
      t_SmaxPR <- as.numeric(unique(thin_df$SmaxPR)) #SmaxPR prior mus
      t_SmaxPR_SD <- thin_df |> #SmaxPR SDs
        select(pop, SmaxPR_SD) |>
        distinct() |>
        pull(SmaxPR_SD) |> #got creative because a SD repeats
        as.numeric() 
      
      thin_list <- list(logRS=t_logRS,
                        N=t_N,
                        pop=t_pop,
                        J=t_J,
                        K=t_K,
                        region=t_region, 
                        reg_pop=t_reg_pop,
                        X=t_X,
                        int_X=t_int_X,
                        SmaxPR=t_SmaxPR,
                        SmaxPR_SD=t_SmaxPR_SD)
      
      thinning_data[[k]] <- thin_list
      
      k <- k + 1
    }
  }
  
  rm(thin_df, thinned_data, n_data, other_data, thin_list, t_logRS, t_N, t_pop, t_J, t_K, 
     t_region, t_reg_pop, t_X, t_int_X, t_SmaxPR, t_SmaxPR_SD)
  
  #model fitting----------------------------------------------------------------------------
  
  #swap delimiter in model names so we can write it as file w/o "." in the name
  thinned_models <- purrr::map(thinning_pops, ~{
    paste0(.x, "_", gsub( ".", "-", thin, fixed = TRUE))}) |>
    unlist()
  
  j <- 1 
  thinned_fits <- list()
  
  for(i in thinned_models){
    
    if(paste0("m1_", i) %in% list.files(here("output/model_fits/thinning")) & refit == FALSE){
      m1_thin <- readRDS(here(paste0("output/model_fits/thinning/m1_", i)))
    } else{
      m1_thin <- stan(file=here("Stan/model1.stan"),
                      data=thinning_data[[j]], 
                      iter = 4000)
      saveRDS(m1_thin, file = here(paste0("output/model_fits/thinning/m1_", i)))
    }
    
    if(paste0("m2a_", i) %in% list.files(here("output/model_fits/thinning")) & refit == FALSE){
      m2a_thin <- readRDS(here(paste0("output/model_fits/thinning/m2a_", i)))
    } else{
      m2a_thin <- stan(file=here("Stan/model2a.stan"),
                       data=thinning_data[[j]], 
                       iter = 4000)
      saveRDS(m2a_thin, file = here(paste0("output/model_fits/thinning/m2a_", i)))
    }
    
    if(paste0("m2b_", i) %in% list.files(here("output/model_fits/thinning")) & refit == FALSE){
      m2b_thin <- readRDS(here(paste0("output/model_fits/thinning/m2b_", i)))
    } else{
      m2b_thin <- stan(file=here("Stan/model2b.stan"),
                       data=thinning_data[[j]], 
                       iter = 4000)
      saveRDS(m2b_thin, file = here(paste0("output/model_fits/thinning/m2b_", i)))
    }
    
    if(paste0("m3_", i) %in% list.files(here("output/model_fits/thinning")) & refit == FALSE){
      m3_thin <- readRDS(here(paste0("output/model_fits/thinning/m3_", i)))
    } else{
      m3_thin <- stan(file=here("Stan/model3.stan"),
                      data=thinning_data[[j]], 
                      iter = 4000)
      saveRDS(m3_thin, file = here(paste0("output/model_fits/thinning/m3_", i)))
    }
    
    if(paste0("m4a_", i) %in% list.files(here("output/model_fits/thinning")) & refit == FALSE){
      m4a_thin <- readRDS(here(paste0("output/model_fits/thinning/m4a_", i)))
    } else{
      m4a_thin <- stan(file=here("Stan/model4a.stan"),
                       data=thinning_data[[j]], 
                       iter = 4000)
      saveRDS(m4a_thin, file = here(paste0("output/model_fits/thinning/m4a_", i)))
    }
    
    if(paste0("m4b_", i) %in% list.files(here("output/model_fits/thinning")) & refit == FALSE){
      m4b_thin <- readRDS(here(paste0("output/model_fits/thinning/m4b_", i)))
    } else{
      m4b_thin <- stan(file=here("Stan/model4b.stan"),
                       data=thinning_data[[j]], 
                       iter = 4000)
      saveRDS(m4b_thin, file = here(paste0("output/model_fits/thinning/m4b_", i)))
    }
    
    thinned_fits[[j]] <- m1_thin
    thinned_fits[[j+length(thinned_models)]] <- m2a_thin #add the total # of models for each to index
    thinned_fits[[j+length(thinned_models)*2]] <- m2b_thin
    thinned_fits[[j+length(thinned_models)*3]] <- m3_thin
    thinned_fits[[j+length(thinned_models)*4]] <- m4a_thin
    thinned_fits[[j+length(thinned_models)*5]] <- m4b_thin
    
    
    j <- j + 1
    
  }
  
  #get draws from thinned models------------------------------------------------------------
  
  #computationally expensive: be patient
  thinned_draws <- map2(thinned_fits, 
                        c(paste0("m1_", thinned_models),
                          paste0("m2a_", thinned_models),
                          paste0("m2b_", thinned_models),
                          paste0("m3_", thinned_models), 
                          paste0("m4a_", thinned_models), 
                          paste0("m4b_", thinned_models)),  
                        ~cbind(rstan:::.make_plot_data(.x, 
                                                       pars = c("log_a_pop", "b_pop"))$samp, 
                               model_name = rep(.y)))  |> #instead of get_draws() because of model naming issue*
    bind_rows()|>
    filter(chain == 1) |>
    mutate(pop_index = as.numeric(str_extract(parameter, "[[:digit:]]+"))) |>
    left_join(select(pop_data, pop, pop_index, region), 
              by = "pop_index") |>
    filter(pop %in% str_to_title(thinning_pops)) |> #still fluff here since pops won't match model_name
    mutate(parameter = ifelse(grepl("log_a", parameter), "log_a", "b")) |>
    tidyr::pivot_wider(names_from = parameter, 
                       values_from = value) |>
    separate(model_name, c("model", "pop_reduced", "thin"), sep = "_") |>
    filter(pop_reduced == tolower(pop)) |> #only retain fits from the pop that was reduced
    select(pop, model, thin, region, log_a, b) |>
    mutate(Smsy = (1-lambert_W0(exp(1-log_a)))/b) |>
    left_join(n, by = c("pop", "thin"))
}
  #*model naming issue - once a stan model is compiled it retains the same name; if you loop 
  #a stan model it'll keep the name of the first time the model was compiled - annoying... 
  
  #then get the smsy bias data--------------------------------------------------------------
  m4_median <- thinned_draws |>
    filter(model == "m4a", thin == 1) |>
    group_by(pop) |>
    summarise(Smsy_med = median(Smsy))
  
  Smsy_bias <- thinned_draws |>
    left_join(m4_median, by = "pop") |>
    mutate(perc_diff_Smsy = perc_diff2(log(Smsy_med), log(Smsy)), 
           pop = factor(pop, levels = thinning_pops),
           n_data = factor(n_data, levels = unique(thinned_draws$n_data)),
           thin = gsub("-", ".", thin, fixed = TRUE))
  
  write.csv(Smsy_bias, here("output/data/Smsy_bias.csv"))
