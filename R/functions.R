#get draws and calculate Smsy from a single model or list of models 
get_draws <-function(x, y){
  draws <- map(x, 
      ~cbind(rstan:::.make_plot_data(.x, 
                                     pars = c("log_a_pop", "b_pop"))$samp, 
             model = .x@model_name)) |>
    bind_rows() |>
    mutate(pop_index = as.numeric(str_extract(parameter, "[[:digit:]]+"))) |>
    left_join(select(y, pop, pop_index, region), 
              by = "pop_index") |>
    mutate(parm = case_when(grepl("log_a", parameter) ~ "log_a", 
                            grepl("b", parameter) ~ "b")) |>
    select(-parameter) |>
    pivot_wider(names_from = parm, 
                values_from = value) |>
    select(model, pop, region, log_a, b) |>
    mutate(Smsy = (1-lambert_W0(exp(1-log_a)))/b)
  return(draws)
} 

#get a table of looic and weights for table 1
get_loo_table <- function(x){
  looic <- map(x, ~round(loo(.x)$estimates[3,1], 2)) |>
    unlist() |>
    cbind(map(x, ~.x@model_name) |>
            unlist())
  
  loos <- map(x, ~loo(.x)) #keen to learn a more eloquent way & not run loo 2x.
    
  
  weights <- loo_model_weights(loos) 
  
  loo_table <- cbind(looic, round(as.numeric(weights), 6)) |>
    as.data.frame() |>
    rename(looic = "V1", 
           model = "V2",  
           weight = "V3") |>
    mutate(weight = as.numeric(weight)*100) |>
    relocate(model, .before = looic) |>
    arrange(desc(weight))
  
  rownames(loo_table) <- NULL
  
  return(loo_table)
}

#get percent difference between 2 things 
perc_diff <- function(a, b){
  diff <- round(((a-b)/a)*100, 2)
  return(diff)
}

#brendan wanted this another way for draw diff from draws and median
perc_diff2 <- function(a, b){
  diff <- round(((b-a)/a)*100, 2)
  return(diff)
}
