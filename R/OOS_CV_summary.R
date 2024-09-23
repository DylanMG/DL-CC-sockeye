library(tidyverse)
library(here)
library(rstan)
library(loo)
#set.seed(123) #setting seed affects which pops are randomly selected, used these for analysis...
#set.seed(1)
#set.seed(69)
#set.seed(2)
set.seed(3)

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

### aggregate the loos -------------------------------------------------------------------
# code prior to this has worked with a single run and written LLs to here("ouput/model_fits/CV") 
  # by run, now we want to summarise multiple runs of OOS CV fits to capture the variability 
  # in which data was left out of sample 
# read in stored LLs from fits by model, where yeo is the last pop with all LLs bound to it
pop.index <- read.csv(here("output/CV/pop_index.csv")) |>
  pull() #hack to make as vector
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

write.csv(loos_table, here("output/CV/loos_table.csv"))

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
