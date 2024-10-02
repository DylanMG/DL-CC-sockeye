library(tidyverse)
library(here)
library(rstan)
library(loo)

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

### aggregate the LLs --------------------------------------------------------------------
# read in stored LLs from fits by model, where yeo is the last pop with all LLs bound to it
  #each object is all LLs from all runs of a given model
pop.index <- read.csv(here("output/CV/pop_index.csv")) |>
  pull() #hack to make vector

m1_LLs<- do.call(bind_cols, #expensive, but works; cbind wants unique names..  
                lapply(here("output/model_fits/CV", list.files(here("output/model_fits/CV"), 
                                                               pattern = "m1_LL_yeo.csv", 
                                                               recursive = TRUE)), read.csv)) |>
  select(-c(1, 994, 1988, 2982, 3976)) |> #janky hack to remove cols with rownames. ugh
  as.matrix() |>
  unname() #remove colnames

m2a_LLs<- do.call(bind_cols, 
                 lapply(here("output/model_fits/CV", list.files(here("output/model_fits/CV"), 
                                                                pattern = "m2a_LL_yeo.csv", 
                                                                recursive = TRUE)), read.csv))|>
  select(-c(1, 994, 1988, 2982, 3976)) |>
  as.matrix()|>
  unname()

m2b_LLs<- do.call(bind_cols, 
                 lapply(here("output/model_fits/CV", list.files(here("output/model_fits/CV"), 
                                                                pattern = "m2b_LL_yeo.csv", 
                                                                recursive = TRUE)), read.csv))|>
  select(-c(1, 994, 1988, 2982, 3976)) |>
  as.matrix()|>
  unname()

m3_LLs<- do.call(bind_cols, 
                 lapply(here("output/model_fits/CV", list.files(here("output/model_fits/CV"), 
                                                                pattern = "m3_LL_yeo.csv", 
                                                                recursive = TRUE)), read.csv))|>
  select(-c(1, 994, 1988, 2982, 3976)) |>
  as.matrix()|>
  unname()

m4a_LLs<- do.call(bind_cols, 
                 lapply(here("output/model_fits/CV", list.files(here("output/model_fits/CV"), 
                                                                pattern = "m4a_LL_yeo.csv", 
                                                                recursive = TRUE)), read.csv))|>
  select(-c(1, 994, 1988, 2982, 3976)) |>
  as.matrix()|>
  unname()

m4b_LLs<- do.call(bind_cols, 
                 lapply(here("output/model_fits/CV", list.files(here("output/model_fits/CV"), 
                                                                pattern = "m4b_LL_yeo.csv", 
                                                                recursive = TRUE)), read.csv))|>
  select(-c(1, 994, 1988, 2982, 3976)) |>
  as.matrix()|>
  unname()

# get elpds
elpd_all <- as.data.frame(rbind(elpd(m1_LLs)$estimates,
                                elpd(m2a_LLs)$estimates,
                                elpd(m2b_LLs)$estimates,
                                elpd(m3_LLs)$estimates,
                                elpd(m4a_LLs)$estimates,
                                elpd(m4b_LLs)$estimates))

# build matrix to feed stacking_weights() - see help file for N*K matrix
lpds <- matrix(NA, 993*5, 6)
for(i in 1:(993*5)){ #data.points*runs
  lpds[i, ] <- c(elpd(as.matrix(m1_LLs[,i]))$pointwise[1], 
                 elpd(as.matrix(m2a_LLs[,i]))$pointwise[1],
                 elpd(as.matrix(m2b_LLs[,i]))$pointwise[1],
                 elpd(as.matrix(m3_LLs[,i]))$pointwise[1],
                 elpd(as.matrix(m4a_LLs[,i]))$pointwise[1],
                 elpd(as.matrix(m4b_LLs[,i]))$pointwise[1])
}

#apply rowwise transformation to fix hard to predict data causing errors when trying in raw space
lpds_Z <- lpds - apply(lpds,1,max) # https://discourse.mc-stan.org/t/stacking-weights-optimization-error/28227/5

st_weights_df <- stacking_weights(lpds_Z) |>
  as.matrix() |>
  as.data.frame()|>
  mutate(Model = c("m1", "m2a", "m2b", "m3", "m4a", "m4b")) |>
  rename(Weight = 1)

elpd_table <- elpd_all |>
  mutate(Metric = rep(c("elpd", "IC"), 6), 
         Model = c(rep("m1", 2), 
                   rep("m2a", 2), 
                   rep("m2b", 2), 
                   rep("m3", 2), 
                   rep("m4a", 2), 
                   rep("m4b", 2))) |>
  left_join(st_weights_df, by = "Model") |>
  select(Model, Metric, Estimate, SE, Weight) |>
  mutate(Weight = round(as.numeric(Weight), 4)) |>
  arrange(desc(Weight))
rownames(elpd_table) <- NULL

write.csv(elpd_table, here("output/CV/elpd_table.csv"))

### by population which performed best overall? ------------------------------------------
pop_weights <- NULL #and weights
pop_weights_hier <- NULL #weights between m2a and m4a to see w. and w/o prior
pop_weights_quo <- NULL #weights between top model and status quo (m1)

j <- 1 #column tracker
for(i in unique(SR_data$pop)){
  n.rows <- pop.index[which(unique(SR_data$pop)== i)] #cols of OOS; same among runs
  #cols <- j:(j+n.cols-1) #account for runs being wide
  rows <- c(j:(j+n.rows-1), #subset cols for run 1
            (j+993):(j+993+n.rows-1), #for run 2 scoot it to start at new run
            (j+993*2):(j+993*2+n.rows-1), #repeat till 5...
            (j+993*3):(j+993*3+n.rows-1),
            (j+993*4):(j+993*4+n.rows-1))
  
  #calc pop's weights among all models ---
  pop_weight <- stacking_weights(lpds_Z[rows,]) |>
    as.matrix() |>
    as.data.frame()|>
    mutate(Model = c("m1", "m2a", "m2b", "m3", "m4a", "m4b"), 
           Pop = i) |>
    rename(Weight = 1)
  
  pop_weights <- rbind(pop_weights, pop_weight)
  
  #calc weights between top(4b) and quo(m1) ---
  lpds_Z_quo <- lpds[,c(1,6)] - apply(lpds[,c(1,6)],1,max) #scale within models we compare
    
  pop_weight_quo <- stacking_weights(lpds_Z_quo[rows,]) |>
    as.matrix() |>
    as.data.frame()|>
    mutate(Model = c("m1", "m4b"), 
           Pop = i) |>
    rename(Weight = 1)
  pop_weights_quo <- rbind(pop_weights_quo, pop_weight_quo)
  
  #and weights with(4b) and without(2b) informative prior ---
  lpds_Z_hier <- lpds[,c(4,6)] - apply(lpds[,c(4,6)],1,max) 
  
  pop_weight_hier <- stacking_weights(lpds_Z_hier[rows,]) |>
    as.matrix() |>
    as.data.frame()|>
    mutate(Model = c("m2b", "m4b"), 
           Pop = i) |>
    rename(Weight = 1)
  pop_weights_hier <- rbind(pop_weights_hier, pop_weight_hier)
  
  j <- j+n.rows #advance index to start of next pop 
}

write.csv(pop_weights, here("output/CV/pop_weights.csv"), row.names = FALSE)
write.csv(pop_weights_hier, here("output/CV/pop_weights_hier.csv"), row.names = FALSE)
write.csv(pop_weights_quo, here("output/CV/pop_weights_quo.csv"), row.names = FALSE)
