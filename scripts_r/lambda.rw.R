library(here)
library(rjags)
library(tidybayes)
library(tidyverse)
library(tidylog)


counts <- read.csv(here::here("data/clean/counts.csv")) %>%
  arrange(herd)
trt <- read.csv(here::here("data/clean/treatments.csv")) %>%
  arrange(herd)


covar.names <- readRDS(here::here("data/clean/forIPM/covariate_names.rds"))
covariates <- read_csv(here::here("data/clean/forIPM/covariates_scaled.csv")) %>%
  dplyr::select(herd, year, !!paste0(covar.names,"_scale"))

#### TREATMENT COMBOS
treatment.combos <- trt %>%
  filter(applied == 1) %>%
  group_by(herd, year) %>%
  summarize(
    trt = paste(
      gsub(" ", "", treatment), # Remove spaces
      collapse = "-" # Add - between treatments
    ),
    .groups = "drop"
  ) %>%
  dplyr::select(herd, year, trt)




##prep data
count.ipm <-counts%>%
  left_join(treatment.combos, by=c("herd","year"))%>%
  mutate(trt = replace_na(trt, "Reference")) %>%
  mutate(
    Est_CL_sd=(((Est_CL-Est_CL.min) + (Est_CL.max-Est_CL))/2)/1.96, ## take average CI bound, go to se via /1.96
    tau=1/(Est_CL_sd^2),
    tau2=1/((Est_CL*0.05)^2),
    yr_id=year-min(year)+1)%>%
  filter(Est_CL>10, Est_CL_sd>0, trt=="Reference")


hd <- tibble(herd=unique(count.ipm$herd),
             herd_num=1:length(herd))

count.ipm <-count.ipm%>%
left_join(hd, by="herd")%>%
  arrange(herd_num,yr_id)

yrs_df <- tibble(year=min(count.ipm$year):max(count.ipm$year),
                 yr_id=1:length(year))
nyr <- length(yrs_df$year)
  
obs <-count.ipm%>%
  select(herd_num, yr_id, Est_CL, tau, tau2)

nobs <- nrow(obs)
nherd<- max(count.ipm$herd_num)


covariates.ipm <-covariates%>%
  left_join(yrs_df, by="year")%>%
  left_join(hd, by="herd")%>%
  drop_na(yr_id, herd_num)%>%
  arrange(herd_num,yr_id)%>%
  rename_at(.vars = vars(ends_with("_scale")),
            .funs = funs(sub("[_]scale$", "", .)))

dist <- covariates.ipm%>%
  select(herd_num, yr_id, dist.hum)%>%
  pivot_wider(names_from=yr_id, values_from=dist.hum)
dist <- dist[,-1]

snow.depth <- covariates.ipm%>%
  select(herd_num, yr_id, snow.depth)%>%
  pivot_wider(names_from=yr_id, values_from=snow.depth)
snow.depth <- snow.depth[,-1]

start.spring <- covariates.ipm%>%
  select(herd_num, yr_id, start.spring)%>%
  pivot_wider(names_from=yr_id, values_from=start.spring)
start.spring <- start.spring[,-1]

warm.dry.summer <- covariates.ipm%>%
  select(herd_num, yr_id, warm.dry.summer)%>%
  pivot_wider(names_from=yr_id, values_from=warm.dry.summer)
warm.dry.summer <- warm.dry.summer[,-1]

rain.on.snow <- covariates.ipm%>%
  select(herd_num, yr_id, rain.on.snow)%>%
  pivot_wider(names_from=yr_id, values_from=rain.on.snow)
rain.on.snow <- rain.on.snow[,-1]

n1 <- count.ipm%>%
  group_by(herd)%>%
  arrange(yr_id)%>%
  slice(1)%>%
  ungroup()%>%
  arrange(herd_num)%>%
  pull(Est_CL)

ipm.dat <- list(
  obs=obs,
  nobs=nobs,
  nherd=nherd,
  nyr=nyr,
  dist=dist,
  snow.depth=snow.depth,
  start.spring=start.spring,
  warm.dry.summer=warm.dry.summer,
  rain.on.snow=rain.on.snow,
  n1=n1
)



mod <-
  "model { 
  # Priors
  # Linear predictor on annual shock/growth
  mu_lambda ~ dnorm(1, 0.01)
  dist_eff ~ dnorm(0, 0.01)
  snow.depth_eff ~ dnorm(0, 0.01)
  start.spring_eff ~ dnorm(0, 0.01)
  warm.dry.summer_eff ~ dnorm(0, 0.01)
  rain.on.snow_eff ~ dnorm(0, 0.01)
  
  # Process and observation error
  tau_process ~ dgamma(0.01, 0.01)
  	
  	##starting vals year 1
for(herd in 1:nherd){	
 N[herd, 1] ~ dpois(n1[herd])
}

  # State model, random walk
  for(herd in 1:nherd){
  for(y in 2:nyr) {
    lambda_det[herd, y] <- mu_lambda + 
    dist_eff * dist[herd, y]+
    snow.depth_eff *snow.depth[herd, y]+
    start.spring_eff *start.spring[herd, y]+
    warm.dry.summer_eff *warm.dry.summer[herd, y]+
    rain.on.snow_eff *rain.on.snow[herd, y]
    
    lambda[herd, y] ~ dnorm(lambda_det[herd, y], tau_process)
    N[herd, y] <- N[herd, y-1] * lambda[herd, y]
  }
  }
  
  # Observation model
  for(i in 1:nobs) {
    obs[i,3] ~ dnorm(N[obs[i,1],obs[i,2]], obs[i,4])
    #obs[i,3] ~ dnorm(N[obs[i,1],obs[i,2]], obs[i,5]) #tau2
    #obs[i,3] ~ dnorm(N[obs[i,1],obs[i,2]], 10) #no error
    #N[h,yr]
  }
  
}"


nth <- 50
nbu <- 20000
nch <- 3
nad <- 20000
nit <- 100000

out <- jagsUI::jags(ipm.dat, 
                           inits =NULL,
                           parameters.to.save=c("dist_eff","snow.depth_eff","start.spring_eff","warm.dry.summer_eff","rain.on.snow_eff"),
                           model.file = textConnection(mod),
                           n.chains = nch,
                           n.cores = nch,
                           n.iter = nit,
                           n.burnin = nbu,
                           n.thin = nth,
                           n.adapt = nad)



out %>%
  gather_draws(dist_eff, snow.depth_eff, start.spring_eff, warm.dry.summer_eff, rain.on.snow_eff) %>%
  median_qi()


out$mean$dist_eff
out$sd$dist_eff
