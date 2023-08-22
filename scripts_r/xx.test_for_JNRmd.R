
## Load Packages
library(tidybayes)
library(ggmcmc)
library(hrbrthemes)
library(tidylog)
library(tidyverse)

# Load data ---------------------------------------------------------------
hd <- read.csv(here::here("data/clean/blueprint.csv"))
hn <- hd %>%
  dplyr::select(herd, herd_num)

##demographic data
cdat <- read_csv(here::here("data/clean/forIPM/cdat.csv"))
edat <- read_csv(here::here("data/clean/forIPM/edat.csv"))
sdat <- read_csv(here::here("data/clean/forIPM/surv.yr.ipm.csv"))
rdat <- read_csv(here::here("data/clean/forIPM/recruit.yr.ipm.csv"))
srdat <- read_csv(here::here("data/clean/forIPM/srdat.csv"))
pdat <- read_csv(here::here("data/clean/forIPM/pdat.csv"))
sexratio_summary <- read.csv(here::here("data/clean/sexratio_summary.csv"))

##covariates
dist.all <- readRDS(here::here("data/clean/forIPM/dist.all.rds"))
clim1 <- readRDS(here::here("data/clean/forIPM/clim1.rds"))
clim2 <- readRDS(here::here("data/clean/forIPM/clim2.rds"))

##other
#first_per_herd <- read_csv(here::here("data/clean/forIPM/first_per_herd.csv"))
month_offset <- read_csv(here::here("data/clean/forIPM/month_offset.csv"))
grp_p <- read_csv(here::here("data/clean/forIPM/grp_p.csv"))
single.params <- read_csv(here::here("data/clean/forIPM/single.params.csv"))
meancount<- read_csv(here::here("data/clean/forIPM/meancount.csv"))
count.otc<- read_csv(here::here("data/clean/forIPM/count.otc.csv"))
count.osc<- read_csv(here::here("data/clean/forIPM/count.osc.csv"))
count.mnka<- read_csv(here::here("data/clean/forIPM/count.mnka.csv"))
yr_df <- read_csv(here::here("data/clean/forIPM/yr_df.csv"))

# Estimate Priors ---------------------------------------------------------------

ipm.prior.dat <- list(
  ns=nrow(sdat),
  nr=nrow(rdat),
  sdat=sdat,
  rdat=rdat,
  dist = dist.all, 
  clim1 = clim1,
  clim2= clim2
)



prior.mod <-
  "model{ 
    # priors
    muS ~ dnorm(logit(0.9995),10)
    beta.dist.s ~ dunif(-10,10)
    beta.clim1.s ~ dunif(-10,10)
    beta.clim2.s ~ dunif(-10,10)
    
    muR ~ dnorm(logit(0.15),10)
    beta.dist.r ~ dunif(-10,10)
    beta.clim1.r ~ dunif(-10,10)
    beta.clim2.r ~ dunif(-10,10)
    beta.season.r ~ dunif(-10,10)
    beta.type.otc.r ~ dunif(-10,10)
    beta.type.osc.r ~ dunif(-10,10)
    beta.type.mnka.r ~ dunif(-10,10)
    
    # likelihood
    for (i in 1:ns){
            sdat[i,3] ~ dbin((1-s[i]), sdat[i,4])
            logit(s[i]) <- muS + 
            beta.dist.s*dist[sdat[i,1],sdat[i,2]] + beta.clim1.s*clim1[sdat[i,1],sdat[i,2]] + beta.clim2.s*clim2[sdat[i,1],sdat[i,2]]
    }
    
        for (i in 1:nr){
            rdat[i,3] ~ dbin(r[i], rdat[i,4])
            logit(r[i]) <- muR + 
            beta.dist.r*dist[rdat[i,1],rdat[i,2]] + beta.clim1.r*clim1[rdat[i,1],rdat[i,2]] + beta.clim2.r*clim2[rdat[i,1],rdat[i,2]] +
            beta.season.r*rdat[i,5]+
            beta.type.otc.r*rdat[i,6] + beta.type.osc.r*rdat[i,7] + beta.type.mnka.r*rdat[i,8]
            
    }

}"


nth <- 50
nbu <- 3000
nch <- 3
nad <- 3000
nit <- 50000

out.priors <- jagsUI::jags(ipm.prior.dat, 
                         inits =NULL,
                         parameters.to.save=c("s", "muS", "beta.dist.s", "beta.clim1.s", "beta.clim2.s", 
                                              "r", "muR", "beta.dist.r", "beta.clim1.r", "beta.clim2.r",  "beta.season.r"),
                         model.file = textConnection(prior.mod),
                         n.chains = nch,
                         n.cores = nch,
                         n.iter = nit,
                         n.burnin = nbu,
                         n.thin = nth,
                         n.adapt = nad)

# mcmcplots::mcmcplot(out.priors$samples, par = c("muR", "beta.clim1.r","beta.clim2.r","beta.dist.r",
#                                          "muS", "beta.clim1.s","beta.clim2.s", "beta.dist.s"))
                                         

priors <- out.priors%>%
  gather_draws(muS, beta.dist.s, beta.clim1.s, beta.clim2.s, 
               muR, beta.dist.r, beta.clim1.r, beta.clim2.r,beta.season.r)%>%
  group_by(.variable)%>%
  summarise(mean=mean(.value),
            sd=sd(.value),
            tau=(1/(sd^2)))

# Prep data for IPM ---------------------------------------------------------------

# Herds and herd number
nherd <- single.params$nherd
nsight_grp <- single.params$nsight_grp
nage <- single.params$nage

#  Years of study
nyr <- single.params$nyr

# sightability info
mean_grp_p <- grp_p$mean_grp_p
mean_grp_ptau <- grp_p$mean_grp_ptau

sgt_grp_ind <- matrix(nrow = nherd, ncol = nyr)
for(sg in 1:nrow(sgt_grp_ind)){
  sgt_grp_ind[sg,] <- hd$sight_grp[hd$herd_num == sg]
}


##starting age class distribution (sensitivity tested this in McNay et al. 2022 and it doesnt have a meaningful effect even if substantially changed)
n1s <- meancount$n1_mean
n1 <- matrix(NA, nrow = nherd, ncol = nage)
for(h in 1:nherd){
  n1[h,1] <- (n1s[h]*0.5)*0.15
  n1[h,2] <- (n1s[h]*0.64)*0.85 
}


meansr <- array(NA, c(1,2))
meansr[1,1] <- sexratio_summary[1,1]
meansr[1,2] <- 1/(sexratio_summary[1,2]^2)

##first year of data (legacy...set at 1 so all herds start at year 1 now)
first <- rep(1,nherd)


## Gather data inputs in a list

ipm_dat <- list(
  nherd = nherd,
  nyr = nyr,
  first = first,
  
  month_offset = month_offset,
  count.otc=count.otc,
  count.osc=count.osc,
  count.mnka=count.mnka,
  
  nc = single.params$nc, 
  ne = single.params$ne,
  ns = single.params$ns,
  nr = single.params$nr,
  nsr = single.params$nsr,
  
  nsight_grp = nsight_grp,
  sight_grp = hd$sight_grp,
  
  mean_grp_p = mean_grp_p,
  mean_grp_ptau = mean_grp_ptau,
  meansr = meansr,
  
  n1 = n1, 
  cdat = cdat,
  edat = edat,
  sdat = sdat,
  srdat = srdat,
  pdat = pdat,
  rdat = rdat,
  
  dist = dist.all, 
  clim1 = clim1,
  clim2= clim2,
  
  prior.dist.s.mean=priors%>%filter(.variable=="beta.dist.s")%>%pull(mean),
  prior.dist.s.tau=priors%>%filter(.variable=="beta.dist.s")%>%pull(tau),
  prior.clim1.s.mean=priors%>%filter(.variable=="beta.clim1.s")%>%pull(mean),
  prior.clim1.s.tau=priors%>%filter(.variable=="beta.clim1.s")%>%pull(tau),
  prior.clim2.s.mean=priors%>%filter(.variable=="beta.clim2.s")%>%pull(mean),
  prior.clim2.s.tau=priors%>%filter(.variable=="beta.clim2.s")%>%pull(tau),
  
  prior.dist.r.mean=priors%>%filter(.variable=="beta.dist.r")%>%pull(mean),
  prior.dist.r.tau=priors%>%filter(.variable=="beta.dist.r")%>%pull(tau),
  prior.clim1.r.mean=priors%>%filter(.variable=="beta.clim1.r")%>%pull(mean),
  prior.clim1.r.tau=priors%>%filter(.variable=="beta.clim1.r")%>%pull(tau),
  prior.clim2.r.mean=priors%>%filter(.variable=="beta.clim2.r")%>%pull(mean),
  prior.clim2.r.tau=priors%>%filter(.variable=="beta.clim2.r")%>%pull(tau)
)


#  Initial values for N to avoid parent node erros
#  Indexing was off here
Nst <- array(NA_integer_, c(nherd, nyr, nage))

for(h in 1:nherd) {
  for(a in 1:nage) {
    Nst[h, first[h]:nyr, a] <- as.integer(
      max(
        c(n1[h, a], cdat$mu[cdat$herd == h], edat$mu[edat$herd == h])))  
  }
}

ipm_inits <- function(){ 
  list(
    N = Nst
  )
}


#  Model parameters to monitor
model_parms <- c(
  
  "totNMF", "lambda",
  
  "S", "R",  "R_adj",
  
  "muS", "muR",
  
  "offset",
  "beta.dist.s",
  "beta.clim1.s",
  "beta.clim2.s",
  "beta.dist.r",
  "beta.clim1.r",
  "beta.clim2.r"
)


# Run IPM ---------------------------------------------------------------
nth <- 50
nbu <- 2000 
nch <- 3
nad <- 2000
nit <- 30000 

out <- jagsUI::jags(ipm_dat, 
                    inits = ipm_inits,
                    model_parms,
                    model.file = here::here("jags/Dist_Clim_IPM_rawVR.txt"),
                    n.chains = nch,
                    n.cores = nch,
                    n.iter = nit,
                    n.burnin = nbu,
                    n.thin = nth,
                    n.adapt = nad)


draws <- out%>%
  gather_draws(muS, beta.dist.s, beta.clim1.s, beta.clim2.s, 
               muR, beta.dist.r, beta.clim1.r, beta.clim2.r)%>%
  mutate(type="Posteriors")%>%
  rbind(out.priors%>%
          gather_draws(muS, beta.dist.s, beta.clim1.s, beta.clim2.s, 
                       muR, beta.dist.r, beta.clim1.r, beta.clim2.r)%>%
          mutate(type="Priors"))

draws%>%
  filter(!.variable%in%c("alpha.r","alpha.s","beta.season.r"))%>%
  ggplot(aes(x=.value, fill=type))+
  geom_density(alpha=0.5)+
  facet_wrap(vars(.variable), scales="free")+
  geom_vline(xintercept = 0, linetype="dashed")


# RUN again but with only S and R data ---------------------------------------------------------------
ipm_dat2 <- list(
  nherd = nherd,
  nyr = nyr,
  first = first,
  
  month_offset = month_offset,
  count.otc=count.otc,
  count.osc=count.osc,
  count.mnka=count.mnka,
  
  nc = 1, 
  ne = 1,
  ns = single.params$ns,
  nr = single.params$nr,
  nsr = 1,
  
  nsight_grp = nsight_grp,
  sight_grp = hd$sight_grp,
  
  mean_grp_p = mean_grp_p,
  mean_grp_ptau = mean_grp_ptau,
  meansr = meansr,
  
  n1 = n1, 
  cdat = cdat%>%slice(1),
  edat = edat%>%slice(1),
  sdat = sdat,
  srdat = srdat%>%slice(1),
  pdat = pdat%>%slice(1),
  rdat = rdat,
  
  dist = dist.all, 
  clim1 = clim1,
  clim2= clim2,
  
  prior.dist.s.mean=priors%>%filter(.variable=="beta.dist.s")%>%pull(mean),
  prior.dist.s.tau=priors%>%filter(.variable=="beta.dist.s")%>%pull(tau),
  prior.clim1.s.mean=priors%>%filter(.variable=="beta.clim1.s")%>%pull(mean),
  prior.clim1.s.tau=priors%>%filter(.variable=="beta.clim1.s")%>%pull(tau),
  prior.clim2.s.mean=priors%>%filter(.variable=="beta.clim2.s")%>%pull(mean),
  prior.clim2.s.tau=priors%>%filter(.variable=="beta.clim2.s")%>%pull(tau),
  
  prior.dist.r.mean=priors%>%filter(.variable=="beta.dist.r")%>%pull(mean),
  prior.dist.r.tau=priors%>%filter(.variable=="beta.dist.r")%>%pull(tau),
  prior.clim1.r.mean=priors%>%filter(.variable=="beta.clim1.r")%>%pull(mean),
  prior.clim1.r.tau=priors%>%filter(.variable=="beta.clim1.r")%>%pull(tau),
  prior.clim2.r.mean=priors%>%filter(.variable=="beta.clim2.r")%>%pull(mean),
  prior.clim2.r.tau=priors%>%filter(.variable=="beta.clim2.r")%>%pull(tau)
)

nth <- 50
nbu <- 2000 
nch <- 3
nad <- 2000
nit <- 30000 

out2 <- jagsUI::jags(ipm_dat2, 
                     inits = ipm_inits,
                     model_parms,
                     model.file = here::here("jags/Dist_Clim_IPM_rawVR.txt"),
                     n.chains = nch,
                     n.cores = nch,
                     n.iter = nit,
                     n.burnin = nbu,
                     n.thin = nth,
                     n.adapt = nad)


draws2 <- out2%>%
  gather_draws(muS,beta.dist.s, beta.clim1.s, beta.clim2.s, 
               muR,beta.dist.r, beta.clim1.r, beta.clim2.r)%>%
  mutate(type="Posteriors-S+R")%>%
  rbind(draws)

draws2%>%
  filter(!.variable%in%c("alpha.r","alpha.s","beta.season.r"))%>%
  ggplot(aes(x=.value, fill=type))+
  geom_density(alpha=0.5)+
  facet_wrap(vars(.variable), scales="free")+
  geom_vline(xintercept = 0, linetype="dashed")


# RUN again but without data when herds at small pop sizes ---------------------------------------------------------------

last <- read_csv(here::here("data","clean","demog.csv"))%>%
  filter(totNMF<20)%>%
  group_by(herd,i)%>%
  summarise(yrs=min(yrs))%>%
  left_join(yr_df, by="yrs")%>%
  ungroup%>%
  select(i,yr_idx)%>%
  mutate(yr_idx=replace_na(yr_idx,1))%>%
  rbind(tibble(i=1:nherd,yr_idx=nyr))%>%
  group_by(i)%>%
  summarise(yrs=min(yr_idx))%>%
  ungroup%>%
  arrange(i)%>%
  rename(herd=i)

ipm_dat3 <- list(
  nherd = nherd,
  nyr = nyr,
  first = first,
  last=last$yrs,
  
  month_offset = month_offset,
  count.otc=count.otc,
  count.osc=count.osc,
  count.mnka=count.mnka,
  
  nc = nrow(cdat%>%left_join(last, by="herd")%>%filter(year<=yrs)), 
  ne = nrow(edat%>%left_join(last, by="herd")%>%filter(year<=yrs)),
  ns = nrow(edat%>%left_join(last, by="herd")%>%filter(year<=yrs)),
  nr = nrow(rdat%>%left_join(last, by="herd")%>%filter(year<=yrs)),
  nsr = nrow(srdat%>%left_join(last, by="herd")%>%filter(year<=yrs)),
  
  nsight_grp = nsight_grp,
  sight_grp = hd$sight_grp,
  
  mean_grp_p = mean_grp_p,
  mean_grp_ptau = mean_grp_ptau,
  meansr = meansr,
  
  n1 = n1, 
  cdat = cdat%>%left_join(last, by="herd")%>%filter(year<=yrs),
  edat = edat%>%left_join(last, by="herd")%>%filter(year<=yrs),
  sdat = sdat%>%left_join(last, by="herd")%>%filter(year<=yrs),
  srdat = srdat%>%left_join(last, by="herd")%>%filter(year<=yrs),
  pdat = pdat%>%left_join(last, by="herd")%>%filter(year<=yrs),
  rdat = rdat%>%left_join(last, by="herd")%>%filter(year<=yrs),
  
  dist = dist.all, 
  clim1 = clim1,
  clim2= clim2,
  
  prior.dist.s.mean=priors%>%filter(.variable=="beta.dist.s")%>%pull(mean),
  prior.dist.s.tau=priors%>%filter(.variable=="beta.dist.s")%>%pull(tau),
  prior.clim1.s.mean=priors%>%filter(.variable=="beta.clim1.s")%>%pull(mean),
  prior.clim1.s.tau=priors%>%filter(.variable=="beta.clim1.s")%>%pull(tau),
  prior.clim2.s.mean=priors%>%filter(.variable=="beta.clim2.s")%>%pull(mean),
  prior.clim2.s.tau=priors%>%filter(.variable=="beta.clim2.s")%>%pull(tau),
  
  prior.dist.r.mean=priors%>%filter(.variable=="beta.dist.r")%>%pull(mean),
  prior.dist.r.tau=priors%>%filter(.variable=="beta.dist.r")%>%pull(tau),
  prior.clim1.r.mean=priors%>%filter(.variable=="beta.clim1.r")%>%pull(mean),
  prior.clim1.r.tau=priors%>%filter(.variable=="beta.clim1.r")%>%pull(tau),
  prior.clim2.r.mean=priors%>%filter(.variable=="beta.clim2.r")%>%pull(mean),
  prior.clim2.r.tau=priors%>%filter(.variable=="beta.clim2.r")%>%pull(tau)
)

nth <- 50
nbu <- 2000 
nch <- 3
nad <- 2000
nit <- 30000 

out3 <- jagsUI::jags(ipm_dat3, 
                     inits = ipm_inits,
                     model_parms,
                     model.file = here::here("jags/Dist_Clim_IPM_rawVR.txt"),
                     n.chains = nch,
                     n.cores = nch,
                     n.iter = nit,
                     n.burnin = nbu,
                     n.thin = nth,
                     n.adapt = nad)

#saveRDS(out5, file = here::here("jags/output/dist_clim_tests_nosmall.rds"))

draws3 <- out3%>%
  gather_draws(muS,beta.dist.s, beta.clim1.s, beta.clim2.s, 
               muR,beta.dist.r, beta.clim1.r, beta.clim2.r)%>%
  mutate(type="Posteriors-nosmall")%>%
  rbind(draws2)


draws3%>%
  filter(!.variable%in%c("beta.season.r","muS","muR"))%>%
  group_by(type,.variable)%>%
  summarise(mean=mean(.value),
            sd=sd(.value))%>%
  ggplot(aes(x=mean, xmin=mean-sd, xmax=mean+sd, y=.variable, color=type))+
  geom_point(alpha=0.5, position=position_dodge(0.4))+
  geom_linerange(alpha=0.5, position=position_dodge(0.4))+
  geom_vline(xintercept = 0, linetype="dashed")



# JOSH CAN IGNORE HERE AND BELOW, JUST ME SLEUTHING ---------------------------------------------------------------



# Something is buggered. As soon as I add N data, disturbance is + with S, why? ---------------------------------------------------------------
lambda <- out$mean$lambda%>%
  as_tibble()%>%
  mutate(herd=rownames(.))%>%
  pivot_longer(-herd)%>%
  rename(year=name,
         lambda=value)%>%
  mutate(year=parse_number(year))

abund <- out$mean$totNMF%>%
  as_tibble()%>%
  mutate(herd=rownames(.))%>%
  pivot_longer(-herd)%>%
  rename(year=name,
         totNMF=value)%>%
  mutate(year=parse_number(year))

S <- out$mean$S%>%
  as_tibble()%>%
  mutate(herd=rownames(.))%>%
  pivot_longer(-herd)%>%
  rename(year=name,
         S=value)%>%
  mutate(year=parse_number(year))

dist.long<- dist.all%>%
  as_tibble()%>%
  mutate(herd=rownames(.))%>%
  pivot_longer(-herd)%>%
  rename(year=name,
         dist=value)%>%
  mutate(year=parse_number(year))

lamba.dist.long <- abund%>%
  left_join(lambda, by=c("herd", "year"))%>%
  left_join(S, by=c("herd", "year"))%>%
  left_join(dist.long, by=c("herd", "year"))%>%
  left_join(hd%>%select(name=herd, herd=herd_num)%>%mutate(herd=as.character(herd)), by="herd")%>%
  filter(year!=1)

ggplot(lamba.dist.long,
       aes(x=dist, y=lambda))+
  geom_point()+
  geom_smooth()

ggplot(lamba.dist.long,
       aes(x=dist, y=S, col=name))+
  geom_point()


ggplot(lamba.dist.long, aes(x=year, y=S, color=totNMF, group=herd))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE, alpha=0.2)

##look at raw data..not really seeing that trend
ggplot(sdat, aes(x=year, y=(1-(events/trials))^365,group=herd))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE, alpha=0.2)


# Look at raw data and IPM data without climate and disturbance as covariates (from treatment effects paper we submitted) ---------------------------------------------------------------
dist.long.all <- dist.long%>%
  rename(herd_num=herd,
         yr_idx=year)%>%
  left_join(hd%>%select(herd,herd_num)%>%mutate(herd_num=as.character(herd_num)), by="herd_num")%>%
  left_join(yr_df%>%rename(year=yrs), by="yr_idx")%>%
  select(herd,year,dist)

dist.long.all%>%
  ggplot(aes(x=year,y=dist))+
  geom_path()+
  facet_wrap(vars(herd))

counts <- read_csv(here::here("data/clean/counts.csv"))%>%
  left_join(dist.long.all, by=c("herd","year"))


counts%>%
  ggplot(aes(x=dist,y=Est_CL, color=herd))+
  geom_point()+
  geom_smooth(method="lm")

counts%>%
  filter(!herd%in%c("Itcha-Ilgachuz","Graham"))%>%
  ggplot(aes(x=dist,y=Est_CL, color=herd))+
  geom_point()+
  geom_smooth(method="lm")

counts%>%
  ggplot(aes(x=dist,y=Est_CL, color=herd))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(vars(herd), scales="free_y")

library(broom)
dist.herds.glm <- counts%>%
  group_by(herd)%>%
  add_count()%>%
  filter(n>1)%>%
  nest(data = c(-herd))%>% 
  mutate(
    fit = map(data, ~ lm(Est_CL ~  dist, data = .x)),
    tidied = map(fit, tidy)
  )%>%
  unnest(tidied) %>% 
  select(-data, -fit)%>%
  filter(term=="dist")


positive.herds <- dist.herds.glm%>%
  filter(term=="dist",
         estimate>0,
         p.value<0.1)

dist.herds.glm%>%
  filter(term=="dist")%>%
  left_join(dist.long.all%>%
              group_by(herd)%>%
              summarise(disturbance=median(dist)),
            by="herd")%>%
  ggplot(aes(x=disturbance,y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, label=herd))+
  geom_point()+
  geom_linerange()+
  geom_label()


dist.herds.glm%>%
  filter(term=="dist",
         !herd%in%c("Tonquin", "Brazeau", "Graham", "Maligne", "Banff", "Hart South"))%>%
  left_join(dist.long.all%>%
              group_by(herd)%>%
              summarise(disturbance=median(dist)),
            by="herd")%>%
  ggplot(aes(x=disturbance,y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, label=herd))+
  geom_point()+
  geom_linerange()+
  geom_label()+
  geom_smooth(method="lm", se=FALSE)


library(lme4)
library(broom.mixed)
glmm <- lmer(Est_CL~dist + (1+ dist|herd), data=counts%>%group_by(herd)%>%add_count()%>%filter(n>3))
  
glmm.result <- glmm%>%
  tidy()%>%
  filter(term=="dist")%>%
  mutate(type="GLMM- disturbance random slope, herd random intercept")%>%
  select(type, term, estimate, std.error, statistic)
  
glm <- lm(estimate~dist, data=dist.herds.glm%>%
            left_join(dist.long.all%>%
                        group_by(herd)%>%
                        summarise(dist=median(dist)),
                      by="herd"))%>%
  tidy()%>%
  filter(term=="dist")%>%
  mutate(type="GLM-random slopes versus mean herd-level disturbance")%>%
select(type, term, estimate, std.error, statistic)


glm2 <- lm(estimate~dist, data=dist.herds.glm%>%
             filter(term=="dist",
                    !herd%in%c("Tonquin", "Brazeau", "Graham", "Maligne", "Banff", "Hart South"))%>%
            left_join(dist.long.all%>%
                        group_by(herd)%>%
                        summarise(dist=median(dist)),
                      by="herd")%>%
             left_join(counts%>%
                         group_by(herd)%>%
                         summarise(abund=median(Est_CL)),
                       by="herd")
                       )%>%
  tidy()%>%
  filter(term=="dist")%>%
  mutate(type="GLM-random slopes versus mean herd-level disturbance (remove outliers/extremes)")%>%
  select(type, term, estimate, std.error, statistic)


rbind(glm,glm2,glmm.result)


## what about abundance/density?


lm(Est_CL ~  dist, data=counts)%>%summary()

counts%>%
  ggplot(aes(x=dist, y=Est_CL))+
  geom_point()+
  geom_smooth()

library(sf)
library(units)
herd.ranges <- st_read(here::here("data/spatial/herds/ipm_herds.shp"))%>%
  st_transform(crs=3005)%>%
  mutate(area=st_area(.)%>%set_units(km^2)%>%as.numeric)

density <- counts%>%
  left_join(herd.ranges%>%
              tibble()%>%
              select(herd,area,ECCC),
            by="herd")%>%
  mutate(D=Est_CL/(area/100))


lm(D ~  dist, data=density)%>%summary()

density%>%
  ggplot(aes(x=dist, y=D))+
  geom_point()+
  geom_smooth()


density%>%
  ggplot(aes(x=dist, y=D, group=herd, color=ECCC))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE)


lm(D ~  dist, data=density)%>%tidy()%>%filter(term=="dist")
lmer(D ~  dist + (1|herd), data=density)%>%tidy()%>%filter(term=="dist")
lmer(D ~  dist + (1+dist|herd), data=density)%>%tidy()%>%filter(term=="dist")

dist.herds.glm.d <- density%>%
  group_by(herd)%>%
  add_count()%>%
  filter(n>1)%>%
  drop_na(D)%>%
  nest(data = c(-herd))%>% 
  mutate(
    fit = map(data, ~ lm(D ~  dist, data = .x)),
    tidied = map(fit, tidy)
  )%>%
  unnest(tidied) %>% 
  select(-data, -fit)%>%
  filter(term=="dist")


dist.herds.glm.d %>%
  left_join(dist.long.all%>%
              group_by(herd)%>%
              summarise(disturbance=median(dist)),
            by="herd")%>%
  ggplot(aes(x=disturbance,y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, label=herd))+
  geom_point()+
  geom_linerange()+
  geom_label()+
  geom_smooth(method="lm", se=FALSE)


dist.herds.glm.d %>%
  filter(term=="dist",
         !herd%in%c("Tonquin", "Brazeau", "Graham", "Maligne", "Banff", "Hart South"))%>%
  left_join(dist.long.all%>%
              group_by(herd)%>%
              summarise(disturbance=median(dist)),
            by="herd")%>%
  ggplot(aes(x=disturbance,y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, label=herd))+
  geom_point()+
  geom_linerange()+
  geom_label()+
  geom_smooth(method="lm", se=FALSE)



# RUN again but without data from the 4 + herds and a few others possibly an issue  ---------------------------------------------------------------


ipm_dat4 <- list(
  nherd = nherd,
  nyr = nyr,
  first = first,
  last=last$yrs,
  
  month_offset = month_offset,
  count.otc=count.otc,
  count.osc=count.osc,
  count.mnka=count.mnka,
  
  nc = nrow(cdat%>%filter(!herd%in%c(positive.herds$herd,"Tonquin", "Brazeau", "Graham", "Maligne", "Banff"))), 
  ne = nrow(edat%>%filter(!herd%in%c(positive.herds$herd,"Tonquin", "Brazeau", "Graham", "Maligne", "Banff"))),
  ns = nrow(edat%>%filter(!herd%in%c(positive.herds$herd,"Tonquin", "Brazeau", "Graham", "Maligne", "Banff"))),
  nr = nrow(rdat%>%filter(!herd%in%c(positive.herds$herd,"Tonquin", "Brazeau", "Graham", "Maligne", "Banff"))),
  nsr = nrow(srdat%>%filter(!herd%in%c(positive.herds$herd,"Tonquin", "Brazeau", "Graham", "Maligne", "Banff"))),
  
  nsight_grp = nsight_grp,
  sight_grp = hd$sight_grp,
  
  mean_grp_p = mean_grp_p,
  mean_grp_ptau = mean_grp_ptau,
  meansr = meansr,
  
  n1 = n1, 
  cdat = cdat%>%filter(!herd%in%c(positive.herds$herd,"Tonquin", "Brazeau", "Graham", "Maligne", "Banff")),
  edat = edat%>%filter(!herd%in%c(positive.herds$herd,"Tonquin", "Brazeau", "Graham", "Maligne", "Banff")),
  sdat = sdat%>%filter(!herd%in%c(positive.herds$herd,"Tonquin", "Brazeau", "Graham", "Maligne", "Banff")),
  srdat = srdat%>%filter(!herd%in%c(positive.herds$herd,"Tonquin", "Brazeau", "Graham", "Maligne", "Banff")),
  pdat = pdat%>%filter(!herd%in%c(positive.herds$herd,"Tonquin", "Brazeau", "Graham", "Maligne", "Banff")),
  rdat = rdat%>%filter(!herd%in%c(positive.herds$herd,"Tonquin", "Brazeau", "Graham", "Maligne", "Banff")),
  
  dist = dist.all, 
  clim1 = clim1,
  clim2= clim2,
  
  prior.dist.s.mean=priors%>%filter(.variable=="beta.dist.s")%>%pull(mean),
  prior.dist.s.tau=priors%>%filter(.variable=="beta.dist.s")%>%pull(tau),
  prior.clim1.s.mean=priors%>%filter(.variable=="beta.clim1.s")%>%pull(mean),
  prior.clim1.s.tau=priors%>%filter(.variable=="beta.clim1.s")%>%pull(tau),
  prior.clim2.s.mean=priors%>%filter(.variable=="beta.clim2.s")%>%pull(mean),
  prior.clim2.s.tau=priors%>%filter(.variable=="beta.clim2.s")%>%pull(tau),
  
  prior.dist.r.mean=priors%>%filter(.variable=="beta.dist.r")%>%pull(mean),
  prior.dist.r.tau=priors%>%filter(.variable=="beta.dist.r")%>%pull(tau),
  prior.clim1.r.mean=priors%>%filter(.variable=="beta.clim1.r")%>%pull(mean),
  prior.clim1.r.tau=priors%>%filter(.variable=="beta.clim1.r")%>%pull(tau),
  prior.clim2.r.mean=priors%>%filter(.variable=="beta.clim2.r")%>%pull(mean),
  prior.clim2.r.tau=priors%>%filter(.variable=="beta.clim2.r")%>%pull(tau)
)

nth <- 50
nbu <- 2000 
nch <- 3
nad <- 2000
nit <- 30000 

out4 <- jagsUI::jags(ipm_dat4, 
                     inits = ipm_inits,
                     model_parms,
                     model.file = here::here("jags/Dist_Clim_IPM_rawVR.txt"),
                     n.chains = nch,
                     n.cores = nch,
                     n.iter = nit,
                     n.burnin = nbu,
                     n.thin = nth,
                     n.adapt = nad)


draws4 <- out4%>%
  gather_draws(muS,beta.dist.s, beta.clim1.s, beta.clim2.s, 
               muR,beta.dist.r, beta.clim1.r, beta.clim2.r)%>%
  mutate(type="Posteriors-no +ve herds")%>%
  rbind(draws3)


draws4%>%
  filter(!.variable%in%c("beta.season.r","muS","muR"))%>%
  group_by(type,.variable)%>%
  summarise(mean=mean(.value),
            sd=sd(.value))%>%
  ggplot(aes(x=mean, xmin=mean-sd, xmax=mean+sd, y=.variable, color=type))+
  geom_point(alpha=0.5, position=position_dodge(0.4))+
  geom_linerange(alpha=0.5, position=position_dodge(0.4))+
  geom_vline(xintercept = 0, linetype="dashed")


