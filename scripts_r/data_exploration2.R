
## Load Packages
library(tidybayes)
library(ggmcmc)
library(hrbrthemes)
library(sf)
library(units)
library(lme4)
library(broom.mixed)
library(MuMIn)
library(mgcv)
library(tidylog)
library(tidyverse)

# Load data ---------------------------------------------------------------
hd <- read.csv(here::here("data/clean/blueprint.csv"))
hn <- hd %>%
  dplyr::select(herd, herd_num)

## demographic data
counts <- read_csv(here::here("data/clean/forIPM/counts.csv"))
surv <- read_csv(here::here("data/clean/forIPM/afs.csv"))
recruit <- read_csv(here::here("data/clean/forIPM/afr.csv"))
sexratios <- read_csv(here::here("data/clean/forIPM/sexratios.csv"))

## other
covar.names <- readRDS(here::here("data/clean/forIPM/covariate_names.rds"))
covar.names_scaled <- paste0(covar.names,"_scale")

covariates <- read_csv(here::here("data/clean/forIPM/covariates_scaled.csv")) %>%
  select(herd, year, !!covar.names_scaled) %>%
  mutate(
    year = year,
    year1 = year + 1,
    year2 = year + 2,
    year5 = year + 5,
    year10 = year + 10,
    year15 = year + 15
  )


yr_df <- read_csv(here::here("data/clean/forIPM/yr_df.csv"))

herd.ranges <- st_read(here::here("data/spatial/herds/ipm_herds.shp")) %>%
  st_transform(crs = 3005) %>%
  mutate(
    area = st_area(.) %>% set_units(km^2) %>% as.numeric(),
    herd = case_when(
      herd == "Bearhole Redwillow" ~ "Narraway BC",
      TRUE ~ herd
    )
  )

# Bring together, remove very imprecise counts (error >50% estimate) ---------------------------------------------------------------
df <- covariates%>%
  # counts
  left_join(
    counts %>%
      mutate(D.error.p = (Est_CL.max - Est_CL) / Est_CL)%>%
      #filter(D.error.p < 0.5) %>%
      select(herd, year, Est_CL, D.error.p) %>%
      add_row(herd = "Narrow Lake", year = 2021, Est_CL = 1, D.error.p=0.2),
    by = c("herd", "year")
  ) %>%
  ## add density
  left_join(
    herd.ranges %>%
      tibble() %>%
      select(herd, area, ECCC),
    by = "herd"
  ) %>%
  mutate(D = Est_CL / (area / 100)) %>%
  # S
  left_join(
    surv %>%
      mutate(error.p = sd / est) %>%
      filter(error.p < 0.5) %>%
      select(herd, year, s = est, s.sd = sd),
    by = c("herd", "year")
  ) %>%
  # R
  left_join(
    recruit %>%
      mutate(error.p = sd / est) %>%
      filter(error.p < 0.5) %>%
      select(herd, year, r = est, r.sd = sd, season_int),
    by = c("herd", "year")
  ) %>%
  # Sex ratio
  left_join(
    sexratios %>%
      mutate(error.p = sratio.sd / sratio) %>%
      filter(error.p < 0.5) %>%
      select(herd, year, sratio, sratio.sd),
    by = c("herd", "year")
  ) %>%
  ## add time lagged vars
  left_join(
    covariates %>%
      select(herd, year = year1, !!covar.names_scaled) %>%
      rename_at(vars(-herd, -year), function(x) paste0(x, "_lag1")),
    by = c("herd", "year")
  ) %>%
  left_join(
    covariates %>%
      select(herd, year = year2, !!covar.names_scaled) %>%
      rename_at(vars(-herd, -year), function(x) paste0(x, "_lag2")),
    by = c("herd", "year")
  ) %>%
  left_join(
    covariates %>%
      select(herd, year = year5, !!covar.names_scaled) %>%
      rename_at(vars(-herd, -year), function(x) paste0(x, "_lag5")),
    by = c("herd", "year")
  ) %>%
  left_join(
    covariates %>%
      select(herd, year = year10, !!covar.names_scaled) %>%
      rename_at(vars(-herd, -year), function(x) paste0(x, "_lag10")),
    by = c("herd", "year")
  ) %>%
  left_join(
    covariates %>%
      select(herd, year = year15, !!covar.names_scaled) %>%
      rename_at(vars(-herd, -year), function(x) paste0(x, "_lag15")),
    by = c("herd", "year")
  ) %>%
  rename_at(vars(year, !!covar.names_scaled), .funs = ~ paste0( .x, "_lag0"))

# Plot map ---------------------------------------------------------------

herd.ranges%>%
  left_join(covariates%>%
              filter(year==2020)%>%
              pivot_longer(!!covar.names_scaled)%>%
              select(herd,name, value),
            by="herd")%>%
  drop_na(value)%>%
  ggplot()+
  geom_sf(inherit.aes = FALSE, aes(fill=value))+
  facet_wrap(vars(name))


# Look at relationships ---------------------------------------------------------------

covar.plot <- covariates %>%
  select(herd, year, !!covar.names_scaled) %>%
  pivot_longer(!!covar.names_scaled, names_to = "vars", values_to = "vals")%>%
  left_join(
    herd.ranges %>%
      tibble() %>%
      select(herd, area, ECCC),
    by = "herd"
  )


ggplot() +
  geom_smooth(
    data = covar.plot,
    aes(x = year, y = vals, group = herd, color = ECCC),
    method = "lm", se = FALSE
  ) +
  geom_smooth(
    data = covar.plot,
    aes(x = year, y = vals), color = "black", linetype = "dashed", method = "gam", formula = y ~ s(x, k = 5)
  ) +
  guides() +
  facet_wrap(vars(vars), scales = "free_y")+
  labs(
    y = "Covariate values"
  )

ggplot() +
  geom_path(
    data = covar.plot,
    aes(x = year, y = vals, group = herd, color = ECCC)
  ) +
  geom_smooth(
    data = covar.plot,
    aes(x = year, y = vals), color = "black", linetype = "dashed", method = "gam", formula = y ~ s(x, k = 5)
  ) +
  guides() +
  facet_wrap(vars(vars), scales = "free")+
  labs(
    y = "Covariate values"
  )


df.mod.lag <- df %>%
  select(
    ECCC, herd, year_lag0, D, s, r, D.error.p, season_int,
    dist.all_scale_lag0:start.spring_scale_lag0,
    dist.all_scale_lag1:start.spring_scale_lag15
  ) %>%
  pivot_longer(dist.all_scale_lag0:start.spring_scale_lag15) %>%
  mutate(
    covariate = str_split(name, "_", simplify = TRUE)[, 1],
    lag = str_split(name, "_", simplify = TRUE)[, 3] %>% extract_numeric(),
    inv.error=case_when(!is.na(D) & D.error.p==0~1/0.25,
                        D==0~1/0.25,
                        TRUE~1/D.error.p),
    weight=inv.error/mean(inv.error,na.rm=TRUE)
  )

ggplot() +
  geom_smooth(
    data = df.mod.lag %>% filter(covariate == "dist.all"),
    aes(x = value, y = D, group = herd, color = ECCC),
    method = "lm", se = FALSE
  ) +
  guides(color = "none") +
  facet_wrap(vars(lag), scales = "free")+
  labs(
    y = "Density (animals per 100 sq.km)",
    x= "Disturbance"
  )


ggplot() +
   geom_smooth(data = df.mod.lag %>%
  filter(lag == 0)%>%drop_na(D),
               aes(x=value, y=D, group=herd, color=ECCC),
               method="lm", se=FALSE)+
  # geom_line(
  #   data = df.mod.lag %>%
  #     filter(lag == 0)%>%drop_na(D),
  #   aes(x = value, y = D, group = herd, color = ECCC)
  # ) +
  geom_smooth(
    data = df.mod.lag %>%
      filter(lag == 0)%>%drop_na(D),
    aes(x = value, y = D), color = "black", linetype = "dashed"
  ) +
  facet_wrap(vars(covariate), scales = "free") +
  labs(
    title = "Univariate effects of climate and disturbance on caribou density",
    y = "Density (animals per 100 sq.km)"
  )

## what's changing through time for these herds?
change <- mgcv::gam(
  year_lag0 ~ 
    s(grow.seas.length_scale_lag0) +
    s(snow.depth_scale_lag0) +
    s(cooler.winter_scale_lag0) +
    s(warm.dry.summer_scale_lag0) +
    s(start.spring_scale_lag0) +
    s(dist.all_scale_lag0) +
  
    s(herd, bs = "re")+
  
  s(grow.seas.length_scale_lag0, herd, bs="re") +
    s(snow.depth_scale_lag0, herd, bs="re") +
    s(cooler.winter_scale_lag0, herd, bs="re") +
    s(warm.dry.summer_scale_lag0, herd, bs="re") +
    s(start.spring_scale_lag0, herd, bs="re") +
    s(dist.all_scale_lag0, herd, bs="re"),
  data = df %>% mutate(herd = as.factor(herd)),
  method = "REML"
)

plot(change, page = 1, scheme = 2)


## correlation matrix
corr_matrix <- cor(df %>%
  select(!!paste0(covar.names_scaled,"_lag0")),use="complete")
ggcorrplot(corr_matrix,
  hc.order = TRUE,
  type = "lower",
  lab = TRUE
)


## by herd
library(corrr)

## vs year_lag0
df %>%
  select(herd, year_lag0, !!paste0(covar.names_scaled,"_lag0")) %>%
  pivot_longer(!!paste0(covar.names_scaled,"_lag0")) %>%
  nest(data = c(-herd, -name)) %>%
  mutate(
    cor = map(data, ~ cor(.x$value, .x$year_lag0), use="complete.obs")
  ) %>%
  unnest(cor) %>%
  select(-data) %>%
  group_by(name) %>%
  summarise(mean.corr = mean(abs(cor)))


# create model data
# remove herds that have been decoupled from disturbance-mediated apparent competition due to increases in prey by other means (reintroductions)
## subset of data
df.mod.dat <- df.mod.lag %>%
  filter(!herd %in% c("Maligne", "Tonquin", "Brazeau", "Banff")) %>%
  mutate(
    herd = as.factor(herd),
    year_lag0 = as.factor(year_lag0)
  ) %>%
  select(-name) %>%
  pivot_wider(names_from = covariate, values_from = value)




## models with herd as random intercept and random slope on vars, so testing how, on average, each variable affects within a herd
### D
gam1.within0 <- mgcv::gam(
  D ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  weights=weight,
  data = df.mod.dat %>% filter(lag == 0)
)

gam1.within1 <- mgcv::gam(
  D ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  weights=weight,
  data = df.mod.dat %>% filter(lag == 1)
)

gam1.within2 <- mgcv::gam(
  D ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  weights=weight,
  data = df.mod.dat %>% filter(lag == 2)
)

gam1.within5 <- mgcv::gam(
  D ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  weights=weight,
  data = df.mod.dat %>% filter(lag == 5)
)

gam1.within10 <- mgcv::gam(
  D ~s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  weights=weight,
  data = df.mod.dat %>% filter(lag == 10)
)

gam1.within15 <- mgcv::gam(
  D ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  weights=weight,
  data = df.mod.dat %>% filter(lag == 15)
)

gam1.within10_modified <- mgcv::gam(
  D ~s(grow.seas.length_10, k=5) +
    s(snow.depth_1, k=5) +
    s(cooler.winter_1, k=5) +
    s(warm.dry.summer_1, k=5) +
    s(start.spring_10, k=5) +
    s(dist.all_10, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length_10, k=5, herd, bs="re") +
    s(snow.depth_1, k=5, herd, bs="re") +
    s(cooler.winter_1, k=5, herd, bs="re") +
    s(warm.dry.summer_1, k=5, herd, bs="re") +
    s(start.spring_10, k=5, herd, bs="re") +
    s(dist.all_10, k=5, herd, bs="re"), 
  weights=weight,
  data = df.mod.dat%>%
    select(herd,year_lag0,D,lag,weight:start.spring)%>%
    pivot_longer(cols=dist.all:start.spring)%>%
    pivot_wider(names_from=c(name, lag), values_from=value)
)

summary(gam1.within0)$r.sq
summary(gam1.within1)$r.sq
summary(gam1.within2)$r.sq
summary(gam1.within5)$r.sq
summary(gam1.within10)$r.sq
summary(gam1.within15)$r.sq
summary(gam1.within10_modified)$r.sq


#model.sel(gam1.within0, gam1.within1, gam1.within2, gam1.within5, gam1.within10)

a <- df.mod.dat%>%group_by(herd)%>%drop_na(D)%>%distinct(herd,year_lag0)%>%mutate(year=as.character(year_lag0)%>%as.numeric())%>%summarise(n=n(), min.yr=min(year), max.yr=max(year))

plot(gam1.within10, page = 1, scheme = 2)


mods <- list(gam1.within0, gam1.within1, gam1.within2, gam1.within5, gam1.within10)
names(mods) <- c("gam1.within0", "gam1.within1", "gam1.within2", "gam1.within5", "gam1.within10")
pred.within.d <- tibble()
for (i in 1:length(mods)) {
  newdata.dist <- tibble(
    dist.all = seq(0, 0.8, length.out = 20),
    clim1 = 0,
    clim2 = 0,
    herd = "Chase"
  )

  newdata.clim1 <- tibble(
    dist.all = 0.3,
    clim1 = seq(min(df$clim1_scale), max(df$clim1_scale), length.out = 20),
    clim2 = 0,
    herd = "Chase"
  )

  newdata.clim2 <- tibble(
    dist.all = 0.3,
    clim1 = 0,
    clim2 = seq(min(df$clim1_scale), max(df$clim1_scale), length.out = 20),
    herd = "Chase"
  )

  pred.within.d <- rbind(
    pred.within.d,
    cbind(
      name = "dist.all",
      value = newdata.dist$dist.all,
      D = predict(mods[[i]], newdata = newdata.dist, type = "response"),
      model = names(mods)[i],
      R2=summary(mods[[i]])$r.sq
    ) %>%
      rbind(
        cbind(
          name = "clim1",
          value = newdata.clim1$clim1,
          D = predict(mods[[i]], newdata = newdata.clim1, type = "response"),
          model = names(mods)[i],
          R2=summary(mods[[i]])$r.sq
        )
      ) %>%
      rbind(
        cbind(
          name = "clim2",
          value = newdata.clim2$clim2,
          D = predict(mods[[i]], newdata = newdata.clim2, type = "response"),
          model = names(mods)[i],
          R2=summary(mods[[i]])$r.sq
        )
      )
  )
}

pred.within.d<- pred.within.d%>%
  mutate(lag=factor(str_sub(model,-2,-1)%>%extract_numeric(), levels=c("0","1","2","5","10","15")))

ggplot(pred.within.d,
       aes(x=as.numeric(value), y=as.numeric(D), color=as.numeric(R2), group=lag, linetype=lag))+
  geom_line()+
  facet_wrap(vars(name), scales="free_x")+
  labs(x="Covariate value",
       y="Density",
       color="R2",
       title="Optimizing year lags within a herd")



## again across herds

gam1.across0 <- mgcv::gam(
  D ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(year_lag0, bs = "re"),
  weights=weight,
  data = df.mod.dat %>% filter(lag == 0),
  method = "REML"
)

gam1.across1 <- mgcv::gam(
  D ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(year_lag0, bs = "re"),
  weights=weight,
  data = df.mod.dat %>% filter(lag == 1),
  method = "REML"
)

gam1.across2 <- mgcv::gam(
  D ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(year_lag0, bs = "re"), 
  weights=weight,
  data = df.mod.dat %>% filter(lag == 2),
  method = "REML"
)

gam1.across5 <- mgcv::gam(
  D ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(year_lag0, bs = "re"), 
  weights=weight,
  data = df.mod.dat %>% filter(lag == 5),
  method = "REML"
)

gam1.across10 <- mgcv::gam(
  D ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(year_lag0, bs = "re"), 
  weights=weight,
  data = df.mod.dat %>% filter(lag == 10),
  method = "REML"
)

gam1.across15 <- mgcv::gam(
  D ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(year_lag0, bs = "re"), 
  weights=weight,
  data = df.mod.dat %>% filter(lag == 15),
  method = "REML"
)

summary(gam1.across0)$r.sq
summary(gam1.across1)$r.sq
summary(gam1.across2)$r.sq
summary(gam1.across5)$r.sq
summary(gam1.across10)$r.sq
summary(gam1.across15)$r.sq

model.sel(gam1.across0, gam1.across1, gam1.across2, gam1.across5, gam1.across10, gam1.across15)

plot(gam1.across2, page = 1, scheme = 2)



mods <- list(gam1.across0, gam1.across1, gam1.across2, gam1.across5, gam1.across10, gam1.across15)
names(mods) <- c("gam1.across0", "gam1.across1", "gam1.across2", "gam1.across5", "gam1.across10", "gam1.across15")
pred.across.d <- tibble()
for (i in 1:length(mods)) {
  newdata.dist <- tibble(
    dist.all = seq(0, 0.8, length.out = 20),
    clim1 = 0,
    clim2 = 0,
    year_lag0 = "2000"
  )
  
  newdata.clim1 <- tibble(
    dist.all = 0.3,
    clim1 = seq(min(df$clim1_scale), max(df$clim1_scale), length.out = 20),
    clim2 = 0,
    year_lag0 = "2000"
  )
  
  newdata.clim2 <- tibble(
    dist.all = 0.3,
    clim1 = 0,
    clim2 = seq(min(df$clim1_scale), max(df$clim1_scale), length.out = 20),
    year_lag0 = "2000"
  )
  
  pred.across.d <- rbind(
    pred.across.d,
    cbind(
      name = "dist.all",
      value = newdata.dist$dist.all,
      D = predict(mods[[i]], newdata = newdata.dist, type = "response"),
      model = names(mods)[i],
      R2=summary(mods[[i]])$r.sq
    ) %>%
      rbind(
        cbind(
          name = "clim1",
          value = newdata.clim1$clim1,
          D = predict(mods[[i]], newdata = newdata.clim1, type = "response"),
          model = names(mods)[i],
          R2=summary(mods[[i]])$r.sq
        )
      ) %>%
      rbind(
        cbind(
          name = "clim2",
          value = newdata.clim2$clim2,
          D = predict(mods[[i]], newdata = newdata.clim2, type = "response"),
          model = names(mods)[i],
          R2=summary(mods[[i]])$r.sq
        )
      )
  )
}

pred.across.d<- pred.across.d%>%
  mutate(lag=factor(str_sub(model,-2,-1)%>%extract_numeric(), levels=c("0","1","2","5","10","15")))

ggplot(pred.across.d,
       aes(x=as.numeric(value), y=as.numeric(D), color=as.numeric(R2), group=lag, linetype=lag))+
  geom_line()+
  facet_wrap(vars(name), scales="free_x")+
  labs(x="Covariate value",
       y="Density",
       color="R2",
       title="Optimizing year lags across a herd")


















### S
gam1.within0 <- mgcv::gam(
  s~s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  data = df.mod.dat %>% filter(lag == 0)
  )

gam1.within1 <- mgcv::gam(
  s~s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  data = df.mod.dat %>% filter(lag == 1)
)

gam1.within2 <- mgcv::gam(
  s~s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  data = df.mod.dat %>% filter(lag == 2)
)

gam1.within5 <- mgcv::gam(
  s~s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  data = df.mod.dat %>% filter(lag == 5)
)

gam1.within10 <- mgcv::gam(
  s~s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  data = df.mod.dat %>% filter(lag == 10)
)

gam1.within15 <- mgcv::gam(
  s~s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  data = df.mod.dat %>% filter(lag == 15)
)

summary(gam1.within0)$r.sq
summary(gam1.within1)$r.sq
summary(gam1.within2)$r.sq
summary(gam1.within5)$r.sq
summary(gam1.within10)$r.sq
summary(gam1.within15)$r.sq

model.sel(gam1.within0, gam1.within1, gam1.within2, gam1.within5, gam1.within10, gam1.within15)


plot(gam1.within2, page = 1, scheme = 2)




### R
gam1.within0 <- mgcv::gam(
  r ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  data = df.mod.dat %>% filter(lag == 0)
)

gam1.within1 <- mgcv::gam(
  r ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  data = df.mod.dat %>% filter(lag == 1)
)

gam1.within2 <- mgcv::gam(
  r ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  data = df.mod.dat %>% filter(lag == 2)
)

gam1.within5 <- mgcv::gam(
  r ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  data = df.mod.dat %>% filter(lag == 5)
)

gam1.within10 <- mgcv::gam(
  r ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  data = df.mod.dat %>% filter(lag == 10)
)

gam1.within15 <- mgcv::gam(
  r ~ s(grow.seas.length, k=5) +
    s(snow.depth, k=5) +
    s(cooler.winter, k=5) +
    s(warm.dry.summer, k=5) +
    s(start.spring, k=5) +
    s(dist.all, k=5) +
    
    s(herd, bs = "re")+
    
    s(grow.seas.length, k=5, herd, bs="re") +
    s(snow.depth, k=5, herd, bs="re") +
    s(cooler.winter, k=5, herd, bs="re") +
    s(warm.dry.summer, k=5, herd, bs="re") +
    s(start.spring, k=5, herd, bs="re") +
    s(dist.all, k=5, herd, bs="re"), 
  data = df.mod.dat %>% filter(lag == 15)
)

summary(gam1.within0)$r.sq
summary(gam1.within1)$r.sq
summary(gam1.within2)$r.sq
summary(gam1.within5)$r.sq
summary(gam1.within10)$r.sq
summary(gam1.within15)$r.sq


model.sel(gam1.within0, gam1.within1, gam1.within2, gam1.within5, gam1.within10, gam1.within15)


plot(gam1.within2, page = 1, scheme = 2)
