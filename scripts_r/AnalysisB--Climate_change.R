library(lme4)
library(broom.mixed)
library(here)
library(tidyverse)
covar.names <- readRDS(here::here("data/clean/forIPM/covariate_names.rds"))
covariates <- read_csv(here::here("data/spatial/covariates_simplified.csv")) %>%
  dplyr::select(herd, year,ECCC, !!covar.names)

covariates%>%
  select(herd,year,ECCC,greenness, dist.hum, dist.nat, snow.depth,warm.dry.summer,start.spring,rain.on.snow)%>%
  pivot_longer(c(-herd,-year,-ECCC))%>%
ggplot() +
  geom_path(
    aes(x = year, y = value, group = herd), color="black", alpha=0.1
  ) +
  geom_smooth(
    aes(x = year, y = value, group=ECCC, color=ECCC), linetype = "dashed", method = "gam", formula = y ~ s(x, k = 5)
  ) +
  guides() +
  facet_wrap(vars(name), scales = "free", nrow=4)+
  labs(
    y = "Covariate values"
  )+
  theme(legend.position = "bottom")



covariates%>%
  select(herd,year,ECCC,snow.depth,warm.dry.summer,start.spring,rain.on.snow)%>%
  mutate_if(is.numeric, scale)%>%
  pivot_longer(c(-herd,-year,-ECCC))%>%
  nest(data = c(-name,-ECCC))%>% 
  mutate(
    fit = map(data, ~ lmer(value ~  year + (1 + year|herd), data = .x)),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied) %>% 
  select(-data, -fit)%>%
  filter(term=="year")%>%
  ggplot(aes(x = estimate, xmin = estimate-(1.96*std.error), xmax = estimate+(1.96*std.error),  y = name, color=ECCC)) +
  geom_linerange(alpha = 0.8, position = position_dodge(width = 0.3)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x="Estimated change through time (standardized)", y="Weather covariate")

