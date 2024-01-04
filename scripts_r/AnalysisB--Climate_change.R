## load packages
library(lme4)
library(broom.mixed)
library(here)
library(tidyverse)

## load helpers
source(here::here("scripts_r/helpers/gg_theme_custom.R"))

## load data
covar.names <- readRDS(here::here("data/clean/forIPM/covariate_names.rds"))
covariates <- read_csv(here::here("data/spatial/covariates_simplified.csv")) %>%
  dplyr::select(herd, year, ECCC, !!covar.names)

covariates %>%
  select(herd, year, ECCC,dist.hum, dist.nat,snow.depth, warm.dry.summer, start.spring, rain.on.snow, cooler.winter,greenness) %>%
  pivot_longer(c(-herd, -year, -ECCC)) %>%
  ggplot() +
  geom_path(
    aes(x = year, y = value, group = herd),
    color = "black", alpha = 0.08
  ) +
  geom_smooth(
    aes(x = year, y = value, group = ECCC, color = ECCC),
    linetype = "dashed", method = "gam", formula = y ~ s(x, k = 5)
  ,se=FALSE) +
  guides() +
  facet_wrap(vars(name), scales = "free", nrow = 4) +
  labs(
    y = "Covariate values",
    x="Year",
    color="Group"
  ) +
  theme(legend.position = "bottom")+
  theme.custom

ggsave(here::here("plots", "AnalysisB_AnnualCovars.png"), width = 6, height = 9, bg = "white")




coef.vals <-covariates %>%
  select(herd, year, ECCC, snow.depth, warm.dry.summer, start.spring, rain.on.snow, cooler.winter) %>%
  mutate_if(is.numeric, scale) %>%
  pivot_longer(c(-herd, -year, -ECCC)) %>%
  nest(data = c(-name, -ECCC)) %>%
  mutate(
    fit = map(data, ~ lmer(value ~ year + (1 + year | herd), data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied) %>%
  select(-data, -fit) %>%
  filter(term == "year")
  
ggplot(data=coef.vals , aes(x = estimate, xmin = estimate - (1.96 * std.error), xmax = estimate + (1.96 * std.error), y = name, color = ECCC)) +
  geom_linerange(alpha = 0.8, position = position_dodge(width = 0.3)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = coef.vals %>% filter(abs(statistic) > 1.96), aes(y = name, x = estimate, label = "*", group=ECCC), position = position_dodge(width = 0.3), color = "black", vjust = 0.3, size = 5, show.legend = FALSE) +
  labs(x = "Estimated change through time (standardized)", y = "Weather covariate", color="Group")+
  theme.custom

ggsave(here::here("plots", "AnalysisB_WeatherChange.png"), width = 7, height = 6, bg = "white")
