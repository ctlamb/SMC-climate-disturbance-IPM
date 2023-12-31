# Packages
library(rgee)
library(sf)
library(mapview)
library(mapedit)
library(lubridate)
library(lwgeom)
library(fasterize)
library(tabularaster)
library(googledrive)
library(terra)
library(units)
library(lme4)
library(broom.mixed)
library(tidyverse)

# Helpers
source(here::here("scripts_r/helpers/gg_theme_custom.R"))

# Load data

## herd data
herds <- st_read(here::here("data/spatial/herds/ipm_herds.shp"))

## annualized ndvi from AVHRR 1982-2015 (8km res?)
ndvi.files <- tibble(path = list.files(here::here("data/spatial/ndvi/biweekly/avhrr"), recursive = TRUE, pattern = "*tif")) %>%
  mutate(
    date = str_sub(path, -12, -5),
    year = str_sub(date, 1, 4) %>% as.numeric(),
    month = str_sub(date, 5, 6) %>% as.numeric(),
    week = str_sub(date, 7, 8) %>% as.numeric(),
    monthweek = str_sub(date, 5, 8) %>% as.numeric()
  )

## annualized EVI 2000-2021 (500 m?)
evi.files <- tibble(path = list.files(here::here("data/spatial/EVI/clean/modis_evi"), recursive = TRUE, pattern = "*tif")) %>%
  mutate(year = extract_numeric(path) %>% as.numeric() %>% abs())



########################################################
## Analyse annualized ndvi from AVHRR 1982-2015
########################################################
ndvi.files.summer <- ndvi.files %>%
  filter(monthweek %in% 702:801) ## use last week of july and first week of august

## extract for each herd-year
herds.geo <- herds %>% st_transform(4326)
summer.ndvi <- tibble()
summer.ndvi.herd <- tibble()
for (i in 1:nrow(ndvi.files.summer)) {
  rast.i <- rast(here::here("data/spatial/ndvi/biweekly/avhrr", ndvi.files.summer$path[i])) %>%
    crop(herds.geo)

  remove <- rast.i[[2]] %>%
    subst(0, 1) %>%
    subst(2:3, NA) %>%
    subst(65535, NA)

  rast.clean <- rast.i[[1]] * remove

  mean.i <- rast.clean %>%
    values() %>%
    mean(na.rm = TRUE)

  summer.ndvi <- summer.ndvi %>%
    rbind(tibble(
      year = ndvi.files.summer$year[i],
      mean = mean.i
    ))

  summer.ndvi.herd <- summer.ndvi.herd %>%
    rbind(terra::extract(rast.clean, herds.geo %>% vect(), bind = TRUE, fun = "mean", na.rm = TRUE) %>%
      data.frame() %>%
      mutate(year = ndvi.files.summer$year[i]) %>%
      select(herd, year, ECCC, ndvi = 10))
}


## look at through time
plot(summer.ndvi$year, summer.ndvi$mean)
lm(summer.ndvi$mean ~ summer.ndvi$year) %>% summary()

summer.ndvi2 <- summer.ndvi %>%
  group_by(year) %>%
  summarise(mean = mean(mean))
plot(summer.ndvi2$year, summer.ndvi2$mean)
lm(summer.ndvi2$mean ~ summer.ndvi2$year) %>% summary()


## plot
ggplot(
  summer.ndvi.herd %>%
    group_by(herd, year, ECCC) %>%
    summarise(mean = mean(ndvi)),
  aes(x = year, y = mean, color = ECCC, group = herd)
) +
  geom_path(alpha = 0.4) +
  geom_smooth(aes(group = ECCC), se = FALSE) +
  theme.custom +
  labs(x = "Year", y = "NDVI")


ggplot(
  summer.ndvi.herd %>%
    group_by(herd, year, ECCC) %>%
    summarise(mean = mean(ndvi)),
  aes(x = year, y = mean)
) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  facet_wrap(vars(herd)) +
  theme.custom +
  labs(title = "NDVI 1982-2020", x = "Year", y = "NDVI")


## mixed effect modelling
lmer(mean ~ year + (1 + year | herd), REML = F, data = summer.ndvi.herd %>%
  group_by(herd, year) %>%
  summarise(mean = mean(ndvi))) %>%
  tidy()

lmer(mean ~ year + (1 + year | herd), REML = F, data = summer.ndvi.herd %>%
       group_by(herd, year, ECCC) %>%
       summarise(mean = mean(ndvi))) %>%
  ranef() %>%
  data.frame()

summer.ndvi.herd %>%
  group_by(herd, year, ECCC) %>%
  summarise(mean = mean(ndvi)) %>%
  nest(data = c(-ECCC)) %>%
  mutate(
    fit = map(data, ~ lmer(mean ~ year + (1 + year | herd) + (1 | herd), data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied) %>%
  select(-data, -fit) %>%
  filter(term == "year")


## just post 2000
summer.ndvi.herd %>%
  filter(year > 2000) %>%
  group_by(herd, year, ECCC) %>%
  summarise(mean = mean(ndvi)) %>%
  nest(data = c(-ECCC)) %>%
  mutate(
    fit = map(data, ~ lmer(mean ~ year + (1 + year | herd) + (1 | herd), data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied) %>%
  select(-data, -fit) %>%
  filter(term == "year")

## check variation by decade
summer.ndvi.herd %>%
  group_by(herd, year, ECCC) %>%
  summarise(mean = mean(ndvi)) %>%
  group_by(herd) %>%
  arrange(year) %>%
  mutate(change = abs(mean - lag(mean))) %>%
  ggplot(aes(x = year, y = change, color = ECCC, group = herd)) +
  geom_point() +
  geom_path(color = "black", alpha = 0.2) +
  geom_smooth(aes(group = ECCC), se = FALSE) +
  theme.custom +
  labs(title = "NDVI Change 1982-2020", x = "Year", y = "NDVI Change")

summer.ndvi.herd %>%
  filter(year > 2000) %>%
  group_by(herd, year, ECCC) %>%
  summarise(mean = mean(ndvi)) %>%
  group_by(herd) %>%
  arrange(year) %>%
  mutate(change = abs(mean - lag(mean))) %>%
  ggplot(aes(x = year, y = change, color = ECCC, group = herd)) +
  geom_point() +
  geom_path() +
  geom_smooth(aes(group = ECCC), se = FALSE) +
  theme.custom +
  labs(title = "NDVI Change 2000-2020", x = "Year", y = "NDVI Change")


############################
## Analyze for EVI 2000-2021
############################
## extract for each herd-year
herds.utm <- herds %>% st_transform(26911)
summer.evi <- tibble()
summer.evi.herd <- tibble()

m <- c(
  0, 0.2, NA,
  0.2, 2, 1
)
rclmat <- matrix(m, ncol = 3, byrow = TRUE)
### mask out disturbed pixels
dist <- rast(here::here("data/spatial/disturbance/2023-11-23_unbuffered/dist.all.unbuffered.tif")) %>%
  project(rast(here::here("data/spatial/EVI/clean/modis_evi", evi.files$path[1]))) %>%
  classify(rclmat, include.lowest = TRUE)

m <- c(
  1, 1, NA,
  NA, NA, 1
)
rclmat <- matrix(m, ncol = 3, byrow = TRUE)
no.dist <- dist %>% classify(rclmat, include.lowest = TRUE)


for (i in 1:nrow(evi.files)) {
  rast.clean <- rast(here::here("data/spatial/EVI/clean/modis_evi", evi.files$path[i]))

  mean.i <- rast.clean %>%
    values() %>%
    mean(na.rm = TRUE)

  summer.evi <- summer.evi %>%
    rbind(tibble(
      year = evi.files$year[i],
      mean = mean.i
    ))

  rast.stack <- c(rast.clean, rast.clean * no.dist, rast.clean * dist)
  names(rast.stack) <- c("evi", "evi.nodist", "evi.dist")

  summer.evi.herd <- summer.evi.herd %>%
    rbind(terra::extract(rast.stack, herds.utm %>% vect(), bind = TRUE, fun = "mean", na.rm = TRUE) %>%
      data.frame() %>%
      mutate(year = evi.files$year[i]) %>%
      select(herd, year, ECCC, evi, evi.nodist, evi.dist) %>%
      cbind(terra::extract(rast.stack, herds.utm %>% vect(), bind = TRUE, fun = \(x) length(na.omit(x)) / length(x)) %>%
        data.frame() %>%
        select(evi.prop = evi, evi.nodist.prop = evi.nodist, evi.dist.prop = evi.dist)))
}

## plot change through time

plot(summer.evi$year, summer.evi$mean)
lm(summer.evi$mean ~ summer.evi$year) %>% summary()

# summer.evi.herd <- summer.evi.herd %>%
#   mutate(
#     evi = evi * 1000,
#     evi.nodist = evi.nodist * 1000,
#     evi.dist = evi.dist * 1000
#   )

ggplot(
  summer.evi.herd,
  aes(x = year, y = evi, color = ECCC, group = herd)
) +
  geom_path(alpha = 0.4) +
  geom_smooth(aes(group = ECCC), se = FALSE) +
  theme.custom +
  labs(title = "EVI 2000-2021", x = "Year", y = "NDVI Change")

ggplot(
  summer.evi.herd %>% pivot_longer(evi:evi.dist) %>% filter(name != "evi"),
  aes(x = year, y = value, color = name)
) +
  geom_point() +
  geom_smooth(aes(group = name), se = FALSE, method = "lm") +
  facet_wrap(vars(herd), scales = "free_y") +
  theme.custom +
  labs(title = "EVI 2000-2021", x = "Year", y = "NDVI Change")

## mixed models
lmer(evi ~ year + (1 + year | herd), REML = F, data = summer.evi.herd %>%mutate_at(vars(evi, evi.nodist, evi.dist), scale) %>%
       mutate(year = year - 1999)) %>%
  tidy()

summer.evi.herd %>%
  mutate_at(vars(evi, evi.nodist, evi.dist), scale) %>%
  mutate(year = year - 1999) %>%
  nest(data = c(-ECCC)) %>%
  mutate(
    fit = map(data, ~ lmer((evi) ~ year + (1 + year | herd), data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied) %>%
  select(-data, -fit) %>%
  filter(term == "year")

## Assess changes by disturbed and nondisturbed areas

coefs.evi <- summer.evi.herd %>%
  # mutate_at(vars(evi, evi.nodist,evi.dist), scale)%>%
  mutate(year = year - 1999) %>%
  group_by(ECCC) %>%
  mutate(mean.evi.group = mean(evi)) %>%
  pivot_longer(evi:evi.dist) %>%
  nest(data = c(-ECCC, -name)) %>%
  mutate(
    fit = map(data, ~ lmer(((value / mean.evi.group) * 100) ~ year + (1 + year | herd), data = .x, control = lmerControl(optimizer = "Nelder_Mead"))),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied) %>%
  select(-data, -fit) %>%
  filter(term == "year") %>%
  left_join(
    summer.evi.herd %>%
      pivot_longer(evi.prop:evi.dist.prop) %>%
      group_by(ECCC, name) %>%
      summarise(prop = mean(value)) %>%
      mutate(name = str_sub(name, 1, -6)),
    by = c("ECCC", "name")
  ) %>%
  mutate(
    estimate.prop = estimate * prop,
    std.error.prop = std.error * prop
  )

ggplot(coefs.evi, aes(y = ECCC, x = estimate.prop, xmin = estimate.prop - std.error.prop, xmax = estimate.prop + std.error.prop, fill = name, color = name)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_linerange(position = position_dodge(width = 0.6)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_text(aes(y = ECCC, x = estimate.prop, label = prop %>% round(2)), vjust = 1.4, position = position_dodge(width = 0.6), show.legend = FALSE) +
  geom_text(data = coefs.evi %>% filter(abs(statistic) > 1.96), aes(y = ECCC, x = estimate.prop, label = "*", group=name), position = position_dodge(width = 0.6), color = "black", vjust = 0., size = 5, show.legend = FALSE) +
  theme.custom +
  labs(x = "Annual % change in EVI", color = "Type", fill = "Type", y = "Group")

ggsave(here::here("plots", "AnalysisA_EviChange.png"), width = 6, height = 5, bg = "white")
