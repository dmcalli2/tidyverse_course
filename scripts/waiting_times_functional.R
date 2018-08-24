# Waiting_times.R
library(tidyverse)
library(stringr)
library(broom)

## Read in all the data from the folder
waits <- list.files("data/cancer_waits/", full.names = TRUE)
wait_names <- list.files("data/cancer_waits/")  %>% 
  str_replace_all("\\.csv", "" )
names(waits) <- wait_names
waits
waits <- map(waits, function(x) read_csv(file = x))

## Examine contents of each dataset
map(waits, dim) %>%  unique()
map(waits, names) %>% unique()
map(waits, function(x) map_chr(x,
                           function(y) y[1]))

## Bind all together in a single dataset
wait <- bind_rows(waits, .id = "year")

## Convert years to a number
wait$year <-   parse_number(wait$year)

## Gather regions as before, adding instruction to exclude year
wait <- wait %>%
  gather(key = "region", value = cmbn, -age_sex, -year) %>%  # gather regions
  separate(age_sex, into = c("age", "sex"), sep = "_") %>%  # separate age and sex
  separate(cmbn, into = c("q50", "q25", "q75"), sep = "\\(|to") %>%  # separate estimate into components
  mutate_at(vars(q50, q25, q75), parse_number) # convert estimates to numbers

## Now I can plot this data, but adding dimensions
wait %>% 
  ggplot(aes(x = year, y = q50, colour = age, group = region)) +
  geom_line() +
  facet_wrap(sex ~region)

## Or aggregate it
wait %>%
  group_by(year, region, sex) %>% 
  summarise_at(vars(q50, q25, q75), median)

## Or run regression models
# original model
mod1 <- wait %>%
  glm(q50 ~ age + sex + region, family = "poisson", data = .)


# Model by lists
# create list
wait_by_yrs <- split(wait, wait$year)

# Run model
wait_yr_mods <- map(wait_by_yrs, 
                    function(x) {
                      glm(q50 ~ age + sex + region, family = "poisson", 
                          data = x)
                      }
                    )

# extract model
wait_yr_mods_coef <- map(wait_yr_mods, tidy) %>% 
  bind_rows(.id = "year") %>% 
  mutate(year = parse_number(year)) %>% 
  filter(term != "(Intercept)")
# plot P-values for each term and year
ggplot(wait_yr_mods_coef, aes(x = year, y = p.value)) +
  geom_point() + 
  facet_wrap(~term) +
  scale_y_log10()
