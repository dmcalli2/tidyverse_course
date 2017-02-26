# Waiting_times.R
library(tidyverse)
library(stringr)
library(broom)

## Create messy data
## Ignore this bit, it is just some code creating a set of messay datasets
wait_raw_all <- map(as.list(paste0("yr", 1990:2015)), function (filename){
  wait_raw <- expand.grid(age = seq(20, 40, 10), sex = c("m", "f"), region = c("Fife", "Lothian"))
  wait_raw$age <- paste(wait_raw$age, wait_raw$age + 9, sep = " to " )
  wait_raw$q50 = rpois(nrow(wait_raw), 100)
  wait_raw$q25 <- pmax(wait_raw$q50 - rpois(nrow(wait_raw), 50), 1)
  wait_raw$q75 <- wait_raw$q50 + rpois(nrow(wait_raw), 50)
  wait_raw$est_iqr <- paste0(wait_raw$q50, " (",wait_raw$q25, " to ", wait_raw$q75, ")")
  wait_raw <- select(wait_raw, age, sex, region, est_iqr)
  wait_raw <- reshape2::acast(wait_raw, age + sex ~ region, value.var = "est_iqr" )
  wait_raw <- as.data.frame(wait_raw)
  wait_raw_write <- rownames_to_column(wait_raw, "age_sex") %>% 
    as_tibble()
  write_csv(wait_raw_write, paste0("data/cancer_waits/", filename, ".csv"))
  wait_raw
})
wait_raw <- wait_raw_all[[1]]
rm(wait_raw_all)

## First move rownames to the dataframe itself
wait<- rownames_to_column(wait_raw, "age_sex") %>% 
  as_tibble()

## Gather regions
wait <- wait %>%
  gather(key = "region", value = cmbn, -age_sex) %>%  # gather regions
  separate(age_sex, into = c("age", "sex"), sep = "_") %>%  # separate age and sex
  separate(cmbn, into = c("q50", "q25", "q75"), sep = "\\(|to") %>%  # separate estimate into components
  mutate_at(vars(q50, q25, q75), parse_number) # convert estimates to numbers

## Now I can plot this data
wait %>% 
  ggplot(aes(x = age, y = q50, tmin = q25, ymax = q75, colour = sex)) +
  geom_point() +
  facet_wrap(~region)

## Or aggregate it
wait %>%
  group_by(region, sex) %>% 
  summarise_at(vars(q50, q25, q75), sum)

## Or run regression models
mod1 <- wait %>%
  glm(q50 ~ age + sex + region, family = "poisson", data = .)

## I can also examine the model outputs in a tidy format
glance(mod1)  # Model-level summaries
tidy(mod1)    # Parameter-level summaries
augment(mod1) # Observation-level summaries
