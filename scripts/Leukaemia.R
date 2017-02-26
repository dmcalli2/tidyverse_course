## Leukaemia exercise
library(tidyverse)
library(stringr)
library(lubridate)
library(Hmisc)
library(broom)

## First load the data
load( file = "data/Fake_leukaemia_data.Rdata")

#####
# Identifiy leukaemia codes
# First search for codes pertaining to Acute myeloid leukaemia
adm_codes$codes[str_detect(adm_codes$codes, "C92")] %>%  
  unique() %>% 
  paste(collapse = "', '")

# Next create a vector of codes for searching
codes_want <- c('C92.2', 'C92.9', 'C92.5', 'C92.3', 'C92.1', 'C92.7', 'C92.4', 'C92.0')

# Select admission episodes where one of these codes appear in the first 4 positions
adm_codes_slct <- adm_codes %>% 
  filter(codes %in% codes_want, pos %in% 1:4) %>% 
  select(codes_id) %>% 
  distinct()

# Select deaths where one of these codes appear in any position
death_codes_slct <- death_codes %>% 
  filter(codes %in% codes_want) %>% 
  select(codes_id) %>% 
  distinct()

# Select any admission where one of these codes appeared 
adm_s <- adm %>%
  semi_join(adm_codes_slct, by = c("adm_code_id" = "codes_id")) %>% 
  arrange(chi, cis, episode) %>% 
  group_by(chi, cis) %>%  
  slice(1) %>%                       ## Reduces data to admission-level data
  ungroup() %>% 
  select(chi, cis, doa, dodis)

# Select the first admission with any of these codes (not worrying about lookback)
adm_s <- adm_s %>% 
  group_by(chi) %>% 
  arrange(cis, dodis) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(chi, doa, dodis)

# Similarly for death
death_s <- death %>% 
  semi_join(death_codes_slct, by = c("death_code_id" = "codes_id")) %>% 
  select(chi, dod)

# Combine deaths and admissions into a single incident events table
events <- bind_rows(adm_s, death_s)
events <- events %>% 
  mutate(do_event  = if_else(is.na(dod), dodis, dod)) %>% 
  arrange(chi, do_event) %>%
  group_by(chi) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(chi, do_event)

## Add patient-level data to the incident events data and aggregate by hb, age, sex and year
patient_s <- patient %>% 
  inner_join(events, by = "chi") %>% 
  mutate(age = (do_event - dob) /365.25,
         age = as.numeric(age) %>%  round(),
         year = year(do_event)) %>% 
  filter(age <= 12, year >= 1995) %>% 
  group_by(hb, age, sex, year) %>% 
  summarise(events = length(do_event))

# Add in denominator data and label data nicely for tables and graphs
final <- inner_join(pop, patient_s, by = c("hb", "age", "sex", "year")) %>% 
  mutate(rate = 100000*events/n,
         sex = factor(sex, levels = 0:1, labels = c("Male", "Female")))

# Run regression model
pois_mod <- glm(events ~ I(year/10) +hb + offset(log(n)), data = final, family = "poisson")

# Create prediction dataset for following plots
# OK if don't follow this section
pred_data <- final %>% 
  select(-n, -events) %>% 
  distinct() %>% 
  complete(hb, age, sex, year)
pred_data$n <- 100000

pred <-  predict(pois_mod, newdata = pred_data, se.fit = TRUE)
pred_data$rate <- exp(pred$fit)
pred_data$lci <- exp(pred$fit - 1.96*pred$se.fit)
pred_data$uci <- exp(pred$fit + 1.96*pred$se.fit)

## Make plot - we will discuss ggplot later today, this is for illustration
final %>% 
  ggplot(aes(x = year, y = rate, colour = hb, fill = hb)) + 
  geom_point(position = "jitter", alpha = 0.2) +
  facet_grid(~sex) +
  geom_ribbon(data = pred_data, aes(ymin = lci, ymax = uci), alpha = 0.5) +
  scale_y_continuous("Incidence rate per 100,000 person-years") +
  scale_x_continuous("Year admitted") +
  scale_fill_discrete("Health Board") +
  scale_color_discrete("Health Board")

## Make table of regression results
pois_mod %>% 
  tidy() %>% 
  mutate(est = exp(estimate),
         lci = exp(est -1.96*std.error),
         uci = exp(est +1.96*std.error),
         p_value = round(p.value,3),
         p_value = if_else(p_value == 0, "<0.001", parse_character(p_value))) %>% 
  mutate_at(vars(est, lci, uci), funs(round(.,2))) %>% 
  select(term, est:p_value) %>% 
  slice(-1)

## Compare Lothian to the rest of Scotland
loth_scot <- final %>%
  mutate(hb = if_else(hb == "Lothian", "Lothian", "Rest of Scotland")) %>% 
  group_by(hb, age, sex, year) %>% 
  summarise_at(vars(-hb), funs("sum")) %>% 
  ungroup()

pois_mod_scot <- update(pois_mod, data = loth_scot)
