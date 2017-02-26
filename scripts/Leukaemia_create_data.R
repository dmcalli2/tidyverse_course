library(tidyverse)
library(stringr)

# Weight the events based on the board
boards <- c('GGC', 'Lanarkshire', 'Tayside', 'Fife', 'Lothian')
popwt <- c(10, 5, 6, 5.5, 8)
names(popwt) <-  boards
admwt <- c(1,1.1,1.15,1.05,2)

## First create the data from a list of ICD-120 codes. This is not part of the exercise.
icd10 <- read_csv("data/allvalid2011icd10.csv", skip = 6)
icd10 <- select(icd10, 1:3)
names(icd10) <- names(icd10) %>% 
  tolower() %>% 
  str_replace_all("\\s", "_")
icd10 <- icd10 %>% 
  filter(nchar(code) >= 4, str_sub(code,1,1) == "C" )

## Pt level data
dobs <- as.Date("1980-01-01") + 1:(365*26)
chis <- map(dobs, ~ paste0(.x, "0",sample(100:999, 5)))

patient <- map2(dobs, chis, ~ tibble(chi = .y, dob = .x)) %>% 
  bind_rows()
rm(dobs, chis)

patient <- patient %>% 
  mutate(sex = sample(0:1, nrow(.), replace = TRUE),
         chi = str_replace_all(chi, "-", ""))


## Admision /death level codes
adm <- tibble(chi = sample(patient$chi, nrow(patient) * 0.2, replace = TRUE))

adm <- adm %>% 
  mutate(cis_n = sample(1:10, size = nrow(.), replace = TRUE)) %>% 
  group_by(chi) %>% 
  do(cis = 1:.$cis_n) %>% 
  unnest() %>% 
  mutate(epi_n = sample(1:4, size = nrow(.), replace = TRUE)) %>% 
  group_by(chi, cis) %>% 
  do(episode = 1:.$epi_n) %>% 
  unnest()

adm <- inner_join(adm, select(patient, chi, dob), by = "chi") %>% 
  mutate(doa = sample(1:(365*20), size = nrow(.), replace = TRUE) + dob,
         doa = ifelse(doa > as.Date("2016-01-01"), NA, doa),
         dod = sample(365:(365*20), size = nrow(.), replace = TRUE) + dob,
         dod = ifelse(dod > as.Date("2016-01-01") | dod < doa, NA, dod)) %>% 
  select(-dob)

adm <- adm %>% 
  group_by(!is.na(dod)) %>% 
  mutate(death_code_id = ifelse(!is.na(dod), seq_along(dod), NA)) %>% 
  ungroup() %>% 
  group_by(!is.na(doa)) %>% 
  mutate(adm_code_id = ifelse(!is.na(doa), seq_along(doa), NA)) %>% 
  ungroup()

adm[ , c("doa", "dod")] <- map(adm[ , c("doa", "dod")],
                                        as.Date, origin = "1970-01-01")

adm <- adm %>% 
  mutate(hb = sample (boards,  nrow(.), replace = TRUE, 
                           prob = popwt * admwt))

patient <- patient %>% 
  inner_join(adm %>%
               select(chi, hb),
             by = "chi")

adm_only <- adm %>%
  select(chi, cis, episode, doa, adm_code_id) %>% 
  filter(!is.na(doa)) %>% 
  mutate(dodis = doa + sample(0:5, nrow(.), replace = TRUE)) %>% 
  mutate(locn = sample ((1:10),  nrow(.), replace = TRUE))

death_only <- adm %>%
  select(chi, dod, death_code_id) %>% 
  filter(!is.na(dod)) %>% 
  group_by(chi) %>% 
  slice(1) %>% 
  ungroup()

rm(adm)

### Admission codes
adm_codes <- tibble(codes_id = rep(adm_only$adm_code_id, 10)) %>% 
  group_by(codes_id) %>% 
  mutate(pos = seq_along(codes_id)) %>% 
  ungroup() %>% 
  filter(pos <= 6) %>% 
  mutate(codes = sample(icd10$code, nrow(.), replace = TRUE))

death_codes <- tibble(codes_id = rep(death_only$death_code_id, 10)) %>% 
  group_by(codes_id) %>% 
  mutate(pos = seq_along(codes_id)) %>% 
  ungroup() %>% 
  filter(pos <=4) %>% 
  mutate(codes = sample(icd10$code, nrow(.), replace = TRUE))

## population data
pop <- expand.grid(hb = unique(patient$hb), age = 0:30, sex = 0:1, year = 1995:2016,
                   stringsAsFactors = FALSE)

pop$n <- 10000 * popwt[pop$hb]
pop <- as_tibble(pop)

adm <- adm_only
rm(adm_only)
death <- death_only
rm(death_only)

save(pop, patient, adm, death, adm_codes, death_codes, file = "data/Fake_leukaemia_data.Rdata")


#### Run code
