library(RMark)
library(tidyverse)

sealdat <- read_csv("data/raw/128L pull 2023_12_05.csv", 
                    show_col_types = FALSE) %>% 
  mutate(observed = if_else(observed == "B", "Breeder", "Non-breeder"), 
         observed_int = if_else(observed == "Breeder", 1, 0), 
         pup_survived = pupseeneveragain > 0) %>% 
  filter(age > 3, 
         year > 1987) %>% 
  mutate(animalID = factor(animalID),
         year_fct = factor(year),
         age10 = age / 10) %>% 
  group_by(animalID) %>% 
  mutate(longevity = ifelse(max(year) < 2020,
                            max(age),
                            NA),
         longevity10 = longevity / 10) %>% 
  ungroup()

sealhist <- sealdat %>% 
  select(animalID, year, observed) %>% 
  mutate(observed = substr(observed, 1, 1)) %>% 
  group_by(animalID, year) %>% 
  summarize(observed = ifelse(any(observed == "B"), "B", "N"),
            .groups = "drop") %>% 
  pivot_wider(names_from = year, 
              values_from = observed,
              values_fill = "0") %>% 
  pivot_longer(-animalID, names_to = "year", values_to = "observation") %>%
  group_by(animalID) %>% 
  summarize(ch = paste(observation, collapse = "")) %>% 
  count(ch, name = "freq")

seal_proc <- process.data(sealhist, model = "Multistrata")
seal_ddl <- make.design.data(seal_proc)

S.stratum <- list(formula =~ stratum)
S.stratumxtime <- list(formula =~ stratum * time)

p.stratum <- list(formula =~ stratum)
p.stratumxtime <- list(formula =~ stratum * time)
p.stratumxtime.fixed <- list(formula =~ stratum * time,
                             fixed = list(time = 4, value = 1))

Psi.s <- list(formula =~ -1 + stratum:tostratum)

model.list <- create.model.list("Multistrata")

model.list <- rbind(model.list,
                    c(S = "S.stratumxtime", 
                      p = "p.stratumxtime",
                      Psi = "Psi.s"))

ms_mod <- mark.wrapper(model.list,
                       data = seal_proc,
                       ddl = seal_ddl,
                       delete = TRUE)

saveRDS(ms_mod, "cluster/ms_mod.rds")















