library(tidyverse)

seals <- sealdat %>% 
  select(animalID, year) %>% 
  mutate(in_study = TRUE)

foo <- read_csv(here::here("data/raw/fullresights.csv"), 
                    show_col_types = FALSE) %>% 
  filter(calyear > 1987, 
         tagsex == "F", 
         timeofyear == "Breeding") %>% 
  mutate(animalID = as.factor(animalID), 
         year = calyear) 

longev <- foo %>% 
  select(animalID, year, yearborn) %>%
  filter(year < 2018) %>%
  group_by(animalID) %>% 
  mutate(longevity = max(year) - yearborn) %>% 
  select(animalID, longevity) %>% 
  unique() 

longev %>% 
  group_by(longevity) %>% 
  count(longevity) %>% 
  ggplot(aes(x = longevity, y = n)) + 
  geom_col(fill = "grey50") + 
  labs(x = "Longevity (Age in Years)", y = "Number of individuals") + 
  theme_classic()

mean(longev$longevity)

foo2 <- left_join(foo, seals, by = c("animalID", "year")) %>% 
  mutate(in_study = ifelse(is.na(in_study), yes = "FALSE", no = "TRUE")) %>% 
  group_by(animalID, year) %>% 
  mutate(n_obs = n_distinct(date))

foo2 %>% filter(animalID == 10297) %>% 
  # select(animalID, date, year, in_study, n_obs) %>% 
    view()

obs_by_breed_status <- foo2 %>% 
  filter(!is.na(withpup)) %>% 
  select(animalID, season, date, withpup, tagsex, yearborn, year, in_study, n_obs) %>% 
  mutate(breed_status = if_else(withpup == 1, "Breeder", "Non-breeder")) %>% 
  ungroup()

n_distinct(obs_by_breed_status$animalID)
#What proportion of seals were seen 4 times in a breeding season? 
bred <- obs_by_breed_status %>% filter(breed_status == "Breeder")
bred_used_df <- bred %>% filter(n_obs > 3) 
#I think this part might be wrong? get max to check
total_bred <- n_distinct(bred$animalID)
bred_used <- n_distinct(bred_used$animalID)
  
#Are non-breeders less likely to be seen 4 times?
nonbred <- obs_by_breed_status %>% filter(breed_status == "Non-breeder")
nonbred_used_df <- nonbred %>% filter(n_obs > 3) 
total_nonbred <- n_distinct(nonbred$animalID)
nonbred_used <- n_distinct(nonbred_used$animalID)

obs_1 <- obs_by_breed_status %>% filter(breed_status == "Breeder") %>% 
  group_by(n_obs) %>% 
  count(n_obs) %>% 
  ggplot(aes(x = n_obs, y = n)) + 
  geom_col() +
  labs(x = "Observation frequency", y = "count", title = "Breeder") + 
  theme_classic()

obs_2 <- obs_by_breed_status %>% filter(breed_status == "Non-breeder") %>% 
  group_by(n_obs) %>% 
  count(n_obs) %>% 
  ggplot(aes(x = n_obs, y = n)) + 
  geom_col() +
  ylim(0, 2500) +
  labs(x = "Observation frequency", y = "count", title = "Non-breeder") + 
  
  theme_classic()

(obs_1 | obs_2)
