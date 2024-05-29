library(marked)
library(tidyverse)

sealdat <- read_csv("data/raw/128L pull 2023_12_05.csv", 
                    show_col_types = FALSE) %>% 
  filter(between(age, 3, 18), 
         year > 1987) %>% 
  mutate(observed = 1)

sealhist <- sealdat %>% 
  group_by(animalID, year) %>% 
  summarize(observed = 1,
            .groups = "drop") %>% 
  pivot_wider(names_from = year, 
              values_from = observed,
              values_fill = 0) %>% 
  pivot_longer(-animalID, 
               names_to = "year", 
               values_to = "observation") %>%
  arrange(animalID, year) %>% 
  group_by(animalID) %>% 
  summarize(ch = paste(observation, collapse = "")) 

seal_proc <- process.data(sealhist)
seal_ddl <- make.design.data(seal_proc)

fit.models <- function() {
  Phi.age <- list(formula = ~age)
  p.age <- list(formula = ~age)
  cml <- create.model.list(c("Phi","p"))
  crm.wrapper(cml, 
              data = seal_proc, 
              ddl = seal_ddl,
              external = FALSE,
              accumulate = FALSE)
}

seal_models <- fit.models()

model_summary <- tibble(
  age = seal_models[["Phi.age.p.age"]][["results"]][["reals"]][["Phi"]][["age"]],
  Phi = seal_models[["Phi.age.p.age"]][["results"]][["reals"]][["Phi"]][["estimate"]],
  p = seal_models[["Phi.age.p.age"]][["results"]][["reals"]][["p"]][["estimate"]]
)

model_summary %>% 
  pivot_longer(-age, names_to = "param", values_to = "value") %>% 
  ggplot(aes(age, value, color = param)) +
  geom_line(aes(group = param)) +
  geom_point() +
  theme_classic()
