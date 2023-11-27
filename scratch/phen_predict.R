
expand_grid(
  age = 4:19,
  observed = c("Breeder", "Non-breeder")
) %>% 
  mutate(ba = predict(phenology_models$model[[1]],
                      newdata = .,
                      re.form = NA),
         bd = predict(phenology_models$model[[2]],
                      newdata = .,
                      re.form = NA), 
         ma = predict(phenology_models$model[[3]],
                      newdata = .,
                      re.form = NA), 
         md = predict(phenology_models$model[[4]],
                      newdata = .,
                      re.form = NA)) %>% 
  mutate(short_trip = ma - bd, 
         molt_dur = md - ma) %>% 
  filter(observed == "Breeder", 
         age %in% c(4,19)) %>% 
  view()

