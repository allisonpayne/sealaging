# Is this the paper you do with Millie?

seal_agent <- function(
    # Parameterize mom's life history here
    ) {
  seal <- list(
    age = 3,
    sex = 0, # female
    pups = 0,
    grandpups = 0,
    recruits = 0,
    alive = TRUE
  )
  while (seal$alive) {
    seal <- seal_agent_year(seal)
  }
  seal
}

seal_agent_year <- function(seal) {
  survived <- rbinom(1, 1, some_probability(seal$age)) # from the age paper
  bred <- rbinom(1, 1, some_probability(seal$age)) # from the age paper
  if (bred == 1) {
    pup = list(
      # this is the strategy switching part
      sex = rbinom(1, 1, some_probability(seal$age))
    )
    if (pup$sex == 0) { 
      # if female
      pup$recruited <- rbinom(1, 1, some_probability)
      pup$lifetime_repro <- some_distribution() # parameterize using whole dataset?
    } else {
      # if male
      pup$recruited <- one_in_nine_maybe() # ???
      pup$lifetime_repro <- any_clue() # ???
    }
    seal$pups <- seal$pups + 1
    seal$recruits <- seal$recruits + pup$recruited
    seal$grandpups <- seal$grandpups + pup$lifetime_repro
  }
  if (survived) {
    seal$age <- seal$age + 1
  } else {
    seal$alive <- FALSE
  }
  seal
}

##pup sex ratio change - taking on a riskier strategy as they get older. well established phenomenoon! 
#seems to be a case of late-in-life change towards a riskier strategy, not maternal effect senescence
#captial breeder? 
#fitness benefits depend on the risk reward ratio for male vs female. we dont know this yet because we only track females
#future work on paternity would allow us to quantitatively test optimal life history strategies in this model system