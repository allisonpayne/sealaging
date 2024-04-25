# load packages
library(RMySQL)
library(lubridate)
library(tidyverse)
library(patchwork)
library(marked)
library(RMark)

animal <- read_csv(here::here("data/raw/animal.csv"))
resight <- read_csv(here::here("data/raw/resight.csv"))

# specify which animals (pups) to select
animals_to_select <- animal %>%
  filter(
    # from Año Nuevo
    region == "AN",
    # tagged as pup/wean
    broad == "YOY",
    # with a flipper tag
    !is.na(origtag),
    # with a correct date
    !is.na(year(tagdate))
  ) #both sex

pups=subset(animals_to_select,momID!=0) #subset of pups where mom ID is known
pups=pups[,c("animalID","tagsex","momID","tagdate")] #some columns from the pup dataset
colnames(pups)=c("pupID","pupsex","animalID","tagdate") #rename the columns
pups$year=year(pups$tagdate) #find year that the pup was tagged

#get rid of pups with no tag sex (149 cases out of ~4,000)
pups=pups %>% subset(pupsex=="M"|pupsex=="F")

# final table for moms
outall <- merge(# merge resight
  resight,
  # and animals to select (only tagdate, tagsex, pupID, and animalID columns)
  animals_to_select %>% dplyr::select(tagdate,animalID,momID,tagsex),
  # by animalID, which is the mom
  by = "animalID") %>%
  # create a new column yearborn
  mutate(yearborn = year(tagdate),
         yday=yday(date),
         timeofyear=if_else(yday>274|yday<74,"Breeding", #between Dec 1 and Mar 15
                            if_else(yday>=74&yday<182,"Molt","Other"))) %>%   #between Mar 15 and July 1
  arrange(date)

outall$calyear=year(outall$date) #make a new column for the calendar year of that observation


#correct pup status for weird character values
outall$withpupCor=NA
outall$withpupCor[outall$withpup%in%c(as.character(1:8),"2+","3+",">8","1 or 2",">1")]=1
outall$withpupCor[outall$withpup==0]=0

write.csv(outall,here::here("data/raw/fullresights.csv"),
          row.names=FALSE)

#Encounter History
g1=as.data.frame.matrix(with(outall,table(animalID,calyear))) #table : combinations of observation year and animal ID
#presence/absence with only encounter history (no breeding info, and no multiple sightings per year)
gg=g1
############################################ NEW***
g1[g1<=3]=0 #seen less than or equal to 3 times means "missing" = 0
g1[g1!=0]=1 #seen more than 3 times is "seen" = 1
#in g1, 0 means not seen that year, 1 means seen (MORE THAN 3 TIMES) that year (the >3 part used to be )

#find out how many pup present sightings there are per season per animal
g2=aggregate(withpupCor~animalID+season,data=outall,FUN=sum) #aggregate pup status by year and animal
g2=subset(g2,withpupCor>0) #subset to only seasons with at least one pup observation

g1b=g1 #make a version where there is a B if they bred
for(i in 1:nrow(g2)) {
  g1b[rownames(g1)==g2$animalID[i],colnames(g1b)==g2$season[i]]="B" #1 means seen, B means seen as breeder
}

#make goal dataframe that starts the "final" dataframe...
#one row for each animal-year combo observation.
goal=data.frame(animalID=rep(as.numeric(rownames(g1b)),each=ncol(g1b)),
                year=rep(as.numeric(colnames(g1b)),times=nrow(g1b)),
                observed=as.vector(t(g1b)))

#find first date of haulout for all the b years
dateobsbreed=outall %>%
  group_by(season,animalID) %>% #by season, so that it doesn't pick up end of year
  subset(timeofyear=="Breeding") %>%
  summarise(firstobsbreed=min(date),
            lastobsbreed=max(date),
            numbreed=length(date)) #this year's breeding season

dateobsmolt=outall %>%
  group_by(calyear,animalID) %>% #by year, so it catches the end of the molt
  subset(timeofyear=="Molt") %>%
  summarise(firstobsmolt=min(date),
            lastobsmolt=max(date)) %>%
  mutate(season=calyear+1)

#stitch those molting and breeding dates together and calculate durations (for breeding season, only if there are 3 or more observations)
dateobs=left_join(dateobsmolt,dateobsbreed,
                  by=c("season","animalID")) %>%
  #only last observed molt less than july 1 ish and only first breeding < march 16 ish
  mutate(firstobsbreeddoy=yday(as.Date(firstobsbreed))) %>%
  mutate(firstobsbreeddoy=if_else(firstobsbreeddoy>300,firstobsbreeddoy-365,firstobsbreeddoy)) %>%
  mutate(lastobsmolt=if_else(yday(lastobsmolt)>185,NA,lastobsmolt)) %>%
  mutate(lastobsmoltdoy=yday(as.Date(lastobsmolt))) %>%
  mutate(breeddur=if_else(numbreed>2,
                          difftime(lastobsbreed,firstobsbreed,units="days"),
                          NA),
         tripdur=difftime(firstobsbreed,lastobsmolt,units="days"),
         moltdur=difftime(lastobsmolt,firstobsmolt,units="days")) %>%
  mutate(calyear=calyear+1)

#change column names so the merge will work
colnames(dateobs)[1]="year"

#add in tagsex
tagsexdat=outall %>%
  group_by(animalID) %>%
  count(tagsex) %>%
  top_n(1)

goal=merge(goal,tagsexdat[,c("animalID","tagsex")],by="animalID")

#merge goal matrix (observations) with dates
goal=left_join(goal,dateobs,by=c("year","animalID"))

goal[goal$observed!="B",5:8]=NA #remove breeding first, last, duration dates from non birth years

#warning message for seal 933 in year 2012 because she had two pups

#now, merge pup information
goal=left_join(goal,
               subset(pups[,c("pupID","pupsex","animalID","year")]),
               by=c("year","animalID"))

age=outall %>%
  group_by(animalID) %>%
  select(animalID,yearborn) %>%
  summarize(yearborn=unique(yearborn))
#warning message because animalID 17 has two entries for 1984?

goal=left_join(goal,
               age,
               by=c("animalID"))
goal$age=goal$year-goal$yearborn

final=subset(goal, age>2)
final=subset(final,observed!=0)

#to calculate whether pups were seen again
#this needs to be one row per seal (animals_to_select)

colnames=as.numeric(colnames(g1)) #years of observation

#add column to tell us which year the pup was born (but only if a pup was born that year)
final= final %>%
  mutate(pupyearborn=ifelse(is.na(pupID),"NA",year),
         firstyearindicator=ifelse(is.na(pupID),"NA",year+1))
final$indices=match(final$firstyearindicator,colnames)

#add column - were they observed the year after their wean date?
final$pupseeneveragain=NA
final$puprecruited=NA

rows=rownames(subset(final,!is.na(pupID)))

for(i in rows){ #for each row of those rows
  #make encounter history subset of g1 where the row names of g1 are equal to the pup ID
  enc=subset(gg,row.names(gg)==final$pupID[rownames(final)==i])
  if(nrow(enc)>0&final$firstyearindicator[rownames(final)==i]!=2024){
    final[i,"pupseeneveragain"]=sum(enc[,final$indices[rownames(final)==i]:ncol(g1)])
    encb=subset(g1b,row.names(g1b)==final$pupID[rownames(final)==i])
    final[i,"puprecruited"]=ifelse(sum(encb=="B")>0,1,0) #is there at least one B in the pup
  }else{
    final[i,"pupseeneveragain"]=0
    final[i,"puprecruited"]=0
  }
}

#subset to females only here (moms) for the main dataset
final = final %>% filter(tagsex=="F")

write.csv(final,here::here("data/raw/128L pull 2023_12_05.csv"),
          row.names=FALSE)
