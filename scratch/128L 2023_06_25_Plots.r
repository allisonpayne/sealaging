#FIGURE CHANGES TO MAKE
#Create annual cycle figure

#################################################################################

setwd("C:/Users/roxan/Documents/Teaching/Current Courses/BIOE 128L/128L 2023/128L Projects")
final=read.csv("128L pull 2023_06_29.csv")
final$observed[final$observed=="B"]="Breeder"
final$observed[final$observed==1]="Non-Breeder"

library(tidyverse)
library(ggthemes)
library(emmeans)

##### MOLT DURATION
agg=final %>% group_by(age,observed) %>%
  summarize(mean=mean(moltdur,na.rm=TRUE),
            sd=sd(moltdur,na.rm=TRUE),
            n=n(),
            se=sd/sqrt(n))

ggplot(data=agg,aes(x=age,y=mean,group=observed,col=observed))+
  geom_point(aes(col=observed))+
  geom_smooth(method="lm",mapping=aes(weight=n))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se,col=observed), width=.2)+
  theme_few()+
  ylab("Molt Haulout Duration (Mean +- SE)")+
  geom_text(data=subset(agg,observed=="Breeder"),aes(label=n,vjust=1,
                                                     col=observed),y=.03) + 
  geom_text(data=subset(agg,observed=="Non-Breeder"),aes(label=n,vjust=-1,
                                                         col=observed),y=.03) 
require(lmerTest)
#summary(lm(as.numeric(moltdur)~age+observed,data=final))
summary(lmer(as.numeric(moltdur)~age+observed+(1|animalID)+(1|year),data=final))
#summary(lm(as.numeric(moltdur)~age*observed,data=final))
summary(lmer(as.numeric(moltdur)~age*observed+(1|animalID)+(1|year),data=final))

#Molting Date
agg=final %>% group_by(age,observed) %>%
  subset(!is.na(lastobsmoltdoy)) %>% 
  summarize(mean=mean(lastobsmoltdoy),
            sd=sd(lastobsmoltdoy),
            n=n(),
            se=sd/sqrt(n))

ggplot(data=agg,aes(x=age,y=mean,group=observed,col=observed))+
  geom_jitter(data=final,aes(x=age,y=lastobsmoltdoy),alpha=0.02)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+
  theme_few()+
  stat_smooth(method="lm",mapping=aes(weight=n))+
  ylab("Last Molting Observation (day of year; Mean +- SE)")+
  geom_text(data=subset(agg,observed=="Breeder"),aes(label=n,vjust=1,
                                                     col=observed),y=80) + 
  geom_text(data=subset(agg,observed=="Non-Breeder"),aes(label=n,vjust=-1,
                                                         col=observed),y=80) 

summary(lmer(lastobsmoltdoy~age+observed+(1|animalID)+(1|year),data=final))
#summary(lm(lastobsmoltdoy~age+observed,data=final))

#TRIP DURATION
agg=final %>% 
  group_by(age,observed) %>%
  summarize(mean=mean(tripdur,na.rm=TRUE),
            sd=sd(tripdur,na.rm=TRUE),
            n=n(),
            se=sd/sqrt(n))

ggplot(data=agg,aes(x=age,y=mean,group=observed,col=observed))+
  #geom_point(data=final,
  #           aes(x=age,y=tripdur),alpha=0.1)+
  geom_point(size=2.5)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+
  theme_few()+
  geom_smooth(method="lm",mapping=aes(weight=n))+
  labs(y="Trip Duration (Mean +- SE)",x="Age")+
  geom_text(data=subset(agg,observed=="Breeder"),aes(label=n,vjust=1,
                                                     col=observed),y=220) + 
  geom_text(data=subset(agg,observed=="Non-Breeder"),aes(label=n,vjust=-1,
                                                         col=observed),y=220) 


#summary(lm(tripdur~age+observed,data=final))
summary(lmer(tripdur~age+observed+(1|animalID)+(1|year),data=final))
m3=lmer(tripdur~age*observed+(1|animalID)+(1|year),data=final)
summary(m3)
emmeans::lstrends(m3,"observed",var="age",infer=T)
emmeans(m3,"observed",infer=T)

#BREEDING DATE
agg=final %>% group_by(age,observed) %>%
  subset(!is.na(firstobsbreeddoy)) %>% 
  summarize(mean=mean(firstobsbreeddoy),
            sd=sd(firstobsbreeddoy),
            n=n(),
            se=sd/sqrt(n))

ggplot(data=agg,aes(x=age,y=mean,group=observed,color=observed))+
  geom_jitter(data=final,aes(x=age,y=firstobsbreeddoy),alpha=0.1)+
  geom_point()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+
  theme_few()+ylim(14,30)+
  stat_smooth(method="lm",mapping=aes(weight=n))+
  ylab("First Breeding Observation (day of year; Mean +- SE)")+
  geom_text(data=subset(agg,observed=="Breeder"),aes(label=n,vjust=1,
                                                     col=observed),y=30) + 
  geom_text(data=subset(agg,observed=="Non-Breeder"),aes(label=n,vjust=-1,
                                                         col=observed),y=30) 

summary(lmer(firstobsbreeddoy~age+(1|animalID)+(1|year),data=final))
summary(lmer(firstobsbreeddoy~age+observed+(1|animalID)+(1|year),data=final))

#####
agg=final %>% group_by(age,observed) %>%
  summarize(mean=mean(breeddur,na.rm=TRUE),
            sd=sd(breeddur,na.rm=TRUE),
            n=n(),
            se=sd/sqrt(n))

ggplot(data=agg,aes(x=age,y=mean,group=observed,col=observed))+
  geom_point(aes(col=observed))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se,col=observed), width=.2)+
  theme_few()+stat_smooth(method="lm",mapping=aes(weight=n))+
  ylab("Breeding Haulout Duration (Mean +- SE)")+
  geom_text(data=subset(agg,observed=="Breeder"),aes(label=n,vjust=1,
                                                     col=observed),y=18) + 
  geom_text(data=subset(agg,observed=="Non-Breeder"),aes(label=n,vjust=-1,
                                                         col=observed),y=18) 

#add breeder?
summary(lmer(as.numeric(breeddur)~age+(1|animalID)+(1|year),data=final))
summary(lmer(as.numeric(breeddur)~age+observed+(1|animalID)+(1|year),data=final))
summary(lmer(as.numeric(breeddur)~age*observed+(1|animalID)+(1|year),data=final))

######
agg2=final %>% group_by(age,observed) %>%
  count(observed) %>% 
  pivot_wider(values_from=n,names_from=observed) %>% 
  mutate(num=`Non-Breeder`+`Breeder`,
         percB=`Breeder`/num,
         se=(percB*(1-percB)/num)^.5)

ggplot(data=subset(agg2,age>3),aes(x=age,y=percB,weight=num))+
  geom_point()+
  geom_errorbar(aes(ymin=percB-se, ymax=percB+se), width=.2)+
  theme_few()+labs(x="Maternal Age",y="% Females Breed")+
  stat_smooth(method="lm", formula = y ~ poly(x, 2), se = TRUE)+
  geom_text(aes(label=num,vjust=1,y=.03)) 
#ggsave("PupSexPoints.png",width=5,height=5)

summary(glmer(as.factor(observed)~poly(age,2)+(1|animalID)+(1|year),
              data=subset(final,age>3),
            family="binomial"))

##### PUP SEX
agg=final %>% group_by(age,pupsex) %>%
  subset(pupsex=="M"|pupsex=="F") %>% 
  count(pupsex)

#PUP SEX
agg2=final %>% group_by(age,pupsex) %>%
  subset(pupsex=="M"|pupsex=="F") %>% 
  count(pupsex) %>% 
  pivot_wider(values_from=n,names_from=pupsex) %>% 
  mutate(num=M+F,
         percf=F/num,
         se=(percf*(1-percf)/num)^.5)

ggplot(data=subset(agg2,age<18),aes(x=age,y=percf,weight=num))+
  geom_point()+
  geom_errorbar(aes(ymin=percf-se, ymax=percf+se), width=.2)+
  geom_hline(yintercept=0.5,lty=2)+
  theme_few()+labs(x="Maternal Age",y="% Pups Female")+
  stat_smooth(method="glm", method.args=list(family="binomial"))+
  geom_text(aes(label=num,vjust=1,y=.03)) 
#ggsave("PupSexPoints.png",width=5,height=5)

summary(glmer(as.factor(pupsex)~age+(1|animalID)+(1|year),
         data=subset(final,age<18),family="binomial"))
summary(glmer(as.factor(pupsex)~age+(1|animalID)+(1|year),
              data=final,family="binomial"))

#### pup seen ever again?
final=final %>% 
  mutate(pupseeneveragainYN=ifelse(pupseeneveragain>0,1,0))

agg=final %>% 
  group_by(age,pupseeneveragainYN) %>%
  subset(!is.na(pupseeneveragainYN)) %>% 
  count(pupseeneveragainYN) %>% 
  pivot_wider(names_from=pupseeneveragainYN,values_from=n) %>% 
  mutate(n=(`0`+`1`),
         perc=`0`/n,
         se=(perc*(1-perc)/n)^.5)

ggplot(data=subset(agg,n>5),aes(x=age,y=perc))+
  geom_errorbar(aes(ymin=perc-se,ymax=perc+se),width=.2)+
  geom_point()+
  labs(x="Age",y="Fraction Pups Not Seen Again")+
  theme_few()+
  geom_smooth(method="glm",method.args=list(family=binomial))+
  geom_text(data=agg,aes(label=n,vjust=-1),y=.5) 

agg=final %>% 
  group_by(age,pupseeneveragainYN,pupsex) %>%
  subset(!is.na(pupseeneveragainYN)) %>% 
  count(pupseeneveragainYN) %>% 
  pivot_wider(names_from=pupseeneveragainYN,values_from=n) %>% 
  mutate(n=(`0`+`1`),
         perc=`0`/n,
         se=(perc*(1-perc)/n)^.5)

ggplot(data=subset(agg,n>5),aes(x=age,y=perc,weight=n,col=pupsex))+
  geom_errorbar(aes(ymin=perc-se,ymax=perc+se),width=.2)+
  geom_point()+
  labs(x="Age",y="Fraction Pups Not Seen Again")+
  theme_few()+
  geom_smooth(method="glm",method.args=list(family=binomial))+
  geom_text(data=subset(agg,pupsex=="M"),aes(label=n,vjust=1,
                                                     col=pupsex),y=.5) + 
  geom_text(data=subset(agg,pupsex=="F"),aes(label=n,vjust=-1,
                                                         col=pupsex),y=.5) 

summary(glmer(pupseeneveragainYN~age+pupsex+(1|animalID)+(1|year),
              data=final,family="binomial"))
summary(glmer(pupseeneveragainYN~age+(1|animalID)+(1|year),
              data=final,family="binomial"))
summary(glmer(pupseeneveragainYN~age+(1|animalID),
              data=final,family="binomial"))

#pup recruited
agg2=final %>% 
  filter(pupsex=="F") %>% 
  group_by(age,puprecruited) %>%
  subset(!is.na(puprecruited)) %>% 
  count(puprecruited) %>% 
  pivot_wider(names_from=puprecruited,values_from=n) %>% 
  mutate(n=(`0`+`1`),
         perc=`0`/n,
         se=(perc*(1-perc)/n)^.5)

ggplot(data=subset(agg2,n>5),aes(x=age,y=perc,weight=n))+
  geom_errorbar(aes(ymin=perc-se,ymax=perc+se),width=.2)+
  geom_point()+xlim(2.5,13.5)+
  stat_smooth(method="glm",method.args = list(family="binomial"))+
  labs(x="Age",y="Fraction Pups Never Breed")+
  theme_few()+
  geom_text(aes(label=n,vjust=1,y=1))

summary(glmer(puprecruited~age+(1|animalID)+(1|year),
              data=final,family="binomial"))

####### TEMPORARY for mid April temporary data pull of "old" animals
ID=c(1366, 1464, 1561, 7838, 8070, 11424,16496,16547,366813)
temp=subset(final, animalID %in% ID)
setwd("C:/Users/roxan/Desktop/128L Projects")
write.csv(temp,"128Loldanimals.csv",row.names=FALSE)


#extra
#Check whether pup survival ~ mom timing (birth or breeding). NOPE. 
ggplot(final,aes(x=lastobsmoltdoy,y=pupseeneveragainYN))+geom_jitter(alpha=.03)+
  geom_smooth(method="glm",method.args=list(family="binomial"))
summary(glm(pupseeneveragainYN~lastobsmoltdoy,data=final,family="binomial"))
summary(glm(pupseeneveragainYN~firstobsbreeddoy,data=final,family="binomial"))

#6. Check whether mom repo success ~ mom timing
ggplot(final,aes(x=observed,y=lastobsmoltdoy))+geom_boxplot()+geom_jitter(alpha=.03)
ggplot(final,aes(x=observed,y=firstobsbreeddoy))+geom_boxplot()+geom_jitter(alpha=.03)

summary(lm(lastobsmoltdoy~observed,data=final))
summary(lm(firstobsbreeddoy~observed,data=final))


#EXTRA
f1=glm(perc~age,data=agg[!agg$age%in%c(10,11,20),],family="binomial",weights=n);summary(f1)
nd1=data.frame(age=3:21)
nd1$pred=predict(f1,newdata=nd1,type="response")
nd1$SE=predict(f1,newdata=nd1,type="response",se.fit=T)$se

