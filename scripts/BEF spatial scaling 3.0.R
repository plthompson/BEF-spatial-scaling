gamma<-0.77

for(r in 1:3000){
  Richness<-round(runif(n = 50,1,max=40))
  Biomass<-c(sapply(Richness,FUN = function(x){(1/(x^gamma))*x}))
  for(i in 1:50){
    patches<-sample(1:50,size = i,replace = F)
    if(r==1 & i == 1){
      results<-data.frame(Scale=i,Richness=sum(Richness[patches]),Biomass=sum(Biomass[patches]))
    } else {
      results<-rbind(results,data.frame(Scale=i,Richness=sum(Richness[patches]),Biomass=sum(Biomass[patches])))
    }
  }
}

library(dplyr)
library(broom)
library(ggplot2)
library(vegan)

ggplot(results,aes(x=Richness,y=Biomass,group=Scale,color=Scale))+
  stat_smooth(method = 'lm')+
  scale_x_log10()+
  scale_y_log10()

ggplot(filter(results,Scale==40),aes(x=Richness,y=Biomass,group=Scale,color=Scale))+
  geom_point()

scale_40<-filter(results,Scale==50)

summary(lm(log(scale_40$Biomass)~log(scale_40$Richness)))

BEF_models<-results %>% 
  group_by(Scale) %>% 
  do(BEF=lm(log(Biomass)~log(Richness),data=.))

BEF.df<-tidy(BEF_models,BEF)

BEF.df<-BEF.df %>% 
  filter(term=="log(Richness)")

ggplot(BEF.df,aes(x=Scale,y=estimate))+
  geom_point()+
  theme_bw()
ggsave("./figures/Slope by scale new.pdf")


Richness<-round(runif(n = 50,1,max=40))

reg.richV<-c(40,300,1000,5000)
for(j in 1:4){
  reg.rich<-reg.richV[j]
  for(r in 1:100){
    Richness<-round(runif(n = 50,1,max=40))
    SpeciesID<-matrix(0,50,reg.rich)
    Biomass<-c(sapply(Richness,FUN = function(x){(1/(x^gamma))*x}))
    for(i in 1:50){
      SpeciesID[i,sample(1:reg.rich,Richness[i],replace = F)]<-1
    }
    for(i in 1:50){
      patches<-sample(1:50,size = i,replace = F)
      if(r==1 & i == 1 & j==1){
        results<-data.frame(Regional_pool=reg.richV[j],Scale=i,Richness=sum(SpeciesID[patches,]>0,na.rm=T),Biomass=sum(Biomass[patches]))
      } else {
        if(i==1){
          results<-rbind(results,data.frame(Regional_pool=reg.richV[j],Scale=i,Richness=sum(SpeciesID[patches,]>0,na.rm=T),Biomass=sum(Biomass[patches])))
        } else {
          results<-rbind(results,data.frame(Regional_pool=reg.richV[j],Scale=i,Richness=sum(colSums(SpeciesID[patches,],na.rm=T)>0),Biomass=sum(Biomass[patches])))
        }
      }
    }
  }
}

ggplot(filter(results,Regional_pool==40,Scale==3),aes(x=Richness,y=Biomass))+
  geom_point()

BEF_models<-results %>% 
  group_by(Scale,Regional_pool) %>% 
  do(BEF=lm(log(Biomass)~log(Richness),data=.))

BEF.df<-tidy(BEF_models,BEF)

BEF.df<-BEF.df %>% 
  filter(term=="log(Richness)")

ggplot(BEF.df,aes(x=Scale,y=estimate,color=as.factor(Regional_pool)))+
  geom_point()+
  theme_bw()+
  scale_color_brewer(type="seq",palette="Set1",name="Regional pool")
ggsave("./figures/Slope by scale new 2.pdf")
