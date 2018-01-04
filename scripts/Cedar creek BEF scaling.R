# code to reproduce the Figure 4 in Thompson et al. - The strength of the biodiversity-ecosystem function relationship depends on spatial scale
# all code written by Patrick Thompson

#packages
library(tidyverse)
library(broom)
library(viridis)
library(ggExtra)
library(cowplot)

#data is available from https://www.cedarcreek.umn.edu/research/data
CCe120<-read.csv(file = "./data/e120_Plant aboveground biomass data.csv",skip = 80,header=T)
head(CCe120)
CCe120_spec_list<-read.csv(file = "./data/e120_specieslist.csv",row.names = NULL)
head(CCe120_spec_list)

CCe120 <- CCe120 %>% 
  filter(!is.na(Year))

summary(CCe120)

Bmass<-CCe120 %>% 
  group_by(Year,Month,Plot,Date,NumSp,SpNum) %>% 
  summarise(Biomass=sum(Biomass..g.m2.,na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(Year,Plot,NumSp) %>% 
  summarise(Biomass=mean(Biomass))

Bmass<-Bmass %>% 
  group_by(Plot) %>%
  filter(n()==14)

Bmass_mono<-Bmass %>%
  ungroup() %>% 
  filter(NumSp==1) %>% 
  group_by(Year) %>% 
  summarise(Mean_mono = mean(Biomass))

Bmass<-left_join(Bmass,Bmass_mono)
Bmass<-Bmass %>% 
  mutate(Rel_Biomass = (Biomass/Mean_mono))

ggplot(Bmass,aes(x= Year,y=Biomass, group=Year))+
  geom_violin()+
  facet_grid(~NumSp)+
  scale_y_log10()

Bmass %>% 
  group_by(Year,NumSp) %>% 
  summarise(Biomass_mean=mean(Biomass),Biomass_sd=sd(Biomass)) %>% 
  ggplot(aes(x=Year,y=Biomass_sd,color=NumSp, group=NumSp))+
  geom_point()+
  geom_line()

local.ab<-Bmass %>% 
  group_by(Year) %>% 
  do(year.lm=lm(log(Biomass)~log(NumSp),data=.)) %>% 
  tidy(year.lm) %>% 
  select(Year,term,estimate) %>% 
  spread(key = term,value = estimate)

names(local.ab)[2:3]<-c("a","b")

ggplot(Bmass,aes(x=NumSp,y=Biomass)) +
  facet_grid(~Year)+
  geom_point()+
  geom_smooth(method = lm)+
  #scale_x_log10()+
  #scale_y_log10()+
  theme_bw()+
  removeGrid()
#ggsave("./figures/Cedar creek local data.pdf", width=20,height=4)

ggplot(Bmass,aes(x=NumSp,y=Biomass,color=as.character(Year),group=Year)) +
  #facet_wrap(~Year)+
  #geom_point()+
  geom_smooth(method = lm,se=T)+
  #scale_color_viridis()+
  scale_x_log10()+
  scale_y_log10()

local_slopes<-Bmass %>% 
  group_by(Year) %>% 
  filter(Biomass>0) %>% 
  do(slopes=lm(log(Rel_Biomass)~log(NumSp),data=.)) %>% 
  tidy(slopes) %>% 
  select(Year,term,estimate) %>% 
  spread(key = term,value = estimate)

names(local_slopes)<- c("Year","a", "b")

col_fun<-colorRampPalette(colors = c("red","green","blue"))
colV<-col_fun(14)

plot(Bmass$Rel_Biomass~Bmass$NumSp)
curve(exp(local_slopes$a[1])*x^local_slopes$b[1],from = 1,to = 16,add=T,col=colV[1])
for(i in 2:14){
  curve(exp(local_slopes$a[i])*x^local_slopes$b[i],from = 1,to = 16,add=T, col=colV[i])
}
legend("topleft",col=colV,legend = 2001:2014, lty=1,bty='n')


Com_mat<-CCe120 %>%
  select(Plot,Achmi:Sornu) %>%
  group_by(Plot)%>%
  filter(row_number()==1)

plot_id<-unique(filter(Bmass,NumSp<16)$Plot)

scales<-c(1,seq(2,30,by=2))
reps<-5000
results.df<-data.frame()
results<-data.frame()
pb <- txtProgressBar(min = 0, max = max(scales), style = 3)
for(rep in 1:1){
print(rep)
  for(i in scales){
  hold<-do.call(rbind,lapply(X = 1:reps,FUN = function(X){
    sampled_region<-Bmass %>% 
      filter(NumSp<16) %>% 
      group_by(Year) %>%
      sample_n(size = i,replace = TRUE) %>%
      left_join(local.ab,by="Year") %>% 
      mutate(Biomass_fake=exp(a)*NumSp^b) %>% 
      left_join(Com_mat,by = "Plot") %>%
      select(Year,Biomass:Sornu) %>% 
      summarise_all(funs(sum)) %>% 
      gather(key = Species, value = Abundance, Achmi:Sornu) %>% 
      group_by(Year) %>% 
      summarise(Div_a = sum(Abundance)/i ,Div_g = sum(Abundance>0),Biomass = mean(Biomass),Biomass_fake=mean(Biomass_fake),Richness_g_unique= sum(Abundance)) %>% 
      mutate(Div_b = Div_g/Div_a) %>% 
      mutate(scale=i*9^2)
    return(sampled_region)
  }))
  results<-rbind(results,hold)
  setTxtProgressBar(pb, i)
  close(pb)
}

results2<-results %>% 
  group_by(scale,Year) %>%
  mutate(Gamma_range=length(unique(Div_g))) %>% 
  filter(Gamma_range>2) %>% 
  do(lm.1=lm(log(Biomass)~log(Div_g),data=.)) %>% 
  tidy(lm.1) %>% 
  filter(term=="log(Div_g)")%>% 
  mutate(Type="actual")

results3<-results %>% 
  group_by(scale,Year) %>%
  mutate(Gamma_range=length(unique(Div_g))) %>% 
  filter(Gamma_range>2) %>% 
  do(lm.1=lm(log(Biomass)~log(Richness_g_unique),data=.)) %>% 
  tidy(lm.1) %>% 
  filter(term=="log(Richness_g_unique)")%>% 
  mutate(Type="unique species")


R2<-results %>% 
  group_by(scale,Year) %>% 
  mutate(Gamma_range=length(unique(Div_g))) %>% 
  filter(Gamma_range>2) %>% 
  do(lm.1=lm(log(Biomass)~log(Div_g),data=.)) %>% 
  glance(lm.1) %>% 
  select(r.squared) %>% 
  mutate(Type="actual")


R2.fixed<-results %>% 
  group_by(scale,Year) %>% 
  mutate(Gamma_range=length(unique(Div_g))) %>% 
  filter(Gamma_range>2) %>% 
  do(lm.1=lm(log(Biomass)~log(Richness_g_unique),data=.)) %>% 
  glance(lm.1) %>% 
  select(r.squared) %>% 
  mutate(Type="unique species")

results2<-bind_rows(results2,results3)
R2<-bind_rows(R2,R2.fixed)

results2<-left_join(results2,R2,by=c("Year","scale","Type"))

results.df<-bind_rows(results.df,results2)
}

save(results.df,file = "./data/BEF_cedar_creek.RData")

Fig4.a<-results.df %>%
  group_by(scale,Type) %>%
  mutate(years=length(unique(Year))) %>% 
  filter(years>10) %>% 
  summarise(lower=quantile(estimate,probs = 0.25,na.rm=TRUE),upper=quantile(estimate,probs = 0.75,na.rm=TRUE),b=quantile(estimate,probs = 0.5,na.rm=TRUE)) %>% 
  ggplot(aes(x=scale,y=b,group=Type,fill=Type,color=Type))+
  scale_fill_grey()+
  scale_color_manual(values = c(1,"grey30"))+
  geom_line()+
  theme_bw()+
  removeGrid()+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.5,color=NA)+
  xlab(expression(paste("Scale m"^"2")))+
  ylab(expression(italic(b) [A]))+
  theme(legend.justification=c(0,0.9), legend.position=c(0.05,0.95))

Fig4.b<-results.df %>%
  group_by(scale,Type) %>%
  mutate(years=length(unique(Year))) %>% 
  filter(years>10) %>% 
  summarise(lower=quantile(r.squared,probs = 0.25,na.rm=TRUE),upper=quantile(r.squared,probs = 0.75,na.rm=TRUE),r.squared=quantile(r.squared,probs = 0.5,na.rm=TRUE)) %>% 
  ggplot(aes(x=scale,y=r.squared,group=Type,fill=Type,color=Type))+
  scale_fill_grey(guide=FALSE)+
  scale_color_manual(values = c(1,"grey30"),guide=FALSE)+
  geom_line()+
  theme_bw()+
  removeGrid()+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.5,color=NA)+
  xlab(expression(paste("Scale m"^"2")))+
  ylab(expression(paste("R"^"2")))

plot_grid(Fig4.a,Fig4.b,labels = c("a)","b)"))
ggsave("./figures/Figure 4.pdf", height=3.5,width = 8)
