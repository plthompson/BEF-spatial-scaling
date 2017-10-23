library(tidyverse)
library(broom)
library(viridis)
library(ggExtra)
library(cowplot)

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
  summarise(Biomass=mean(Biomass)) #should I be taking the average biomass across all samples in a year?

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
ggsave("./figures/Cedar creek local data.pdf", width=20,height=4)

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
ggsave("./figures/Cedar creek - BEF_scale.pdf", height=3.5,width = 8)


ggplot(filter(Regional_div,Year==2014),aes(x=Div_g, y=Biomass/scale))+
  geom_point()+
  facet_wrap(~scale)+
  geom_smooth(method = 'lm')+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  xlab("Species richness")
ggsave("./figures/BEF scale CC - 2014.pdf", width = 10,height = 8)


slopes<-Regional_div %>% 
  group_by(scale,Year) %>% 
  do(BEF_slope = lm(log(Biomass) ~ log(Div_g), data = .)) %>% 
  tidy(BEF_slope) %>% 
  filter(term=="log(Div_g)")

ggplot(slopes,aes(x=scale,y = estimate, color=as.character(Year), group=Year))+
  geom_point()+
  #geom_line()+
  geom_smooth(method = 'lm', formula = y~poly(x,3), se = F)+
  scale_color_discrete(name="Year")+
  theme_bw()+
  #scale_y_continuous(breaks = seq(0,1.5,by=0.25),minor_breaks = seq(0.25,1.25,by=0.25))+
  ylab("b")+
  xlab("Scale m^2")+
  removeGrid()
ggsave("./figures/Cedar creek - BEF_scale.pdf", height=6,width = 8)

#plot beta diversity by space relationship
Beta_mean<-Regional_div %>% 
  filter(Year==2005) %>%
  group_by(scale) %>% 
  summarise(mean_beta = mean(Div_b))


ggplot(filter(Regional_div, Year==2005),aes(x=scale,y = Div_b))+
  geom_point()+
  geom_line(data=Beta_mean,aes(x=scale,y=mean_beta))+
  #geom_line()+
  #geom_smooth(method = 'lm', se = F)+
  theme_bw()+
  xlim(0,1700)+
  ylab("beta diversity")+
  xlab("Scale m^2")+
  removeGrid()


Regional_div %>% 
  group_by(scale,Year) %>% 
  summarise(Biomass=max(Biomass)/min(Biomass), Richness=max(Div_g)/min(Div_g)) %>% 
  gather(key = Type,value = Ratio,Biomass:Richness) %>% 
  ggplot(aes(x=scale,y=Ratio,color=Type))+
  geom_point()+
  facet_wrap(~Year)+
  scale_y_log10()+
  theme_bw()+
  removeGrid()+
  ylab("max/min")
ggsave("./figures/Ratios.pdf",width = 8,height=6)

Regional_div %>% 
  filter(scale==81*20) %>%
  ggplot(aes(x=Div_g,y=Biomass))+
  geom_point()+
  facet_grid(~Year)+
  stat_smooth(method="lm")

Regional_div %>% 
  filter(scale==81) %>%
  ggplot(aes(x=Div_g,y=Biomass))+
  geom_point()+
  facet_wrap(~Year)+
  stat_smooth(method="lm")+
  scale_x_log10(breaks=c(1,2,4,8))+
  scale_y_log10()+
  xlab("Species richness")+
  theme_bw()+
  removeGrid()
ggsave("./figures/Local BEF.pdf",width = 8,height=8)

BEF_scale_effect<-slopes %>% 
  group_by(Year) %>% 
  do(scale_slope=lm(estimate~scale,data=.)) %>% 
  tidy(scale_slope) %>% 
  filter(term=="scale") %>% 
  select(Year,estimate) 

names(BEF_scale_effect)<-c("Year","Scale_slope")

ggplot(BEF_scale_effect,aes(x=Year,y=estimate))+
  geom_point()+
  theme_bw()

local.ab<-Bmass %>% 
  group_by(Year) %>% 
  do(year.lm=lm(log(Biomass)~log(NumSp),data=.)) %>% 
  tidy(year.lm) %>% 
  select(Year,term,estimate) %>% 
  spread(key = term,value = estimate)
names(local.ab)<-c("Year","a","b")

local.ab %>% 
  gather(key = parameter,value=value,a:b) %>% 
  ggplot(aes(x=Year,y=value))+
  geom_point()+
  facet_wrap(~parameter, scale="free")+
  theme_bw()+
  removeGrid()
ggsave("./figures/Local_ab.pdf", width = 7,height = 4)


left_join(BEF_scale_effect,local.ab) %>% 
  ggplot(aes(x=b,y=Scale_slope))+
  geom_point()



#test with no fixed species pool####
plot_id<-unique(filter(Bmass,NumSp<16)$Plot)
reps<-999
scales<-c(1,2,4,6,8,10,12,14,16,18,20)
Regional_div_infSp<-data.frame()
pb <- txtProgressBar(min = 0, max = reps, style = 3)
for(i in scales){
  print(paste("Scale = ", i, sep=""))
  for(r in 1:reps){
    sampled_plots<-sample(plot_id,size = i,replace=F)
    
    sampled_region<-Bmass %>%
      filter(NumSp!=16) %>%
      group_by(Year) %>% 
      filter(Plot %in% sampled_plots) %>% 
      summarise(Biomass=sum(Biomass),NumSp=sum(NumSp))
    
    sampled_region$scale <- i*9^2
    sampled_region$rep <- r
    Regional_div_infSp<-rbind(Regional_div_infSp,sampled_region)
    setTxtProgressBar(pb, r)
  }
  close(pb)
}

ggplot(filter(Regional_div_infSp,Year==2014),aes(x=NumSp, y=Biomass))+
  geom_point()+
  facet_wrap(~scale)+
  geom_smooth(method = 'lm')+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  xlab("Species richness")
ggsave("./figures/BEF scale CC - InfSp - 2014.pdf", width = 10,height = 8)

slopes_infSp<-Regional_div_infSp %>% 
  group_by(scale,Year) %>% 
  do(BEF_slope = lm(log(Biomass) ~ log(NumSp), data = .)) %>% 
  tidy(BEF_slope) %>% 
  filter(term=="log(NumSp)")

ggplot(slopes_infSp,aes(x=scale,y = estimate, color=as.character(Year), group=Year))+
  geom_point()+
  #geom_line()+
  geom_smooth(method = 'lm', formula = y~poly(x,1), se = F)+
  scale_color_discrete(name="Year")+
  theme_bw()+
  xlim(0,1700)+
  scale_y_continuous(breaks = seq(0,1.5,by=0.5),minor_breaks = seq(0.25,1.25,by=0.25))+
  ylab("b")+
  xlab("Scale m^2")+
  removeGrid()


#test with fixed ab####
local.ab<-Bmass %>% 
  group_by(Year) %>% 
  do(year.lm=lm(log(Biomass)~log(NumSp),data=.)) %>% 
  tidy(year.lm) %>% 
  select(Year,term,estimate) %>% 
  spread(key = term,value = estimate)

names(local.ab)<-c("Year","a","b")

reps<-500
scales<-c(1,seq(5,50,by=5))
Regional_div.fixed.ab<-data.frame()
pb <- txtProgressBar(min = 0, max = reps, style = 3)
for(i in scales){
  print(paste("Scale = ", i, sep=""))
  for(r in 1:reps){
    sampled_plots<-sample(plot_id,size = i,replace=TRUE)
    
    sampled_region<-Bmass %>% 
      filter(Plot %in% sampled_plots) %>% 
      left_join(local.ab,by="Year") %>%
      group_by(Year) %>% 
      mutate(Biomass_est=exp(a)*NumSp^b) %>% 
      left_join(Com_mat,by = "Plot") %>%
      select(Year,Biomass_est:Sornu) %>%
      summarise_all(funs(sum)) %>% 
      gather(key = Species, value = Abundance, Achmi:Sornu) %>%
      group_by(Year) %>% 
      summarise(Div_a = sum(Abundance)/i ,Div_g = sum(Abundance>0), Div_b = 1-Div_a/Div_g,Biomass = mean(Biomass_est))
    
    sampled_region$scale <- i
    sampled_region$rep <- r
    Regional_div.fixed.ab<-rbind(Regional_div.fixed.ab,sampled_region)
    setTxtProgressBar(pb, r)
  }
  close(pb)
}

save(Regional_div.fixed.ab,file = "./data/fixed_ab.RData")

ggplot(filter(Regional_div.fixed.ab,Year==2014),aes(x=Div_g,y=Biomass))+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(~scale)+
  scale_x_log10()+
  scale_y_log10()

slopes.fixed.ab<-Regional_div.fixed.ab %>% 
  group_by(scale,Year) %>% 
  do(BEF_slope = lm(log(Biomass) ~ log(Div_g), data = .)) %>% 
  tidy(BEF_slope) %>% 
  filter(term=="log(Div_g)")

R2.fixed.ab<-Regional_div.fixed.ab %>%
  group_by(scale,Year) %>% 
  do(BEF_slope = lm(log(Biomass) ~ log(Div_g), data = .)) %>% 
  glance(BEF_slope) %>% 
  select(r.squared)

ggplot(slopes.fixed.ab,aes(x=scale,y = estimate, color=as.character(Year), group=Year))+
  geom_point()+
  #geom_line()+
  stat_smooth(method = "gam", formula = y ~ s(x, k = 5), size = 1,se=F)+
  #geom_smooth(method = 'lm', formula = y~poly(x,3), se = F)+
  scale_color_discrete(name="Year")+
  theme_bw()+
  scale_y_continuous(breaks = seq(0,1.5,by=0.5),minor_breaks = seq(0.25,1.25,by=0.25))+
  ylab("b")+
  xlab("Scale m^2")+
  removeGrid()
ggsave("./figures/Cedar creek - BEF_scale - fixed ab.pdf", height=6,width = 8)

ggplot(R2.fixed.ab,aes(x=scale,y = r.squared, color=as.character(Year), group=Year))+
  geom_point()+
  #geom_line()+
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 5), size = 1,se=F)+
  #geom_smooth(method = 'lm', formula = y~poly(x,3), se = F)+
  scale_color_discrete(name="Year")+
  theme_bw()+
  scale_y_continuous(breaks = seq(0,1.5,by=0.5),minor_breaks = seq(0.25,1.25,by=0.25))+
  ylab("b")+
  xlab("Scale m^2")+
  removeGrid()


#bootstrap local b####
boot.b<-data.frame()
reps<-5
pb <- txtProgressBar(min = 0, max = reps, style = 3)
for(r in 1:reps){
  sampled_local<-Bmass %>%
    filter(NumSp!=16) %>%
    group_by(Year) %>% 
    sample_n(size =length(plot_id),replace=T) %>% 
    do(local.lm=lm(log(Biomass)~log(NumSp),data=.)) %>% 
    tidy(local.lm) %>% 
    filter(term=="log(NumSp)") %>% 
    select(Year,estimate) %>% 
    mutate(Scale=1)
  
  sampled.region.df<-data.frame()
  for(i in 1:999){
    sampled_plots<-sample(plot_id,size = 10,replace=F)
    sampled_region_10<-Bmass %>%
      filter(NumSp!=16) %>%
      group_by(Year) %>% 
      filter(Plot %in% sampled_plots) %>% 
      left_join(Com_mat,by = "Plot") %>%
      select(Year,Biomass:Sornu) %>%
      summarise_each(funs(sum)) %>% 
      gather(key = Species, value = Abundance, Achmi:Sornu) %>%
      group_by(Year) %>% 
      summarise(Div_a = sum(Abundance)/i ,Div_g = sum(Abundance>0), Div_b = 1-Div_a/Div_g,Biomass = mean(Biomass)) %>% 
      mutate(Scale=10)
    
    sampled_plots<-sample(plot_id,size = 20,replace=F)
    sampled_region_20<-Bmass %>%
      filter(NumSp!=16) %>%
      group_by(Year) %>% 
      filter(Plot %in% sampled_plots) %>% 
      left_join(Com_mat,by = "Plot") %>%
      select(Year,Biomass:Sornu) %>%
      summarise_each(funs(sum)) %>% 
      gather(key = Species, value = Abundance, Achmi:Sornu) %>%
      group_by(Year) %>% 
      summarise(Div_a = sum(Abundance)/i ,Div_g = sum(Abundance>0), Div_b = 1-Div_a/Div_g,Biomass = mean(Biomass)) %>% 
      mutate(Scale=20)
    sampled.region.df<-bind_rows(sampled.region.df,sampled_region_10,sampled_region_20)
  }
  
  sampled_region<-sampled.region.df %>% 
    group_by(Year,Scale) %>% 
    do(reg.lm=lm(log(Biomass)~log(Div_g),data=.)) %>% 
    tidy(reg.lm) %>% 
    filter(term=="log(Div_g)") %>% 
    select(Year,estimate,Scale)
  
  boot.b<-bind_rows(boot.b,sampled_local,sampled_region)
  setTxtProgressBar(pb, r)
}
close(pb)

ggplot(boot.b,aes(x=as.character(Year),y=estimate,fill=Scale,group=interaction(Year,Scale)))+
  geom_violin(draw_quantiles = c(0.025,0.5,0.975))+
  theme_bw()+
  removeGrid()+
  scale_fill_viridis()
ggsave("./figures/Bootstrapped local b.pdf", height=5, width=6)



ggplot(Regional_div,aes(x=Div_g, y=Biomass))+
  geom_point()+
  facet_wrap(~scale, scales="free")+
  geom_smooth(method = 'lm',formula = y~poly(x,2))



ggplot(Regional_div,aes(x=Div_g, y=Biomass, group=scale,color=scale))+
  #geom_point(pch=1)+
  facet_wrap(~Year)+
  geom_smooth(method = 'lm', se=F)+
  scale_color_viridis()+
  theme_bw()+
  removeGrid()+
  scale_x_log10()+
  scale_y_log10()
ggsave("./figures/Cedar creek - BEF_scale 2.pdf", height=6,width = 8)




slopes2<-Regional_div %>% 
  group_by(scale,Year) %>% 
  filter(Biomass!=0) %>% 
  do(BEF_slope = lm(log(Biomass) ~ log(Div_g), data = .)) %>% 
  tidy(BEF_slope) %>% 
  select(scale,Year,term,estimate) %>% 
  spread(key = term,value=estimate)

names(slopes2)<-c("scale","Year","a","b")


plot(Biomass~Div_g,data=filter(Regional_div,scale==20,year==2012), type='n', ylim=c(0,7000), xlim=c(1,18))
for(s in 1:20){
  sl<-slopes2 %>% 
    filter(Year==2012,scale==s)
  curve(exp(sl$a[1])*x^sl$b[1],from = 1,to = 18,add=T, col=col_fun(20)[s])
}



ggplot(slopes,aes(x=scale,y = estimate, group=Year))+
  geom_point()+
  #geom_line()+
  facet_grid(~Year)+
  geom_smooth(method = 'lm', formula = y~poly(x,1), se = F, color="forestgreen")+
  #scale_color_viridis(breaks=c(2001,2005,2010,2014))+
  theme_bw()+
  xlim(0,1700)+
  scale_y_continuous(breaks = seq(0,1.5,by=0.5),minor_breaks = seq(0.25,1.25,by=0.25))+
  ylab("BEF slope")+
  xlab("Scale m^2")+
  removeGrid()
ggsave("./figures/CC BEF scale - Years - all equal b.pdf",width = 20,height=4)


save(slopes,Regional_div,file = "./data/BEF_cedar_creek.RData")

