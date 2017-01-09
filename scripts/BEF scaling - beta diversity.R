library(ggplot2)
library(dplyr)
library(broom)
library(ggExtra)

all.results<-data.frame()
for(S_total in c(10,40,400,4000)){
  species_pool<-1:S_total
  
  species_comp<-matrix(NA,2000,10)
  species_comp[1,]<-sapply(1:10,function(x){sample(species_pool,1,replace = T)})
  S<-sapply(1:10,FUN = function(x){length(unique(species_comp[1,1:x]))})
  
  results<-data.frame(Scale=1:10,S=S,EF=1:10*(1^0.23),S_mean=1,beta_div="yes")
  
  for(i in 2:2000){
    species_comp[i,sample(1:10,1,replace=F)]<-sample(species_pool,1,replace = F)
    local_S<-apply(species_comp,2,function(x) {length(unique(x[!is.na(x)]))} )
    local_EF<-local_S^0.23
    reg_EF<-sapply(1:10,FUN = function(x){sum(local_EF[1:x])})
    reg_S<-sapply(1:10,FUN= function(x){
      region<-species_comp[,1:x]
      length(unique(c(region[!is.na(region)])))})
    results<-rbind(results, data.frame(Scale=1:10,S=reg_S,EF=reg_EF,S_mean=mean(local_S),beta_div="yes"))
  }
  
  local_S<-1:round(max(results$S))
  local_EF<-local_S^0.23
  reg_EF<-local_EF*10
  reg_S<-local_S
  results<-rbind(results, data.frame(Scale=10,S=reg_S,EF=reg_EF,S_mean=reg_S,beta_div="no"))
  results$species_pool<-S_total
  
  all.results<-rbind(all.results,results)
}

#all.results$Scale<-factor(all.results$Scale,levels = c(1:10,"10 equal"),ordered=T)

ggplot(all.results,aes(x=S,y=EF,color=as.factor(Scale),group=interaction(Scale,beta_div),linetype=beta_div))+
  geom_line()+
  facet_wrap(~species_pool,scale='free')+
  scale_color_grey(start=0.8,end=0,name="Scale")+
  xlab("Species richness")+
  ylab("Ecosystem function")+
  scale_linetype(name="Beta\ndiversity")+
  theme_bw()+
  removeGrid()
ggsave("./figures/Beta div EF scaling.png",width = 7,height=6)

ggplot(all.results,aes(x=S,y=EF,color=as.factor(Scale),group=interaction(Scale,beta_div),linetype=beta_div))+
  geom_line()+
  facet_wrap(~species_pool,scale='free')+
  scale_color_grey(start=0.8,end=0,name="Scale")+
  xlab("Species richness")+
  ylab("Ecosystem function")+
  scale_linetype(name="Beta\ndiversity")+
  theme_bw()+
  removeGrid()+
  scale_x_log10(breaks=c(1,10,100,1000))+
  scale_y_log10(breaks=c(1,10))
ggsave("./figures/Beta div EF scaling (log).png",width = 7,height=6)
