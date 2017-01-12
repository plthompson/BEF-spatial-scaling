#packages####
library(dplyr)
library(broom)
library(ggplot2)
library(ggExtra)
library(cowplot)

#Figure S1####
spec_acc<-c(0,0.5,1)
results<-data.frame()
for(s in spec_acc){
  accum<-(1:10)^s
  SR<-sapply(1:20,function(x){x*accum})
  EF<-t(sapply(1:10,function(x){x*SR[1,]^0.26}))
  
  results<-rbind(results,data.frame(SR=c(SR),EF=c(EF),Scale=1:10,spec_acc=s))
}

results<-results %>% 
  mutate(slope=log(EF)/log(SR))

results$spec_acc<-paste("Sp. accum. rate = ",results$spec_acc,sep="")

S1a<-ggplot(results,aes(x=SR,y=EF,group=as.factor(Scale),color=as.factor(Scale)))+
  facet_wrap(~spec_acc)+
  geom_line()+
  theme_bw()+
  removeGrid()+
  scale_color_grey(start=0.8,end=0,name="Scale")+
  xlab("Species richness")+
  ylab("Ecosystem functioning")
  
S1b<-ggplot(results,aes(x=SR,y=EF,group=as.factor(Scale),color=as.factor(Scale)))+
  facet_wrap(~spec_acc)+
  geom_line()+
  theme_bw()+
  removeGrid()+
  scale_color_grey(start=0.8,end=0,name="Scale")+
  xlab("Species richness")+
  ylab("Ecosystem functioning")+
  scale_x_log10()+
  scale_y_log10(breaks=c(1,5,10,20))

plot_grid(S1a,S1b,ncol = 1)
ggsave("./figures/Figure S1.pdf",width = 8,height=6)

#Figure S2####


#Figure S3####
#variation in slopes
sites<-50
species<-40

slopes.df<-data.frame()
for(i in 1:1000){
  slopes<-rnorm(n = sites,mean = 0.26, sd = 0.11)
  L_SR<-seq(1:species)
  L_bmass<-sapply(L_SR,FUN = function(x){x^slopes})
  
  Reg_bmass<-apply(L_bmass,2,function(x){sapply(1:sites,FUN = function(y) {sum(x[1:y])})})
  slopes<-apply(Reg_bmass,1,function(x){coef(lm(log(x)~log(c(1:species))))[2]})
  
  slopes.df<-rbind(slopes.df,data.frame(Slopes=slopes,Scale=1:sites))
}

slopes.summary<-slopes.df %>% 
  group_by(Scale) %>% 
  summarise(lower=quantile(Slopes,probs = 0.25),upper=quantile(Slopes,probs=0.75),Slopes=mean(Slopes))

S3a<-ggplot(slopes.summary,aes(x=Scale,y=Slopes))+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3,color=NA)+
  geom_line()+
  theme_bw()+
  removeGrid()+
  ylab("Biodiversity ecosystem function slope")

#variation in local richness
slopes2.df<-data.frame()
for(i in 1:1000){
  slopes<-0.26
  L_SR<-round(runif(sites,min = 1,max=20))
  L_SR<-sapply(1:20, FUN = function(x){(L_SR+x)})
  
  L_EF<-L_SR^0.26
  
  Reg_EF<-apply(L_EF,2,function(x){sapply(1:sites,FUN = function(y) {sum(x[1:y])})})
  Reg_S<-apply(L_SR,2,function(x){sapply(1:sites,FUN = function(y) {sum(x[1:y])})})
  slopes<-sapply(1:sites,function(x){coef(lm(log(Reg_EF[x,])~log(Reg_S[x,])))[2]})
  
  slopes2.df<-rbind(slopes2.df,data.frame(Slopes=slopes,Scale=1:sites))
}

slopes.summary2<-slopes2.df %>% 
  group_by(Scale) %>% 
  summarise(lower=quantile(Slopes,probs = 0.25),upper=quantile(Slopes,probs=0.75),slope_mean=mean(Slopes))

S3b<-ggplot(slopes.summary2,aes(x=Scale,y=slope_mean))+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3,color=NA)+
  geom_line()+
  theme_bw()+
  removeGrid()+
  ylab("Biodiversity ecosystem function slope")

plot_grid(S3a,S3b,ncol = 2,labels = c("a)","b)"))
ggsave("./figures/Figure S3.pdf",width = 8,height=4)

#Figure S4####
all.results<-data.frame()
for(S_total in c(10,40,400,4000)){
  species_pool<-1:S_total
  
  species_comp<-matrix(NA,2000,10)
  species_comp[1,]<-sapply(1:10,function(x){sample(species_pool,1,replace = T)})
  S<-sapply(1:10,FUN = function(x){length(unique(species_comp[1,1:x]))})
  
  results<-data.frame(Scale=1:10,S=S,EF=1:10*(1^0.26),S_mean=1,beta_div="yes")
  
  for(i in 2:2000){
    species_comp[i,sample(1:10,1,replace=F)]<-sample(species_pool,1,replace = F)
    local_S<-apply(species_comp,2,function(x) {length(unique(x[!is.na(x)]))} )
    local_EF<-local_S^0.26
    reg_EF<-sapply(1:10,FUN = function(x){sum(local_EF[1:x])})
    reg_S<-sapply(1:10,FUN= function(x){
      region<-species_comp[,1:x]
      length(unique(c(region[!is.na(region)])))})
    results<-rbind(results, data.frame(Scale=1:10,S=reg_S,EF=reg_EF,S_mean=mean(local_S),beta_div="yes"))
  }
  
  local_S<-1:round(max(results$S))
  local_EF<-local_S^0.26
  reg_EF<-local_EF*10
  reg_S<-local_S
  results<-rbind(results, data.frame(Scale=10,S=reg_S,EF=reg_EF,S_mean=reg_S,beta_div="no"))
  results$species_pool<-S_total
  
  all.results<-rbind(all.results,results)
}

S4a<-ggplot(all.results,aes(x=S,y=EF,color=as.factor(Scale),group=interaction(Scale,beta_div),linetype=beta_div))+
  geom_line()+
  facet_wrap(~species_pool,nrow = 1,scale="free")+
  scale_color_grey(start=0.8,end=0,name="Scale")+
  xlab("Species richness")+
  ylab("Ecosystem function")+
  scale_linetype(name="Beta\ndiversity")+
  theme_bw()+
  removeGrid()

S4b<-ggplot(all.results,aes(x=S,y=EF,color=as.factor(Scale),group=interaction(Scale,beta_div),linetype=beta_div))+
  geom_line()+
  facet_wrap(~species_pool,nrow = 1,scale="free")+
  scale_color_grey(start=0.8,end=0,name="Scale")+
  xlab("Species richness")+
  ylab("Ecosystem function")+
  scale_linetype(name="Beta\ndiversity")+
  theme_bw()+
  removeGrid()+
  scale_x_log10()+
  scale_y_log10()

plot_grid(S4a,S4b,ncol = 1)
ggsave("./figures/Figure S4.pdf",width = 10,height=7)
