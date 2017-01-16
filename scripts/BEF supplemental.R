#Code by Patrick Thompson
#In collaboration with Andrew Gonzalez and Forest Isbell 
#January 2017

#packages####
library(dplyr)
library(broom)
library(ggplot2)
library(ggExtra)
library(cowplot)

#Figure S1####
spec_acc<-c(0,0.5,1) #three rates of species accumulation
results<-data.frame()
for(s in spec_acc){
  accum<-(1:10)^s #species accumulation across 10 spatial scales
  SR<-sapply(1:20,function(x){x*accum}) #simulate local richness gradients from 1:20 and correspond richness at 10 spatial scales
  EF<-t(sapply(1:10,function(x){x*SR[1,]^0.26})) #calculate ecosystem functioning at all scales and levels of local richness
  
  results<-rbind(results,data.frame(SR=c(SR),EF=c(EF),Scale=1:10,spec_acc=s))
}

results<-results %>% 
  mutate(slope=log(EF)/log(SR)) #calculate BEF slope

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
Local_BEF_slope<-seq(-0.2,0.7,length=1000) #produce range of local b  
function_5<-15^Local_BEF_slope #caculate the functioning of 15 species at each level of b
function_4<-14^Local_BEF_slope #calculate the functioning of 14 species at each level of b
function_change<-function_5-function_4 #calculate the difference
inequality.df<-data.frame(BEF_slope=Local_BEF_slope,Function_change=function_change)

Local_SR<-seq(1,40,length=1000) #produce range of intial local richness
function_1<-Local_SR^0.26 #calculate corresponding EF
function_2<-(Local_SR+1)^0.26 #calculate EF with one additional species
function_change<-function_2-function_1 #calculate the difference
inequality.df2<-data.frame(Local_SR=Local_SR,Function_change=function_change)

S2a<-ggplot(inequality.df,aes(x=BEF_slope,y=Function_change))+
  geom_line()+
  theme_bw()+
  removeGrid()+
  xlab(expression(paste("Local ", italic("b"))))+
  ylab("Change in function")
  

S2b<-ggplot(inequality.df2,aes(x=Local_SR,y=Function_change))+
  geom_line()+
  theme_bw()+
  removeGrid()+
  xlab("Initial local species richness")+
  ylab("Change in function")

plot_grid(S2a,S2b,ncol = 2,labels = c("a)","b)"))
ggsave("./figures/Figure S2.pdf",width = 8,height=4)

#Figure S3####
#variation in slopes
sites<-50 #number of local patches
species<-40 #number of species

slopes.df<-data.frame()
for(i in 1:1000){
  slopes<-rnorm(n = sites,mean = 0.26, sd = 0.11) #draw local b from normal distribution
  L_SR<-seq(1:species) #gradient of local richness
  L_EF<-sapply(L_SR,FUN = function(x){x^slopes}) #calculate local EF for each level of richness with each local b
  
  Reg_EF<-apply(L_EF,2,function(x){sapply(1:sites,FUN = function(y) {sum(x[1:y])})}) #calculate EF at larger spatial scales
  slopes<-apply(Reg_EF,1,function(x){coef(lm(log(x)~log(c(1:species))))[2]}) #calculate b at larger spatial scales
  
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
  ylab(expression(italic("b")))+
  xlab("Spatial scale")
#variation in local richness
slopes2.df<-data.frame()
for(i in 1:1000){
  slopes<-0.26 #set local b
  L_SR<-round(runif(sites,min = 1,max=20)) #generate random local SR
  L_SR<-sapply(1:20, FUN = function(x){(L_SR+x)}) #diversity change
  
  L_EF<-L_SR^0.26 #calculate local functioning
  
  Reg_EF<-apply(L_EF,2,function(x){sapply(1:sites,FUN = function(y) {sum(x[1:y])})}) #calculate regional functioning
  Reg_S<-apply(L_SR,2,function(x){sapply(1:sites,FUN = function(y) {sum(x[1:y])})}) #calculate regional SR
  slopes<-sapply(1:sites,function(x){coef(lm(log(Reg_EF[x,])~log(Reg_S[x,])))[2]}) #calculate b at each spatial scale
  
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
  ylab(expression(italic("b")))+
  xlab("Spatial scale")

plot_grid(S3a,S3b,ncol = 2,labels = c("a)","b)"))
ggsave("./figures/Figure S3.pdf",width = 8,height=4)
