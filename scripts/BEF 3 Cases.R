library(tidyverse)
library(broom)
library(viridis)
library(ggExtra)
library(cowplot)

#simulate BEF slope relationships####
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}


Scale_fun<-function(r,scale = 5,bV = bV) {
  regional_rich.df<-do.call(rbind,lapply(1:3,function(case) {
    logit_newV<-c(0,0.01,0.025,0.05,0.1,0.25,0.5,0.75,1)
    #simulate local richness
    repeat {
      local_richness<-rep(round(rnorm(n = 1,mean=10,sd=3)),scale)
      if(mean(local_richness>0)==1) break
    }
    if(case==1){
      repeat {
        local_richness<-round(rnorm(n = scale,mean=10,sd=3))
        if(mean(local_richness>0)==1) break
      }
    }
    
    #estimate local biomass
    biomass<-sapply(X = bV, FUN = function(b){
      if(case == 2){
        if(scale==1){
          biomass<-local_richness^rnorm(scale,mean = b,sd = 0.11)
        } else{
          biomass<-sum(local_richness^rnorm(scale,mean = b,sd = 0.11))
        }
      } else {
        if(scale==1){
          biomass<-local_richness^b
        } else{
          biomass<-sum(local_richness^b)
        }
      }
    })
    
    #estimate regional richness
    if(case==3){
      regional_richness<-rep(local_richness[1],length(logit_newV))
      if(scale>1){
        for(j in 2:scale){
          if(local_richness[j]==1){
            regional_richness_hold<-regional_richness+sapply(logit_newV,function(k){rbinom(n = local_richness[j],size = 1,prob = logit2prob(5-k*regional_richness))})
            regional_richness<-apply(rbind(regional_richness,regional_richness_hold),2, max)
          } else {
            regional_richness_hold<-regional_richness+colSums(sapply(logit_newV,function(k){rbinom(n = local_richness[j],size = 1,prob = logit2prob(5-k*regional_richness))}))
            regional_richness<-apply(rbind(regional_richness,regional_richness_hold),2, max)
          }
        }
      }
      return(data.frame(Div_a = mean(local_richness), Div_g = regional_richness, Div_b = 1-(mean(local_richness)/regional_richness),biomass = rep(biomass,each = length(regional_richness)),local_b=rep(bV,each = length(regional_richness)), logit_new = logit_newV,scale = scale,case = case))
      
    } else {
      regional_richness<-sum(local_richness)
      return(data.frame(Div_a = mean(local_richness), Div_g = regional_richness, Div_b = 1-(mean(local_richness)/regional_richness),biomass = biomass,local_b=bV, logit_new = NA,scale = scale,case = case))
    }
  }))
  
  
  regional_rich.df$rep <- r
  return(regional_rich.df)
}

repsV<-100
reps1<-1:100
scales<-1:50
bV<-seq(-0.25,1.25, by = 0.25)
pb <- txtProgressBar(min = 0, max = max(scales), style = 3)
slopes.df<-data.frame()
for(r in 1:repsV){
  print(r)
  for(i in scales){
    setTxtProgressBar(pb, i)
    
    scale_BEF<-do.call(rbind, lapply(reps1, FUN = Scale_fun,scale=i,bV=bV))
    
    slopes.df.temp<-scale_BEF %>%
      filter(Div_g>0) %>% 
      group_by(scale, logit_new,case,local_b) %>% 
      do(lm.scale=lm(log(biomass)~log(Div_g),data=.)) %>%
      tidy(lm.scale) %>% 
      select(term,estimate,scale, logit_new,local_b,case) %>% 
      spread(key = term,value = estimate)
    
    names(slopes.df.temp)[5:6]<-c("a","b")
    
    R2<-scale_BEF %>%
      filter(Div_g>0) %>% 
      group_by(scale, logit_new,local_b,case) %>% 
      do(lm.scale=lm(log(biomass)~log(Div_g),data=.)) %>%
      glance(lm.scale) %>% 
      select(r.squared,scale, logit_new,local_b,case)
    
    slopes.df.temp$r.squared<-R2$r.squared
    
    slopes.df<-bind_rows(slopes.df,slopes.df.temp)
  }
}

save(slopes.df,file = "./data/Simulated scale BEF.RData")

slopes_mean<-slopes.df %>% 
  group_by(scale,logit_new,local_b,case) %>% 
  summarise(lower = quantile(b,probs = 0.25), upper = quantile(b,probs = 0.75),b = mean(b))

ggplot(filter(slopes_mean, case!=3),aes(x=scale,y=b, group = local_b, fill=as.character(local_b)))+
  geom_hline(yintercept = bV, linetype=2)+
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.5)+
  facet_wrap(~case)+
  #geom_hline(yintercept = 0, linetype=1,size=0.2)+
  geom_line()+
  #geom_smooth(method = "gam",formula = y ~ s(x,k=4),se=F)+
  #geom_point()+
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 3), size = 1)+
  #geom_smooth(method="lm", formula = y~poly(x,3),se=F)+
  scale_fill_brewer(palette = "Set1",name="local b")+
  theme_bw()+
  removeGrid()
ggsave("./figures/Case I & II.pdf",width = 8,height = 6)


ggplot(filter(slopes_mean,local_b==0.25, case==3),aes(x=scale,y=b,group=as.factor(logit_new),fill=as.factor(logit_new),color=as.factor(logit_new)))+
  geom_hline(yintercept = 0.25, linetype=2)+
  geom_hline(yintercept = 0, linetype=1,size=0.2)+
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.5,color=NA)+
  geom_line(size=1)+
  #geom_smooth(method = "gam",formula = y ~ s(x,k=4),se=F)+
  #geom_point()+
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 3), size = 1)+
  #geom_smooth(method="lm", formula = y~poly(x,3),se=F)+
  scale_color_hue(name=expression(italic(p)),c = 100,l = 60)+
  theme_bw()+
  removeGrid()
ggsave("./figures/Case III.pdf",width = 6,height = 4)

ggplot(filter(slopes_mean, case==3),aes(x=scale,y=b,group=as.factor(logit_new),fill=as.factor(logit_new),color=as.factor(logit_new)))+
  geom_hline(aes(yintercept = local_b), linetype=2)+
  geom_hline(yintercept = 0, linetype=1,size=0.2)+
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.5,color=NA)+
  geom_line()+
  facet_wrap(~local_b)+
  scale_color_hue(name=expression(italic(p)),c = 100,l = 60)+
  scale_fill_hue(guide=F)+
  theme_bw()+
  removeGrid()
ggsave("./figures/Case III - all.pdf",width = 8,height = 6)


ggplot(filter(slopes.df,local_b==0.5),aes(x=Scale,y=b,group=beta_dist,color=beta_dist))+
  geom_hline(yintercept = 0.5, linetype=2)+
  geom_hline(yintercept = 0, linetype=1)+
  geom_line()+
  #geom_smooth(method = "gam",formula = y ~ s(x,k=4),se=F)+
  #geom_point()+
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 3), size = 1)+
  #geom_smooth(method="lm", formula = y~poly(x,3),se=F)+
  scale_color_hue(name="Beta distribution")+
  facet_wrap(~logit_new)+
  theme_bw()+
  removeGrid()
ggsave("./figures/Simulated slope by scale.pdf",width = 10,height = 6.5)

ggplot(slopes.df,aes(x=Scale,y=r.squared,group=beta_dist,color=beta_dist))+
  geom_line()+
  #geom_point()+
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 3), size = 1)+
  #geom_smooth(method="lm", formula = y~poly(x,3),se=F)+
  scale_color_hue(name="Beta distribution")+
  facet_wrap(~logit_new)+
  theme_bw()+
  removeGrid()
ggsave("./figures/Simulated R2 by scale.pdf",width = 10,height = 6.5)




logits<-data.frame(sapply(logit_newV,function(k) logit2prob(5+k*1:200)))
names(logits)<-logit_newV
logits$Richness<-1:200

logits<-logits %>% 
  gather(key = Parameter,value = Probability,-Richness)

ggplot(logits, aes(x=Richness,y=Probability, color=as.factor(as.numeric(Parameter))))+
  geom_line()+
  theme_bw()+
  removeGrid()+
  scale_x_log10()
ggsave("./figures/Richness accumulation prob.pdf",height = 3,width=5)


#Jensen's inequality case I and II####
func_change<-data.frame(sapply(bV, FUN = function(x) { return((seq(2,41)^x)-(seq(1,40)^x))}))

names(func_change)<-bV
func_change$initial_S<-1:40

fig2a<-func_change %>% 
  gather(key = b, value = function_change,-initial_S) %>% 
  ggplot(aes(x=initial_S,y=function_change,color=b))+
  geom_line(size=1)+
  theme_bw()+
  removeGrid()+
  scale_color_hue(guide=F)+
  ylab("Change in functioning")+
  xlab("Initial richness")+
  ylim(-0.2,3.2)


func_change<-data.frame(sapply(bV, FUN = function(x) { return((5^seq(from = x-0.1,x+0.1, length=40))
                                                              -(4^seq(from = x-0.1,x+0.1, length=40)))}))
b.df<-data.frame(sapply(bV, FUN = function(x) { return(seq(from = x-0.1,x+0.1, length=40))}))

names(func_change)<-bV
names(b.df)<-bV


func_change<-func_change %>% 
  gather(key = mean_b, value = function_change) 
func_change$b<-unlist(c(b.df))

fig2b<-ggplot(func_change,aes(x=b,y=function_change,color=mean_b))+
  geom_line(size=1)+
  theme_bw()+
  removeGrid()+
  ylab("Change in functioning")+
  scale_color_hue(name="mean b")+
  xlab("Local b")+
  ylim(-0.2,3.2)+
  theme(legend.justification=c(-0.1,1.05), legend.position=c(0,1))

plot_grid(fig2a,fig2b,ncol = 2,labels = c("a)","b)"))
ggsave("./figures/Figure 2.pdf",width = 8,height=4)


#estimate probability of that new species is unique####
reps<-99
pb <- txtProgressBar(min = 0, max = reps, style = 3)
NewSp.df<-data.frame()
for(r in 1:reps){
  setTxtProgressBar(pb, r)
  sampled_plots<-sample(unique(filter(Bmass,NumSp<16)$Plot),size = 20,replace=F)
  
  sampled_region<-Bmass %>%
    filter(NumSp!=16) %>%
    filter(Year==2014) %>% 
    filter(Plot %in% sampled_plots) %>%
    mutate(Sim_biomass = exp(5.5)*NumSp^b) %>% 
    left_join(Com_mat,by = "Plot") %>% 
    ungroup() %>% 
    select(Achmi:Sornu)
  
  
  for(j in 2:20){
    sampled_region_sub<-sampled_region %>% 
      select(which(sampled_region[j,]==1))
    New<-c(sampled_region_sub[j,]>0 & colSums(sampled_region_sub[1:(j-1),])==0)
    NewSp.df<-rbind(NewSp.df,data.frame(New=New,Richness=sum(colSums(sampled_region[1:(j-1),])>0),Scale=j,rep=r))
  }
}

model<-glm(New~Richness,family = binomial(link = "logit"),data=NewSp.df)
summary(model)


