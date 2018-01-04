# code to reproduce the results in Thompson et al. - The strength of the biodiversity-ecosystem function relationship depends on spatial scale
# all code written by Patrick Thompson
# this script will produce all figures except for Figure 4 - this is produced by the scrip called "Cedar creek BEF scaling.R"

#packages####
library(tidyverse)
library(broom)
library(viridis)
library(ggExtra)
library(cowplot)

#logit function
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#simulation function####
Scale_fun<-function(r,scale = 5,bV = bV) {
  regional_rich.df<-do.call(rbind,lapply(1:5,function(case) {
    logit_newV<-c(0,0.01,0.025,0.05,0.1,0.25,0.5,0.75)
    #simulate local richness
    repeat {
      local_richness<-rep(round(rnorm(n = 1,mean=10,sd=3)),scale)
      if(mean(local_richness>0)==1) break
    }
    if(case == 3){
      repeat {
        local_richness<-round(rnorm(n = scale,mean=10,sd=3))
        if(mean(local_richness>0)==1) break
      }
    }
    
    #estimate local biomass
    biomass<-sapply(X = bV, FUN = function(b){
      if(case == 4){
        if(scale==1){
          biomass<-local_richness^rnorm(scale,mean = b,sd = 0.11)
        } else{
          biomass<-sum(local_richness^rnorm(scale,mean = b,sd = 0.11))
        }
      } else{
        if(case == 2){
          repeat {
          a<-rnorm(n = scale,mean = 5, sd = 2)
          if (sum(a>0) ==scale) break
          }
          if(scale==1){
            biomass<-a*local_richness^b
          } else{
            biomass<-sum(a*local_richness^b)
          }
        } else {
          if(scale==1){
            biomass<-local_richness^b
          } else{
            biomass<-sum(local_richness^b)
          }
        }
      }
    }
    )
    
    #estimate regional richness
    if(case == 5){
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

#run simulation####
repsV<-100
reps1<-1:2000
scales<-c(1:5,seq(10,50,by=5))
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

save(slopes.df,bV,file = "./data/Simulated scale BEF.RData")

slopes_mean<-slopes.df %>% 
  group_by(scale,logit_new,local_b,case) %>% 
  summarise(lower = quantile(b,probs = 0.25), upper = quantile(b,probs = 0.75),b = median(b),
            lower.r2 = quantile(r.squared,probs = 0.25,na.rm=TRUE), upper.r2 = quantile(r.squared,probs = 0.75,na.rm=TRUE),r2 = median(r.squared,na.rm=TRUE)) %>% 
  mutate(casetext = paste("Case",as.roman(case)))
  

bV.df<-data.frame(sapply(bV, FUN = function(b) {
  seq(1,30)^b
}))
names(bV.df)<-bV
bV.df$species<-1:30

bV.df<-bV.df %>% 
  gather(key = b, value = Y, -species)

#Figure 1####
Fig.1a<-ggplot(filter(bV.df, b>0.75), aes(x = species, y = Y, group = b))+
  geom_line(size = 1)+
  #facet_wrap(~b, scales = "free")+
  ylab("Ecosystem functioning")+
  xlab("Species richness")+
  theme_bw()+
  removeGrid()

Fig.1b<-ggplot(filter(bV.df, b<1), aes(x = species, y = Y, group = b))+
  geom_line(size = 1)+
  #facet_wrap(~b, scales = "free")+
  ylab("Ecosystem functioning")+
  xlab("Species richness")+
  theme_bw()+
  removeGrid()

plot_grid(Fig.1a,Fig.1b, labels = c("a)", "b)"))
ggsave("./figures/Fig. 1.pdf", width = 9, height = 4)


#Figure S1####
ggplot(filter(slopes_mean, case != 5),aes(x=scale,y=b, group = local_b))+
  #geom_hline(yintercept = bV, linetype=3)+
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.5)+
  #facet_grid(.~case, scales = "free")+
  facet_grid(local_b~casetext, scales = "free")+
  #geom_hline(yintercept = 0, linetype=1,size=0.2)+
  geom_line()+
  #geom_smooth(method = "gam",formula = y ~ s(x,k=4),se=F)+
  #geom_point()+
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 3), size = 1)+
  #geom_smooth(method="lm", formula = y~poly(x,3),se=F)+
  #scale_fill_brewer(palette = "Set1",name= expression(paste("mean ", italic(b) [i]), sep=""))+
  theme_bw()+
  removeGrid()+
  xlab("Scale")+
  scale_y_continuous(breaks = seq(-0.5,1.5,by=0.01))+
  ylab(expression(italic(b) [A]))+
  theme(strip.background = element_rect(colour="NA", fill="NA"),strip.text.x = element_text(size = 10))
ggsave("./figures/Fig. S1.pdf",width = 7,height = 5)

#Figure 2####
Fig2a<-ggplot(filter(slopes_mean, case != 5),aes(x=scale,y=b-local_b, group = local_b, color = as.factor(local_b), fill = as.factor(local_b)))+
  #geom_ribbon(aes(ymin = lower-local_b,ymax = upper-local_b), alpha = 0.5, color= NA)+
  facet_grid(.~casetext, scales = "free")+
  geom_line(size = 1)+
  #scale_fill_brewer(palette = "Set1",name= expression(paste("mean ", italic(b) [i]), sep=""))+
  #scale_color_brewer(palette = "Set1",name= expression(paste("mean ", italic(b) [i]), sep=""))+
  scale_color_viridis(option = "D",discrete = TRUE,name= expression(paste("mean ", italic(b) [i]), sep=""))+
  theme_bw()+
  removeGrid()+
  xlab("Scale")+
  scale_y_continuous()+
  ylab(expression(italic(b) [A] - bar(italic(b) [i])))+
  theme(strip.background = element_rect(colour="NA", fill="NA"),strip.text.x = element_text(size = 10))


Fig2b<-ggplot(filter(slopes_mean, case != 5),aes(x=scale,y=r2,color = as.factor(local_b), fill = as.factor(local_b)))+
  geom_ribbon(aes(ymin = lower.r2,ymax = upper.r2), alpha = 0.5, color = NA)+
  facet_grid(.~casetext)+
  #geom_hline(yintercept = 0, linetype=1,size=0.2)+
  geom_line(size = 1)+
  scale_color_viridis(option = "D",discrete = TRUE,name= expression(paste("mean ", italic(b) [i]), sep=""))+
  scale_fill_viridis(option = "D",discrete = TRUE,name= expression(paste("mean ", italic(b) [i]), sep=""), guide = F)+
  theme_bw()+
  removeGrid()+
  xlab("Scale")+
  ylab(expression(paste("R"^"2")))+
  theme(strip.background = element_rect(colour="NA", fill="NA"),strip.text.x = element_text(size = 10))
plot_grid(Fig2a, Fig2b,nrow = 2, labels = c("a)", "b)"))
ggsave("./figures/Fig. 2.pdf",width = 8,height = 6)

  
#Figure 3####
Fig3.a<-ggplot(filter(slopes_mean,local_b==0.25, case == 5,logit_new<1),aes(x=scale,y=b,group=as.factor(logit_new),fill=as.factor(logit_new),color=as.factor(logit_new)))+
  geom_hline(yintercept = 0.25, linetype=2)+
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.5,color=NA)+
  geom_line(size=1)+
  #scale_color_hue(name=expression(italic(p)),c = 100,l = 60,h.start = 50)+
  #scale_fill_hue(name=expression(italic(p)),c = 100,l = 60,guide = FALSE,h.start = 50,)+
  scale_color_viridis(discrete = TRUE, option = "A", name=expression(italic(p)), guide=F, end = 0.95)+
  scale_fill_viridis(discrete = TRUE, option = "A", name=expression(italic(p)), guide = F,end = 0.95)+
  theme_bw()+
  removeGrid()+
  xlab("Scale")+
  ylab(expression(italic(b) [A]))


Fig3.b<-ggplot(filter(slopes_mean,local_b==0.25, case == 5, logit_new<1),aes(x=scale,y=r2,group=as.factor(logit_new),fill=as.factor(logit_new),color=as.factor(logit_new)))+
  geom_ribbon(aes(ymin = lower.r2,ymax = upper.r2), alpha = 0.5,color=NA)+
  geom_line(size=1)+
  #scale_color_hue(name=expression(italic(p)),c = 100,l = 60,h.start = 50)+
  #scale_fill_hue(name=expression(italic(p)),c = 100,l = 60,guide = FALSE,h.start = 50,)+
  scale_color_viridis(discrete = TRUE, option = "A", name=expression(italic(p)), end = 0.95)+
  scale_fill_viridis(discrete = TRUE, option = "A", name=expression(italic(p)), end = 0.95)+
  theme_bw()+
  removeGrid()+
  scale_y_continuous(limits = c(0,1))+
  xlab("Scale")+
  ylab(expression(paste("R"^"2")))+ 
  theme(legend.justification=c(0,0), legend.position=c(0,0),legend.background = element_rect(fill="NA", size=0.5, linetype="solid"),
        legend.key.size =  unit(0.15, "in"))

plot_grid(Fig3.a,Fig3.b,labels = c("a)","b)"))
ggsave("./figures/Fig. 3.pdf",width = 8,height = 3.5)


#Figure S3####
ggplot(filter(slopes_mean, case == 5),aes(x=scale,y=b,group=as.factor(logit_new),fill=as.factor(logit_new),color=as.factor(logit_new)))+
  geom_hline(yintercept = 0,size=0.5)+
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.5,color=NA)+
  geom_line(size=1)+
  facet_wrap(~local_b)+
  #scale_color_hue(name=expression(italic(p)),c = 100,l = 60,h.start = 50)+
  #scale_fill_hue(name=expression(italic(p)),c = 100,l = 60,guide = FALSE,h.start = 50,)+
  scale_color_viridis(discrete = TRUE, option = "A", name=expression(italic(p)))+
  scale_fill_viridis(discrete = TRUE, option = "A", name=expression(italic(p)))+
  theme_bw()+
  removeGrid()+
  xlab("Scale")+
  ylab(expression(italic(b)[A]))+
  theme(strip.background = element_rect(colour="NA", fill="NA"),strip.text.x = element_text(size = 10))
ggsave("./figures/Fig S3.pdf",width = 8,height = 6)


#logit exploration####
logit_newV<-c(0,0.01,0.025,0.05,0.1,0.25,0.5,0.75)
logits<-data.frame(sapply(logit_newV,function(k) logit2prob(5-k*1:200)))

names(logits)<-logit_newV
logits$Richness<-1:200

logits<-logits %>% 
  gather(key = Parameter,value = Probability,-Richness)

ggplot(logits, aes(x=Richness,y=Probability, color=as.factor(as.numeric(Parameter))))+
  geom_line(size=1)+
  theme_bw()+
  scale_color_viridis(name=expression(italic(p)), discrete = TRUE)+
  removeGrid()
ggsave("./figures/Richness accumulation prob.pdf",height = 3,width=5)


#Figure S2####
func_change<-data.frame(sapply(bV, FUN = function(x) { return((seq(2,41)^x)-(seq(1,40)^x))}))

names(func_change)<-bV
func_change$initial_S<-1:40

figS2a<-func_change %>% 
  gather(key = b, value = function_change,-initial_S) %>% 
  ggplot(aes(x=initial_S,y=function_change,color=b))+
  geom_line(size=1.2)+
  theme_bw()+
  removeGrid()+
  scale_color_viridis(discrete = TRUE,name=expression(paste("mean ",italic(b) [i])), guide = FALSE)+
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

figS2b<-ggplot(func_change,aes(x=b,y=function_change,color=mean_b))+
  geom_line(size=1.2)+
  theme_bw()+
  removeGrid()+
  scale_color_viridis(discrete = TRUE,name=expression(paste("mean ",italic(b) [i])))+
  ylab("Change in functioning")+
  xlab(expression(italic(b) [i]))+
  ylim(-0.2,3.2)+
  theme(legend.justification=c(-0.1,1.05), legend.position=c(0,1))

plot_grid(figS2a,figS2b,ncol = 2,labels = c("a)","b)"))
ggsave("./figures/Figure S2.pdf",width = 8,height=4)

#Figure S4####
logit_newV<-c(0,0.01,0.025,0.05,0.1,0.25,0.5,0.75)
SAR.df<-data.frame()
for(r in 1:25){
  for(scale in 1:50){
    repeat {
      local_richness<-rep(round(rnorm(n = 1,mean=10,sd=3)),scale)
      if(mean(local_richness>0)==1) break
    }
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
    SAR.df<-bind_rows(SAR.df,data.frame(Div_g = regional_richness, B1 = logit_newV, scale = scale, rep=r))
  }
}


SAR.df %>% 
  group_by(B1) %>% 
  do(lm1=lm(log(Div_g) ~ log(scale), data = .)) %>% 
  tidy(lm1) %>% 
  filter(term == "log(scale)") %>% 
  ggplot(aes(x=B1, y = estimate))+
  geom_point()+
  theme_bw()+
  removeGrid()+
  ylab("z")+
  xlab(expression(B [1]))
ggsave("./figures/Figure S4.pdf", width = 6, height = 4)


