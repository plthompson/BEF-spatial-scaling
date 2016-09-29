library(ggplot2)
library(ggExtra)
library(cowplot)
library(viridis)
library(dplyr)

scales=c(1,2,5,10,25,50,100,200,400)
richness_scaling<-seq(0,1,by=0.1)#c(0,0.1,0.25,0.5,0.75,1)
gamma_dist<-c("Fixed", "Beta 5,5", "Beta 2,5")
for(g in 1:length(gamma_dist)){
  reps<-500
  for(r in 1:reps){
    print(r)
    for(i in scales){
      if(i==1){
        if(g==1){
          gamma<-0.77
        } 
        if(g==2){
          gamma<-rbeta(n=10000,shape1 = 5,shape2=5)
          gamma<-c(gamma-mean(gamma)+0.77)[1]
        }
        if(g==3){
          gamma<-1-rbeta(n=10000,shape1 = 2,shape2=5)
          gamma<-c(gamma-mean(gamma)+0.77)[1]
        }
        Bmass<-sapply(1:40,FUN = function(x){(1/(x^gamma))*x})
        hold.df<-data.frame(Production=Bmass*i,Richness=(1:40)*i,Scale=i,Richness_scaling=0,g_dist=gamma_dist[g])
        results.df<-hold.df
      } else {
        for(j in richness_scaling){
          if(g==1){
            gamma<-rep(0.77,i)
          } 
          if(g==2){
            gamma<-rbeta(n=10000,shape1 = 5,shape2=5)
            gamma<-c(gamma-mean(gamma)+0.77)[1:i]
          }
          if(g==3){
            gamma<-1-rbeta(n=10000,shape1 = 2,shape2=5)
            gamma<-c(gamma-mean(gamma)+0.77)[1:i]
          }
          Bmass<-sapply(1:40,FUN = function(x){(1/(x^gamma))*x})
          Bmass<-colSums(Bmass)
          results.df<-rbind(results.df,data.frame(Production=Bmass,Richness=1:40*i^j,Scale=i,Richness_scaling=j,g_dist=gamma_dist[g]))
        }
      }
    }
    Slope.df.temp<-results.df %>%
      group_by(Scale,Richness_scaling,g_dist) %>%
      do(mod = lm(log(Production) ~ log(Richness), data = .)) %>%
      mutate(Slope = summary(mod)$coeff[2]) %>%
      select(-mod) %>% 
      ungroup() %>% 
      group_by(Scale,Richness_scaling,g_dist) %>% 
      summarise(Slope=mean(Slope))
    if(r==1 & g==1){
      Slope.df<-Slope.df.temp
    } else {
      Slope.df<-rbind(Slope.df, Slope.df.temp)
    }
  }
}

ggplot(Slope.df,aes(x=Scale,y=Slope,group=Richness_scaling,color=Richness_scaling)) +
  geom_point()+
  facet_wrap(~g_dist)+
  geom_smooth(method='lm',formula=y~x,se = F)+
  scale_x_log10()+
  scale_color_viridis()+
  geom_abline(slope = 0,intercept=0.23, lty=2)
ggsave("./figures/BEF by scale and gamma distribution.pdf", width = 10,height=5)

ggplot(Slope.df,aes(x=Scale,y=Slope,group=Richness_scaling,color=Richness_scaling)) +
  facet_wrap(~g_dist)+
  geom_smooth(method='lm',formula=y~x,se = F)+
  scale_x_log10()+
  scale_color_viridis()+
  geom_abline(slope = 0,intercept=0.23, lty=2)

Slope_means<-Slope.df%>%
  group_by(Scale,Richness_scaling,g_dist) %>% 
  summarise(Slope=mean(Slope))

ggplot(Slope_means, aes(x=Scale,y=Slope,group=Richness_scaling,color=Richness_scaling))+
  geom_point()+
  facet_grid(~g_dist)+
  scale_color_viridis()+
  scale_x_log10()
ggsave("./figures/BEF by scale and gamma distribution.pdf", width = 10,height=5)



ggplot(Slope.df,aes(x=Scale,y=Slope,group=Richness_scaling,color=Richness_scaling)) +
  facet_wrap(~g_dist)+
  geom_smooth(method='lm',formula=y~x,se = F)+
  scale_x_log10()+
  scale_color_viridis()+
  geom_abline(slope = 0,intercept=0.23, lty=2)



slopes.df<-results.df %>%
  group_by(Scale,Richness_scaling,rep) %>%
  do(mod = lm(log(Production) ~ log(Richness), data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2]) %>%
  select(-mod)

ggplot(slopes.df,aes(x=Scale,y=Slope,group=Richness_scaling,color=Richness_scaling))+
  geom_point()+
  geom_smooth(method='lm',formula=y~poly(x,2))+
  scale_color_viridis()+
  scale_fill_viridis()+
  scale_x_log10()+
  geom_abline(intercept=0.3,slope=0)

slope_means<-slopes.df %>% 
  ungroup() %>% 
  group_by(Scale,Richness_scaling) %>% 
  summarise(Slope=mean(Slope))

ggplot(slope_means,aes(x=Scale,y=Slope,group=Richness_scaling,color=Richness_scaling)) +
  geom_point()+
  geom_smooth(method='lm',formula=y~poly(x,2))+
  scale_x_log10()+
  scale_color_viridis()


g1<-ggplot(results.df,aes(x=Richness,y=Production,group=interaction(Scale,Richness_scaling),color=Richness_scaling))+
  geom_line(size=1)+
  scale_color_viridis()+
  #scale_color_brewer(palette = "YlGnBu")+
  #removeGrid()+
  theme(legend.justification=c(1,0), legend.position=c(1,0))

g2<-ggplot(results.df,aes(x=Richness,y=Production,group=interaction(Scale,Richness_scaling),color=Richness_scaling))+
  geom_line(size=1)+
  scale_color_viridis(guide=F)+
  removeGrid()+
  scale_x_log10(breaks=c(1,10,100,1000,10000))+
  scale_y_log10(breaks=c(1,10,100,1000))

plot_grid(g1, g2, labels=c("a", "b"), ncol = 2, nrow = 1) 
ggsave(filename = "./figures/BEF slopes - scaling SR.pdf",height = 4,width=8.5)
