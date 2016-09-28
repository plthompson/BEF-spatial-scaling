library(ggplot2)
library(ggExtra)
library(cowplot)

Bmass<-sapply(1:40,FUN = function(x){(1/(x^0.77))*x})
summary(lm(log(Bmass)~log(1:40)))
summary(lm(log(2*Bmass)~log(1:40)))
summary(lm(log(2*Bmass)~log(seq(from=2,to=80,by=2))))

scales=c(1,2,5,10,25,50,100,200,400)
turnoverV<-c(0,0.01,0.1,0.33,0.66,1)
for(i in scales){
  if(i==1){
    hold.df<-data.frame(Production=Bmass*i,Richness=(1:40)*i,Scale=i,Turnover=0)
    results.df<-hold.df
  } else {
    for(turn in turnoverV){
      results.df<-rbind(results.df,data.frame(Production=Bmass*i,Richness=seq(from=(1+turn*1*(i-1)),to=(40+40*turn*(i-1)),length=40),Scale=i,Turnover=turn))
    }
  }
}
results.df$Turnover<-as.factor(results.df$Turnover)

g1<-ggplot(results.df,aes(x=Richness,y=Production,group=interaction(Scale,Turnover),color=Turnover))+
  geom_line(size=1)+
  scale_color_brewer(palette = "YlGnBu")+
  removeGrid()+
  theme(legend.justification=c(1,0), legend.position=c(1,0))

g2<-ggplot(results.df,aes(x=Richness,y=Production,group=interaction(Scale,Turnover),color=Turnover))+
  geom_line(size=1)+
  scale_color_brewer(palette = "YlGnBu",guide=F)+
  removeGrid()+
  scale_x_log10(breaks=c(1,10,100,1000,10000))+
  scale_y_log10(breaks=c(1,10,100,1000))

plot_grid(g1, g2, labels=c("a", "b"), ncol = 2, nrow = 1) 
ggsave(filename = "./figures/BEF slopes - new method.pdf",height = 4,width=8.5)
  


plot(1:40,Bmass,type='l',ylim=c(1,120),xlim=c(1,40*50),xlab="Species richness")
lines(1:40,Bmass*2)
lines(c(1:40)*2,Bmass*2)
lines(c(1:40)*10,Bmass*10)
lines(c(1:40),Bmass*10)
lines(c(1:40)*25,Bmass*25)
lines(c(1:40),Bmass*25)
lines(c(1:40)*50,Bmass*50)
lines(c(1:40),Bmass*50)

plot(1:40,Bmass,type='l',ylim=c(1,120),xlim=c(1,40*50),xlab="Species richness",log='xy')
lines(1:40,Bmass*2)
lines(c(1:40)*2,Bmass*2)
lines(c(1:40)*10,Bmass*10)
lines(c(1:40),Bmass*10)
lines(c(1:40)*25,Bmass*25)
lines(c(1:40),Bmass*25)
lines(c(1:40)*50,Bmass*50)
lines(c(1:40),Bmass*50)
