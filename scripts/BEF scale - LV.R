library(tidyverse)
library(broom)
library(ggExtra)
library(drc)

meta_dyn_model_loss_gain<-function(disp = 0.1, type = "competitive", spatial_env = TRUE, temporal_env = FALSE,loss=TRUE, global_disp=FALSE,spatial_period=2){
  patches<-60
  species<-100
  interval<-150
  Tmax<-interval*species
  
  Env_perform<-function(env,z,zmax=z+5,sig_p){
    Tmat<-matrix(rep(env,each=species),species,length(env)) 
    wT<-exp(-((Tmat-z)/2*rep(sig_p,each=patches))^2)
    #wT2<-1-((Tmat-z)/(z-zmax))^2
    #wT[Tmat>=z]<-wT2[Tmat>=z]
    wT[wT<0]<-0
    wT<-wT-1
    return(wT)
  }
  
  #environmental trait
  z<-seq(-0.5,1.5,length=species)#runif(n = species, min = -0.5,max = 1.5)
  sig_p<-0.05 #rate of performance decay as z differs from local environment
  envPeriod<-5000 # temporal period of environmental sin wave
  
  if(spatial_env == TRUE){
    env <- 0.5 * (sin(seq(0, spatial_period*pi, length.out = patches) + (2 * pi * 1) / envPeriod) + 1) # local environmental conditions at time t
    plot(env, type="l")
    A<-t(Env_perform(env,z,sig_p = sig_p)*2000)
    
  } else {
    A<-matrix(0,nrow=60, ncol=40)
  }
  
  r<-1 #intrinsic rate of increase
  
  bdiag1<- -.2*(z*0.2+0.8)
  if(type == "neutral"){
    B=matrix(-0.2,nrow=species,ncol=species)
    N<-sapply(X = 1:species,FUN = function(X){
      hold<-rep(0,patches)
      hold[sample(1:patches,size = 5,replace=F)]<-5
      return(hold)})
  }
  
  if(type == "competitive"){
    B<-matrix(runif(species*patches)*-0.05,nrow=species,ncol=species)#matrix(runif(species*patches)*-0.15,nrow=species,ncol=species)
    N<-matrix(1,ncol=species,nrow=patches)
  }
  
  if(type == "priority"){
    B<-matrix(rnorm(species*species,mean = bdiag1+0.01,sd=0.005),nrow=species,ncol=species)
    N<-sapply(X = 1:species,FUN = function(X){
      hold<-rep(0,patches)
      hold[sample(1:patches,size = 5,replace=F)]<-5
      return(hold)})
  }
  diag(B)<-bdiag1
  
  colonizeOrder<-sample(1:species,size = species,replace = F)
  
  if(loss==TRUE){
    N<-matrix(1,ncol=species,nrow=patches)
  } else{
    N<-matrix(0,ncol=species,nrow=patches)
    N[,colonizeOrder[1]]<-1
  }
  
  Nsave<-array(NA,dim=c(patches,species,Tmax))
  Esave<-matrix(NA,nrow = Tmax, ncol = patches)
  
  colV<-seq(from = interval,to=Tmax,by=interval)
  
  disp_matrix<-matrix(0,nrow = patches,ncol = patches)
  for(p in 1:(patches-1)){
    disp_matrix[p,p+1]<-0.5
    disp_matrix[p+1,p]<-0.5
  }
  disp_matrix[1,patches]<-0.5
  disp_matrix[patches,1]<-0.5
  
  if(global_disp==TRUE){
    disp_matrix<-matrix(1/(patches-1),nrow = patches,ncol = patches)
    diag(disp_matrix)<-0
  }
  
  pb <- txtProgressBar(min = 0, max = Tmax, style = 3)
  for(i in 1:Tmax){
    if(sum(i==colV)==1){
      if(loss==TRUE){
        N[,colonizeOrder[which(colV==i)+1]]<-0
      } else{
        N[,colonizeOrder[which(colV==i)+1]]<-1
      }
    }
    if(temporal_env==TRUE){
      if(spatial_env==TRUE){
        env <- 0.5 * (sin(seq(0, 2*pi, length.out = patches) + (2 * pi * i) / envPeriod) + 1) # local environmental conditions at time t
      } else{
        env <- 0.5 * (sin(rep(0,patches) + (2 * pi * i) / envPeriod) + 1) # local environmental conditions at time t
        
      }
      
      A<-t(Env_perform(env,z,sig_p = sig_p)*2000)
    }
    
    if(spatial_env==T | temporal_env==T){
      Esave[i,]<-env
    }
    
    Nsave[,,i]<-N
    
    
    
    Nt<-N*exp(r+N%*%B+A)
    Nt<-Nt+disp*disp_matrix%*%Nt-Nt*disp
    
    #Nt<-Nt+rnorm(n = patches*species,mean = 0,sd=Nt*0.01)*(Nt>0)
    Nt[Nt<0.05]<-0
    N<-Nt
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(list(Nsave,Esave))  
}


# note - play with the strength of interspecific competition coefficients
# play with dispersal rate to alter species sorting and mass effects
# inter = intra gives local neutrality, 

reps<-30
results.df<-data.frame()
dispV<-c(0,0.5)
spatial_periodV<-c(2,4,8,16)
for(disp in dispV){
  for(spP in spatial_periodV){
    for(r in 1:reps){
      print(paste("rep = ",r," dispersal = ", disp, " period = ",spP,sep=""))
      output<-meta_dyn_model_loss_gain(disp = 0.1,loss = TRUE, spatial_period = spP)
      Nsave<-output[[1]]
      Esave<-output[[2]]
      
      matplot(t(Nsave[1,,seq(1,dim(Nsave)[3],by=10)]), type='l')
      
      plot(rowSums(Nsave[,,149])~Esave[dim(Nsave)[3],], ylab="Biomass", xlab="Environment")
      plot(rowSums(Nsave[,,149]>0)~Esave[dim(Nsave)[3],], ylab="Species richness", xlab="Environment")
      plot(rowSums(Nsave[,,149]>0)~rowSums(Nsave[,,149]), xlab="Biomass", ylab="Species richness")
      
      plot(rowSums(Nsave[,,149]>0),rowSums(Nsave[,,149]), ylab="Biomass", xlab="Species richness")
      
      
      plot(apply(Nsave[,,seq(299,dim(Nsave)[3]-1,by=300)],3,sum), type="b")
      plot(colSums(Nsave[1,,seq(299,dim(Nsave)[3]-1,by=300)]), type="b")
      
      results<-data.frame()
      for(x in 1:dim(Nsave)[1]){
        if(x==1){
          com_data<-Nsave[x,,seq(299,dim(Nsave)[3]-1,by=300)]
          EF<-colSums(com_data)
          S<-colSums(com_data>0)
        } else{
          com_data<-Nsave[1:x,,seq(299,dim(Nsave)[3]-1,by=300)]
          EF<-apply(com_data,3,sum)
          S<-colSums(apply(com_data,3,colSums)>0)
        }
        results<-rbind(results,data.frame(Rpool = 1:dim(Nsave)[2],EF=EF,S=S,scale=x,rep=r,env_period=spP,dispersal= disp))
      }
      results.df<-rbind(results.df,results)
    }
  }
}



ggplot(filter(results.df,rep==1, scale %in% c(1,5,10,15,20,30,40,50),S>0),aes(x=S,y=EF))+
  facet_wrap(~scale)+
  geom_point()+
  geom_smooth(method="lm")+
  scale_x_log10()+
  scale_y_log10()

slopes<-results.df %>%
  filter(S>0) %>% 
  group_by(scale,rep,env_period,dispersal) %>% 
  do(lm1=lm(log(EF)~log(S),data=.)) %>% 
  tidy(lm1) %>% 
  filter(term=="log(S)")

R2<-results.df %>% 
  group_by(scale,rep,env) %>% 
  filter(S>0) %>% 
  do(lm1=lm(log(EF)~log(S),data=.)) %>% 
  glance(lm1)

slopes_mean <-slopes %>% 
  group_by(scale,env_period,dispersal) %>% 
  summarise(estimate=mean(estimate,na.rm=T))

ggplot(slopes,aes(x=scale,y=estimate,group=rep))+
  facet_grid(dispersal~env_period)+
  geom_line()+
  geom_line(data=slopes_mean,aes(x=scale,y=estimate,group=interaction(dispersal,env_period)),col=2,size=1)+
  scale_color_hue(guide=F)+
  theme_bw()+
  removeGrid()+
  ylab("b")
ggsave("./figures/model2.pdf",width = 10,height = 8)

ggplot(R2,aes(x=scale,y=r.squared,group=rep))+
  facet_wrap(~env)+
  geom_line()+
  scale_color_hue(guide=F)+
  theme_bw()+
  removeGrid()


mm.slopes<-results.df %>% 
  group_by(scale,rep,env_period,dispersal) %>% 
  filter(S>0) %>% 
  do(mm=drm (EF ~ S, fct = MM.2(),data=.))

mm.slopes.coef<-data.frame()
for(i in 1:nrow(mm.slopes)){
  hold<-mm.slopes[i,]
  mm.slopes.coef<-rbind(mm.slopes.coef,data.frame(scale=hold$scale,rep=hold$rep,env_period=hold$env_period,dispersal=hold$dispersal ,d=coef(hold$mm[[1]])[1],e=coef(hold$mm[[1]])[2]))
}

mm.slopes.coef.mean<-mm.slopes.coef %>% 
  filter(e<100) %>% 
  group_by(scale,env_period,dispersal) %>% 
  summarise(e=mean(e,na.rm=T))

ggplot(filter(mm.slopes.coef,e<100),aes(x=scale,y=e,group=rep))+
  facet_grid(dispersal~env_period)+
  geom_line()+
  geom_line(data=mm.slopes.coef.mean,aes(x=scale,y=e,group=interaction(dispersal,env_period)),col=2,size=1)+
  theme_bw()+
  removeGrid()+
  ylab("Richness at half saturation")
ggsave("./figures/model2 mm.pdf",width = 10,height = 8)

#within a rep####
results<-data.frame()
for(x in 1:dim(Nsave)[1]){
  for(y in 1:(dim(Nsave)[1]-x+1)){
    if(x==1){
      com_data<-Nsave[y,,seq(299,dim(Nsave)[3]-1,by=300)]
      EF<-colSums(com_data)
      S<-colSums(com_data>0)
    } else{
      com_data<-Nsave[y:(y+x-1),,seq(299,dim(Nsave)[3]-1,by=300)]
      EF<-apply(com_data,3,sum)
      S<-colSums(apply(com_data,3,colSums)>0)
    }
    results<-rbind(results,data.frame(Rpool = 1:dim(Nsave)[2],EF=EF,S=S,scale=x,rep=1,subunit=y,env="autocorrelated"))
  }
}


results %>% 
  group_by(scale,subunit,env) %>% 
  do(lm1=lm(log(EF)~log(S),data=.)) %>% 
  tidy(lm1) %>% 
  filter(term=="log(S)") %>% 
  ggplot(aes(x=scale,y=estimate,color=subunit))+
  facet_wrap(~env)+
  geom_point()+
  #scale_color_hue(guide=F)+
  theme_bw()+
  removeGrid()+
  ylab("b")
ggsave("./figures/model2 - 1 rep.pdf",width = 7,height = 5)

