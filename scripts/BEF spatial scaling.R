#Lotka-Volterra simulation model 
#simulates a 20x20 patch metacommunity
#simulates BEF relationship at multiple spatial scales

#code written by Patrick Thompson in collaboration with Andrew Gonzalez
#September 2016

#packages####
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggExtra)
library(mvtnorm)
library(ncf)
library(viridis)

#functions####
matsplitter<-function(M, r1, c1) {
  r<-sqrt(r1)
  c<-sqrt(c1)
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  cv
}

dist.matrix <- function(side)
{
  row.coords <- rep(1:side, times=side)
  col.coords <- rep(1:side, each=side)
  row.col <- data.frame(row.coords, col.coords)
  D <- dist(row.col, method="euclidean", diag=TRUE, upper=TRUE)
  D <- as.matrix(D)
  return(D)
}

# function that simulates the autocorrelated 2D array with a given side,
# and with exponential decay given by lambda
# (the mean mu is constant over the array, it equals to global.mu)
cor.surface <- function(side, global.mu, lambda)
{
  D <- dist.matrix(side)
  # scaling the distance matrix by the exponential decay
  SIGMA <- exp(-lambda*D)
  mu <- rep(global.mu, times=side*side)
  # sampling from the multivariate normal distribution
  M <- matrix(nrow=side, ncol=side)
  M[] <- rmvnorm(1, mu, SIGMA)
  return(M)
}

#model function####
BEF_spat_scale<-function(disp=0,Int_type="Comp",nSpecies=40){
  burn<-500
  Tmax<-burn #length of simulation 1000 is enough to reach equilibrium
  
  nCom_rows<-20 #number of rows in landscape matrix
  nCom_cols<-20 #number of columns in landscape matrix
  nCom<-nCom_cols*nCom_rows #total number of patches in landscape matrix
  
  #number of species in each trophic level
  if(Int_type=="Trophic"){
    nPrey<-nSpecies*0.5
    nHerb<-nSpecies*0.3
    nPred<-nSpecies*0.2
  } else{
    nPrey<-nSpecies
    nHerb<-0
    nPred<-0
  }
  
  #Identity vectors for each trophic level
  preyV<-1:nPrey
  if(nHerb>0){
    herbV<-(nPrey+1):(nPrey+nHerb)
    pred<-(nSpecies-nPred+1):(nSpecies)
  }
  
  #species environmental niches
  #Opt1<-c(seq(1-10,nCom_rows+10,length=nPrey), seq(1,nCom_rows,length=nHerb),seq(1,nCom_rows,length=nPred))
  Opt1<-c(runif(n = nPrey,min = 1-10,max=nCom_rows+10),runif(n = nHerb,min = 1-10,max=nCom_rows+10),runif(n = nPred,min = 1-10,max=nCom_rows+10))
  T_Norm<-apply(t(Opt1),2,dnorm,sd=50,x=seq(1,nCom_rows))*300
  A1<-(T_Norm-max(T_Norm))
  
  #Opt2<-c(seq(1-10,nCom_rows+10,length=nPrey)[sample(nPrey,replace = F)], seq(1,nCom_rows,length=nHerb)[sample(nHerb,replace = F)],seq(1,nCom_rows,length=nPred)[sample(nPred,replace = F)])
  Opt2<-c(runif(n = nPrey,min = 1-10,max=nCom_rows+10),runif(n = nHerb,min = 1-10,max=nCom_rows+10),runif(n = nPred,min = 1-10,max=nCom_rows+10))
  T_Norm<-apply(t(Opt2),2,dnorm,sd=50,x=seq(1,nCom_rows))*300
  A2<-(T_Norm-max(T_Norm))
  
  traits<-cbind(Opt1,Opt2)
  row.names(traits)<-1:nSpecies
  
  #interaction matricies####
  #competitive
  weight=0.125 #weight interaction strength
  
  if(Int_type=="Comp" | Int_type=="NoInt"){
    b11=-.15 #mean interspecific competition strength
    bdiag1=-.02 #intraspecific competition strength
    BB=b11*matrix(runif(nPrey*nPrey),nPrey,nPrey)
    BB=weight*BB
    diag(BB)<-bdiag1
    B<-BB
    if(Int_type=="NoInt"){B=diag(diag(BB))}
  } else {if(Int_type=="Mixed"){
    b11=-.15 #mean interspecific competition strength
    bdiag1=-.2 #intraspecific competition strength
    BB=b11*matrix(runif(nPrey*nPrey),nPrey,nPrey)
    BB=weight*BB
    diag(BB)<-bdiag1
    BI<-BB
    
    BB<-matrix(-1,nSpecies,nSpecies)
    int.n<-sum(BB[upper.tri(BB)])*-1
    BB[upper.tri(BB)][sample(int.n, replace=F)<=(0.35*int.n)]<-0.5
    BB[lower.tri(BB)][t(BB)[lower.tri(BB)]>0][sample(0.35*int.n, replace=F)<(0.10*int.n)]<-0.5
    B<-BB*-BI
  } else {if(Int_type == "Trophic"){
    b11=-0.1
    b12=-0.3
    b21=0.1
    b23=-.1
    b32=.08
    bdiag1=-.02
    bdiag2=-.015
    
    B11=b11*matrix(runif(nPrey*nPrey),nPrey,nPrey)
    B12=b12*matrix(runif(nPrey*nHerb),nPrey,nHerb)
    B13=matrix(0,nPrey,nPred)
    B21=b21*matrix(runif(nHerb*nPrey),nHerb,nPrey)
    B22=matrix(0,nHerb,nHerb)
    B23=b23*matrix(runif(nHerb*nPred),nHerb,nPred)
    B31=matrix(0,nPred,nPrey)
    B32=b32*matrix(runif(nPred*nHerb),nPred,nHerb)
    B33=matrix(0,nPred,nPred)
    BB=rbind(cbind(B11 ,B12, B13),cbind(B21,B22, B23),cbind(B31, B32, B33))
    BB=weight*BB
    diag(BB)<-bdiag1
    diag(BB[(nPrey+nHerb+1):nSpecies,(nPrey+nHerb+1):nSpecies])<-bdiag2
    B<-BB
  }}}
  
  C1<-c(rep(0.05,nPrey),rep(0,nSpecies-nPrey))
  
  #dispersal####
  disp_mat<-matrix(1/(nCom-1),nCom,nCom)
  diag(disp_mat)<-0

  
  autocorrV<-c(1,0.95,0.9,0.75,0.6,0.5,0)
  Morans<-data.frame(MoransI=NA,autocorr=rep(autocorrV,each=14),Distance=NA)
  for(ac in autocorrV){
    
    #environmental gradients
    Env<-data.frame(Env1=seq(1:nCom_rows),Env2=rep(1:nCom_cols,each=nCom_rows))
    Env_spatial<-Env
    
    rand_env<-(1-ac)*nCom
    rand_patches<-sample(x = nCom,size = rand_env,replace = F)
    re_random_patches<-sample(x = rand_env,size = rand_env,replace=F)
    
    Env[rand_patches,]<-Env[rand_patches,][re_random_patches,]
    
    x.coord <- rep(1:20, each=20)
    y.coord <- rep(1:20, times=20)
    xy <- data.frame(x.coord, y.coord)
    
    # all paiwise euclidean distances between the cells
    xy.dist <- dist(xy)
    
    # PCNM axes of the dist. matrix (from 'vegan' package)
    
    # using 8th PCNM axis as my atificial z variable
    z.value <- matrix(Env$Env1,20,20)
    
    ncf.cor <- correlog(x.coord, y.coord, c(z.value),increment = 2,resamp = 1)
    Morans$MoransI[Morans$autocorr==ac]<-ncf.cor$correlation
    Morans$Distance[Morans$autocorr==ac]<-ncf.cor$mean.of.class

    
    # if(lambda=="random"){
    #   Env <- data.frame(Env1=c(decostand(matrix(rnorm(nCom), nrow=nCom_rows, ncol=nCom_rows),method="range")),
    #                     Env2=c(decostand(matrix(rnorm(nCom), nrow=nCom_rows, ncol=nCom_rows),method="range")))
    # } else {  Env <- data.frame(Env1=c(decostand(cor.surface(side=nCom_rows, lambda=as.numeric(lambda), global.mu=0),method = "range")),
    #                             Env2=c(decostand(cor.surface(side=nCom_rows, lambda=as.numeric(lambda), global.mu=0),method = "range")))
    # }
    # Env<-1+Env*19
    
    #model####
    X=array(NA,dim=c(nCom,nSpecies,Tmax))
    X[,,1]<-10
    X_outside<-rep(0,nCom*nSpecies)
    X_outside[1:round(nSpecies*nCom*0.0005)]<-0.005
    X_inside<-rep(0,nCom*nSpecies)
    X_inside[1:round(nSpecies*nCom*0.1)]<-1
    hold<-X[,,1]
    
    for(l in 1:(Tmax-1)){
      X[,,l+1]<-X[,,l]*exp(rep(C1,nCom)+X[,,l]%*%B+A1[Env$Env1,]+A2[Env$Env2,])+(disp_mat%*%X[,,l])*disp-X[,,l]*disp
      X[,,l+1][(X[,,l+1]<10^-2.5)]<-0
      #print(l)
    }
    X_final<-X[,,Tmax]
    
    sizes<-c(1,4,16,25,100,400)
    for(i in 1:length(sizes)){
      mats<-matsplitter(matrix(1:400,20,20),r1=sizes[i],c1 =sizes[i])
      for(j in 1:dim(mats)[3]){
        X_sub<-X_final[c(mats[,,j]),]
        if(i==1){
          data.sub<-data.frame(scale=sizes[i],SR=sum(X_sub>0),Biomass=sum(X_sub))
        } else{
          data.sub<-data.frame(scale=sizes[i],SR=sum(colSums(X_sub)>0),Biomass=sum(X_sub))
        }
        data.sub$autocorr<-ac
        if(i==1 & j==1 & ac == autocorrV[1]){
          results.df<-data.sub
        } else {results.df<-rbind(results.df,data.sub)}
      }
    }
  }
  return(results.df)
}

#run model over range of gamma diversity####
G_div<-c(5,10,20,40,80)
reps<-10
for(r in 1:reps){
  print(r)
  for(g in G_div){
    hold<-BEF_spat_scale(nSpecies = g)
    hold$Gamma<-g
    hold$rep<-r
    if(g==G_div[1] & r==1){
      results.df<-hold
    } else {results.df<-rbind(results.df,hold)}
  }
}

results.df$autocorr_f<-as.factor(results.df$autocorr)

ggplot(Morans,aes(x=Distance,y=MoransI,group=autocorr,color=factor(autocorr)))+
  geom_point()+
  geom_line()+
  scale_color_brewer(palette = "YlGnBu",name="Env.\nautocorrelation")+
  #theme_bw()+
  removeGrid()
ggsave("./figures/MoransI.pdf",width = 6,height = 4)


ggplot(results.df,aes(x=SR,y=Biomass,group=autocorr_f,color=autocorr_f))+
  geom_point()+
  facet_wrap(~scale,scales="free")+
  #geom_smooth(method = 'lm',formula = y ~ poly(x,3))+
  theme_bw()+
  scale_color_brewer(palette = "YlGnBu",name="Env.\nautocorrelation")+
  removeGrid()+
  xlab("Species richness")
ggsave(filename = "./figures/BEF curves - raw.pdf",width = 13,height = 4)

ggplot(filter(results.df,SR!=0),aes(x=SR,y=Biomass,group=autocorr_f,color=autocorr_f))+
  geom_point()+
  facet_grid(~scale)+
  geom_smooth(method = 'lm',formula = y ~ x)+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  removeGrid()+
  scale_color_brewer(palette = "YlGnBu",name="Env.\nautocorrelation")+
  xlab("Species richness")
ggsave(filename = "./figures/BEF curves.pdf",width = 13,height = 4)

ggplot(filter(results.df,autocorr_f==0),aes(x=SR,y=Biomass,group=scale,color=scale))+
         geom_smooth(method = 'lm',formula = y ~ x)+
         theme_bw()+
         scale_x_log10()+
         scale_y_log10()+
         removeGrid()+
         xlab("Species richness")

slopes.df<-results.df %>%
  filter(SR!=0) %>% 
  group_by(scale,autocorr_f,rep) %>% # You can add here additional grouping variables if your real data set enables it
  do(mod = lm(log(Biomass) ~ log(SR), data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2]) %>%
  select(-mod)

ggplot(slopes.df,aes(x=scale,y=Slope,group=autocorr_f,color=autocorr_f, fill=autocorr_f))+
  #geom_point(pch=21,aes(fill=autocorr_f),color="black")+
  #theme_bw()+
  scale_x_log10()+
  geom_smooth(method = 'lm',formula = y~poly(x,2))+
  removeGrid()+
  scale_color_brewer(palette = "YlGnBu",name="Env.\nautocorrelation")+
  scale_fill_brewer(palette = "YlGnBu",name="Env.\nautocorrelation")+
  ylab("BEF slope")+
  xlab("Spatial scale")
ggsave(filename = "./figures/BEF slope by scale.pdf",width = 6,height = 4)

slopes.mean<-slopes.df %>% 
  ungroup() %>% 
  group_by(scale,autocorr_f) %>% 
  summarise(slope=mean(Slope))

ggplot(slopes.mean,aes(x=scale,y=autocorr_f,fill=slope))+
  geom_tile()+
  scale_x_log10()

pdf("./figures/Slope 3D surface.pdf",width = 8,height = 8)
wireframe(slope ~ log(scale)*autocorr_f, data = slopes.mean,
          #xlab = "X Coordinate (feet)", ylab = "Y Coordinate (feet)",
          #main = "Surface elevation data",
          drape = TRUE,
          colorkey = TRUE,
          col.regions = viridis(n = 100),
          screen = list(z = -20, x = -70, y = 5))
dev.off()


#at just one level of gamma####
ggplot(filter(results.df,Gamma==40),aes(x=SR,y=Biomass,group=autocorr_f,color=autocorr_f))+
  geom_point()+
  facet_wrap(~scale,scales = "free")+
  geom_smooth(method = 'lm',formula = y ~ poly(x,2))+
  theme_bw()+
  removeGrid()+
  xlab("Species richness")

ggplot(filter(results.df,Gamma==80),aes(x=SR,y=Biomass,group=autocorr_f,color=autocorr_f))+
  geom_point()+
  facet_wrap(~scale,scales = "free_y")+
  geom_smooth(method = 'lm',formula = y ~ x)+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  removeGrid()+
  xlab("Species richness")

slopes.df<-results.df %>%
  filter(Gamma==80,SR!=0) %>% 
  group_by(scale,autocorr_f) %>% # You can add here additional grouping variables if your real data set enables it
  do(mod = lm(log(Biomass) ~ log(SR), data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2]) %>%
  select(-mod)

ggplot(slopes.df,aes(x=scale,y=Slope,group=autocorr_f,color=autocorr_f))+
  geom_point()+
  theme_bw()+
  scale_x_log10()+
  geom_smooth(method = 'lm',formula = y~poly(x,2))+
  removeGrid()+
  xlab("Spatial scale")

ggplot(results.df,aes(x=autocorr,y=Biomass))+
  geom_point()+
  facet_wrap(~scale,scales = "free")

