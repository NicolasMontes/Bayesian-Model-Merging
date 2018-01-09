# Auxiliar functions

vero=function(hmm,f){
  
  g=list()
  for (i in 1:nrow(f)){
    g[[i]]=unlist(strsplit(f[i,],split = ""))
  }
  
  v=1
  j=NULL
  for (i in 1:length(g)) {v=v*sum((FW1=exp(forward(hmm = hmm,(unlist(g[i])))))[,length((unlist(g[i])))]) }
  
  return(v)
}

priori=function(hmm){
  p_mod=((nrow(hmm$transProbs))^(-sum((hmm$transProbs[1,])>0)))*((nrow(hmm$emissionProbs))^(-sum((hmm$emissionProbs[1,])>0)))
  
  for (z in 2:nrow(hmm$transProbs)){
    p_mod=p_mod*((nrow(hmm$transProbs))^(-sum((hmm$transProbs[z,])>0)))*((ncol(hmm$emissionProbs))^(-sum((hmm$emissionProbs[z,])>0)))
  }
  return(p_mod)
}

reesViterbi=function(hmm,f){
  
  g=list()
  for (i in 1:nrow(f)){
    g[[i]]=unlist(strsplit(f[i,],split = ""))
  }
  
  j=NULL
  for (i in 1:length(g)) {j=c(j,unlist(g[i])) }
  
  j1=NULL
  for (i in 1:length(g)) {j1=c(j1,viterbi(hmm,unlist(g[i]))) }
  j2=cbind(j,j1)
  
  #vector con todos los estados generadores
  x=as.numeric(j1[-length(j1)])
  # vector de los estados sin el primer elemento y con una F al final
  y=as.numeric(j1[-1])
  #junto los dos vectores para tener las transiciones
  #para tener las transiciones , no tomo en cuenta cuando arranca de F
  e=which(x==999)
  #probabilidades iniciales
  sp=c(x[1],y[e])
  g1=as.data.frame(table(sp),stringsAsFactors = F)
  sp1=rep(0,length(unique(c(x,y))))
  
  
  
  x=x[-e];y=y[-e]
  
  lev=unique(sort(c(x,y)))
  y<-factor(y,levels =lev)
  x<-factor(x,levels =lev)
  
  
  t=prop.table(table(x,y),margin = 1)
  t[t=="NaN"]=0
  #matriz de emisiones
  e1=as.matrix(prop.table(table(j2[,2],j2[,1]),1))
  e1=e1[mixedorder(rownames(e1)),]
  #estados
  s=rownames(e1)
  #alfabeto
  al=colnames(e1)
  
  #probabilidades iniciales
  sp2=data.frame(cbind(s,sp1),stringsAsFactors = F)
  sp2[sp2$s %in% (g1$sp),]$sp1=as.numeric(g1$Freq)
  sp2$sp1=as.numeric(sp2$sp1)
  sp2$sp1=sp2$sp1*(1/length(sp))
  
  
  #primero se identifica el modelo inicial
  hmm=initHMM(States = s,Symbols = al,transProbs=t,emissionProbs=e1,startProbs =sp2$sp1)
  return(hmm)
  
}

modmeg=function(hmm,f) {
  #en u voy a ir guardando las posterioris
  u=NULL
  #en ff vo a ir guardando los indices que se van a juntar
  ff=NULL
  #pruebo mergear cada uno con cada uno
  for (i in 1:(length(hmm$States)-1)){
    j=i+1
    while (   j<length(hmm$States) )
    {
      TP=hmm$transProbs
      #elimino fila y columna j
      TP[,i]=TP[,i]+TP[,j]
      TP=TP[,-j]
      TP[i,]=TP[i,]+TP[j,]
      TP=TP[-j,]
      
      
      TP[i,]=TP[i,]/sum(!TP[i,]==0)
      #TP[i,]=TP[i,]/sum(TP[i,])
      
      
      
      
      ##reestimo emisiones
      EP=hmm$emissionProbs
      EP[i,]=(EP[i,]+EP[j,])/sum((EP[i,]+EP[j,]))
      EP=EP[-j,]
      
      #restimo probabilidades iniciales
      SP=hmm$startProbs
      SP[i]=(SP[i]+SP[j])
      SP=SP[-j]
      
      
      #corro el modelo con un nuevo merge
      hmm2 = initHMM(hmm$States[-j],hmm$Symbols,
                     transProbs=TP,
                     emissionProbs=EP,
                     startProbs =SP)
      
      #reestimo probabilidades
      hmm2=reesViterbi(hmm2,f)
      
      #print(c(i,j))
      ff=c(ff,paste(i,j,sep="-"))
      j=j+1
      print(tail(ff,1))
      
      v=vero(hmm2,f)
      
      #prob a priori
      #p_mod=priori(hmm2)
      #voy guardando las prob a posteriori
      #u=c(u,v*p_mod)
      u=c(u,v)
      
      
      #m=list(m,hmm2)
      
    }
    
  }
  #identifico el candidato a mergiar
  e=(which.max(u))
  
  #i y j son los candidatos a mergiar
  #cuantos digitos tienen los dos estados
  estados=as.vector(unlist(strsplit(ff[e],"-")[1]))
  
  i=as.numeric(estados[1]);j=as.numeric(estados[2])
  
  TP=hmm$transProbs
  TP[,i]=TP[,i]+TP[,j]
  TP=TP[,-j]
  TP[i,]=TP[i,]+TP[j,]
  TP[i,]=TP[i,]/sum(!TP[i,]==0)
  TP[i,]=TP[i,]/sum(TP[i,])
  TP=TP[-j,]
  
  
  
  ##reestimo emisiones
  EP=hmm$emissionProbs
  EP[i,]=(EP[i,]+EP[j,])/sum((EP[i,]+EP[j,]))
  EP=EP[-j,]
  
  #restimo probabilidades iniciales
  SP=hmm$startProbs
  SP[i]=(SP[i]+SP[j])
  SP=SP[-j]
  
  
  
  #corro el modelo con un nuevo merge
  hmm3 = initHMM(hmm$States[-j],hmm$Symbols,
                 transProbs=TP,
                 emissionProbs=EP,
                 startProbs =SP)
  #reestimo probabilidades
  hmm3=reesViterbi(hmm3,f)
  
  return(hmm3)
}