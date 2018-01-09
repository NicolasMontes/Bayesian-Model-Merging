
m0=function(f){
  #primero hago una lista en donde pongo cada nombre como un vector y cada letra un caracter
  g=list()
  for (i in 1:nrow(f)){
    g[[i]]=unlist(strsplit(f[i,],split = ""))
  }
  
  #Pmodelo M0 a partir de los datos
  #creo un vector con todos los simbolos observados concatenados
  j=NULL
  for (i in 1:length(g)) {j=c(j,unlist(g[i])) }
  #creo un vector con todos los simbolos concatenados y al lado el estado
  #que lo genero
  j2=as.data.frame(cbind(j,cumsum(!j==".")),stringsAsFactors = F)
  j2[j2$j==".",2]="999"
  
  #vector con todos los estados generadores
  x=as.numeric(j2[-nrow(j2),2])
  # vector de los estados sin el primer elemento y con una F al final
  y=as.numeric(j2[-1,2])
  #junto los dos vectores para tener las transiciones
  #para tener las transiciones , no tomo en cuenta cuando arranca de F
  e=which(x==999)
  #probabilidades iniciales
  sp=c(1,y[e])
  sp1=rep(0,length(unique(c(x,y))))
  sp1[as.numeric(sp)]=1
  sp1=sp1*(1/length(sp))
  x=x[-e];y=y[-e]
  
  
  lev=unique(sort(c(x,y)))
  y<-factor(y,levels =lev)
  x<-factor(x,levels =lev)
  
  
  t=prop.table(table(x,y),margin = 1)
  t[t=="NaN"]=0
  #matriz de emisiones
  e1=as.matrix(prop.table(table(j2[,2],j2[,1]),1))
  e1=e1[mixedorder(rownames(e1)), ]
  #estados
  s=rownames(e1)
  #alfabeto
  al=colnames(e1)
  #primero se identifica el modelo inicial
  hmm = initHMM(States = s,Symbols = al,transProbs=t,emissionProbs=e1,startProbs =sp1)
  return(hmm)
  
}