# Bayesian Model Merging
source(file = "../functions/bmm_aux.R")

bmm=function(hmm,f,lam){
  states=(length(hmm$States)-1)
  i=1
  s=0
  v0=vero(hmm,f)
  p0=priori(hmm)
  
  sink("hist.txt")
  print(paste("step=;",i,"; length=;",length(hmm$States), "; Log lik=;", log(v0),"; Log prior=;", log(p0), "; time=;",Sys.time()))
  sink()
  
  
  
  while  (i<states & s==0){
    #busca el mejor candidato
    hmm1=modmeg(hmm,f)
    v1=vero(hmm1,f)
    p1=priori(hmm1)
    
    {
      if ( (log(v1)+lam*log(p1))>=(log(v0)+lam*log(p0))) {
        i=i+1
        hmm=hmm1
        v0=v1
        p0=p1
      }
      else {s=1}
    }
    
    sink(file = "hist.txt", append = TRUE)
    print(paste("step=;",i,"; length=;",length(hmm$States), "; Log lik=;", log(v1),"; Log prior=;", log(p1), "; time=;",Sys.time()))
    sink()
    
    b=paste("models/hmm_",i,".RData",sep = "")
    save(hmm,file =b)
  }
  return(hmm)
}