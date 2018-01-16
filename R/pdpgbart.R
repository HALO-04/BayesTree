pdpgbart = function (
   x.train, y.train, x.test=matrix(0.0,0,0),
   model = "./pgbart.model",
   if_trained = FALSE,train_result=NULL,
   xind=1:2, levs=NULL, levquants=c(.05,(1:9)/10,.95),
   pl=TRUE,  plquants=c(.05,.95),
   sigest=NA, sigdf=3, sigquant=.90,
   k=2.0,
   power=NA, base=.95,
   binaryOffset=0,
   ntree=200,
   ndpost=1000, nskip=100, usepg=TRUE, numparticles = 10,
   printevery=100, keepevery=1, keeptrainfits=TRUE,
   usequants=FALSE, numcut=100, printcutoffs=0,
   verbose=TRUE,
   ...
)
{  xind=1:ncol(x.train)
   n = nrow(x.train)
   nvar = length(xind)
   nlevels = rep(0,nvar)
   if(is.null(levs)) {
      levs = list()
      for(i in 1:nvar) {
         ux = unique(x.train[,xind[i]])
	 if(length(ux) < length(levquants)) levs[[i]] = sort(ux)
	 else levs[[i]] = unique(quantile(x.train[,xind[i]],probs=levquants))
      }
   }
   nlevels = unlist(lapply(levs,length))
   x.test=NULL
   for(i in 1:nvar) {
      for(v in levs[[i]]) {
         temp = x.train
         temp[,xind[i]] = v
         x.test = rbind(x.test,temp)
      }
   }
   if (if_trained==TRUE){
     pdbrt=train_result
   }else{
     pdbrt = pgbart_train(x.train,y.train,x.test,model= "./pgbart.model",
                          sigest, sigdf, sigquant,k,power, base,
                          binaryOffset, ntree,ndpost, nskip, usepg, 
                          numparticles ,printevery, keepevery, 
                          keeptrainfits,usequants, numcut, printcutoffs)
   }
   fdr = list()
   cnt=0
   for(j in 1:nvar) {
      fdrtemp=NULL
      for(i in 1:nlevels[j]) {
         cind = cnt + ((i-1)*n+1):(i*n)
         fdrtemp = cbind(fdrtemp,(apply(pdbrt$yhat.test[,cind],1,mean)))
      }
      fdr[[j]] = fdrtemp
      cnt = cnt + n*nlevels[j]
   }
   if(is.null(colnames(x.train))) xlbs = paste('x',xind,sep='')
   else xlbs = colnames(x.train)[xind]
   if('sigma' %in% names(pdbrt)) {
   retval = list(fd = fdr,levs = levs,xlbs=xlbs,
      bartcall=pdbrt$call,yhat.train=pdbrt$yhat.train,
      first.sigma=pdbrt$first.sigma,sigma=pdbrt$sigma,
      yhat.train.mean=pdbrt$yhat.train.mean,sigest=pdbrt$sigest,y=pdbrt$y)
   } else {
   retval = list(fd = fdr,levs = levs,xlbs=xlbs,
      bartcall=pdbrt$call,yhat.train=pdbrt$yhat.train,
      y=pdbrt$y)
   }
   class(retval) = 'pdpgbart'
   if(pl) plot.pgbart(retval,plquants=plquants)
   return(retval)
}
