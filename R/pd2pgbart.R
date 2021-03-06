pd2pgbart = function (
   x.train, y.train,x.test=matrix(0.0,0,0),
   model = "./pgbart.model",
   if_trained = FALSE,train_result=NULL,
   xind=1:2,levs=NULL, levquants=c(.05,(1:9)/10,.95),
   pl=TRUE, plquants=c(.05,.95), 
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
{
   n = nrow(x.train)
   nlevels = rep(0,2)
   if(is.null(levs)) {
      levs = list()
      for(i in 1:2) {
         ux = unique(x.train[,xind[i]])
	 if(length(ux) <= length(levquants)) levs[[i]] = sort(ux)
	 else levs[[i]] = unique(quantile(x.train[,xind[i]],probs=levquants))
      }
   } 
   nlevels = unlist(lapply(levs,length))
   xvals <- as.matrix(expand.grid(levs[[1]],levs[[2]]))
   nxvals <- nrow(xvals)
   if (ncol(x.train)==2){
      cat('special case: only 2 xs\n')
      x.test = xvals
   } else {
      x.test=NULL
      for(v in 1:nxvals) {
         temp = x.train
         temp[,xind[1]] = xvals[v,1]
         temp[,xind[2]] = xvals[v,2]
         x.test = rbind(x.test,temp)
      }
   }
   if (if_trained==TRUE){
     pdbrt=train_result
   }else{
     pdbrt = pgbart_train(x.train,y.train,x.test,
                          model= "./pgbart.model",
                          sigest, sigdf, sigquant,k,power, base,
                          binaryOffset, ntree,ndpost, nskip, usepg, 
                          numparticles ,printevery, keepevery, 
                          keeptrainfits,usequants, numcut, printcutoffs,
                          verbose)
   }
   if (ncol(x.train)==2) {
      fdr = pdbrt$yhat.test
   } else {
      fdr = NULL 
      for(i in 1:nxvals) {
         cind =  ((i-1)*n+1):(i*n)
         fdr = cbind(fdr,(apply(pdbrt$yhat.test[,cind],1,mean)))
      }
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
   class(retval) = 'pd2pgbart'
   if(pl) plot.pgbart(retval,plquants=plquants)
   return(retval)
}

