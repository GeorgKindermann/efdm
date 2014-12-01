setupstatespace<-function(filename)
  #This is a function for initializing the state space.
  #(Does the same as core has to do when it starts, 
  #should ultimately be written just once somewhere!
{
  factors.tmp<-readLines(filename)
  nfact<<-length(factors.tmp)
  factors.tmp<-strsplit(factors.tmp," ")
  factnames<<-character(nfact)
  for(i in 1:nfact) {
    factnames[i]<<-factors.tmp[[i]][1]
    factors.tmp[[i]]<-factors.tmp[[i]][-1]
  }
  names(factors.tmp)<-factnames
  factlvls<<-factors.tmp
  factdims<<-sapply(factlvls,length)
}


pre.estimate<-function(inputfilename,factorsfilename,changing,resultfilename)
  #This is a function for performing an estimation of a transition matrix
  #"beforehand", and not during the running of the whole model/simulation
  #Arguments:
  #inputfilename: 
  #-the name of the estimation control file (that would be given 
  #in the activities control file)
  #factorsfilename:
  #-the name of the state space defining control file (factors.txt by default)
  #changing: 
  #-names of the factors that "have" the transitions (that would be given in
  #the activities control file)
  #resultfilename:
  #-name of the file where the result is to be saved (do use the extension .RData
  #if you plan to read this in during model/simulation run!)


  #the call of this function would then look something like this:
  #pre.estimate("estiminput.txt","factors.txt","vol age","nomgmtP.RData")
  #(with the needed data and possible prior available of course)
{
  setupstatespace(factorsfilename)
  basedimnames<<-strsplit("vol age",split=" ")[[1]]
  otherdimnames<<-setdiff(factnames,basedimnames)
  basedim<<-prod(factdims[basedimnames])
  otherdim<<-factdims[otherdimnames]
  if(!exists("estimatetransprobs")) source("efdmestim.r")
  #create a name for the object based on the time and date
  #(I want to be fairly certain it is unique...)
  oname<-paste0(strsplit(date()," ")[[1]],collapse="")
  oname<-paste0(strsplit(oname,":")[[1]],collapse="")
  oname<-paste0("estimate",oname,collapse="")
  #note, this should not be relevant to the user, the file name
  #in which this object will be saved is a different thing
  assign(x=oname,value=estimatetransprobs(inputfilename))
  save(list=oname,file=resultfilename)
}

makeathinP<-function(inputfilename,voldrops,resultfilename)
  #Special purpose function for making a transition matrix out of
  #a pre-made "no management" transition matrix, in case of a
  #"vol-age" model. The idea is that thinning means "drop so and
  #so many volume classes, then grow as usual", leading to moving 
  #columns right.  
  
  #In the following nvol refers to the number of volume classes
  #and nage to the number of age classes. These can be deduced
  #from the inputs: the length of voldrops = nvol, nage can 
  #then be calculated from the dimensions of the input matrix.
  
  #Assumption 1: the pre-made matrix in the input file is a saved
  #R object of type array of dimensions
  #(nvol*nage) x (nvol*nage) x something, 
  #where something is either nothing (single forest type), or 
  #a single number (that many forest types defined by a single factor),
  #or several numbers (that many forest type defining factors). This
  #will be the case if it was created by "estimatetransprobs" function.

  #Assumption 2: The "first" 2-dimensional array is made of blocks
  #of size nvol x nvol arranged in nage x nage array. That is,
  #the first nvol columns refer to the first age class (and not the
  #other way around, first nage columns referring to the first vol
  #class.) This will be the case if the matrix in the input file
  #was created by the aforementioned function and specifying the 
  #changing factors to be "vol age", in that order (the order is
  #what matters, not the names used)

  #Arguments:
  #inputfilename
  #-name of the save file where the R object for "no management" matrix is
  #voldrops
  #-a vector of integers, how many levels to drop at each level (length 
  #equal to number of volume levels)
  #-resultfilename
  #-name of the file where the result is to be saved (do use the extension .RData
  #if you plan to read this in during model/simulation run!)

{
  oname<-load(inputfilename)
  tmp<-get(oname)
  nvol<-length(voldrops)
  dims<-dim(tmp)
  nage<-dims[1]/nvol
  dim(tmp)<-c(dims[1:2],prod(dims[-(1:2)]))
  ind<-1:nvol-voldrops
  ind<-rep((0:(nage-1))*nvol,each=nvol)+ind
  tmp<-tmp[,ind,]
  dim(tmp)<-dims
  #create a name for the object based on the time and date
  #(I want to be fairly certain it is unique...)
  oname<-paste0(strsplit(date()," ")[[1]],collapse="")
  oname<-paste0(strsplit(oname,":")[[1]],collapse="")
  oname<-paste0("thin",oname,collapse="")
  #note, this should not be relevant to the user, the file name
  #in which this object will be saved is a different thing
  assign(x=oname,value=tmp)
  save(list=oname,file=resultfilename)
}

neutralprior<-function(nvol,nage,resultfilename)
  #This is a function for creating a simple prior for the
  #estimation procedure, for a "vol-age" model, implementing
  #the idea "grow one age class older, stay in the same vol
  #class, unless already in the highest age class, in which
  #case stay where you are"
  #Arguments:
  #nvol and nage
  #-number of vol and age classes
  #resultfilename
  #-name of the (text) file the output will be written in
{
  tmp<-diag(nage)
  tmp<-tmp[,c(2:nage,nage)]
  diag(nvol)%x%tmp
  write.table(tmp%x%diag(nvol),file=resultfilename,row.names=FALSE,col.names=FALSE)  
}

