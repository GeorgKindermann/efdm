dir.create("res")
setwd("./res/")

#Start efdm
tmp<-readLines("../dat/efdminput.txt")
tmp<-matrix(unlist(strsplit(tmp,split=" ")),ncol=2,byrow=TRUE)
VARS<-tmp[,1]
apply(tmp,1,FUN=function(x) assign(x[1],x[2],pos=.GlobalEnv))
source("../../../EFDMcode/v2.0/efdmcore.r")
if(OUTPUTREQUEST_FILENAME!=0) {
  source("../../../EFDMcode/v2.0/efdmoutput.r")
  outputcalc(OUTPUTREQUEST_FILENAME)
}

q()


#INSTEAD OF:
#source("../../../EFDMcode/v2.0/efdmcore.r")

#describing the factors:
factlvls <- strsplit(readLines(STATESPACE_FILENAME), " ")
factlvls <- lapply(setNames(factlvls, lapply(factlvls, "[", 1)), "[", -1)
factdims <- vapply(factlvls, length, 1)

#describing the activities
tmp <- strsplit(readLines(ACTIVITIES_FILENAME), " ")
actnames <- vapply(tmp, "[", character(1), 1)
actmethod <- vapply(tmp, "[", character(1), 2)
actfiles <- vapply(tmp, "[", character(1), 3)
activities <- lapply(setNames(tmp, actnames), "[", -1:-3)
rm(tmp)

#setting up the transition matrices
transmats <- lapply(setNames(actfiles, actnames), function(i) get(load(i)))

#setting up the initial state
tmp <- read.table(INITSTATE_FILENAME,header=TRUE)
initstate <- array(0, factdims, factlvls)
initstate[as.matrix(tmp[names(dimnames(initstate))])]  <- tmp$area

#setting up the probabilities of activities
tmp <- read.table(ACTPROBS_FILENAME,header=TRUE)
actproblist <- lapply(setNames(actnames, actnames), function(i) {
  "[<-"(array(0, factdims, factlvls), as.matrix(tmp[names(dimnames(initstate))]), tmp[,i])
})

#Divide by activities
state <- lapply(names(actproblist), function(i) {
  initstate * actproblist[[i]]
})

resultstates <- list(step0=state)

length(activities)
names(factlvls)
tmp <- aperm(state[[1]], 1:3)
dim(tmp) <- c(prod(dim(tmp)[1:2]),1,prod(dim(tmp)[-(1:2)]))
transmati<-transmats[[1]]
apply(abind(transmati,tmp,along=2),3,mmf)

nstate<-newstate(state)  
state<-rowSums(do.call(cbind,lapply(nstate,as.vector)))
dim(state)<-factdims
state<-dividebyA(state)
resultstates[[paste0("step",i)]]<-state



mmf<-function(M)
  #simply for doing the necessary "matrix times a vector" multiplications
  #somewhat efficiently
{
  m<-dim(M)[1]
  M[,1:m]%*%M[,m+1]
}

newstate<-function(oldstatediv) 
  #for producing the next state from the current one, using the obtained
  #activity and transition probabilities
{
  result<-oldstatediv
  for(i in 1:nact) {
    nr<-1:nfact
    names(nr)<-factnames
    acti<-activities[[i]]
    factnames.here<-c(acti,setdiff(factnames,acti))
    tmp<-aperm(oldstatediv[[actnames[i]]],nr[factnames.here])
    #tmp is the subpopulation receiving the activity i, with appropriate dims
    dim(tmp)<-c(prod(dim(tmp)[1:length(acti)]),1,prod(dim(tmp)[-(1:length(acti))]))
    transmati<-transmats[[i]]
    #transmati is the transition matrix corresponding to activity i
    dim(transmati)<-c(dim(transmati)[1:2],prod(dim(transmati)[-(1:2)]))
    tmp<-apply(abind(transmati,tmp,along=2),3,mmf)
    #tmp now has the new states of this subpopulation
    dim(tmp)<-factdims[factnames.here]
    names(nr)<-factnames.here
    tmp<-aperm(tmp,nr[factnames])
    dimnames(tmp)<-dimnames(tmp)
    result[[actnames[i]]]<-tmp
  }

  result
}
