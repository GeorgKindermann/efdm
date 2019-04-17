###################################################################
# EFDM core
# January 25th 2013
# Author: Seija Sirkiä
#   
#   
#   Copyright 2013 European Union
#   
#   Licensed under the EUPL, Version 1.1 or – as soon they
#   will be approved by the European Commission - subsequent
#   versions of the EUPL (the "Licence");
#   You may not use this work except in compliance with the
#   Licence.
#   You may obtain a copy of the Licence at:
#   
#   http://joinup.ec.europa.eu/software/page/eupl/licence-eupl
#   
#   Unless required by applicable law or agreed to in
#   writing, software distributed under the Licence is
#   distributed on an "AS IS" basis,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
#   express or implied.
#   See the Licence for the specific language governing
#   permissions and limitations under the Licence.
#   
###################################################################

#Disclaimer: this collection of code should be seen merely as an implementation
#of an idea, a prototype, not an actual finished program. It relies heavily
#on use of global variables, in other words, functions use variables outside 
#their own environment. Assumptions are generally not checked and if they do
#not hold, there is no telling what will happen.

#General comment: the most delicate thing here is keeping all the multiway
#arrays in a knowable format, while allowing for differing amounts of
#factors and their levels etc. This is easily broken if the code is edited
#carelessly.

#Input files: 
# factors.txt (or any other name, remember to edit the name in the code if
# some other name is used)
# -one line per a factor: first the name of the factor, then the names of 
#  the levels, separated by whitespaces
# -order is irrelevant, but preserved as a common default everywhere
#
# activities.txt (or any other name, see above)
# -one line per an activity:
#  *first, the name of the activity
#  *second, the method of obtaining the transition probabilities involved;
#   read: read them from a file; estimate: estimate based on data from 
#   a file; custom: using a custom function provided by the user
#  *third, the name of the file for probabilities/data, or the name of the
#   custom function (see below for format)
#  *fourth, names of the factors that are prone to change under this activity
#  everything separated by whitespaces
#
# files for ready to use transition probabilities:
# -the name of the file must be given in activities.txt as the third entry
#  for the activity
# -option one: one matrix, given as a simple square table of numbers with 
#  an equal number of rows an columns, separated by whitespaces
# -option two: a sequence of such matrices, for different levels or level
#  combinations of one or more factors. In this case there must be lines 
#  in the file declaring for which factor level combination the transition 
#  matrix following it is, as in 
#  factor1=level1 factor2=level1 on the first line and
#  factor1=level2 factor2=level1 on the line after the first matrix
#  The levels for the "leftmost" factor name "run the fastest", similar 
#  to array indices in R, as shown above
#
# files for data for the estimation
# -the name of the file must be given in activities.txt as the third entry
# -"ordinary" whitespace separated table format: names of the columns on the 
#  first line, one line per plot or similar, observed at two different time
#  points
# -there must be one column per each factor that does not change, named by
#  its name
# -there must be two columns per each factor that does change, named 
#  factorname0 and factorname1, for the "before" observation and "after"
#  observation, respectively, and where factorname is the name of the factor
#
# priors for Bayes-like estimation: 
# -the name of the file must be the same as for the data file, given in 
#  activities.txt as the third entry, superceded by "prior_"
# -one square table of prior frequencies, similar to the file of a single
#  ready to use transition probability matrix
#
# statespace.txt (or any other name, again)
# -"ordinary" whitespace separated table format: names of the columns on the 
#  first line, one line per state
# -there must be one column per each factor, named by their names
# -there must be one column for the amount/area of units in the state, named
#  freely
# -there must be one column per activity, named by the name of the activity
#  giving the probability that a unit in the state receives this activtiy
#
# number of simulation rounds is given as a variable near the end of the file,
# just before the simulation loop







#--------------------------------
# functions needed:

library(abind) #for function abind

readtransprobs<-function(filename)
  #for reading a ready made transition matrix, or several
{
  rawdata<-readLines(filename)
  if(length(rawdata)==basedim) # if so, 
    #there's just one transition matrix used everywhere
  {
    tmp<-matrix(as.numeric(unlist(strsplit(rawdata,split=" "))),
                ncol=basedim,byrow=TRUE)
    transmat<-rep(tmp,times=prod(otherdim))
    dim(transmat)<-c(basedim,basedim,otherdim)
    return(transmat)
  }
  #otherwise, there are several matrices, for different factor levels
  
  #Which factors are involved? This should be mentioned on the first line
  tmp<-unlist(strsplit(rawdata[1],split=" "))
  tmp<-unlist(strsplit(tmp,split="="))
  fs<-tmp[seq(1,length(tmp),by=2)] #first, third etc are the factor names
  fn<-factdims[fs]
  #rest of the file is the matrices, separated by "header" lines
  frows<-seq(from=1,length=prod(fn),by=basedim+1) #these are the header lines
  tmp<-matrix(as.numeric(unlist(strsplit(rawdata[-frows],split=" "))),
              ncol=basedim,byrow=TRUE)
  #Actually, the header lines are thrown away! It doesn't matter what's on them,
  #but they probably help keeping the file organized as assumed
  dim(tmp)<-c(basedim,fn,basedim)
  tmp<-aperm(tmp,c(1,length(fn)+2,1:length(fn)+1))
  os<-setdiff(otherdimnames,fs)
  transmat<-rep(tmp,times=prod(factdims[os]))
  dim(transmat)<-c(basedim,basedim,fn,factdims[os])
  nr<-1:length(otherdimnames)
  names(nr)<-c(fs,os)
  transmat<-aperm(transmat,c(1,2,nr[otherdimnames]+2))
  return(transmat)
}

estimatetransprobs<-function(filename)
  #for estimating a transition matrix from data using the Bayes-like 
  #algorithm explained elsewhere
{
  pairdata<-read.table(filename,header=TRUE)
  
  #the following ensures that factors are represented in the objects correctly
  #and consistently
  for(j in 1:length(otherdimnames))
    pairdata[[otherdimnames[j]]]<-factor(pairdata[[otherdimnames[j]]],
                                         levels=factlvls[[otherdimnames[j]]],ordered=TRUE)
  
  name0<-paste0(basedimnames,0)
  name1<-paste0(basedimnames,1)
  for(j in 1:length(basedimnames))
  {
    pairdata[[name0[j]]]<-factor(pairdata[[name0[j]]],
                                 levels=factlvls[[basedimnames[j]]],ordered=TRUE)
    pairdata[[name1[j]]]<-factor(pairdata[[name1[j]]],
                                 levels=factlvls[[basedimnames[j]]],ordered=TRUE)
  }
  #prior is read from another file:
  tmp<-readLines(paste0("prior_",filename))
  #the order of factors, "from most influential to least influential", is
  #given on the first line, together with corresponding weights
  factorder<-unlist(strsplit(unlist(strsplit(tmp[1],split=" ")),split="="))
  dim(factorder)<-c(2,length(factorder)/2)
  pairdata<-pairdata[,c(name1,name0,factorder[1,])]
  
  #start of estimations algorithm
  freq<-table(pairdata)
  dim(freq)<-c(basedim,basedim,factdims[factorder[1,]])
  s<-matrix(as.numeric(unlist(strsplit(tmp[-1],split=" "))),ncol=basedim,byrow=TRUE)
  priorweights<-c(1,as.numeric(factorder[2,]))
  #To try and explain what is going on in the following for loop:
  
    #1: On the 1st round, tmp is 2-dimensional, the frequency table of
  #the changing factor levels, pooled over every non-changing factor.
  #The 2-dimensional table of prior "frequencies" in s is summed to 
  #to tmp, and saved again to s.
  
  #2: On the 2nd round, tmp is 3-dimensional, standing for the 2-dimensional
  #frequency tables of the changing factor levels, one for each level of the
  #first non-changing factor.
  #The 2-dimensional table of prior "frequencies" in s obtained on the first
  #round is summed to each 2-dimensional factor1 specific subtable 
  #of tmp, and the 3-dimensional result is saved again to s.
  
  #3: On the 3rd round, tmp is 4-dimensional, standing for the 2-dimensional
  #frequency tables of the changing factor levels, one for each combination 
  #of the levels of the first 2 non-changing factors.
  #The 2-dimensional tables of prior "frequencies" in s obtained on the second
  #round for each level of the first non-changing factor, are summed to each 
  #2-dimensional factor1 and factor2 specific subtable of tmp, and the 
  #4-dimensional result is saved again to s.
  
  #And so on.
  for (i in 1:length(priorweights))
  {
    tmp<-apply(freq,1:(1+i),sum)
    s<-sweep(tmp,1:max(2,i),s*priorweights[i],'+')
  }
  #in the end, s contains all the numerator values for the singular 
  #probability estimates, and the denominators are the column sums of
  #all those 2-dimensional tables, obtained here:
  N<-apply(s,2:length(dim(s)),sum)
  #then, numbers in s are divided by the corresponding number in N:
  transmat<-sweep(s,2:length(dim(s)),N,'/')
  #end of estimation algorithm
  
  nr<-1:length(otherdimnames)
  names(nr)<-factorder[1,]
  transmat<-aperm(transmat,c(1,2,nr[otherdimnames]+2))
  return(transmat)
}

mmf<-function(M)
  #simply for doing the necessary "matrix times a vector" multiplications
  #somewhat efficiently
{
  m<-dim(M)[1]
  M[,1:m]%*%M[,m+1]
}

newstate<-function(oldstate) 
  #for producing the next state from the current one, using the obtained
  #activity and transition probabilities
{
  result<-0*oldstate
  for(i in 1:nact) {
    nr<-1:nfact
    names(nr)<-factnames
    acti<-activities[[i]]
    factnames.here<-c(acti,setdiff(factnames,acti))
    tmp<-aperm(oldstate*actproblist[[i]],nr[factnames.here])
    #tmp is the subpopulation receiving the activity i
    dim(tmp)<-c(prod(dim(tmp)[1:length(acti)]),1,prod(dim(tmp)[-(1:length(acti))]))
    transmati<-transmats[[i]]
    #transmati is the transition matrix corresponding to activity i
    dim(transmati)<-c(dim(transmati)[1:2],prod(dim(transmati)[-(1:2)]))
    tmp<-apply(abind(transmati,tmp,along=2),3,mmf)
    #tmp now has the new states of this subpopulation
    dim(tmp)<-factdims[factnames.here]
    names(nr)<-factnames.here
    tmp<-aperm(tmp,nr[factnames])
    result<-result+tmp
  }
  dimnames(result)<-dimnames(oldstate)
  result
}

#-----------
#initializations:

#describing the factors:
factors.tmp<-readLines("factors.txt")
nfact<-length(factors.tmp)
factors.tmp<-strsplit(factors.tmp," ")
factnames<-character(nfact)
for(i in 1:nfact) {
  factnames[i]<-factors.tmp[[i]][1]
  factors.tmp[[i]]<-factors.tmp[[i]][-1]
}
factlvls<-factors.tmp
names(factlvls)<-factnames
factdims<-sapply(factlvls,length)
rm(factors.tmp)

#describing the activities
activities.tmp<-readLines("activities.txt")
nact<-length(activities.tmp)
activities<-strsplit(activities.tmp," ")
actnames<-character(nact)
actmethod<-character(nact)
actfiles<-character(nact)
rm(activities.tmp)

for(i in 1:nact) {
  actnames[i]<-activities[[i]][1]
  actmethod[i]<-activities[[i]][2]
  actfiles[i]<-activities[[i]][3]
  activities[[i]]<-activities[[i]][-(1:3)]
}
names(activities)<-actnames

#setting up the transition matrices
transmats<-list()
for(i in 1:length(activities))
{
 #the levels of these factors may change under this activity:
 basedimnames<-activities[[i]] 
 #these are fixed under this activity:
 otherdimnames<-setdiff(factnames,basedimnames)
 basedim<-prod(factdims[basedimnames])
 otherdim<-factdims[otherdimnames]

 if(actmethod[i]=="read") 
   transmat.forthisround <- readtransprobs(actfiles[i]) 
 if(actmethod[i]=="estimate") 
   transmat.forthisround <- estimatetransprobs(actfiles[i])
 if(actmethod[i]=="custom") 
   transmat.forthisround<- eval(call(actfiles[i]))
   #the functions that are called here have to be written and loaded by
   #the user!

 #the rest is making sure dimensions have correct and sensible names
 lvlnames<-factlvls[[basedimnames[1]]]
 if(length(basedimnames>1))
 for(j in 2:length(basedimnames))
 {
   lvlnames.more<-factlvls[[basedimnames[j]]]
   lvlnames<-paste(rep(lvlnames,times=length(lvlnames.more)),
         rep(lvlnames.more,each=length(lvlnames)),sep=".")
 }
 dimnames(transmat.forthisround)<-c(list(lvlnames,lvlnames),factlvls[otherdimnames])
 transmats<-c(transmats,list(transmat.forthisround))
 names(transmats)[i]<-actnames[i]
}

rm(basedim,basedimnames,otherdim,otherdimnames,transmat.forthisround,lvlnames,
   lvlnames.more)

#setting up the initial state and the probabilities of activities
statespace<-read.table("statespace.txt",header=TRUE)
for(fname in factnames)
  statespace[[fname]]<-factor(statespace[[fname]],
                      levels=factlvls[[fname]],ordered=TRUE)
statespace<-statespace[do.call(order,statespace[rev(factnames)]),]

initstate<-statespace[[setdiff(names(statespace),union(factnames,actnames))]]
dim(initstate)<-factdims
dimnames(initstate)<-factlvls

actproblist<-list()
for(i in 1:nact) {
  actproblist[[i]]<-array(statespace[[actnames[i]]],dim=factdims)
  dimnames(actproblist[[i]])<-factlvls
}

#-----------------------------
#finally, the actual simulation

resultstates<-statespace[factnames]
resultstates[["step0"]]<-as.vector(initstate)

nrofsteps<-10 #how many steps should be taken?

state<-initstate
for (i in 1:nrofsteps) {
  state<-newstate(state)
  resultstates[[paste0("step",i)]]<-as.vector(state)
}


#if all went well, the object resultstates now has the state of the population
#at the initial step and further timesteps asked

#it can be written in to a file and read in to e.g. a spreadsheet in another 
#program
write.table(resultstates,"resultstates.txt")

#or used within R