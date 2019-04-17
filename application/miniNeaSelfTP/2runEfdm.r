dir.create("res")
setwd("./res/")

#Start efdm
tmp<-readLines("../dat/efdminput.txt")
tmp<-matrix(unlist(strsplit(tmp,split=" ")),ncol=2,byrow=TRUE)
VARS<-tmp[,1]
apply(tmp,1,FUN=function(x) assign(x[1],x[2],pos=.GlobalEnv))
source("../../../EFDMcode/nea/efdmcore.r")
if(OUTPUTREQUEST_FILENAME!=0) {
  source("../../../EFDMcode/v2.0/efdmoutput.r")
  outputcalc(OUTPUTREQUEST_FILENAME)
}
