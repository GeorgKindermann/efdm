##################################################################
# EFDMoutput
# December 12th, 2014
# Authors: Maria Eronen, Seija Sirkiä, Tuula Packalen
# Original Author: Seija Sirkiä, January 25th 2013
#   
#   
#   Copyright 2014 European Union
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
outputcalc <- function(outputrequests){
  
requests<-readLines(outputrequests)
#Reads the file which includes output requests

nline<-length(requests)
#Number of lines in outputrequests file

for(i in seq(1,nline-1,by=2)){ #1 3 5 ...
#tmp<-readLines(outputrequests)

outputrequest<-requests[i]
#Reads the output request from a textfile
  
strsplit(outputrequest, " ") -> strlist
#Splits the string from the textfile to pieces
  
outputvar <-strlist[[1]][1]
byvar <- strlist[[1]][3:length(strlist[[1]])]
#The 1st and 3rd word from the string

coeffile<-read.table(requests[i+1],header=TRUE) 
#Reads the coefficient file

var.acts<-paste(outputvar,actnames,sep=".")
acts <- intersect(names(coeffile),var.acts)
#Number of activities

if(length(acts) != 0){
  #If there are some activities
  
    for(j in setdiff(var.acts,acts)){          
      coeffile[[j]] <- 0
    }
    #Activities which are not active are 0.
  
  output<-merge(rawoutput,coeffile)
    #Combines two data together
  
  tmp <- matrix(0,nrow=nrow(output), ncol=nrofsteps+1)
    #Matrix which has the same number of columns as there 
    #are steps (or the size of variable nrofsteps)
  
  for(k in 1:nact){ 
    table <-as.matrix(subset(output,select=paste0("step", 0:nrofsteps,".", actnames[k])))
    #Chooses the steps to another matrix
    coef <- output[[var.acts[k]]]
    #Chooses the coefficients of activities
    tmp <- sweep(table,1,coef,FUN="*") + tmp    
  }
    
  out <- cbind(output[factnames], tmp)
    #Combines two data together
  names(out) <- c(factnames,paste0("step",0:nrofsteps))
    #Gives names to some columns
    
}else{
  #If there are not any activities
  
  output<-merge(resultstates,coeffile)
  #Combines two data together
  
  table <-as.matrix(subset(output,select=paste0("step",0:nrofsteps)))
  #Chooses the steps to another matrix
  coef <- output[[outputvar]]  
  tmp <- sweep(table,1,coef,FUN="*")
  
  out <- cbind(output[factnames], tmp)
  #Combines two data together
  names(out) <- c(factnames,paste0("step",0:nrofsteps))
  #Gives names to some columns
}

output<-aggregate(x=subset(out,select=paste0("step",0:nrofsteps)),by=subset(out,select=byvar),FUN=sum) 
#Splits data to subsets by variable byvar and sums them 
#output <- colSums(output[,paste0("step", 0:nrofsteps)])
#Sums the columns of steps

write.table(output, file=paste0(outputvar, "_by_",paste(byvar,collapse="_"),".txt"))
#Saves the results to a text file
efdmplot(output,outputvar,byvar)

} #For-statement ends
} #Function ends

efdmplot <- function(dname,outputvar,byvar){
  #Draws barplots of data frames given as argument or in a file
  if(is.character(dname)) 
    #dname is a name of a file
  {
    if(missing(outputvar)|missing(byvar))
      #the names of the variables have to be deduced from the 
      #file name. THIS ASSUMES WINDOWS STYLE FILE NAMES!!      
    {
      out <- strsplit(dname, "[\\]") #Splits the path to elements
      out <- tail(out[[1]], n=1) #Chooses the last one
      out <- strsplit(out, "_") 
      #Splits the last to elements (again)
      
      outputvar <- out[[1]][1] #Chooses the first one
      byvar <- out[[1]][3:length(out[[1]])] 
      #Chooses elements from 3rd to last
      
      byvarlast <- strsplit(tail(byvar,n=1), "[.]") 
      #Splits the last of byvar
      byvar <- c(byvar[-length(byvar)],byvarlast[[1]][1]) 
      #Combines byvar together
      byvar <- as.character(byvar)
    }
    output <- read.table(dname) #Reads the data for barplot    
  } #if statement (dname is a file name) ends
  
  else #dname is presumably a data frame
    output<-dname

  
  #M<-as.matrix(output[-(1:length(byvar))])
  M<-as.matrix(output[setdiff(names(output),byvar)])
  #Data to matrix without the names of byvars
  
  #rownames(M) <- as.character(paste(output[[1]],output[[length(byvar)]]))
  rownames(M) <- do.call(paste,output[byvar])
  #Rownames to data
  
  jpeg(paste0("barplot_of_totals_ of_",outputvar,"_by_",paste(byvar,collapse="_"),".jpg"))
  barplot(M, legend.text=TRUE, col=c(1:nrow(M)), ylab=outputvar)
  dev.off()
  #Barplot with totals
  
  jpeg(paste0("barplot_of_",outputvar,"_by_",paste(byvar,collapse="_"),".jpg"))
  barplot(M,beside=TRUE, legend.text=TRUE, col=c(1:nrow(M)), ylab=outputvar)
  dev.off()
  #Barplot without totals
  
} #Function ends

