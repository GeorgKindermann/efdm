Minimal documentation for the European Forest Dynamic Model

Copyright 2013 European Union


Disclaimer: this collection of code should be seen merely as an implementation
of an idea, a prototype, not an actual finished program. It relies heavily
on use of global variables, in other words, functions use variables outside 
their own environment. Assumptions are generally not checked and if they do
not hold, there is no telling what will happen.

General comment: the most delicate thing here is keeping all the multiway
arrays in a knowable format, while allowing for differing amounts of
factors and their levels etc. This is easily broken if the code is edited
carelessly.

Input files: 
 factors.txt (or any other name, remember to edit the name in the code if
 some other name is used)
 -one line per a factor: first the name of the factor, then the names of 
  the levels, separated by whitespaces
 -order is irrelevant, but preserved as a common default everywhere

 activities.txt (or any other name, see above)
 -one line per an activity:
  *first, the name of the activity
  *second, the method of obtaining the transition probabilities involved;
   read: read them from a file; estimate: estimate based on data from 
   a file; custom: using a custom function provided by the user
  *third, the name of the file for probabilities/data, or the name of the
   custom function (see below for format)
  *fourth, names of the factors that are prone to change under this activity
  everything separated by whitespaces

 files for ready to use transition probabilities:
 -the name of the file must be given in activities.txt as the third entry
  for the activity
 -option one: one matrix, given as a simple square table of numbers with 
  an equal number of rows an columns, separated by whitespaces
 -option two: a sequence of such matrices, for different levels or level
  combinations of one or more factors. In this case there must be lines 
  in the file declaring for which factor level combination the transition 
  matrix following it is, as in 
  factor1=level1 factor2=level1 on the first line and
  factor1=level2 factor2=level1 on the line after the first matrix
  The levels for the "leftmost" factor name "run the fastest", similar 
  to array indices in R, as shown above

 files for data for the estimation
 -the name of the file must be given in activities.txt as the third entry
 -"ordinary" whitespace separated table format: names of the columns on the 
  first line, one line per plot or similar, observed at two different time
  points
 -there must be one column per each factor that does not change, named by
  its name
 -there must be two columns per each factor that does change, named 
  factorname0 and factorname1, for the "before" observation and "after"
  observation, respectively, and where factorname is the name of the factor

 priors for Bayes-like estimation: 
 -the name of the file must be the same as for the data file, given in 
  activities.txt as the third entry, superceded by "prior_"
 -one square table of prior frequencies, similar to the file of a single
  ready to use transition probability matrix

 statespace.txt (or any other name, again)
 -"ordinary" whitespace separated table format: names of the columns on the 
  first line, one line per state
 -there must be one column per each factor, named by their names
 -there must be one column for the amount/area of units in the state, named
  freely
 -there must be one column per activity, named by the name of the activity
  giving the probability that a unit in the state receives this activtiy

 number of simulation rounds is given as a variable near the end of the file,
 just before the simulation loop
