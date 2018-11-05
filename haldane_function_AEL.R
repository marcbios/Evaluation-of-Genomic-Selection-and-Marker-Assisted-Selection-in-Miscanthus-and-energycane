

#haldanes function
#Calculate likelihood of recombination based on genetic map
#X = genetic map file
#cM-col = column of genetic map file containing cM position
# m1 and m2 correspond to the markers that function is applied to

haldane <- function(x,cM_col,m1,m2){
  
  # c = probability of an odd number of recombination events between two markers
  
  # c = (1-e^(-2m))/2  
  
M <- x[,cM_col]/100
a <- M[m1]
b <- M[m2]

# where m = map distance between 2 markers
m <- b-a

c = (1-exp((-2)*m))/2
return(c)
}


