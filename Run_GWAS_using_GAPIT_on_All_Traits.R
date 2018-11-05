
run.GAPIT.on.the.simulated.traits <- function(the.phenotypes = NULL, the.genotypes = NULL, the.output.dir = NULL)
{
  
  setwd(the.output.dir)
  myGD <- data.frame(colnames(the.genotypes)[-c(1:5)],t(the.genotypes[,-c(1:5)]))
  myGM <- data.frame(the.genotypes[,c(1,3,4)])
  #Read in a phenotype
  for (p in 2:ncol(the.phenotypes)){
    print(paste("===================== Now working on trait ", (p-1), " !!!!!!! ============================================="))
      myY <- the.phenotypes[,c(1,p)]
  
      #Step 2: Run GAPIT
      myGAPIT <- try(GAPIT(
        Y=myY, 
        PCA.total=0,
        GD = myGD,
        GM = myGM,
        Geno.View.output=FALSE,
        group.from = nrow(myY),
        group.to = nrow(myY)
      ))
      
      if(inherits(myGAPIT , "try-error")){
        print(paste("GAPIT most likely stopped working because you have only one marker anchored a chromosome and/or linakge group. Please double check the input genotypic data", sep = ""))
        break
      }#end if(inherits(myGAPIT , "try-error"))
     
  }#for (p in 1:100){

}#end run.GAPIT.on.the.simulated.traits
