

GWAS.current.generation.validate.MAS.in.BC.generation <- function(the.phenotypes.current.gen = NA, the.genotypes.current.gen = NA,
                                                                  the.genotypes.BC = NA,  the.phenotypes.BC = NA,
                                                                  the.output.dir = NA, GS.or.MAS = NA,
                                                                  number.of.markers.to.use.for.MAS = NA){
  
  the.genotypes <- data.frame(the.genotypes.current.gen, the.genotypes.BC[,-c(1:5)])
  the.phenotypes <- rbind(the.phenotypes.current.gen[,1:2], the.phenotypes.BC)
  pred <- (nrow(the.phenotypes.current.gen)+1):nrow(the.phenotypes)
  
  
  setwd(the.output.dir)
  myGD <- data.frame(colnames(the.genotypes.current.gen)[-c(1:5)],t(the.genotypes.current.gen[,-c(1:5)]))
  myGM <- data.frame(the.genotypes.current.gen[,c(1,3,4)])
  
  
  
#Read in a phenotype
#k <- number.of.folds - 1
r.gy <- NULL
the.coefficients <- NULL
count <- 0
for (p in 2:ncol(the.phenotypes)){
  print(paste("===================== Now working on trait ", (p-1), " !!!!!!! ============================================="))
  myY <- the.phenotypes.current.gen[,c(1,p)]

  #for (i in 0:k){ 
    #pick out of he training set (t.s.) (k-1) folds
    
    #print(paste("----------------Starting ", (i+1), " out of ", number.of.folds, " folds-----------------", sep = ""))
    
    #pred <- random.sample[((increment*i)+1):min(((increment*i)+increment) , sample.size)]
    #train <- random.sample[-(((increment*i)+1):min(((increment*i)+increment) , sample.size))] 
    
    #Perform GWAS on T.S.
    #myY.train <- myY[train,]
  
 
      
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
      
      
      #Read in the GWAS output file
      
      #Note to Brian and Alex: we need to get this to work for multiple traits
      Gwas.output<-read.csv(paste("GAPIT..",colnames(myY)[2],".GWAS.Results.csv", sep = ""))

      #Extract the SNP name; also extract the chromosome and bp information
      print(paste("-------Now fitting the peak markers from the current generation into the BC generation ----------------------", sep = ""))
      
      these.markers <- as.character(Gwas.output[1:number.of.markers.to.use.for.MAS,1])
      
      the.genotypic.data.on.these.markers <- the.genotypes[which(as.character(the.genotypes$Snp) %in% these.markers),]
      the.genotypic.data.on.these.markers.formatted.for.lm <- data.frame(colnames(the.genotypic.data.on.these.markers)[-c(1:5)],t(the.genotypic.data.on.these.markers[,-c(1:5)]))
      colnames(the.genotypic.data.on.these.markers.formatted.for.lm)[1] <- "Taxa.Names"
      myY.plus.order <- data.frame(the.phenotypes[,c(1,p)],1:nrow(the.phenotypes))
      colnames(myY.plus.order)[1] <- "<Trait>"
      
      data.for.lm.almost <- data.frame(myY.plus.order, the.genotypic.data.on.these.markers.formatted.for.lm[,-1])
      #data.for.lm.almost <- data.for.lm.almost[order(data.for.lm.almost[,3]),]
    
      
   
      data.for.lm.train <- data.for.lm.almost[-pred,]
      data.for.lm.pred <- data.for.lm.almost[pred,]
      equation.for.lm <- paste(colnames(data.for.lm.train)[2], "~", colnames(data.for.lm.train)[4],sep = "")
      for(index in 5:ncol(data.for.lm.train)) equation.for.lm <- paste(equation.for.lm,colnames(data.for.lm.train)[index],sep = "+")
      
      
      lm.model.fitted.in.train <- lm(equation.for.lm, data = data.for.lm.train)
      the.predicted.MAS.values <- predict(lm.model.fitted.in.train, newdata = data.for.lm.pred)
      
      these.observed.and.predicted.phenotypic.values <- data.frame(data.for.lm.pred[,1:2], the.predicted.MAS.values)
      colnames(these.observed.and.predicted.phenotypic.values) <- c("Taxa", "Observed.Value", "Predicted.Value")
      #if(i == 0){
      the.observed.and.predicted.phenotypic.values <- these.observed.and.predicted.phenotypic.values
      #}else{
       # the.observed.and.predicted.phenotypic.values <- rbind(the.observed.and.predicted.phenotypic.values,these.observed.and.predicted.phenotypic.values)        
      #}#end if(i == 0)
    
    #Measure correclation between OBS and Pred in validation set (V.S.)
    r.gy <- c(r.gy, cor(these.observed.and.predicted.phenotypic.values[,3], these.observed.and.predicted.phenotypic.values[,2]))
    
    #Fit a linear regression model where the observed values are the response variable and the predicted values is the explanatory variable 
    the.fitted.regression.model <- lm(these.observed.and.predicted.phenotypic.values[,2] ~ these.observed.and.predicted.phenotypic.values[,3])
    
    the.coefficients <- c(the.coefficients, the.fitted.regression.model$coefficients[1], the.fitted.regression.model$coefficients[2])
  
  #}#end for loop through the folds (12/14/2017) 
  #calcualte the average and std. deviation
  r.gy <- c(r.gy,(p-1))
  the.coefficients <- c(the.coefficients, (p-1))
  
  r.gy.output <- t(as.matrix(r.gy))
  the.coefficients.output <- t(as.matrix(the.coefficients))
  
  #Sort the table of observed and predicted phenotypic values from smallest to largest
  the.observed.and.predicted.phenotypic.values <- the.observed.and.predicted.phenotypic.values[order(the.observed.and.predicted.phenotypic.values$Predicted.Value, decreasing = TRUE),]
  
  #Output the table of observed and predicted phenotypic values
  write.table(the.observed.and.predicted.phenotypic.values ,paste( "MAS.Predicted.Values.trait.",(p-1),"trained.in.current.gen.validated.in.",GS.or.MAS,".BC.Gen.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
  
  
  
  
  if(count == 0){
    r.gy.m <- r.gy.output
    the.coefficients.m <- the.coefficients.output
    
  }else{
    r.gy.m <- rbind(r.gy.m,r.gy.output) 
    the.coefficients.m <- rbind(the.coefficients.m, the.coefficients.output)
   
  }#end if(count == 0)
  
  
  
  
  #Create a list of the top X% lines with optimal predicted phenotypes
  
  count <- count+1
  r.gy <- NULL
  the.coefficients <- NULL
  #the.observed.and.predicted.phenotypic.values <- NULL
  
  #Reset the "the.observed.and.predicted.phenotypic.values" to NULL
}#for (p in 1:100){

colnames(r.gy.m)<-c("Correlation","Trait")
colnames(the.coefficients.m) <- c("Intercept", "Slope", "trait")

write.table(r.gy.m ,paste( "MAS.Correl.obs.and.pred.trained.in.current.gen.validated.in.",GS.or.MAS,".BC.Gen.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)

write.table(the.coefficients.m, paste( "MAS.SLR.intercept.and.slope.trained.in.current.gen.validated.in.",GS.or.MAS,".BC.Gen.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)


return(list(r.gy.m = r.gy.m, the.coefficients.m = the.coefficients.m, the.observed.and.predicted.phenotypic.values = the.observed.and.predicted.phenotypic.values))

}#end the function

