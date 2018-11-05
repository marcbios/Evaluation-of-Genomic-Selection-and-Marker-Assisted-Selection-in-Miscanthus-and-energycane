
identify.top.performing.lines.via.MAS <- function(the.phenotypes = NULL, the.genotypes = NULL, the.output.dir = NULL,
                                                  number.of.folds = NULL, proportion.to.select = NULL, number.of.markers.to.use.for.MAS = NULL,
                                                  user.input.seed.number = FALSE, seed.number = NULL){
  
  
  setwd(the.output.dir)
  s=proportion.to.select
  myGD <- data.frame(colnames(the.genotypes)[-c(1:5)],t(the.genotypes[,-c(1:5)]))
  myGM <- data.frame(the.genotypes[,c(1,3,4)])
  
  sample.size <- nrow(the.phenotypes)
  if(!user.input.seed.number) seed.number <- sample(-1000000:1000000,1, replace = FALSE)
  set.seed(seed.number)
  sequence.sample <- rep(1:sample.size)
  random.sample <- sample(1:sample.size, replace = FALSE)
  increment <- ceiling(length(random.sample)/number.of.folds) 
  
#Read in a phenotype
k <- number.of.folds - 1
r.gy <- NULL
myCI <- NULL
the.coefficients <- NULL
count <- 0
for (p in 2:ncol(the.phenotypes)){
  print(paste("===================== Now working on trait ", (p-1), " !!!!!!! ============================================="))
  myY <- the.phenotypes[,c(1,p)]

  for (i in 0:k){ 
    #pick out of he training set (t.s.) (k-1) folds
    
    print(paste("----------------Starting ", (i+1), " out of ", number.of.folds, " folds-----------------", sep = ""))
    
    pred <- random.sample[((increment*i)+1):min(((increment*i)+increment) , sample.size)]
    train <- random.sample[-(((increment*i)+1):min(((increment*i)+increment) , sample.size))] 
    
    #Perform GWAS on T.S.
    myY.train <- myY[train,]
  
 
      
      #Step 2: Run GAPIT
      myGAPIT <- try(GAPIT(
        Y=myY.train, 
        PCA.total=0,
        GD = myGD,
        GM = myGM,
        Geno.View.output=FALSE,
        group.from = nrow(myY.train),
        group.to = nrow(myY.train)
      ))
      
      if(inherits(myGAPIT , "try-error")){
        print(paste("GAPIT most likely stopped working because you have only one marker anchored a chromosome and/or linakge group. Please double check the input genotypic data", sep = ""))
        break
      }#end if(inherits(myGAPIT , "try-error"))
      
      
      #Read in the GWAS output file
      
      #Note to Brian and Alex: we need to get this to work for multiple traits
      Gwas.output<-read.csv(paste("GAPIT..",colnames(myY.train)[2],".GWAS.Results.csv", sep = ""))

      #Extract the SNP name; also extract the chromosome and bp information
      print(paste("-------Now fitting the peak marker from training set into validation set model for fold ", (i+1), " -----------", sep = ""))
      
      these.markers <- as.character(Gwas.output[1:number.of.markers.to.use.for.MAS,1])
      
      the.genotypic.data.on.these.markers <- the.genotypes[which(as.character(the.genotypes$Snp) %in% these.markers),]
      the.genotypic.data.on.these.markers.formatted.for.lm <- data.frame(colnames(the.genotypic.data.on.these.markers)[-c(1:5)],t(the.genotypic.data.on.these.markers[,-c(1:5)]))
      colnames(the.genotypic.data.on.these.markers.formatted.for.lm)[1] <- "Taxa.Names"
      myY.plus.order <- data.frame(myY,1:nrow(myY))
      colnames(myY.plus.order)[1] <- "<Trait>"
      
      data.for.lm.almost <- merge(myY.plus.order, the.genotypic.data.on.these.markers.formatted.for.lm, by.x = "<Trait>", by.y = "Taxa.Names")
      data.for.lm.almost <- data.for.lm.almost[order(data.for.lm.almost[,3]),]
    
      
   
      data.for.lm.train <- data.for.lm.almost[train,]
      data.for.lm.pred <- data.for.lm.almost[pred,]
      equation.for.lm <- paste(colnames(data.for.lm.train)[2], "~", colnames(data.for.lm.train)[4],sep = "")
      for(index in 5:ncol(data.for.lm.train)) equation.for.lm <- paste(equation.for.lm,colnames(data.for.lm.train)[index],sep = "+")
      
      
      lm.model.fitted.in.train <- lm(equation.for.lm, data = data.for.lm.train)
      the.predicted.MAS.values <- predict(lm.model.fitted.in.train, newdata = data.for.lm.pred)
      
      these.observed.and.predicted.phenotypic.values <- data.frame(data.for.lm.pred[,1:2], the.predicted.MAS.values)
      colnames(these.observed.and.predicted.phenotypic.values) <- c("Taxa", "Observed.Value", "Predicted.Value")
      
      
      if(i == 0){
        the.observed.and.predicted.phenotypic.values <- these.observed.and.predicted.phenotypic.values
      }else{
        the.observed.and.predicted.phenotypic.values <- rbind(the.observed.and.predicted.phenotypic.values,these.observed.and.predicted.phenotypic.values)        
      }#end if(i == 0)
    
    #Measure correclation between OBS and Pred in validation set (V.S.)
    r.gy <- c(r.gy, cor(these.observed.and.predicted.phenotypic.values[,3], these.observed.and.predicted.phenotypic.values[,2]))
    x.p=these.observed.and.predicted.phenotypic.values[,c(1,3)]
    y.o=these.observed.and.predicted.phenotypic.values[,c(1,2)]
    myCI <- c(myCI, round(CI(x.p,y.o,s=s,top=T),2))
    #Fit a linear regression model where the observed values are the response variable and the predicted values is the explanatory variable 
    the.fitted.regression.model <- lm(these.observed.and.predicted.phenotypic.values[,2] ~ these.observed.and.predicted.phenotypic.values[,3])
    
    the.coefficients <- c(the.coefficients, the.fitted.regression.model$coefficients[1], the.fitted.regression.model$coefficients[2])
  
  } 
  #calcualte the average and std. deviation
  r.gy <- c(r.gy, mean(r.gy), sd(r.gy),(p-1))
  the.coefficients <- c(the.coefficients, (p-1))
  
  myCI <- c(myCI, mean(myCI), sd(myCI), (p-1))
  myCI.output <- t(as.matrix(myCI))
  
  r.gy.output <- t(as.matrix(r.gy))
  the.coefficients.output <- t(as.matrix(the.coefficients))
  
  #Sort the table of observed and predicted phenotypic values from smallest to largest
  the.observed.and.predicted.phenotypic.values <- the.observed.and.predicted.phenotypic.values[order(the.observed.and.predicted.phenotypic.values$Predicted.Value, decreasing = TRUE),]
  
  #Output the table of observed and predicted phenotypic values
  write.table(the.observed.and.predicted.phenotypic.values ,paste( "Predicted.Values.MAS.trait.",(p-1),".txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
  
  number.of.taxa.to.select <- round((proportion.to.select*nrow(myY)),0)
  
  
  if(count == 0){
    r.gy.m <- r.gy.output
    myCI.m <- myCI.output
    the.coefficients.m <- the.coefficients.output
    the.selected.taxa <- as.character(the.observed.and.predicted.phenotypic.values[1:number.of.taxa.to.select,1])
  }else{
    r.gy.m <- rbind(r.gy.m,r.gy.output)
    myCI.m <- rbind(myCI.m, myCI.output)
    the.coefficients.m <- rbind(the.coefficients.m, the.coefficients.output)
    the.selected.taxa <- c(the.selected.taxa, as.character(the.observed.and.predicted.phenotypic.values[1:number.of.taxa.to.select,1]))
  }#end if(count == 0)
  
  
  
  
  #Create a list of the top X% lines with optimal predicted phenotypes
  
  count <- count+1
  r.gy <- NULL
  myCI <- NULL
  the.coefficients <- NULL
  the.observed.and.predicted.phenotypic.values <- NULL
  
  #Reset the "the.observed.and.predicted.phenotypic.values" to NULL
}#for (p in 1:100){

colnames(r.gy.m)<-c("fold.1","fold.2","fold.3","fold.4","fold.5","mean","sd","trait")
colnames(the.coefficients.m) <- c("Intercpet.fold.1", "Slope.fold.1", "Intercpet.fold.2", "Slope.fold.2",
                                  "Intercpet.fold.3", "Slope.fold.3", "Intercpet.fold.4", "Slope.fold.4",
                                  "Intercpet.fold.5", "Slope.fold.5","trait")

write.table(r.gy.m ,paste( "MAS.Correlation.between.obs.and.pred.with.",number.of.folds,".folds.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)

write.table(the.coefficients.m, paste( "MAS.SLR.intercept.and.slope.Y.observed.X.pred.with.",number.of.folds,".folds.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)

the.selected.taxa <- the.selected.taxa[!the.selected.taxa%in%parents.vec]
the.unique.selected.taxa <- unique(the.selected.taxa)
write.table(the.unique.selected.taxa, paste( "MAS.the.unique.selected.taxa.",number.of.folds,".folds.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)

colnames(myCI.m)<-c("fold.1","fold.2","fold.3","fold.4","fold.5","mean","sd","trait")
write.table(myCI.m,paste( "CI.betw.obs.pred.",number.of.folds,".folds.MAS.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)

return(the.unique.selected.taxa)
  
}#end the function

