



rrBLUP.current.generation.validate.in.BC.generation <- function(the.phenotypes.current.gen = NULL, the.genotypes.current.gen = NULL,
                                                  the.genotypes.BC = NULL, the.phenotypes.BC = NULL,
                                                  the.output.dir = NULL, GS.or.MAS = NULL){
  
  the.genotypes <- data.frame(the.genotypes.current.gen, the.genotypes.BC[,-c(1:5)])
  
  the.phenotypes <- rbind(the.phenotypes.current.gen[,c(1:2)], the.phenotypes.BC)
  
  setwd(the.output.dir)
  
  y <- as.matrix(the.phenotypes[,2])
  
  G <- as.numeric(t(the.genotypes[,-c(1:5)]))
  
  G <- matrix(G, nrow(y), nrow(the.genotypes)) #Note: nrow(the genotypes) is specifying the number of columsns for G, whcih we ant
         #to be equal to the number of markers
  
  G <- G - 1
  
  cv.for.rrBLUP <- (as.matrix(rep(1, length(y))))
  
  #Calculate the kinship matrix in rrBLUP
  A1 <- A.mat(G,shrink=TRUE)
  
  library(sommer)
  A2 <- A.mat(G,shrink=TRUE) # additive relationship matrix 
  D1 <- D.mat(G,shrink=TRUE) # dominance relationship matrix 
  E1 <- E.mat(G,shrink=TRUE) # epistatic relationship matrix
  
  M <-tcrossprod(G)/ncol(G)
  X <- G
  
  #burn=10
  #Iter=600
  
  #Save all of the above work into an object
  #save.image("Workspace_20170815.Rdata")
  
  sample.size <- length(y)
  pred <- (nrow(the.phenotypes.current.gen)+1):nrow(the.phenotypes) #Obtain the indicies of the indivdiuals that are in the backcross generation
  
  

  count <- 0
  r.gy <- NULL
  r.gy.E <- NULL
  
  r.gy.BA <- NULL
  r.gy.RK <- NULL
  r.gy.SV <- NULL
  
  
  the.coefficients <- NULL
  the.coefficients.E <- NULL
  
  the.coefficients.BA <- NULL
  the.coefficients.RK <- NULL
  the.coefficients.SV <- NULL
  
  for (p in 2:ncol(the.phenotypes)){
    if((p >2)&(floor(p/10)==p/10)) {print(paste("--------------------------Working on the ", (p-1), "'th trait--------------------------------",sep=""))}
    
    y <- the.phenotypes[,p]
    
    Za <- diag(length(y)) 
    Zd <- diag(length(y)) 
    Ze <- diag(length(y))
    #for (i in 0:k){
    #print(paste("-------Now fitting the RR-BLUP model for fold ", (i+1), " -----------", sep = ""))
    #pred <- random.sample[((increment*i)+1):min(((increment*i)+increment) , sample.size)]
    #train <- random.sample[-(((increment*i)+1):min(((increment*i)+increment) , sample.size))] 
    
  
    yNA <- y
    yNA <- as.vector(yNA)
    yNA[pred] <- NA
    
    data1 <- data.frame(y=yNA,gid=1:length(y), cv = cv.for.rrBLUP)
    the.cv.names <- NULL
    for(j in 1:ncol(cv.for.rrBLUP)) the.cv.names <- c(the.cv.names, paste("CV_",j,sep = ""))
    
    colnames(data1) <- c("y","gid", the.cv.names)
    
    
    rownames(A1) <- 1:nrow(A1) #A1 is created on line 114
    ans1 <- kin.blup(data1,K=A1,geno="gid",pheno="y", covariate = the.cv.names)
    #Measure correclation between OBS and Pred in validation set (V.S.)
    r.gy <- c(r.gy, cor(ans1$g[pred], y[pred]) )
    #Fit a linear regression model, where the Y variable is the observed value and the x variabls is the predicted value
    the.fitted.model <- lm(y[pred] ~ ans1$g[pred])
    the.coefficients <- c(the.coefficients, the.fitted.model$coefficients[1], the.fitted.model$coefficients[2])
    
    
    # Fit sommer models
    #### ADDITIVE MODEL #### 
    rownames(A2) <- 1:nrow(A2)
    ETA.A <- list(add=list(Z=Za,K=A2)) 
    ans.A <- mmer(Y=data1$y, Z=ETA.A) 
    
    
    #### ADDITIVE-DOMINANCE MODEL ####
    rownames(D1) <- 1:nrow(D1)
    ETA.D <- list(add=list(Z=Za,K=A2), dom=list(Z=Zd,K=D1)) 
    ans.D <- mmer(Y=data1$y, Z=ETA.D) 
    
    #### ADDITIVE-DOMINANCE-EPISTASIS MODEL ####
    rownames(E1) <- 1:nrow(E1)
    ETA.E <- list(add=list(Z=Za,K=A2), dom=list(Z=Zd,K=D1), epi=list(Z=Ze,K=E1)) 
    ans.E <- mmer(Y=data1$y, Z=ETA.E) 
    r.gy.E <- c(r.gy.E, cor(ans.E$fitted.y[pred], y[pred]) )
    the.fitted.model.E <- lm(y[pred] ~ ans.E$fitted.y[pred])
    the.coefficients.E <- c(the.coefficients.E, the.fitted.model.E$coefficients[1], the.fitted.model.E$coefficients[2])
    
    
    #### BayesA Parametric MODEL #### 
    ETA<-list(list(X=X,model='BayesA')) 
    fm.BA<-BGLR(y=yNA,ETA=ETA,response_type="gaussian" ,nIter=Iter, burnIn=burn)
    r.gy.BA <- c(r.gy.BA, cor(fm.BA$yHat[pred], y[pred]) )
    the.fitted.model.BA <- lm(y[pred] ~ fm.BA$yHat[pred])
    the.coefficients.BA <- c(the.coefficients.BA, the.fitted.model.BA$coefficients[1], the.fitted.model.BA$coefficients[2])
    
    #### Bayes RKHS Semi Parametric MODEL #### 
    ETA<-list(list(K=M,model='RKHS')) 
    fm.RK<-BGLR(y=yNA,ETA=ETA,response_type="gaussian" ,nIter=Iter, burnIn=burn)
    r.gy.RK <- c(r.gy.RK, cor(fm.RK$yHat[pred], y[pred]) )
    the.fitted.model.RK <- lm(y[pred] ~ fm.RK$yHat[pred])
    the.coefficients.RK <- c(the.coefficients.RK, the.fitted.model.RK$coefficients[1], the.fitted.model.RK$coefficients[2])
    
    #### SVM Non Parametric Model ########
    Gtrain <- G[-pred,]
    Gtest <- G[pred,]
    ySVM <- y
    names(ySVM) <- 1:length(y)
    ytrain <- ySVM[-pred]
    ytest <- ySVM[pred]
    
    svp_w <- ksvm(Gtrain, ytrain, type="eps-svr", kernel = "rbfdot")
    yhat <- predict(svp_w, Gtest)
    rownames(yhat) <- names(ytest)
    
    r.gy.SV <- c(r.gy.SV, cor(yhat, ytest) )
    the.fitted.model.SV <- lm(ytest ~ yhat)
    the.coefficients.SV <- c(the.coefficients.SV, the.fitted.model.SV$coefficients[1], the.fitted.model.SV$coefficients[2])
    
    
    #Obtain an object that has three columns. The first column is the taxa names in the validation population, the second column is the observed
    # phenotypic value, and the third column is the predicted phenotypic value
    the.taxa.in.the.validation.population <- as.character(the.phenotypes[as.numeric(rownames(ans1$g[pred])),1]) 
    the.observed.and.predicted.phenotypic.values <- data.frame(the.taxa.in.the.validation.population, y[pred], ans1$g[pred])
    colnames(the.observed.and.predicted.phenotypic.values) <- c("Taxa", "Observed.Value", "Predicted.Value")
    
    
    Emat <- as.vector(ans.E$fitted.y)
    names(Emat) <- 1:length(Emat)
    the.taxa.in.the.validation.population.E <- as.character(the.phenotypes[as.numeric(names(Emat[pred])),1]) 
    the.observed.and.predicted.phenotypic.values.E <- data.frame(the.taxa.in.the.validation.population.E, y[pred], Emat[pred])
    colnames(the.observed.and.predicted.phenotypic.values.E) <- c("Taxa", "Observed.Value", "Predicted.Value")
    
    
    BAmat <- as.vector(fm.BA$yHat)
    names(BAmat) <- 1:length(BAmat)
    the.taxa.in.the.validation.population.BA <- as.character(the.phenotypes[as.numeric(names(BAmat[pred])),1]) 
    the.observed.and.predicted.phenotypic.values.BA <- data.frame(the.taxa.in.the.validation.population.BA, y[pred], BAmat[pred])
    colnames(the.observed.and.predicted.phenotypic.values.BA) <- c("Taxa", "Observed.Value", "Predicted.Value")
    
    
    RKmat <- as.vector(fm.RK$yHat)
    names(RKmat) <- 1:length(RKmat)
    the.taxa.in.the.validation.population.RK <- as.character(the.phenotypes[as.numeric(names(RKmat[pred])),1]) 
    the.observed.and.predicted.phenotypic.values.RK <- data.frame(the.taxa.in.the.validation.population.RK, y[pred], RKmat[pred])
    colnames(the.observed.and.predicted.phenotypic.values.RK) <- c("Taxa", "Observed.Value", "Predicted.Value")
    
    #SVmat <- as.vector(yhat)
    #names(SVmat) <- 1:length(SVmat)
    the.taxa.in.the.validation.population.SV <- as.character(the.phenotypes[as.numeric(rownames(yhat)),1]) 
    the.observed.and.predicted.phenotypic.values.SV <- data.frame(the.taxa.in.the.validation.population.SV, ytest, yhat)
    colnames(the.observed.and.predicted.phenotypic.values.SV) <- c("Taxa", "Observed.Value", "Predicted.Value")
    

    #}#end for (i in 0:k)
    #calcualte the average and std. deviation
    r.gy <- c(r.gy, (p-1))
    the.coefficients <- c(the.coefficients, (p-1))
    r.gy.output <- t(as.matrix(r.gy))
    the.coefficients.output <- t(as.matrix(the.coefficients))
    
    # sommer Epistasis
    r.gy.E <- c(r.gy.E, (p-1))
    the.coefficients.E <- c(the.coefficients.E, (p-1))
    r.gy.output.E <- t(as.matrix(r.gy.E))
    the.coefficients.output.E <- t(as.matrix(the.coefficients.E))
    
    
    # BayesA
    r.gy.BA <- c(r.gy.BA, (p-1))
    the.coefficients.BA <- c(the.coefficients.BA, (p-1))
    r.gy.output.BA <- t(as.matrix(r.gy.BA))
    the.coefficients.output.BA <- t(as.matrix(the.coefficients.BA))
    
   #BayesRK
    r.gy.RK <- c(r.gy.RK, (p-1))
    the.coefficients.RK <- c(the.coefficients.RK, (p-1))
    r.gy.output.RK <- t(as.matrix(r.gy.RK))
    the.coefficients.output.RK <- t(as.matrix(the.coefficients.RK))
    
    #SV
    r.gy.SV <- c(r.gy.SV, (p-1))
    the.coefficients.SV <- c(the.coefficients.SV, (p-1))
    r.gy.output.SV <- t(as.matrix(r.gy.SV))
    the.coefficients.output.SV <- t(as.matrix(the.coefficients.SV))
    
    
    #Sort the table of observed and predicted phenotypic values from smallest to largest
    the.observed.and.predicted.phenotypic.values <- the.observed.and.predicted.phenotypic.values[order(the.observed.and.predicted.phenotypic.values$Predicted.Value, decreasing = TRUE),]
    the.observed.and.predicted.phenotypic.values.E <- the.observed.and.predicted.phenotypic.values.E[order(the.observed.and.predicted.phenotypic.values.E$Predicted.Value, decreasing = TRUE),]
    the.observed.and.predicted.phenotypic.values.BA <- the.observed.and.predicted.phenotypic.values.BA[order(the.observed.and.predicted.phenotypic.values.BA$Predicted.Value, decreasing = TRUE),]
    the.observed.and.predicted.phenotypic.values.RK <- the.observed.and.predicted.phenotypic.values.RK[order(the.observed.and.predicted.phenotypic.values.RK$Predicted.Value, decreasing = TRUE),]
    the.observed.and.predicted.phenotypic.values.SV <- the.observed.and.predicted.phenotypic.values.SV[order(the.observed.and.predicted.phenotypic.values.SV$Predicted.Value, decreasing = TRUE),]
    
    
    #Output the table of observed and predicted phenotypic values
    write.table(the.observed.and.predicted.phenotypic.values ,paste("rrB.Pred.Val.trt.",(p-1),"trnd.in.current.gen.vald.in.",GS.or.MAS,".BC.Gen.txt",".txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
    write.table(the.observed.and.predicted.phenotypic.values.E ,paste( "Pred.Val.trt.smADE.trait",(p-1),"trnd.in.current.gen.vald.in.",GS.or.MAS,".BC.Gen.txt",".txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
    write.table(the.observed.and.predicted.phenotypic.values.BA ,paste( "Pred.Val.trt.BA.trait",(p-1),"trnd.in.current.gen.vald.in.",GS.or.MAS,".BC.Gen.txt",".txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
    write.table(the.observed.and.predicted.phenotypic.values.RK ,paste( "Pred.Val.trt.RK.trait",(p-1),"trnd.in.current.gen.vald.in.",GS.or.MAS,".BC.Gen.txt",".txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
    write.table(the.observed.and.predicted.phenotypic.values.SV ,paste( "Pred.Val.trt.SV.trait",(p-1),"trnd.in.current.gen.vald.in.",GS.or.MAS,".BC.Gen.txt",".txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
    
    
    if(count == 0){
      r.gy.m <- r.gy.output
      r.gy.m.E <- r.gy.output.E
      r.gy.m.BA <- r.gy.output.BA
      r.gy.m.RK <- r.gy.output.RK
      r.gy.m.SV <- r.gy.output.SV
      
      the.coefficients.m <- the.coefficients.output
      the.coefficients.m.E <- the.coefficients.output.E
      the.coefficients.m.BA <- the.coefficients.output.BA
      the.coefficients.m.RK <- the.coefficients.output.RK
      the.coefficients.m.SV <- the.coefficients.output.SV
      
      
      
      
      
    }else{
      r.gy.m <- rbind(r.gy.m,r.gy.output)
      the.coefficients.m <- rbind(the.coefficients.m, the.coefficients.output)
      
      r.gy.m.E <- rbind(r.gy.m.E,r.gy.output.E) 
      the.coefficients.m.E <- rbind(the.coefficients.m.E, the.coefficients.output.E)
      
      r.gy.m.BA <- rbind(r.gy.m.BA,r.gy.output.BA) 
      the.coefficients.m.BA <- rbind(the.coefficients.m.BA, the.coefficients.output.BA)
      
      r.gy.m.RK <- rbind(r.gy.m.RK,r.gy.output.RK) 
      the.coefficients.m.RK <- rbind(the.coefficients.m.RK, the.coefficients.output.RK)
      
      r.gy.m.SV <- rbind(r.gy.m.SV,r.gy.output.SV) 
      the.coefficients.m.SV <- rbind(the.coefficients.m.SV, the.coefficients.output.SV)
      
      
    }#end if(count == 0)
    
  
   
    
    #Create a list of the top X% lines with optimal predicted phenotypes
    
    count <- count+1
    r.gy <- NULL
    the.coefficients <- NULL
    
    r.gy.E <- NULL
    the.coefficients.E <- NULL
    
    r.gy.BA <- NULL
    the.coefficients.BA <- NULL
    
    r.gy.RK <- NULL
    the.coefficients.RK <- NULL
    
    r.gy.SV <- NULL
    the.coefficients.SV <- NULL
        #the.observed.and.predicted.phenotypic.values <- NULL
    
    #Reset the "the.observed.and.predicted.phenotypic.values" to NULL
    
  }#for (p in 1:100){
  
  
  colnames(r.gy.m)<-c("Correlation","Trait")
  colnames(the.coefficients.m) <- c("Intercept", "Slope", "trait")
  
  write.table(r.gy.m ,paste( "Cor.obs.and.pred.trnd.in.curr.gen.vald.in.",GS.or.MAS,".BC.Gen.rrB.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
  
  write.table(the.coefficients.m, paste( "SLR.intpt.and.slp.trnd.in.curr.gen.vald.in.",GS.or.MAS,".BC.Gen.rrB.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)

  ls.rrBLUP <- list(r.gy.m = r.gy.m, the.coefficients.m = the.coefficients.m, the.observed.and.predicted.phenotypic.values = the.observed.and.predicted.phenotypic.values)
  
  
  colnames(r.gy.m.E)<-c("Correlation","Trait")
  colnames(the.coefficients.m.E) <- c("Intercept", "Slope", "trait")
  
  write.table(r.gy.m.E ,paste( "Cor.obs.and.pred.trnd.in.curr.gen.vald.in.",GS.or.MAS,".BC.Gen.E.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
  
  write.table(the.coefficients.m.E, paste( "SLR.intpt.and.slp.trnd.in.curr.gen.vald.in.",GS.or.MAS,".BC.Gen.E.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
  
  ls.sommerE <- list(r.gy.m.E = r.gy.m.E, the.coefficients.m.E = the.coefficients.m.E, the.observed.and.predicted.phenotypic.values.E = the.observed.and.predicted.phenotypic.values.E)
  
  
  colnames(r.gy.m.BA)<-c("Correlation","Trait")
  colnames(the.coefficients.m.BA) <- c("Intercept", "Slope", "trait")
  
  write.table(r.gy.m.BA ,paste( "Cor.obs.and.pred.trnd.in.curr.gen.vald.in.",GS.or.MAS,".BC.Gen.BA.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
  
  write.table(the.coefficients.m.BA, paste( "SLR.intpt.and.slp.trnd.in.curr.gen.vald.in.",GS.or.MAS,".BC.Gen.BA.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
  
  ls.BA <- list(r.gy.m.BA = r.gy.m.BA, the.coefficients.m.BA = the.coefficients.m.BA, the.observed.and.predicted.phenotypic.values.BA = the.observed.and.predicted.phenotypic.values.BA)
  
  
  colnames(r.gy.m.RK)<-c("Correlation","Trait")
  colnames(the.coefficients.m.RK) <- c("Intercept", "Slope", "trait")
  
  write.table(r.gy.m.RK ,paste( "Cor.obs.and.pred.trnd.in.curr.gen.vald.in.",GS.or.MAS,".BC.Gen.RK.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
  
  write.table(the.coefficients.m.RK, paste( "SLR.intpt.and.slp.trnd.in.curr.gen.vald.in.",GS.or.MAS,".BC.Gen.RK.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
  
  ls.RK <- list(r.gy.m.RK = r.gy.m.RK, the.coefficients.m.RK = the.coefficients.m.RK, the.observed.and.predicted.phenotypic.values.RK = the.observed.and.predicted.phenotypic.values.RK)
  
  
  
  colnames(r.gy.m.SV)<-c("Correlation","Trait")
  colnames(the.coefficients.m.SV) <- c("Intercept", "Slope", "trait")
  
  write.table(r.gy.m.SV ,paste( "Cor.obs.and.pred.trnd.in.curr.gen.vald.in.",GS.or.MAS,".BC.Gen.SV.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
  
  write.table(the.coefficients.m.SV, paste( "SLR.intpt.and.slp.trnd.in.curr.gen.vald.in.",GS.or.MAS,".BC.Gen.SV.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
  
  ls.SV <- list(r.gy.m.SV = r.gy.m.SV, the.coefficients.m.SV = the.coefficients.m.SV, the.observed.and.predicted.phenotypic.values.SV = the.observed.and.predicted.phenotypic.values.SV)
  
  return(list(ls.rrBLUP, ls.sommerE, ls.BA, ls.RK, ls.SV))
}#end rrBLUP.current.generation.validate.in.BC.generation





