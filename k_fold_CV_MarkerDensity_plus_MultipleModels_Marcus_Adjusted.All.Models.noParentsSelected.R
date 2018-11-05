run.k.fold.CV.using.rrBLUP <- function(the.phenotypes = NULL, the.genotypes = NULL, the.output.dir = NULL, number.of.folds = NULL,
                                       proportion.to.select = NULL, user.input.seed.number = FALSE, seed.number = NULL, geno.proportion=NULL, ncycles=NULL){
  
    setwd(the.output.dir)
    
    seed.number = seednum.specific
    
    y <- as.matrix(the.phenotypes[,2])
    G <- as.numeric(t(the.genotypes[,-c(1:5)]))
    
    G <- matrix(G, nrow(y), nrow(the.genotypes)) #Note: nrow(the genotypes) is specifying the number of columsns for G, whcih we ant
      #to be equal to the number of markers
 
  
  #for(props in 1:length(props.vec)){
   #geno.proportion <- props.vec[props]
    
    for(ncycle in 1:ncycles){#
      
          colnames(G) <- 1:ncol(G)
          test_proptn <- round(ncol(G)*geno.proportion)
          set.seed(seed.numer.cycle[ncycle])
          train.G = as.matrix(sample(1:ncol(G), test_proptn, replace=F))
            
          Gt <- G[,train.G]
            
            
          G2 <- Gt
          Gt <- Gt - 1
          
          cv.for.rrBLUP <- (as.matrix(rep(1, length(y))))
          #cv.for.sommer <- (as.matrix(rep(1, length(y))))
          
          #Calculate the kinship matrix in rrBLUP
          Gt <- as.matrix(Gt)
          A1 <- A.mat(Gt,shrink=TRUE)
          library(sommer)
          A2 <- A.mat(Gt,shrink=TRUE) # additive relationship matrix 
          D1 <- D.mat(Gt,shrink=TRUE) # dominance relationship matrix 
          E1 <- E.mat(Gt,shrink=TRUE) # epistatic relationship matrix
         
          M <-tcrossprod(G2)/ncol(G2)
          X <- G2
          
          #burn=10
          #Iter=600
          
          #Save all of the above work into an object
          #save.image("Workspace_20170815.Rdata")
          
          sample.size <- length(y)
          #if(!user.input.seed.number) seed.number <- sample(-1000000:1000000,1, replace = FALSE)
          seed.number = seednum.specific
          set.seed(seed.number)
          sequence.sample <- rep(1:sample.size)
          random.sample <- sample(1:sample.size, replace = FALSE)
          increment <- ceiling(length(random.sample)/number.of.folds) 
          
          write.table(paste("Here_is_seed_number:_",seed.number, sep = ""), paste("Seed.number.for.",number.of.folds,".fold.CV.txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
          #have a "for" loop, start it at 0, and end it at 4
          #I am setting up "k" to denote the nubmer of folds - 1. This is done
          # so that the for loop will work correctly.
          r.gy.m<-matrix(NA,nrow=1,ncol=(number.of.folds+4))
          r.gy.m.E<-matrix(NA,nrow=1,ncol=(number.of.folds+4))
          r.gy.m.BA<-matrix(NA,nrow=1,ncol=(number.of.folds+4))
          r.gy.m.BC<-matrix(NA,nrow=1,ncol=(number.of.folds+4))
          r.gy.m.RK<-matrix(NA,nrow=1,ncol=(number.of.folds+4))
          
          #setwd(wd)
          
          
          
          #Read in a phenotype
          count <- 0
          k <- number.of.folds-1
          r.gy <- NULL
          r.gy.E <- NULL
          r.gy.BA <- NULL
          #r.gy.BB <- NULL
          r.gy.BC <- NULL
          #r.gy.BL <- NULL
          r.gy.RK <- NULL
          #r.gy.SV <- NULL
          
          the.coefficients <- NULL
          the.coefficients.E <- NULL
          the.coefficients.BA <- NULL
          #the.coefficients.BB <- NULL
          the.coefficients.BC <- NULL
         # the.coefficients.BL <- NULL
          the.coefficients.RK <- NULL
          #the.coefficients.SV <- NULL
          

          #for (p in 2:ncol(the.phenotypes)){
          for (p in 2:ncol(the.phenotypes)){
            if((p >2)&(floor(p/10)==p/10)) {print(paste("--------------------------Working on the ", (p-1), "'th trait--------------------------------",sep=""))}
            
            y <- as.matrix(the.phenotypes[,p])
            Za <- diag(length(y)) 
            Zd <- diag(length(y)) 
            Ze <- diag(length(y))
            
            
            for (i in 0:k){
              print(paste("-------Now fitting the RR-BLUP model for fold ", (i+1), " -----------", sep = ""))
              pred <- random.sample[((increment*i)+1):min(((increment*i)+increment) , sample.size)]
              train <- random.sample[-(((increment*i)+1):min(((increment*i)+increment) , sample.size))] 
              
            
              yNA <- y
              yNA <- as.vector(yNA)
              yNA[pred] <- NA
              
              data1 <- data.frame(y=yNA,gid=1:length(y), cv = cv.for.rrBLUP)
            
              the.cv.names <- NULL
              for(j in 1:ncol(cv.for.rrBLUP)) the.cv.names <- c(the.cv.names, paste("CV_",j,sep = ""))
              
              # Fit rrBLP
              colnames(data1) <- c("y","gid", the.cv.names)
              rownames(A1) <- 1:nrow(A1) #A1 is created on line 114
              ans1 <- kin.blup(data1,K=A1,geno="gid",pheno="y", covariate = the.cv.names)
              #Measure correclation between OBS and Pred in validation set (V.S.)
              r.gy <- c(r.gy, cor(ans1$g[pred], y[pred]))
              #Fit a linear regression model, where the Y variable is the observed value and the x variabls is the predicted value
              the.fitted.model <- lm(y[pred] ~ ans1$g[pred])
              the.coefficients <- c(the.coefficients, the.fitted.model$coefficients[1], the.fitted.model$coefficients[2])
              
              # Fit sommer models
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
          
              
              #### BayesC Parametric MODEL #### 
              ETA<-list(list(X=X,model='BayesC')) 
              fm.BC<-BGLR(y=yNA,ETA=ETA,response_type="gaussian" ,nIter=Iter, burnIn=burn)
              r.gy.BC <- c(r.gy.BC, cor(fm.BC$yHat[pred], y[pred]) )
              the.fitted.model.BC <- lm(y[pred] ~ fm.BC$yHat[pred])
              the.coefficients.BC <- c(the.coefficients.BC, the.fitted.model.BC$coefficients[1], the.fitted.model.BC$coefficients[2])
              
              
              #### Bayes RKHS Semi Parametric MODEL #### 
              ETA<-list(list(K=M,model='RKHS')) 
              fm.RK<-BGLR(y=yNA,ETA=ETA,response_type="gaussian" ,nIter=Iter, burnIn=burn)
              r.gy.RK <- c(r.gy.RK, cor(fm.RK$yHat[pred], y[pred]) )
              the.fitted.model.RK <- lm(y[pred] ~ fm.RK$yHat[pred])
              the.coefficients.RK <- c(the.coefficients.RK, the.fitted.model.RK$coefficients[1], the.fitted.model.RK$coefficients[2])
              
              
              #Obtain an object that has three columns. The first column is the taxa names in the validation population, the second column is the observed
              # phenotypic value, and the third column is the predicted phenotypic value
              the.taxa.in.the.validation.population <- as.character(the.phenotypes[as.numeric(rownames(ans1$g[pred])),1]) 
              these.observed.and.predicted.phenotypic.values <- data.frame(the.taxa.in.the.validation.population, y[pred], ans1$g[pred])
              colnames(these.observed.and.predicted.phenotypic.values) <- c("Taxa", "Observed.Value", "Predicted.Value")
              
              Emat <- as.vector(ans.E$fitted.y)
              names(Emat) <- 1:length(Emat)
              the.taxa.in.the.validation.population.E <- as.character(the.phenotypes[as.numeric(names(Emat[pred])),1]) 
              these.observed.and.predicted.phenotypic.values.E <- data.frame(the.taxa.in.the.validation.population.E, y[pred], Emat[pred])
              colnames(these.observed.and.predicted.phenotypic.values.E) <- c("Taxa", "Observed.Value", "Predicted.Value")
              
              BAmat <- as.vector(fm.BA$yHat)
              names(BAmat) <- 1:length(BAmat)
              the.taxa.in.the.validation.population.BA <- as.character(the.phenotypes[as.numeric(names(BAmat[pred])),1]) 
              these.observed.and.predicted.phenotypic.values.BA <- data.frame(the.taxa.in.the.validation.population.BA, y[pred], BAmat[pred])
              colnames(these.observed.and.predicted.phenotypic.values.BA) <- c("Taxa", "Observed.Value", "Predicted.Value")
              
                   
              BCmat <- as.vector(fm.BC$yHat)
              names(BCmat) <- 1:length(BCmat)
              the.taxa.in.the.validation.population.BC <- as.character(the.phenotypes[as.numeric(names(BCmat[pred])),1]) 
              these.observed.and.predicted.phenotypic.values.BC <- data.frame(the.taxa.in.the.validation.population.BC, y[pred], BCmat[pred])
              colnames(these.observed.and.predicted.phenotypic.values.BC) <- c("Taxa", "Observed.Value", "Predicted.Value")
              
              
              RKmat <- as.vector(fm.RK$yHat)
              names(RKmat) <- 1:length(RKmat)
              the.taxa.in.the.validation.population.RK <- as.character(the.phenotypes[as.numeric(names(RKmat[pred])),1]) 
              these.observed.and.predicted.phenotypic.values.RK <- data.frame(the.taxa.in.the.validation.population.RK, y[pred], RKmat[pred])
              colnames(these.observed.and.predicted.phenotypic.values.RK) <- c("Taxa", "Observed.Value", "Predicted.Value")
              
             
              if(i == 0){
                the.observed.and.predicted.phenotypic.values <- these.observed.and.predicted.phenotypic.values
                the.observed.and.predicted.phenotypic.values.E <- these.observed.and.predicted.phenotypic.values.E
                the.observed.and.predicted.phenotypic.values.BA <- these.observed.and.predicted.phenotypic.values.BA
                the.observed.and.predicted.phenotypic.values.BC <- these.observed.and.predicted.phenotypic.values.BC
                the.observed.and.predicted.phenotypic.values.RK <- these.observed.and.predicted.phenotypic.values.RK
        
                
              }else{
                the.observed.and.predicted.phenotypic.values <- rbind(the.observed.and.predicted.phenotypic.values,these.observed.and.predicted.phenotypic.values)
                the.observed.and.predicted.phenotypic.values.E <- rbind(the.observed.and.predicted.phenotypic.values.E,these.observed.and.predicted.phenotypic.values.E)
                the.observed.and.predicted.phenotypic.values.BA <- rbind(the.observed.and.predicted.phenotypic.values.BA,the.observed.and.predicted.phenotypic.values.BA)
                the.observed.and.predicted.phenotypic.values.BC <- rbind(the.observed.and.predicted.phenotypic.values.BC,these.observed.and.predicted.phenotypic.values.BC)
                the.observed.and.predicted.phenotypic.values.RK <- rbind(the.observed.and.predicted.phenotypic.values.RK,these.observed.and.predicted.phenotypic.values.RK)
                
              }#end if(i == 0)
          
            }#end for (i in 0:k)
            
            #calcualte the average and std. deviation
            # rrBLUP
            r.gy <- c(r.gy, mean(r.gy), sd(r.gy),(p-1))
            the.coefficients <- c(the.coefficients, (p-1))
            r.gy.output <- t(as.matrix(r.gy))
            the.coefficients.output <- t(as.matrix(the.coefficients))
            
            # sommer Epistasis
            r.gy.E <- c(r.gy.E, mean(r.gy.E), sd(r.gy.E),(p-1))
            the.coefficients.E <- c(the.coefficients.E, (p-1))
            r.gy.output.E <- t(as.matrix(r.gy.E))
            the.coefficients.output.E <- t(as.matrix(the.coefficients.E))
            
            # BayesA
            r.gy.BA <- c(r.gy.BA, mean(r.gy.BA), sd(r.gy.BA),(p-1))
            the.coefficients.BA <- c(the.coefficients.BA, (p-1))
            r.gy.output.BA <- t(as.matrix(r.gy.BA))
            the.coefficients.output.BA <- t(as.matrix(the.coefficients.BA))
            
            #BayesC
            r.gy.BC <- c(r.gy.BC, mean(r.gy.BC), sd(r.gy.BC),(p-1))
            the.coefficients.BC <- c(the.coefficients.BC, (p-1))
            r.gy.output.BC <- t(as.matrix(r.gy.BC))
            the.coefficients.output.BC <- t(as.matrix(the.coefficients.BC))
            
            #BayesRK
            r.gy.RK <- c(r.gy.RK, mean(r.gy.RK), sd(r.gy.RK),(p-1))
            the.coefficients.RK <- c(the.coefficients.RK, (p-1))
            r.gy.output.RK <- t(as.matrix(r.gy.RK))
            the.coefficients.output.RK <- t(as.matrix(the.coefficients.RK))
            
            #Sort the table of observed and predicted phenotypic values from smallest to largest
            the.observed.and.predicted.phenotypic.values <- the.observed.and.predicted.phenotypic.values[order(the.observed.and.predicted.phenotypic.values$Predicted.Value, decreasing = TRUE),]
            the.observed.and.predicted.phenotypic.values.E <- the.observed.and.predicted.phenotypic.values.E[order(the.observed.and.predicted.phenotypic.values.E$Predicted.Value, decreasing = TRUE),]
            the.observed.and.predicted.phenotypic.values.BA <- the.observed.and.predicted.phenotypic.values.BA[order(the.observed.and.predicted.phenotypic.values.BA$Predicted.Value, decreasing = TRUE),]
            the.observed.and.predicted.phenotypic.values.BC <- the.observed.and.predicted.phenotypic.values.BC[order(the.observed.and.predicted.phenotypic.values.BC$Predicted.Value, decreasing = TRUE),]
            the.observed.and.predicted.phenotypic.values.RK <- the.observed.and.predicted.phenotypic.values.RK[order(the.observed.and.predicted.phenotypic.values.RK$Predicted.Value, decreasing = TRUE),]
             
            
            #Output the table of observed and predicted phenotypic values
            #write.table(the.observed.and.predicted.phenotypic.values ,paste( "Pred.Vals.rrBLUP.trait.",(p-1),".txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
            #write.table(the.observed.and.predicted.phenotypic.values.E ,paste( "Pred.Vals.sommer.ADE.trait",(p-1),".txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
            #write.table(the.observed.and.predicted.phenotypic.values.BA ,paste( "Pred.Vals.BGLR.BA.trait",(p-1),".txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
            #write.table(the.observed.and.predicted.phenotypic.values.BC ,paste( "Pred.Vals.BGLR.BC.trait",(p-1),".txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
            #write.table(the.observed.and.predicted.phenotypic.values.RK ,paste( "Pred.Vals.BGLR.RK.trait",(p-1),".txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
             
            number.of.taxa.to.select <- round((proportion.to.select*length(y)),0)
            
            
            if(count == 0){
              r.gy.m <- r.gy.output
              r.gy.m.E <- r.gy.output.E
              r.gy.m.BA <- r.gy.output.BA
              r.gy.m.BC <- r.gy.output.BC
              r.gy.m.RK <- r.gy.output.RK
               
              the.coefficients.m <- the.coefficients.output
              the.selected.taxa <- as.character(the.observed.and.predicted.phenotypic.values[1:number.of.taxa.to.select,1])
              
              the.coefficients.m.E <- the.coefficients.output.E
              the.selected.taxa.E <- as.character(the.observed.and.predicted.phenotypic.values.E[1:number.of.taxa.to.select,1])
              
              the.coefficients.m.BA <- the.coefficients.output.BA
              the.selected.taxa.BA <- as.character(the.observed.and.predicted.phenotypic.values.BA[1:number.of.taxa.to.select,1])
              
              the.coefficients.m.BC <- the.coefficients.output.BC
              the.selected.taxa.BC<- as.character(the.observed.and.predicted.phenotypic.values.BC[1:number.of.taxa.to.select,1])
              
              the.coefficients.m.RK <- the.coefficients.output.RK
              the.selected.taxa.RK <- as.character(the.observed.and.predicted.phenotypic.values.RK[1:number.of.taxa.to.select,1])
              
              
            }else{
              r.gy.m <- rbind(r.gy.m,r.gy.output) 
              the.coefficients.m <- rbind(the.coefficients.m, the.coefficients.output)
              the.selected.taxa <- c(the.selected.taxa, as.character(the.observed.and.predicted.phenotypic.values[1:number.of.taxa.to.select,1]))
              
              r.gy.m.E <- rbind(r.gy.m.E,r.gy.output.E) 
              the.coefficients.m.E <- rbind(the.coefficients.m.E, the.coefficients.output.E)
              the.selected.taxa.E <- c(the.selected.taxa.E, as.character(the.observed.and.predicted.phenotypic.values.E[1:number.of.taxa.to.select,1]))
              
              r.gy.m.BA <- rbind(r.gy.m.BA,r.gy.output.BA) 
              the.coefficients.m.BA <- rbind(the.coefficients.m.BA, the.coefficients.output.BA)
              the.selected.taxa.BA <- c(the.selected.taxa.BA, as.character(the.observed.and.predicted.phenotypic.values.BA[1:number.of.taxa.to.select,1]))
              
              r.gy.m.BC <- rbind(r.gy.m.BC,r.gy.output.BC) 
              the.coefficients.m.BC <- rbind(the.coefficients.m.BC, the.coefficients.output.BC)
              the.selected.taxa.BC <- c(the.selected.taxa.BC, as.character(the.observed.and.predicted.phenotypic.values.BC[1:number.of.taxa.to.select,1]))
              
              r.gy.m.RK <- rbind(r.gy.m.RK,r.gy.output.RK) 
              the.coefficients.m.RK <- rbind(the.coefficients.m.RK, the.coefficients.output.RK)
              the.selected.taxa.RK <- c(the.selected.taxa.RK, as.character(the.observed.and.predicted.phenotypic.values.RK[1:number.of.taxa.to.select,1]))
              
              
            }#endif(count == 0)
        
            
            #Create a list of the top X% lines with optimal predicted phenotypes
            
            count <- count+1
            
            #Reset the "the.observed.and.predicted.phenotypic.values" to NULL
            r.gy <- NULL
            the.coefficients <- NULL
            the.observed.and.predicted.phenotypic.values <- NULL
            
            r.gy.E <- NULL
            the.coefficients.E <- NULL
            the.observed.and.predicted.phenotypic.values.E <- NULL
            
            r.gy.BA <- NULL
            the.coefficients.BA <- NULL
            the.observed.and.predicted.phenotypic.values.BA <- NULL
            
            r.gy.BC <- NULL
            the.coefficients.BC <- NULL
            the.observed.and.predicted.phenotypic.values.BC <- NULL
            
            r.gy.RK <- NULL
            the.coefficients.RK <- NULL
            the.observed.and.predicted.phenotypic.values.RK <- NULL
            
            dev.off()
            
          }#for (p in 1:100){
          
          colnames(r.gy.m)<-c("fold.1","fold.2","fold.3","fold.4","fold.5","mean","sd","trait")
          
          colnames(the.coefficients.m) <- c("Intercpet.fold.1", "Slope.fold.1", "Intercpet.fold.2", "Slope.fold.2",
                                            "Intercpet.fold.3", "Slope.fold.3", "Intercpet.fold.4", "Slope.fold.4",
                                            "Intercpet.fold.5", "Slope.fold.5","trait")
          
          write.table(r.gy.m ,paste( "Cor.betw.obs.pred.",ncycle,"_Cycle_",number.of.folds,".folds.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
          
          write.table(the.coefficients.m, paste( "SLR.int.slope.Y.obs.X.pred",ncycle,"_Cycle_",number.of.folds,".folds.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
          
          the.selected.taxa <- the.selected.taxa[!the.selected.taxa%in%parents.vec]
          the.unique.selected.taxa <- unique(the.selected.taxa)
          write.table(the.selected.taxa, paste("the.sel.taxa.RR_BLUP",ncycle,"_Cycle",".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
          write.table(the.unique.selected.taxa, paste("the.uniq.sel.taxa.RR_BLUP",ncycle,"_Cycle",".txt"), sep="\t", quote=F, row.names=F, col.names=F)
          
          
          colnames(r.gy.m.E)<-c("fold.1","fold.2","fold.3","fold.4","fold.5","mean","sd","trait")
          
          colnames(the.coefficients.m.E) <- c("Intercpet.fold.1", "Slope.fold.1", "Intercpet.fold.2", "Slope.fold.2",
                                              "Intercpet.fold.3", "Slope.fold.3", "Intercpet.fold.4", "Slope.fold.4",
                                              "Intercpet.fold.5", "Slope.fold.5","trait")
          
          write.table(r.gy.m.E ,paste( "Cor.betw.obs.pred.",ncycle,"_Cycle_",number.of.folds,".folds.sommerADE.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
          
          write.table(the.coefficients.m.E, paste( "SLR.int.slope.Y.obs.X.pred",ncycle,"_Cycle_",number.of.folds,".folds.sommerADE.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
          
          the.selected.taxa.E <- the.selected.taxa.E[!the.selected.taxa.E%in%parents.vec]
          the.unique.selected.taxa.E <- unique(the.selected.taxa.E)
          
          write.table(the.selected.taxa.E, paste("the.sel.taxa.E",ncycle,"_Cycle",".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
          write.table(the.unique.selected.taxa.E, paste("the.uniq.sel.taxa.sommer.E",ncycle,"_Cycle",".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
          
          colnames(r.gy.m.BA)<-c("fold.1","fold.2","fold.3","fold.4","fold.5","mean","sd","trait")
          
          colnames(the.coefficients.m.BA) <- c("Intercpet.fold.1", "Slope.fold.1", "Intercpet.fold.2", "Slope.fold.2",
                                               "Intercpet.fold.3", "Slope.fold.3", "Intercpet.fold.4", "Slope.fold.4",
                                               "Intercpet.fold.5", "Slope.fold.5","trait")
          
          write.table(r.gy.m.BA ,paste( "Cor.betw.obs.pred.",ncycle,"_Cycle_",number.of.folds,".folds.BGLR.BA.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
          
          write.table(the.coefficients.m.BA, paste( "SLR.int.slope.Y.obs.X.pred",ncycle,"_Cycle_",number.of.folds,".folds.BGLR.BA.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
          
          the.selected.taxa.BA <- the.selected.taxa.BA[!the.selected.taxa.BA%in%parents.vec]
          the.unique.selected.taxa.BA <- unique(the.selected.taxa.BA)
          
          write.table(the.selected.taxa.BA, paste("the.sel.taxa.BA",ncycle,"_Cycle_",".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
          write.table(the.unique.selected.taxa.BA, paste("the.uniq.sel.taxa.BGLR.BA",ncycle,"_Cycle_",".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
          
          
          colnames(r.gy.m.BC)<-c("fold.1","fold.2","fold.3","fold.4","fold.5","mean","sd","trait")
          
          colnames(the.coefficients.m.BC) <- c("Intercpet.fold.1", "Slope.fold.1", "Intercpet.fold.2", "Slope.fold.2",
                                               "Intercpet.fold.3", "Slope.fold.3", "Intercpet.fold.4", "Slope.fold.4",
                                               "Intercpet.fold.5", "Slope.fold.5","trait")
          
          write.table(r.gy.m.BC ,paste( "Cor.betw.obs.pred.",ncycle,"_Cycle_",number.of.folds,".folds.BGLR.BC.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
          
          write.table(the.coefficients.m.BC, paste( "SLR.int.slope.Y.obs.X.pred",ncycle,"_Cycle_",number.of.folds,".folds.BGLR.BC.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
          
          the.selected.taxa.BC <- the.selected.taxa.BC[!the.selected.taxa.BC%in%parents.vec]
          the.unique.selected.taxa.BC <- unique(the.selected.taxa.BC)
          
          write.table(the.selected.taxa.BC, paste("the.sel.taxa.BC",ncycle,"_Cycle_",".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
          write.table(the.unique.selected.taxa.BC, paste("the.uniq.sel.taxa.BGLR.BC",ncycle,"_Cycle_",".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
        
          
          colnames(r.gy.m.RK)<-c("fold.1","fold.2","fold.3","fold.4","fold.5","mean","sd","trait")
          
          colnames(the.coefficients.m.RK) <- c("Intercpet.fold.1", "Slope.fold.1", "Intercpet.fold.2", "Slope.fold.2",
                                               "Intercpet.fold.3", "Slope.fold.3", "Intercpet.fold.4", "Slope.fold.4",
                                               "Intercpet.fold.5", "Slope.fold.5","trait")
          
          write.table(r.gy.m.RK ,paste( "Cor.betw.obs.pred.",ncycle,"_Cycle_",number.of.folds,".folds.BGLR.RK.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
          
          write.table(the.coefficients.m.RK, paste( "SLR.int.slope.Y.obs.X.pred",ncycle,"_Cycle_",number.of.folds,".folds.BGLR.RK.txt",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE,col.names = TRUE)
          
          the.selected.taxa.RK <- the.selected.taxa.RK[!the.selected.taxa.RK%in%parents.vec]
          the.unique.selected.taxa.RK <- unique(the.selected.taxa.RK)
          
          write.table(the.selected.taxa.RK, paste("the.sel.taxa.BGLR.RK",ncycle,"_Cycle",".txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)
          write.table(the.unique.selected.taxa.RK, paste("the.uniq.sel.taxa.BGLR.RK",ncycle,"_Cycle",".txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)
        
          
          #return(list(the.unique.selected.taxa=the.unique.selected.taxa, the.unique.selected.taxa.A=the.unique.selected.taxa.A, the.unique.selected.taxa.D=the.unique.selected.taxa.D, 
                      #the.unique.selected.taxa.E=the.unique.selected.taxa.E, the.unique.selected.taxa.BA=the.unique.selected.taxa.BA, the.unique.selected.taxa.BB=the.unique.selected.taxa.BB,
                      #the.unique.selected.taxa.BC=the.unique.selected.taxa.BC, the.unique.selected.taxa.BL=the.unique.selected.taxa.BL, the.unique.selected.taxa.BR=the.unique.selected.taxa.BR,
                      #the.unique.selected.taxa.RK=the.unique.selected.taxa.RK, the.unique.selected.taxa.SV=the.unique.selected.taxa.SV))
        
          
              }# Cycles per density
            # Marker density
        }#end run.k.fold.CV.using.rrBLUP
        
        

