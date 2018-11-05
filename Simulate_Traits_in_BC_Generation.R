


#This function will create simulated data for a trait with genes that following a geometric series.
# The user will input the number of QTL, the heritabilities to be explored, the additive effect,
# and the output directory

create.simluated.data.in.backcross.generation <- function(Additive.QTN.number = NULL, Dominance.QTN.number = NULL, Epistatic.QTN.number = NULL,
                                                          heritabilities.vector = NULL, additive.effect = NULL, dominance.effect = NULL, 
                                                          amount.of.partial.or.overdominance = NULL, epistatic.effect = NULL, replicates = NULL,
                                                          genotypes = NULL, output.dir = NULL, GS.or.MAS = NULL){
 
  
  #Set the working directory
  setwd(paste(home.dir,"/", output.dir, sep = ""))
  
  
  Add.QTN.genotypic.information.from.parental.generation <- read.table(paste("Genotypic.information.for.", Additive.QTN.number,".Additive.QTN",
                                                    ".txt", sep = ""), head = TRUE)
  
  index.of.additive.markers.in.order <- NULL
  for(i in 1:nrow(Add.QTN.genotypic.information.from.parental.generation)){
    index.of.additive.markers.in.order <- c(index.of.additive.markers.in.order, which(as.character(genotypes[,1]) == as.character(Add.QTN.genotypic.information.from.parental.generation[i,1])))
  }
  Add.QTN.genotypic.information <- genotypes[index.of.additive.markers.in.order,]
  
  
  Dom.QTN.genotypic.information.from.parental.generation <- read.table(paste("Genotypic.information.for.", Dominance.QTN.number,".Dominance.QTN",
                                                    ".txt", sep = ""), head = TRUE)
  
  index.of.dominance.markers.in.order <- NULL
  for(i in 1:nrow(Dom.QTN.genotypic.information.from.parental.generation)){
    index.of.dominance.markers.in.order <- c(index.of.dominance.markers.in.order, which(as.character(genotypes[,1]) == as.character(Dom.QTN.genotypic.information.from.parental.generation[i,1])))
  }
  Dom.QTN.genotypic.information <- genotypes[index.of.dominance.markers.in.order,]
  
  
  
  Epi.QTN.genotypic.information.from.parental.generation <- read.table(paste("Genotypic.information.for.", Epistatic.QTN.number,".Epistatic.QTN",
                                                    ".txt", sep = ""), head = TRUE)
  
  index.of.epistatic.markers.in.order <- NULL
  for(i in 1:nrow(Epi.QTN.genotypic.information.from.parental.generation)){
    index.of.epistatic.markers.in.order <- c(index.of.epistatic.markers.in.order, which(as.character(genotypes[,1]) == as.character(Epi.QTN.genotypic.information.from.parental.generation[i,1])))
  }
  Epi.QTN.genotypic.information <- genotypes[index.of.epistatic.markers.in.order,]
  
  
  
  #Create a "base line" trait, which is basically just the additive effects; this is what we would see if 
  # the heritability of the simulated trait were 1
  additive.effect.trait.object <- t(Add.QTN.genotypic.information[,-c(1:5)]) #this was originally the base.line.trait.object
  
  dominance.effect.trait.object.temp <- t(Dom.QTN.genotypic.information[,-c(1:5)])
  #Convert the numeric genotypes from additrive to 1 if heterozygous and 0 if homozygous
  dominance.effect.trait.object <- NULL
  for(i in 1:ncol(dominance.effect.trait.object.temp)){
    this.add.to.dom.coded.vector <- dominance.effect.trait.object.temp[,i]
    this.add.to.dom.coded.vector[which(this.add.to.dom.coded.vector == 2)] <- 0
    this.add.to.dom.coded.vector[which(this.add.to.dom.coded.vector == 1)] <- amount.of.partial.or.overdominance
    dominance.effect.trait.object <- cbind(dominance.effect.trait.object, this.add.to.dom.coded.vector)
  }#end for(i in 1:ncol(dominance.effect.trait.object.temp))
  colnames(dominance.effect.trait.object) <- colnames(dominance.effect.trait.object.temp)
  
  
  epistatic.effect.trait.object <-t(Epi.QTN.genotypic.information[,-c(1:5)])
  #AEL Changed: - epistatic.effect.trait.object<- epistatic.effect.trait.object[,number.of.epistasis]
  
  colnames(additive.effect.trait.object) <- paste("Chr_",Add.QTN.genotypic.information[,3], "_",Add.QTN.genotypic.information[,4],sep = "")
  
  colnames(dominance.effect.trait.object) <- paste("Chr_",Dom.QTN.genotypic.information[,3], "_",Dom.QTN.genotypic.information[,4],sep = "")
  
  colnames(epistatic.effect.trait.object) <- paste("Chr_",Epi.QTN.genotypic.information[,3], "_",Epi.QTN.genotypic.information[,4],sep = "")
  
  #make base.line.trait additive.component and epistatic.component
  additive.component<- as.data.frame(matrix(0, nrow = nrow(additive.effect.trait.object), ncol = 1))
  dominance.component <- as.data.frame(matrix(0, nrow = nrow(dominance.effect.trait.object), ncol = 1))
  epistatic.component<- as.data.frame(matrix(0, nrow = nrow(epistatic.effect.trait.object), ncol = 1))
  #base.line.trait <- as.data.frame(matrix(0, nrow = nrow(base.line.trait.object), ncol = 1)) 
  #Simulate each of the the additive, dominance, and epistatic effects. For example, the effect size of the
  # j^th largest additive QTN is "additive.effect^i"
  for(i in 1:Additive.QTN.number) additive.component <- additive.component + (additive.effect.trait.object[,i]*(additive.effect^i))
  rownames(additive.component) <- rownames(additive.effect.trait.object)
  colnames(additive.component) <- "Additive.effect"
  additive.genetic.variance <- var(additive.component, na.rm = TRUE)
  
  for(i in 1:Dominance.QTN.number) dominance.component <- dominance.component + (dominance.effect.trait.object[,i]*(dominance.effect^i))
  rownames(dominance.component) <- rownames(dominance.effect.trait.object)
  colnames(dominance.component) <- "Dominance.effect"
  dominance.genetic.variance <- var(dominance.component, na.rm = TRUE)
  
  last.number.of.this.loop <- Epistatic.QTN.number - 1
  for(i in 0:last.number.of.this.loop) epistatic.component <- epistatic.component + ((epistatic.effect.trait.object[,((2*i)+1)]*epistatic.effect.trait.object[,((2*i)+2)])*(epistatic.effect^(i+1)))
  rownames(epistatic.component) <- rownames(epistatic.effect.trait.object)
  colnames(epistatic.component) <- "Epistatic.effect"
  epistatic.genetic.variance<- var(epistatic.component, na.rm = TRUE)
  
  #Set the row names of the base.line.trait object to the new names
  base.line.trait <- additive.component+dominance.component+epistatic.component
  #base.line.trait.with.new.taxa <- merge(base.line.trait, taxa.name.converter, by.x = "row.names", 
  #                                       by.y = "Old_Taxa_ID")
  
  #the.new.taxa.ids <- as.character(base.line.trait.with.new.taxa[,2])
  #base.line.trait <- as.matrix(base.line.trait.with.new.taxa[,2], nrow = nrow(base.line.trait.with.new.taxa))
 # rownames(base.line.trait) <- as.character(base.line.trait[,3])
  
  
  #For loop through the vector of heritabilities
for(i in heritabilities.vector){
    #If heritability is zero
    if(i == 0){ 
      #Simulate m replicates of n N(0,b) random variables, where b = additive.genetic.variance
      the.seed.number.vector <- NULL
      for(j in 1:replicates){
        seed.number <- sample(-1000000:1000000, 1)
        set.seed(seed.number)
        the.normal.random.variables <- rnorm(nrow(base.line.trait), mean = 0, sd = 1)
        if(j == 1){
          simulated.data <- the.normal.random.variables
        }else{
          simulated.data <- cbind(simulated.data, the.normal.random.variables)
          colnames(simulated.data)[j] <- paste(colnames(simulated.data)[j],".",j,sep = "")
        }
        the.seed.number.vector <- c(the.seed.number.vector, seed.number)
      }
      
      #Format the output file for the simulated phenotypes
      simulated.data <- cbind(rownames(base.line.trait),simulated.data)
      colnames(simulated.data)[1] <- "<Trait>"
      colnames(simulated.data)[2] <- "the.normal.random.variables.1"
      
      
      #Output the m replicates and the seed numbers, formatted for TASSEL
      write.table(the.seed.number.vector, paste("Seed.number.for.", replicates,".Reps",
                                                ".Herit.",i,".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
      write.table(simulated.data, paste("Simulated.Data.", replicates,".Reps",
                                        ".Herit.",i,".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
      
      
    }else{
      #Calcualte V_e, the residual variance
      #residual.variance <- (additive.genetic.variance*(1-i))/i
      residual.variance <-( ((additive.genetic.variance + dominance.genetic.variance + epistatic.genetic.variance)/i) - additive.genetic.variance - dominance.genetic.variance - epistatic.genetic.variance)
      #Use V_e to generate n replicates of N(0, Ve) random variables
      the.seed.number.vector <- NULL
      col.name.vector <- NULL
      for(j in 1:replicates){
        seed.number <- sample(-1000000:1000000, 1)
        set.seed(seed.number)
        the.normal.random.variables <- rnorm(nrow(base.line.trait), mean = 0, sd = sqrt(residual.variance))
        the.base.line.trait.plus.normal.random.variables <- base.line.trait+the.normal.random.variables
        if(j == 1){
          simulated.data <- the.base.line.trait.plus.normal.random.variables
        }else{
          simulated.data <- cbind(simulated.data, the.base.line.trait.plus.normal.random.variables)
          colnames(simulated.data)[j] <- paste(colnames(simulated.data)[j],".",j,sep = "")
        }
        the.seed.number.vector <- c(the.seed.number.vector, seed.number)
        col.name.vector <- c(col.name.vector, paste("Heritability_",i, "_Rep_", j, sep = ""))
      }
      
      colnames(simulated.data)  <- col.name.vector
      
      #Format the output file for the simulated phenotypes
      simulated.data <- cbind(rownames(base.line.trait),simulated.data)
      colnames(simulated.data)[1] <- "<Trait>"
      
      #Output the m replicates and the seed numbers, formatted for TASSEL
      write.table(the.seed.number.vector, paste("Seed.number.for.",GS.or.MAS,".BC1.Generation.", replicates,".Reps",
                                                ".Herit.",i,".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
      write.table(simulated.data, paste("Simulated.Data.in.",GS.or.MAS,".BC1.Generation.", replicates,".Reps",
                                        ".Herit.",i,".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
      
    }   
  }#End for(i in heritabilities.vector)
 return(simulated.data)
  
}#end "create.simluated.data.in.backcross.generation()"




