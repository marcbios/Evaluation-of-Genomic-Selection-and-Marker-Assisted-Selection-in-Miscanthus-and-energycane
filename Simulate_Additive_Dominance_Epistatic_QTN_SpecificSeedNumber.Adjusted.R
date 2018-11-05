
#This function will create simulated data for a trait with genes that following a geometric series.
# The user will input the number of QTL, the heritabilities to be explored, the additive effect,
# and the output directory

create.simluated.data <- function(Additive.QTN.number = NULL, Dominance.QTN.number = NULL, Epistatic.QTN.number = NULL,heritabilities.vector = NULL, additive.effect = NULL, 
                                  dominance.effect = NULL, amount.of.partial.or.overdominance = 1, epistatic.effect = NULL, replicates = NULL, 
                                  data.in.hapmap.format = TRUE, file.G = NULL,  file.Ext.G = NULL, file.from = NULL, 
                                  file.to = NULL, file.fragment = NULL, output.dir = NULL, numeric.genotype.filename = NULL,
                                  impute.missing.data.with.hets = NULL){
  #Experimental code
  setwd(home.dir)
  
  if(data.in.hapmap.format){
    count.filenum <- 0
    for(filenum in file.from:file.to){
      myFRG=GAPIT.Fragment(file.path=NULL,file.from=NULL, file.to=NULL,file.total=1,file.G=file.G,
                           file.Ext.G=file.Ext.G,seed=123,SNP.effect="Add",SNP.impute="Middle",
                           genoFormat=NULL, file.GD=NULL, file.Ext.GD=NULL, file.GM=NULL, file.Ext.GM=NULL, file = filenum, file.fragment=file.fragment,
                           LD.chromosome=NULL,LD.location=NULL,LD.range=NULL, Create.indicator = FALSE, Major.allele.zero = FALSE)
      
      if(count.filenum == 0){
        all.FRGGs <- myFRG$G
      }else{
        all.FRGGs <- rbind(all.FRGGs, myFRG$G[-1,])
      }#End if(count.filenum == 0)
      count.filenum <- count.filenum+1
    }#end for(filenum in file.from:file.to) 
    hm=GAPIT.HapMap(G = all.FRGGs,SNP.effect="Add",SNP.impute="Major")
    #####################################
    #Obtain the mafs of all SNPs
    
    #Total number of lines
    ns <- nrow(hm$GD)
    
    #Sum of the allele scores for each SNP
    ss <- apply(hm$GD, 2, sum)
    
    #Combine two situations: one where the allele coded as "2" is major; one where "0" is coded as major.
    maf.matrix <- rbind((.5*ss/ns), (1-(0.5*ss/ns)))
    
    #Copy the minor allele frequencies for all SNPs
    maf <- apply(maf.matrix, 2, min)
    
    #Find out which SNPs have MAF < 0.05
    snps.below.0.05.maf <- which(maf < 0.05)
    
    # Remove these SNPs from hm$GD
    
    #hm.GD.without.snps.below.0.05.maf <- hm$GD[,-snps.below.0.05.maf]
    
    ###############################
    #temp code#
    hm.GD.without.snps.below.0.05.maf <- hm$GD
    
    genotypes <- data.frame(hm$GI[,1], rep(NA,nrow(hm$GI)),hm$GI[,2:3],rep(NA,nrow(hm$GI)),
                            t(hm.GD.without.snps.below.0.05.maf))
    
    colnames(genotypes) <- c("Snp", "allele", "chr", "pos", "cm", t(as.character(hm$GT)))
    
  }else{
    the.genotypic.data.directly.read.in <- read.table(this.numeric.genotype.filename, head = TRUE)
    the.genotypic.data.directly.read.in.transposed <- t(the.genotypic.data.directly.read.in[,-1])
    colnames(the.genotypic.data.directly.read.in.transposed) <- t(the.genotypic.data.directly.read.in[,1])
    if(impute.missing.data.with.hets){
     data.temp <-  the.genotypic.data.directly.read.in.transposed[,-c(1,2)]
     data.temp[which(is.na(data.temp))] <- 1
     the.genotypic.data.directly.read.in.transposed[,-c(1,2)] <- data.temp
    }#end if(impute.missing.data.with.hets)
    genotypes <- data.frame(rownames(the.genotypic.data.directly.read.in.transposed), rep(NA,nrow(the.genotypic.data.directly.read.in.transposed)), 
                            the.genotypic.data.directly.read.in.transposed[,1], the.genotypic.data.directly.read.in.transposed[,2],
                            rep(NA,nrow(the.genotypic.data.directly.read.in.transposed)), the.genotypic.data.directly.read.in.transposed[,-c(1,2)])
    
    colnames(genotypes)[1:5] <- c("Snp", "allele", "chr", "pos", "cm")
  }#end if(data.in.hapmap.format)
 
  #End experimental code
  
  #Create a working directory for the output results:
  dir.create(paste(home.dir,"/", output.dir, sep = ""))
  
  #Set the working directory
  setwd(paste(home.dir,"/", output.dir, sep = ""))
  
  #Randomly select (without replacement) k additive QTN, and assign an effect size
  #seed.number <- sample(-1000000:1000000, 1)
  #seed.number <- seednum.specific
  seed.number=seednum.specific.adt
  # Use the seed number to generate uniform random variables
  set.seed(seed.number)
  
  vector.of.add.QTN <- sample(1:nrow(genotypes), Additive.QTN.number, replace = FALSE)
  Add.QTN.genotypic.information <- genotypes[vector.of.add.QTN,]
  
  
  #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
  # number used
  write.table(paste("Here_is_the_seed_number:_",seed.number, sep = ""), paste("Seed.number.for.", Additive.QTN.number,"Add.QTN",
                                                                              ".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
  write.table(Add.QTN.genotypic.information, paste("Genotypic.information.for.", Additive.QTN.number,".Additive.QTN",
                                                   ".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
  
  
  
 
  #Randomly select (without replacement) k dominance QTN, and assign an effect size
  #seed.number <- sample(-1000000:1000000, 1)
  #seed.number <- seednum.specific
  seed.number=seednum.specific.dom
  # Use the seed number to generate uniform random variables
  set.seed(seed.number)
  
  vector.of.dom.QTN <- sample(1:nrow(genotypes), Dominance.QTN.number, replace = FALSE)
  Dom.QTN.genotypic.information <- genotypes[vector.of.dom.QTN,]
  
  
  #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
  # number used
  write.table(paste("Here_is_the_seed_number:_",seed.number, sep = ""), paste("Seed.number.for.", Dominance.QTN.number,"Add.QTN",
                                                                              ".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
  write.table(Dom.QTN.genotypic.information, paste("Genotypic.information.for.", Dominance.QTN.number,".Dominance.QTN",
                                                   ".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
  
  
  #Randomly select (without replacement) 2*k epistatic QTN, and assign an effect size
  #seed.number <- sample(-1000000:1000000, 1)
  #seed.number <- seednum.specific
  seed.number=seednum.specific.epi
  # Use the seed number to generate uniform random variables
  set.seed(seed.number)
  
  vector.of.epi.QTN <- sample(1:nrow(genotypes), (2*Epistatic.QTN.number), replace = FALSE)
  Epi.QTN.genotypic.information <- genotypes[vector.of.epi.QTN,]
  
  
  #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
  # number used
  write.table(paste("Here_is_the_seed_number:_",seed.number, sep = ""), paste("Seed.number.for.", Epistatic.QTN.number,"Epi.QTN",
                                                                              ".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
  write.table(Epi.QTN.genotypic.information, paste("Genotypic.information.for.", Epistatic.QTN.number,".Epistatic.QTN",
                                                   ".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
  
  
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
  
  #seednum.from <- 130001
  #seednum.to <- 130100
  #seednum.by <- 1
  #seednum.seq <- seq(from=seednum.from, to=seednum.to, by=seednum.by)
  #seednum.seq has been replaced by the seednum.seq from replicates.
  #Therefore each replicate is simulated by a particular seed number constant for coparison across polyrad and clare
  
  
  #For loop through the vector of heritabilities
  for(i in heritabilities.vector){
    #If heritability is zero
    if(i == 0){ 
      #Simulate m replicates of n N(0,b) random variables, where b = additive.genetic.variance
      the.seed.number.vector <- NULL
      for(j in 1:replicates){
        #seed.number <- sample(-1000000:1000000, 1)
        seed.number <- seednum.seq[j]
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
        #seed.number <- sample(-1000000:1000000, 1)
        seed.number <- seednum.seq[j]
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
      write.table(the.seed.number.vector, paste("Seed.number.for.", replicates,".Reps",
                                                ".Herit.",i,".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
      write.table(simulated.data, paste("Simulated.Data.", replicates,".Reps",
                                        ".Herit.",i,".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
      
    }   
  }#End for(i in heritabilities.vector)
 return( list(genotypes = genotypes, simulated.data = simulated.data))
  
}#end "create.simluated.data()"




