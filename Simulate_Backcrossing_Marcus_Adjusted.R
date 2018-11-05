
create.backcross.population <- function(the.genotypic.data = NULL, the.individuals.of.interest = NULL, the.parent = NULL){
 
  #nam_ind <- unique(names(the.individuals.of.interest[[1]]))
  #the.individuals.of.interest <- as.vector(the.individuals.of.interest[[1]])
  
   the.genotypic.data.for.making.crosses <- the.genotypic.data[,c(1:5, which(colnames(the.genotypic.data) == the.parent),
                                                                 which(colnames(the.genotypic.data) %in% the.individuals.of.interest) )  ]
  
  #For loop through the chromosomes
  for(chr in unique(the.genotypic.data.for.making.crosses$chr)){
    the.genotypic.data.for.making.crosses.this.chr <- the.genotypic.data.for.making.crosses[which(the.genotypic.data.for.making.crosses$chr == chr),]
    the.genotypic.data.progeny.this.chr <- the.genotypic.data.for.making.crosses.this.chr[,1:5]
    #For individual i
    for(i in 7:ncol(the.genotypic.data.for.making.crosses.this.chr)){
      #Subdivide individual i and the parent into two strands
      ind.i.strand.1 <- NULL
      ind.i.strand.2 <- NULL
      parent.strand.1 <- NULL
      parent.strand.2 <- NULL
      for(genotype in 1:nrow(the.genotypic.data.for.making.crosses.this.chr)){
        if(the.genotypic.data.for.making.crosses.this.chr[genotype,i] == 2){
          ind.i.strand.1 <- c(ind.i.strand.1, "A")
          ind.i.strand.2 <- c(ind.i.strand.2, "A")
        }
        if(the.genotypic.data.for.making.crosses.this.chr[genotype,i] == 1){
          ind.i.het.unif.rv <- runif(1, min = 0, max = 1) 
          if(ind.i.het.unif.rv <= 0.5){
            ind.i.strand.1 <- c(ind.i.strand.1, "a")
            ind.i.strand.2 <- c(ind.i.strand.2, "A")
          }else{
            ind.i.strand.1 <- c(ind.i.strand.1, "A")
            ind.i.strand.2 <- c(ind.i.strand.2, "a")
          }#if(ind.i.het.unif.rv <= 0.5)
        }
        if(the.genotypic.data.for.making.crosses.this.chr[genotype,i] == 0){
          ind.i.strand.1 <- c(ind.i.strand.1, "a")
          ind.i.strand.2 <- c(ind.i.strand.2, "a")
        }
        
        
        if(the.genotypic.data.for.making.crosses.this.chr[genotype,6] == 2){
          parent.strand.1 <- c(parent.strand.1, "A")
          parent.strand.2 <- c(parent.strand.2, "A")
        }
        if(the.genotypic.data.for.making.crosses.this.chr[genotype,6] == 1){
          parent.strand.het.unif.rv <- runif(1, min = 0, max = 1) 
          if(parent.strand.het.unif.rv <= 0.5){
            parent.strand.1 <- c(parent.strand.1, "a")
            parent.strand.2 <- c(parent.strand.2, "A")
          }else{
            parent.strand.1 <- c(parent.strand.1, "A")
            parent.strand.2 <- c(parent.strand.2, "a")
          }#if(ind.i.het.unif.rv <= 0.5)
        }
        if(the.genotypic.data.for.making.crosses.this.chr[genotype,6] == 0){
          parent.strand.1 <- c(parent.strand.1, "a")
          parent.strand.2 <- c(parent.strand.2, "a")
        }
        
      }#end for(genotype in 1:nrow(the.genotypic.data.for.making.crosses.this.chr))
       
  
      #For the i'th individual and for the parent, randomly select one of the two strands
      ind.i.strand.uniform.RV <- runif(1, min = 0, max = 1) 
      if(ind.i.strand.uniform.RV <= 0.5){
        ind.i.strand.in.progeny <- ind.i.strand.1
        the.other.ind.i.strand <- ind.i.strand.2
      }else{
        ind.i.strand.in.progeny <- ind.i.strand.2
        the.other.ind.i.strand <- ind.i.strand.1
      }#end if(ind.i.strand.uniform.RV <= 0.5)
      
      parent.strand.uniform.RV <- runif(1, min = 0, max = 1) 
      if(parent.strand.uniform.RV <= 0.5){
        parent.strand.in.progeny <- parent.strand.1
        the.other.parent.strand <- parent.strand.2
      }else{
        parent.strand.in.progeny <- parent.strand.2
        the.other.parent.strand <- parent.strand.1
      }# end if(parent.strand.uniform.RV <= 0.5)
      
      
      #For the j'th interval between two markers
      for(j in 2:nrow(the.genotypic.data.for.making.crosses.this.chr)){
        #Use Haldane's mapping function to calculate the probability of an odd number of crossover events occurring
        c <- haldane(the.genotypic.data.for.making.crosses.this.chr, 4,(j-1), j)
        #Simulate whether or not  a crossover event occurs for the strand from the i'th individual
        ind.i.unif.rv <- runif(1, min = 0, max = 1) 
        if(ind.i.unif.rv < c){
          #Change the observed strand into the other strand
          ind.i.strand.in.progeny.temp <- c(ind.i.strand.in.progeny[1:(j-1)], the.other.ind.i.strand[j:nrow(the.genotypic.data.for.making.crosses.this.chr)])
          
          #Change the other strand into the observed strand
          the.other.ind.i.strand.temp <- c(the.other.ind.i.strand[1:(j-1)], ind.i.strand.in.progeny[j:nrow(the.genotypic.data.for.making.crosses.this.chr)])
         
          #Reset the names of all of the strands
          ind.i.strand.in.progeny <- ind.i.strand.in.progeny.temp
          the.other.ind.i.strand <- the.other.ind.i.strand.temp
  
        }# end if(this.unif.RV < c)
  
        #Simulate whether or not  a crossover event occurs for the strand from the parent
        parent.unif.rv <- runif(1, min = 0, max = 1) 
        if(parent.unif.rv < c){
          #Change the observed strand into the other strand
          parent.strand.in.progeny.temp <- c(parent.strand.in.progeny[1:(j-1)], the.other.parent.strand[j:nrow(the.genotypic.data.for.making.crosses.this.chr)])
          
          #Change the other strand into the observed strand
          the.other.parent.strand.temp <- c(the.other.parent.strand[1:(j-1)], parent.strand.in.progeny[j:nrow(the.genotypic.data.for.making.crosses.this.chr)])
          
          #Reset the names of all of the strands
          parent.strand.in.progeny <- parent.strand.in.progeny.temp
          the.other.parent.strand <- the.other.parent.strand.temp
          
        }# end if(this.unif.RV < c)
        
        
      }#End for the j'th interval between two markers for(j in 2:nrow(the.genotypic.data.for.making.crosses.this.chr))
      
      #Append this genotypic information into a genotype file that is formatted in the same way as the genotype file
      the.genotypic.information.this.progeny <- NULL
      for(marker.index in 1:length(parent.strand.in.progeny)){
        allele.from.ind.i <- ind.i.strand.in.progeny[marker.index]
        allele.from.parent <- parent.strand.in.progeny[marker.index]
        if((allele.from.ind.i == "A")  & (allele.from.parent == "A")  ){
          the.genotypic.information.this.progeny <- c(the.genotypic.information.this.progeny, 2) 
        }
        if( ((allele.from.ind.i == "A")  & (allele.from.parent == "a")) | ((allele.from.ind.i == "a")  & (allele.from.parent == "A"))  ){
          the.genotypic.information.this.progeny <- c(the.genotypic.information.this.progeny, 1) 
        }
        if((allele.from.ind.i == "a")  & (allele.from.parent == "a")  ){
          the.genotypic.information.this.progeny <- c(the.genotypic.information.this.progeny, 0) 
        }
        
        the.name.of.this.progeny <- paste(colnames(the.genotypic.data.for.making.crosses.this.chr)[i], "_X_", the.parent)
        
      }#end for(marker.index in 1:nrow(parent.strand.in.progeny))
      
      the.genotypic.data.progeny.this.chr <- data.frame(the.genotypic.data.progeny.this.chr, the.genotypic.information.this.progeny)
      colnames(the.genotypic.data.progeny.this.chr)[ncol(the.genotypic.data.progeny.this.chr)] <- the.name.of.this.progeny
  
    } #End for individual i (  for(i in 7:nrow(the.genotypic.data.for.making.crosses.this.chr))  )
    #End product: genotypic information on the backcrossed individual; this should be appended to a data set
    if(chr ==1){
      the.genotypic.data.progeny <- the.genotypic.data.progeny.this.chr
    }else{
      the.genotypic.data.progeny <- rbind(the.genotypic.data.progeny, the.genotypic.data.progeny.this.chr)
    }#end if(chr == 1)
  
    print(paste("------------------------------------ Just finished simulating progeny in Chromosome ", chr, "!!!!!!!!!!!--------------------------------------------------", sep = ""))
  }#End for(chr in unique(the.genotypic.data.for.making.crosses$chr))
  return(the.genotypic.data.progeny)
}#end create.backcross.population









