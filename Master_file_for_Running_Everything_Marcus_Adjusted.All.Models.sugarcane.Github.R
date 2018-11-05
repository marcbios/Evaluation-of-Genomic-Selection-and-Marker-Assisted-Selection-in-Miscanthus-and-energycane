rm(list=ls())
####################

setwd("C:/Users/omo/Desktop/PostDoc/Manuscript1/Sugarcane/NewAnalysis.09202018/")
home.dir <- getwd()
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("multtest")
library(sommer)
library(BGLR)
library(kernlab)
library(rrBLUP)
library('MASS')
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")


source("http://zzlab.net/GAPIT/previous/gapit_functions_20160408.txt")

source("http://zzlab.net/GAPIT/emma.txt")

# load functions
source("C:/Users/omo/Desktop/PostDoc/Manuscript1/Miscanthus/Functions/Simulate_Additive_Dominance_Epistatic_QTN_Marcus_Adjusted_05.08.2018.R")
source("C:/Users/omo/Desktop/PostDoc/Manuscript1/Miscanthus/Functions/CoincidenceIndex.R")
source("C:/Users/omo/Desktop/PostDoc/Manuscript1/Miscanthus/Functions/k_fold_CV_using_rrBLUP.m.s_plus_MultipleModels_Marcus_Adjusted.5GS.Models.09.18.2018.noParentsSelected.R")

source("C:/Users/omo/Desktop/PostDoc/Manuscript1/Miscanthus/Functions/Run_GWAS_using_GAPIT_on_All_Traits.R")
source("C:/Users/omo/Desktop/PostDoc/Manuscript1/Miscanthus/Functions/k-fold.CV.GWAS.in.Training.MAS.in.pred.09182018.noParentsSelected.CI.R")
source("C:/Users/omo/Desktop/PostDoc/Manuscript1/Miscanthus/Functions/haldane_function_AEL.R")
source("C:/Users/omo/Desktop/PostDoc/Manuscript1/Miscanthus/Functions/Simulate_Backcrossing_Marcus_Adjusted.04.19.2018.roundAllele.R")
source("C:/Users/omo/Desktop/PostDoc/Manuscript1/Miscanthus/Functions/Simulate_Traits_in_BC_Generation_20171214.R")
source("C:/Users/omo/Desktop/PostDoc/Manuscript1/Miscanthus/Functions/Fit_rrBLUP_plus_Multimodels_in_Current_Generation_Validate_in_BC_Generation_Marcus_Adjusted.5GS.Models.08.22.2018.R")
source("C:/Users/omo/Desktop/PostDoc/Manuscript1/Miscanthus/Functions/Do_GWAS_in_Current_Gen_Validated_MAS_in_BC_Gen.R")



###############
#User input

#Number of additive QTN (k)
this.Additive.QTN.number <- 1

#Number of dominance QTN (l)
this.Dominance.QTN.number <- 1

#Number of epistatic QTN (m)
this.Epistatic.QTN.number <- 2


#Vector of heritabilities to investigate
this.heritabilities.vector <- 0.3

#Size of the additive effect of the largest QTL (must be (-1,1) but preferably (0,1))

#psi <- (this.Additive.QTN.number - 1)/(this.Additive.QTN.number + 1)
#psi <- (this.Dominance.QTN.number - 1)/(this.Dominance.QTN.number + 1)
psi <- ((this.Epistatic.QTN.number*2) - 1)/((this.Epistatic.QTN.number*2) + 1)
this.additive.effect <- 0.0

#Size of the dominance effect of the largest QTL (must be (-1,1) but preferably (0,1))
this.dominance.effect <- 0.0
this.amount.of.partial.or.overdominance <- 1.0 #If there is no partial or overdominance, then set this to 1

#Consider partial dominance

#Size of the epistatic effect of the largest QTL (must be (-1,1) but preferably (0,1)|
this.epistatic.effect <- psi


#Generate seed numbers to use for simulating phenotypes
seednum.from <- 130001
seednum.to <- 130100
seednum.by <- 1
seednum.seq <- seq(from=seednum.from, to=seednum.to, by=seednum.by)
seednum.seq
length(seednum.seq)

seednum.specific=130164
seed.number=seednum.specific

seednum.specific.adt=130164
seed.number=seednum.specific.adt

seednum.specific.dom=140160
seed.number=seednum.specific.dom

seednum.specific.epi=150160
seed.number=seednum.specific.epi

burn=1000
Iter=2500

parents.vec <- c("CP95", "CP88") # this is to ensure that parents are not selected by MAS and GS in the F1
#Number of replicates of simulated phenotypes for each heritability (m)
#this.replicates <- 20 # change replicates to the length of seed numbers
this.replicates=length(seednum.seq)


#IMPORTANT!!!!!!! Pelase make sure that your genotypic data are sorted by chromosomes and by ascending cM positions within
# each chromosome.
these.data.in.hapmap.format <- FALSE

this.numeric.genotype.filename <- "C:/Users/omo/Desktop/PostDoc/Manuscript1/Sugarcane/data/Scane.GenomicData.175obs.4606markers.txt"


this.file.from = 1 # what are these values for? Chromosome or what?
this.file.to = 1
this.file.fragment = 1000000

this.impute.missing.data.with.hets <- TRUE

model_names <- c("rrBLUP", "sommerE", "BA", "RK", "SV")

#Output directory
this.output.dir <- paste(this.Additive.QTN.number,"_Add_QTN_", this.Dominance.QTN.number,"_Dom_QTN_",this.Epistatic.QTN.number,"_Epi_QTN_h.2_",
                    this.heritabilities.vector,"_add.eff_", this.additive.effect,"", "_dom.eff_", this.dominance.effect, "_epi.eff_", 
                    this.epistatic.effect,"_reps_", this.replicates, sep = "")

#this.output.dir <- "/Users/omo/Google Drive/PostDoc_Scripts/Results"
################
#Create the simulated data
these.pheno.and.geno.data <- create.simluated.data(Additive.QTN.number = this.Additive.QTN.number, Dominance.QTN.number = this.Dominance.QTN.number, Epistatic.QTN.number = this.Epistatic.QTN.number,
                      heritabilities.vector = this.heritabilities.vector, additive.effect = this.additive.effect, dominance.effect = this.dominance.effect, 
                      amount.of.partial.or.overdominance = this.amount.of.partial.or.overdominance, epistatic.effect = this.epistatic.effect, 
                      replicates = this.replicates, data.in.hapmap.format = these.data.in.hapmap.format, file.G = this.file.G, file.Ext.G = this.file.Ext.G, 
                      file.from = this.file.from, file.to = this.file.to, 
                      file.fragment = this.file.fragment, output.dir = this.output.dir, numeric.genotype.filename = this.numeric.genotype.filename,
                      impute.missing.data.with.hets = this.impute.missing.data.with.hets)



setwd(home.dir)

#seed.number <- sample(-1000000:1000000,1, replace = FALSE)
seednum.specific=130164
seed.number=seednum.specific



#these.pheno.and.geno.data$genotypes # these.pheno.and.geno.data$simulated.data
#Run k-fold cross validation on all 100 traits. The same folds are being used for each trait. #these.pheno.and.geno.data$genotypes #these.pheno.and.geno.data$genotypes #these.pheno.and.geno.data$genotypes
#dev.off()
the.taxa.selected.from.GS <- run.k.fold.CV.using.rrBLUP(the.phenotypes = these.pheno.and.geno.data$simulated.data, the.genotypes = these.pheno.and.geno.data$genotypes,
                           the.output.dir = this.output.dir, number.of.folds = 5, proportion.to.select = 0.10, user.input.seed.number = TRUE, seed.number = seed.number)
save(the.taxa.selected.from.GS, file="the.taxa.selected.from.GS.rda")
save.image("F1.GS.rda")
setwd(home.dir)


#MAS-based prediction; in each training set of each trait; run GWAS; pick the top markers, fit a model with only the top markers in the training set; 
# use this fitted model to predict which lines have optimal phenotypic performance
the.taxa.selected.from.MAS <- identify.top.performing.lines.via.MAS(the.phenotypes = these.pheno.and.geno.data$simulated.data, the.genotypes = these.pheno.and.geno.data$genotypes, 
                                      the.output.dir = this.output.dir, number.of.folds = 5, proportion.to.select = 0.10, number.of.markers.to.use.for.MAS = 4,
                                                  user.input.seed.number = TRUE, seed.number = seed.number)
save.image("F1.MAS.rda")
setwd(home.dir)
save.image()
#Backcross the individuals selected from GS to one of the two parents

#geoBC.GS <- list(length(unique(names(the.taxa.selected.from.GS))))
#names(geoBC.GS) <- names(the.taxa.selected.from.GS)

for(z in 1:length(unique(names(the.taxa.selected.from.GS)))){
  the.taxa.selected.from.GS_z <- the.taxa.selected.from.GS[z]
  the.taxa.selected.from.GS_z <- the.taxa.selected.from.GS_z[[1]]
  the.taxa.selected.from.GS_z <- as.vector(the.taxa.selected.from.GS_z)
  #genotypes.combine.genotypes
  the.taxa.selected.from.GS_z <- the.taxa.selected.from.GS_z[!the.taxa.selected.from.GS_z%in%parents.vec]
  the.backcross.population.inds.selected.from.GS <- create.backcross.population(the.genotypic.data = these.pheno.and.geno.data$genotypes, the.individuals.of.interest = the.taxa.selected.from.GS_z, 
                                                                              the.parent = "CP95") #Kaskade.Justin these.pheno.and.geno.data$genotypes
  
  #geoBC.GS[z] <- the.backcross.population.inds.selected.from.GS
  assign(paste("the.backcross.population.inds.selected.from.GS", model_names[z],sep="."), the.backcross.population.inds.selected.from.GS)
  
  setwd(home.dir)
}


#Backcross the individuals selected from MAS to one of the two parents
the.taxa.selected.from.MAS <- the.taxa.selected.from.MAS[!the.taxa.selected.from.MAS%in%parents.vec]
the.backcross.population.inds.selected.from.MAS <- create.backcross.population(the.genotypic.data = these.pheno.and.geno.data$genotypes, the.individuals.of.interest = the.taxa.selected.from.MAS, 
                                                                              the.parent = "CP95") #Kaskade.Justin

#Simulate the same phenotype with the exact same genetic architecture in the backcross population from the GS selected individuals
for(zz in 1:length(model_names)){

  the.backcross.population.inds.selected.from.GS <- get(paste("the.backcross.population.inds.selected.from.GS", model_names[zz],sep="."))
  #select only Goldengate markers for simulations
  
phenotypes.in.GS.backcross.population <- create.simluated.data.in.backcross.generation(Additive.QTN.number = this.Additive.QTN.number, Dominance.QTN.number = this.Dominance.QTN.number, 
                                              Epistatic.QTN.number = this.Epistatic.QTN.number,
                                              heritabilities.vector = this.heritabilities.vector, additive.effect = this.additive.effect, dominance.effect = this.dominance.effect, 
                                              amount.of.partial.or.overdominance = this.amount.of.partial.or.overdominance, epistatic.effect = this.epistatic.effect, replicates = 1,
                                              genotypes = the.backcross.population.inds.selected.from.GS, output.dir = this.output.dir, GS.or.MAS = "GS")
                                                          #new gg +POLYRAD COMBINED GENOMIC DATA
assign(paste("phenotypes.in.GS.backcross.population", model_names[zz],sep="."), phenotypes.in.GS.backcross.population)
phenotypes.in.GS.backcross.population=NULL

setwd(home.dir)
}

#Simulate the same phenotype with the exact same genetic architecture in the backcross population from the MAS selected individuals
# Select only Goldengate markers for phenotype simulations

phenotypes.in.MAS.backcross.population <- create.simluated.data.in.backcross.generation(Additive.QTN.number = this.Additive.QTN.number, Dominance.QTN.number = this.Dominance.QTN.number, 
                                               Epistatic.QTN.number = this.Epistatic.QTN.number,
                                               heritabilities.vector = this.heritabilities.vector, additive.effect = this.additive.effect, dominance.effect = this.dominance.effect, 
                                               amount.of.partial.or.overdominance = this.amount.of.partial.or.overdominance, epistatic.effect = this.epistatic.effect, replicates = 1,
                                               genotypes = the.backcross.population.inds.selected.from.MAS, output.dir = this.output.dir, GS.or.MAS = "MAS")

setwd(home.dir)


for(zzz in 1:length(model_names)){
  

phenotypes.in.GS.backcross.population_zzz <- get(paste("phenotypes.in.GS.backcross.population", model_names[zzz],sep="."))


the.backcross.population.inds.selected.from.GS_zzz <- get(paste("the.backcross.population.inds.selected.from.GS", model_names[zzz],sep="."))

#Fit rrBLUP model in the current generation and validate in the GS selected backcross population 
performance.of.GS.in.GS.BC1 <- rrBLUP.current.generation.validate.in.BC.generation(the.phenotypes.current.gen = these.pheno.and.geno.data$simulated.data, the.genotypes.current.gen = these.pheno.and.geno.data$genotypes,
                                                    the.genotypes.BC = the.backcross.population.inds.selected.from.GS_zzz, the.phenotypes.BC = phenotypes.in.GS.backcross.population_zzz,
                                                    the.output.dir = this.output.dir, GS.or.MAS = "GS")
assign(paste("performance.of.GS.in.GS.BC1", model_names[zzz],sep="."), performance.of.GS.in.GS.BC1)

save(performance.of.GS.in.GS.BC1, file=paste("performance.of.GS.in.GS.BC1",model_names[zzz], "rda",sep="."))
performance.of.GS.in.GS.BC1=NULL

setwd(home.dir)
}



#Fit GWAS model in the current generation and validate in the MAS selected backcross population
performance.of.MAS.in.MAS.BC1 <- GWAS.current.generation.validate.MAS.in.BC.generation(the.phenotypes.current.gen = these.pheno.and.geno.data$simulated.data,
                                                                                      the.genotypes.current.gen = these.pheno.and.geno.data$genotypes,
                                                                                      the.genotypes.BC = the.backcross.population.inds.selected.from.MAS,
                                                                                      the.phenotypes.BC = phenotypes.in.MAS.backcross.population,
                                                                                      the.output.dir = this.output.dir, GS.or.MAS = "MAS",
                                                                                      number.of.markers.to.use.for.MAS = 4)


setwd(home.dir)









