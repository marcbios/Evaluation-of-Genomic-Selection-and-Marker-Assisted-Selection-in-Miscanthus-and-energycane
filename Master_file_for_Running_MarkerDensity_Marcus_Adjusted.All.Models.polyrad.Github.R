rm(list=ls())
####################

setwd("C:/Users/omo/Desktop/PostDoc/Manuscript1/Marker.Density/Polyrad.Discret")
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
#source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/Simulate_Additive_Dominance_Epistatic_QTN_Marcus_Adjusted_03.23.2018.R")
source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/Simulate_Additive_Dominance_Epistatic_QTN_Marcus_Adjusted_05.08.2018.R")
#source("C:/Users/olatoye-temp/Desktop/Marcus/SourceCodesForPipeline/Simulate_Additive_Dominance_Epistatic_With_EffectDirection_QTN_20180413.R")
#source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/k_fold_CV_using_rrBLUP_plus_MultipleModels_Marcus_Adjusted.All.Models.03.23.2018.R")
source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/k_fold_CV_MarkerDensity_plus_MultipleModels_Marcus_Adjusted.All.Models.07.11.2018.noParentsSelected.R")

source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/Run_GWAS_using_GAPIT_on_All_Traits.R")
#source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/k-fold.CV.GWAS.in.Training.MAS.in.pred.R")
source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/k-fold.CV.GWAS.in.Training.MAS.in.pred.04192018.noParentsSelected.R")
source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/haldane_function_AEL.R")
#source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/Simulate_Backcrossing_Marcus_Adjusted.03.23.2018.R")
source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/Simulate_Backcrossing_Marcus_Adjusted.04.19.2018.roundAllele.R")
source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/Simulate_Traits_in_BC_Generation_20171214.R")
source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/Fit_rrBLUP_plus_Multimodels_in_Current_Generation_Validate_in_BC_Generation_Marcus_Adjusted.03.13.2018.R")
source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/Do_GWAS_in_Current_Gen_Validated_MAS_in_BC_Gen.R")
#library(rrBLUP)

gg <- read.delim("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/GoldenGateSNPs_85obs_241markers_edited_03292018.txt", header=T)
dim(gg)
ggnames <- colnames(gg[,-1])
###############
#User input

#Proportion of Marker density
geno.prop=0.25
#Percent
number.of.cycles=1

seed.numer.cycle <- seq(from=101, to=200, by=20)

#Number of additive QTN (k)
this.Additive.QTN.number <- 1

#Number of dominance QTN (l)
this.Dominance.QTN.number <- 4

#Number of epistatic QTN (m)
this.Epistatic.QTN.number <- 1


#Vector of heritabilities to investigate
this.heritabilities.vector <- 0.3

#Size of the additive effect of the largest QTL (must be (-1,1) but preferably (0,1))
#psi <- (this.Additive.QTN.number - 1)/(this.Additive.QTN.number + 1)
psi <- (this.Dominance.QTN.number - 1)/(this.Dominance.QTN.number + 1)
#psi <- ((this.Epistatic.QTN.number*2) - 1)/((this.Epistatic.QTN.number*2) + 1)

this.additive.effect <- 0.0

#Size of the dominance effect of the largest QTL (must be (-1,1) but preferably (0,1))
this.dominance.effect <- psi
this.amount.of.partial.or.overdominance <- 1.0 #If there is no partial or overdominance, then set this to 1

#Consider partial dominance

#Size of the epistatic effect of the largest QTL (must be (-1,1) but preferably (0,1)|
this.epistatic.effect <- 0.0


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

parents.vec <- c("Kaskade.Justin", "Zebrinus.Justin") # this is to ensure that parents are not selected by MAS and GS in the F1
#Number of replicates of simulated phenotypes for each heritability (m)
#this.replicates <- 20 # change replicates to the length of seed numbers
this.replicates=length(seednum.seq)


#IMPORTANT!!!!!!! Pelase make sure that your genotypic data are sorted by chromosomes and by ascending cM positions within
# each chromosome.
these.data.in.hapmap.format <- FALSE

this.numeric.genotype.filename <- "http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/Functions/GoldenGateSNPs_85obs_241markers_edited_03292018.txt"
#this.numeric.genotype.filename <- "~/Desktop/PostDoc/Manuscript1/Simulation.Analysis/clare5/Clare05.GenomicData.GG.85obs.3044snps.txt"
#this.numeric.genotype.filename <- "/Users/omo/Desktop/PostDoc/Manuscript1/Simulation.Analysis/Polyrad.noLinkage/PolyRad.noLinkage.GG.85obs.3044snps.txt"

#backcross.genotype.gg.poly <- read.delim("GoldenGate.Polyrad1.Use.txt", header=T)

#genotypes.impute <- read.delim("/Users/omo/Desktop/PostDoc/Manuscript1/Data/Miscanthus/Polyrad/"180510ClareMap1Genotypes_polyRAD_NoLinkage_continuous.csv", header=T)
genotypes.impute <- read.delim("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/ImputationStrategy/polyRAD/data/PolyRad.nolinkage.discreet.GG.85obs.3044snps.txt", header=T)
#genotypes.impute <- read.delim("/Users/omo/Desktop/PostDoc/Manuscript1/Data/Miscanthus/Polyrad/PolyRad.nolinkage.continuous.GG.85obs.3044snps.txt", header=T)
#genotypes.impute <- read.delim("~/Desktop/PostDoc/Manuscript1/Simulation.Analysis/clare5/Clare05.GenomicData.GG.85obs.3044snps.txt", header=T)
genotypes.impute.directly.read.in <- genotypes.impute
genotypes.impute.directly.read.in.transposed <- t(genotypes.impute.directly.read.in[,-1])
colnames(genotypes.impute.directly.read.in.transposed) <- t(genotypes.impute.directly.read.in[,1])

genotypes.impute.temp <-  genotypes.impute.directly.read.in.transposed[,-c(1,2)]
genotypes.impute.temp[which(is.na(genotypes.impute.temp))] <- 1
genotypes.impute.directly.read.in.transposed[,-c(1,2)] <- genotypes.impute.temp

genotypes.impute.genotypes <- data.frame(rownames(genotypes.impute.directly.read.in.transposed), rep(NA,nrow(genotypes.impute.directly.read.in.transposed)), 
                                         genotypes.impute.directly.read.in.transposed[,1], genotypes.impute.directly.read.in.transposed[,2],
                        rep(NA,nrow(genotypes.impute.directly.read.in.transposed)), genotypes.impute.directly.read.in.transposed[,-c(1,2)])

colnames(genotypes.impute.genotypes)[1:5] <- c("Snp", "allele", "chr", "pos", "cm")
ppnames=rownames(genotypes.impute.genotypes)
length(ppnames)

#this.file.G="SNP55K_maize282_AGPv2_20100513_" 
#this.file.Ext.G = "hmp.txt"

genotypes.combine <- read.delim("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript1/Miscanthus/ImputationStrategy/polyRAD/data/GoldenGate.PolyradnoLinkageD.Combine.SNPs_85obs_3182markers_04302018.txt", header=T)
genotypes.combine.directly.read.in <- genotypes.combine
genotypes.combine.directly.read.in.transposed <- t(genotypes.combine.directly.read.in[,-1])
colnames(genotypes.combine.directly.read.in.transposed) <- t(genotypes.combine.directly.read.in[,1])

genotypes.combine.temp <-  genotypes.combine.directly.read.in.transposed[,-c(1,2)]
genotypes.combine.temp[which(is.na(genotypes.combine.temp))] <- 1
genotypes.combine.directly.read.in.transposed[,-c(1,2)] <- genotypes.combine.temp

genotypes.combine.genotypes <- data.frame(rownames(genotypes.combine.directly.read.in.transposed), rep(NA,nrow(genotypes.combine.directly.read.in.transposed)), 
                                          genotypes.combine.directly.read.in.transposed[,1], genotypes.combine.directly.read.in.transposed[,2],
                                          rep(NA,nrow(genotypes.combine.directly.read.in.transposed)), genotypes.combine.directly.read.in.transposed[,-c(1,2)])

colnames(genotypes.combine.genotypes)[1:5] <- c("Snp", "allele", "chr", "pos", "cm")
genotypes.combine.genotypes <- genotypes.combine.genotypes[order(genotypes.combine.genotypes$chr, genotypes.combine.genotypes$pos),]
genotypes.combine.genotypes[1:6,1:12]

this.file.from = 1 # what are these values for? Chromosome or what?
this.file.to = 1
this.file.fragment = 1000000

this.impute.missing.data.with.hets <- TRUE

model_names <- c("rrBLUP", "sommerE", "BA", "BC", "RK")

#Output directory
this.output.dir <- paste((geno.prop*100),"%_MarkerDensity_",this.Additive.QTN.number,"_Add_QTN_", this.Dominance.QTN.number,"_Dom_QTN_",this.Epistatic.QTN.number,"_Epi_QTN_h.2_",
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


pheno_sim <- these.pheno.and.geno.data$simulated.data
dim(pheno_sim)
pheno.names <- data.frame(rownames(pheno_sim))
names(pheno.names)[1] <- "Taxa"


geno.imp.names <- data.frame(colnames(genotypes.impute.genotypes[,-c(1:5)]))
length(geno.imp.names)
names(geno.imp.names)[1] <- "Taxa"

QC.names <- merge(pheno.names, geno.imp.names, by="Taxa")
pheno_simulated <- pheno_sim[match(QC.names$Taxa, rownames(pheno_sim)),]

geno_imp.non.imp <- genotypes.impute.genotypes[,c(1:5, match(QC.names$Taxa, colnames(genotypes.impute.genotypes)))]

dim(pheno_simulated)
head(pheno_simulated)

dim(geno_imp.non.imp)
head(colnames(geno_imp.non.imp[,-c(1:5)]))

#geno_imp.non.imp # pheno_simulated
#Run k-fold cross validation on all 100 traits. The same folds are being used for each trait. #these.pheno.and.geno.data$genotypes #geno_imp.non.imp #these.pheno.and.geno.data$genotypes
#dev.off()

the.taxa.selected.from.GS <- run.k.fold.CV.using.rrBLUP(the.phenotypes = these.pheno.and.geno.data$simulated.data, the.genotypes = geno_imp.non.imp,
                           the.output.dir = this.output.dir, number.of.folds = 5, proportion.to.select = 0.10, user.input.seed.number = TRUE, seed.number = seed.number, geno.proportion=geno.prop, ncycles=number.of.cycles)
#save(the.taxa.selected.from.GS, file="the.taxa.selected.from.GS.rda")  
save.image()
setwd(home.dir)

