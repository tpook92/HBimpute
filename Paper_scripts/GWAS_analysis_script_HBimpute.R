## GWAS Comparison HBimpute:

# All Data sets are available upon request
# Mail to torsten.pook@uni-goettingen.de


# To run the analysis with multiple random seeds
# Paper includes nr 1 - 10'000
args <- commandArgs(TRUE)
nr <- as.numeric(args[1])

# R package use to perform the GWAS analysis
library(statgenGWAS)

# Read in genotype data
{

  # All Data sets are available upon request
  # Mail to torsten.pook@uni-goettingen.de

  # load("HBimpute_HBseq.RData)
  # finally end up with following datasets:

  Z_hbseq # HBseq
  Z_hbseq2 # HBseq (only markers overlapping with array data)
  Z_hbseq3 # HBseq (large ((weaker filtering)))
  Z_array  # HBarray
  Z_beagle5 # Dataset imputed with BEAGLE 5.0
  Z_stitch # Dataset imputed with STITCH
  Z_beagle41 # Dataset imputed with BEAGLE 4.1

  # All datasets contain 321 columns and different number of rows depending on the number of markers

}


trait <- nr

# Generation of a Trait
set.seed(trait)
print(trait)

# Only place QTLs on markers with MAF > 0.1
e1 <- sample((1:nrow(Z_chip))[rowMeans(Z_chip)>0.2 & rowMeans(Z_chip)<1.8], 5)
e2 <- sample((1:nrow(Z_hbseq))[rowMeans(Z_hbseq)>0.2 & rowMeans(Z_hbseq)<1.8], 5)

# Equal size QTLs
eff1 <- rbinom(5,1,0.5)*2-1
eff2 <- rbinom(5,1,0.5)*2-1

# half of the effects based on array data, half based on HBseq
y_real <- colSums(Z_chip[e1,]*eff1)+colSums(Z_hbseq[e2,]*eff2)
y <- y_real + rnorm(length(y_real), sd=sd(y_real))

# GWAS analysis
{
  Z_chip <- t(Z_chip)
  map <- data.frame(colnames(Z_chip), 1, 1:ncol(Z_chip), "A", "B")
  colnames(map) <- c("SNP.names", "chr", "pos", "allele1", "allele2")
  rownames(map) <- map[,1]
  gData <- createGData(geno = Z_chip, map=map)
  pheno <- data.frame("all", row.names(Z_chip), y)
  colnames(pheno) <- c("Experiment", "genotype", "trait")
  phenolist <- split(x = pheno[c("genotype", "trait")],
                     f = pheno[["Experiment"]])
  gData<- createGData(gData = gData, pheno = phenolist)
  GWASresult<- runSingleTraitGwas(gData = gData,
                                  trials = "all",
                                  traits = c("trait"))
  a <- cbind(GWASresult$GWAResult$all$pos, -log(GWASresult$GWAResult$all$pValue, base=10))

  print("chip-gwas finished")


  Z_beagle5 <- t(Z_beagle5)
  map <- data.frame(colnames(Z_beagle5), 1, 1:ncol(Z_beagle5), "A", "B")
  colnames(map) <- c("SNP.names", "chr", "pos", "allele1", "allele2")
  rownames(map) <- map[,1]
  gData <- createGData(geno = Z_beagle5, map=map)
  pheno <- data.frame("all", row.names(Z_beagle5), y)
  colnames(pheno) <- c("Experiment", "genotype", "trait")
  phenolist <- split(x = pheno[c("genotype", "trait")],
                     f = pheno[["Experiment"]])
  gData<- createGData(gData = gData, pheno = phenolist)
  GWASresult<- runSingleTraitGwas(gData = gData,
                                  trials = "all",
                                  traits = c("trait"))
  b <- cbind(GWASresult$GWAResult$all$pos, -log(GWASresult$GWAResult$all$pValue, base=10))

  print("BEAGLE-seq finished")

  if(length(rownames(Z_array))==0){
    rownames(Z_array) <- paste0("SNP", 1:nrow(Z_array))
  }
  Z_array <- t(Z_array)
  map <- data.frame(colnames(Z_array), 1, 1:ncol(Z_array), "A", "B")
  colnames(map) <- c("SNP.names", "chr", "pos", "allele1", "allele2")
  rownames(map) <- map[,1]
  gData <- createGData(geno = Z_array, map=map)
  pheno <- data.frame("all", row.names(Z_array), y)
  colnames(pheno) <- c("Experiment", "genotype", "trait")
  phenolist <- split(x = pheno[c("genotype", "trait")],
                     f = pheno[["Experiment"]])
  gData<- createGData(gData = gData, pheno = phenolist)
  GWASresult<- runSingleTraitGwas(gData = gData,
                                  trials = "all",
                                  traits = c("trait"))
  c <- cbind(GWASresult$GWAResult$all$pos, -log(GWASresult$GWAResult$all$pValue, base=10))
  print("HB-array finished")

  if(length(rownames(Z_hbseq))==0){
    rownames(Z_hbseq) <- paste0("SNP", 1:nrow(Z_hbseq))
  }
  Z_hbseq <- t(Z_hbseq)
  map <- data.frame(colnames(Z_hbseq), 1, 1:ncol(Z_hbseq), "A", "B")
  colnames(map) <- c("SNP.names", "chr", "pos", "allele1", "allele2")
  rownames(map) <- map[,1]
  gData <- createGData(geno = Z_hbseq, map=map)
  pheno <- data.frame("all", row.names(Z_hbseq), y)
  colnames(pheno) <- c("Experiment", "genotype", "trait")
  phenolist <- split(x = pheno[c("genotype", "trait")],
                     f = pheno[["Experiment"]])
  gData<- createGData(gData = gData, pheno = phenolist)
  GWASresult<- runSingleTraitGwas(gData = gData,
                                  trials = "all",
                                  traits = c("trait"))
  d <- cbind(GWASresult$GWAResult$all$pos, -log(GWASresult$GWAResult$all$pValue, base=10))
  print("HB-seq finished")

  if(length(rownames(Z_hbseq2))==0){
    rownames(Z_hbseq2) <- paste0("SNP", 1:nrow(Z_hbseq2))
  }
  Z_hbseq2 <- t(Z_hbseq2)
  map <- data.frame(colnames(Z_hbseq2), 1, 1:ncol(Z_hbseq2), "A", "B")
  colnames(map) <- c("SNP.names", "chr", "pos", "allele1", "allele2")
  rownames(map) <- map[,1]
  gData <- createGData(geno = Z_hbseq2, map=map)
  pheno <- data.frame("all", row.names(Z_hbseq2), y)
  colnames(pheno) <- c("Experiment", "genotype", "trait")
  phenolist <- split(x = pheno[c("genotype", "trait")],
                     f = pheno[["Experiment"]])
  gData<- createGData(gData = gData, pheno = phenolist)
  GWASresult<- runSingleTraitGwas(gData = gData,
                                  trials = "all",
                                  traits = c("trait"))
  e <- cbind(GWASresult$GWAResult$all$pos, -log(GWASresult$GWAResult$all$pValue, base=10))
  print("HB-seq filter finished")

  if(length(rownames(Z_hbseq3))==0){
    rownames(Z_hbseq3) <- paste0("SNP", 1:nrow(Z_hbseq3))
  }
  Z_hbseq3 <- t(Z_hbseq3)
  map <- data.frame(colnames(Z_hbseq3), 1, 1:ncol(Z_hbseq3), "A", "B")
  colnames(map) <- c("SNP.names", "chr", "pos", "allele1", "allele2")
  rownames(map) <- map[,1]
  gData <- createGData(geno = Z_hbseq3, map=map)
  pheno <- data.frame("all", row.names(Z_hbseq3), y)
  colnames(pheno) <- c("Experiment", "genotype", "trait")
  phenolist <- split(x = pheno[c("genotype", "trait")],
                     f = pheno[["Experiment"]])
  gData<- createGData(gData = gData, pheno = phenolist)
  GWASresult<- runSingleTraitGwas(gData = gData,
                                  trials = "all",
                                  traits = c("trait"))
  f <- cbind(GWASresult$GWAResult$all$pos, -log(GWASresult$GWAResult$all$pValue, base=10))
  print("HB-seq-large finished")

  if(length(rownames(Z_beagle41))==0){
    rownames(Z_beagle41) <- paste0("SNP", 1:nrow(Z_beagle41))
  }
  Z_beagle41 <- t(Z_beagle41)
  map <- data.frame(colnames(Z_beagle41), 1, 1:ncol(Z_beagle41), "A", "B")
  colnames(map) <- c("SNP.names", "chr", "pos", "allele1", "allele2")
  rownames(map) <- map[,1]
  gData <- createGData(geno = Z_beagle41, map=map)
  pheno <- data.frame("all", row.names(Z_beagle41), y)
  colnames(pheno) <- c("Experiment", "genotype", "trait")
  phenolist <- split(x = pheno[c("genotype", "trait")],
                     f = pheno[["Experiment"]])
  gData<- createGData(gData = gData, pheno = phenolist)
  GWASresult<- runSingleTraitGwas(gData = gData,
                                  trials = "all",
                                  traits = c("trait"))
  g <- cbind(GWASresult$GWAResult$all$pos, -log(GWASresult$GWAResult$all$pValue, base=10))
  print("BEAGLE 4.1 finished")

  if(length(rownames(Z_stitch))==0){
    rownames(Z_stitch) <- paste0("SNP", 1:nrow(Z_stitch))
  }

  Z_stitch <- t(Z_stitch)
  map <- data.frame(colnames(Z_stitch), 1, 1:ncol(Z_stitch), "A", "B")
  colnames(map) <- c("SNP.names", "chr", "pos", "allele1", "allele2")
  rownames(map) <- map[,1]
  gData <- createGData(geno = Z_stitch, map=map)
  pheno <- data.frame("all", row.names(Z_stitch), y)
  colnames(pheno) <- c("Experiment", "genotype", "trait")
  phenolist <- split(x = pheno[c("genotype", "trait")],
                     f = pheno[["Experiment"]])
  gData<- createGData(gData = gData, pheno = phenolist)
  GWASresult<- runSingleTraitGwas(gData = gData,
                                  trials = "all",
                                  traits = c("trait"))
  h <- cbind(GWASresult$GWAResult$all$pos, -log(GWASresult$GWAResult$all$pValue, base=10))
  print("STITCH finished")

  Z_chip1 <- Z_chip[,(1:(ncol(Z_chip)/10)) * 10]
  if(length(rownames(Z_chip1))==0){
    rownames(Z_chip1) <- paste0("SNP", 1:nrow(Z_chip1))
  }
  map <- data.frame(colnames(Z_chip1), 1, 1:ncol(Z_chip1), "A", "B")
  colnames(map) <- c("SNP.names", "chr", "pos", "allele1", "allele2")
  rownames(map) <- map[,1]
  gData <- createGData(geno = Z_chip1, map=map)
  pheno <- data.frame("all", row.names(Z_chip1), y)
  colnames(pheno) <- c("Experiment", "genotype", "trait")
  phenolist <- split(x = pheno[c("genotype", "trait")],
                     f = pheno[["Experiment"]])
  gData<- createGData(gData = gData, pheno = phenolist)
  GWASresult<- runSingleTraitGwas(gData = gData,
                                  trials = "all",
                                  traits = c("trait"))
  a1 <- cbind(GWASresult$GWAResult$all$pos, -log(GWASresult$GWAResult$all$pValue, base=10))

  print("chip-gwas 50k finished")

  Z_chip1 <- Z_chip[,(1:(ncol(Z_chip)/50)) * 50]
  if(length(rownames(Z_chip1))==0){
    rownames(Z_chip1) <- paste0("SNP", 1:nrow(Z_chip1))
  }
  map <- data.frame(colnames(Z_chip1), 1, 1:ncol(Z_chip1), "A", "B")
  colnames(map) <- c("SNP.names", "chr", "pos", "allele1", "allele2")
  rownames(map) <- map[,1]
  gData <- createGData(geno = Z_chip1, map=map)
  pheno <- data.frame("all", row.names(Z_chip1), y)
  colnames(pheno) <- c("Experiment", "genotype", "trait")
  phenolist <- split(x = pheno[c("genotype", "trait")],
                     f = pheno[["Experiment"]])
  gData<- createGData(gData = gData, pheno = phenolist)
  GWASresult<- runSingleTraitGwas(gData = gData,
                                  trials = "all",
                                  traits = c("trait"))
  a2 <- cbind(GWASresult$GWAResult$all$pos, -log(GWASresult$GWAResult$all$pValue, base=10))


  print("chip-gwas 10k finished")

}

# Store results

save(file=paste0("gwas_rerun5/GWAS_results_rerun5_", trait, ".RData"), list=c("a", "a1", "a2", "b", "c", "d", "e","f", "g", "h",
                                                                                "e1", "e2", "eff1", "eff2",
                                                                                "y", "y_real"))

# To have a smaller datasets to read faster also just store GWAS highs with -log p-value > 1.5
a <- a[a[,2]>1.5,]
a1 <- a1[a1[,2]>1.5,]
a2 <- a2[a2[,2]>1.5,]
b <- b[b[,2]>1.5,]
c <- c[c[,2]>1.5,]
d <- d[d[,2]>1.5,]
e <- e[e[,2]>1.5,]
f <- f[f[,2]>1.5,]
g <- g[g[,2]>1.5,]
h <- h[h[,2]>1.5,]
rownames(a) <- NULL
rownames(a1) <- NULL
rownames(a2) <- NULL

save(file=paste0("gwas_rerun5/GWAS_results_rerun5_cop_", trait, ".RData"), list=c("a", "a1", "a2", "b", "c", "d", "e","f", "g","h",
                                                                                    "e1", "e2", "eff1", "eff2",
                                                                                    "y", "y_real"))

