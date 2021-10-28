## GWAS Comparison HBimpute:

# All Data sets are available upon request
# Mail to torsten.pook@uni-goettingen.de


# Read in genotypic and phenotypic datasets:


pheno <- read.table("Genetic_Datasets/MAZE_BLUEs_acrossLocations_DHperse2017_v1.csv", sep=";")

# Various different datasets were use. This script only includes three exemplary datasets!
Z_chip # 600k array data
Z_hbseq # HBimpute HB-seq
# Exemplary down-sampling of the array data
Z_chip1 <- Z_chip[round(seq(1, nrow(Z_chip), length.out=10000)),] # 10k array data

# If you have the miraculix package you can use a faster algorithm for calculation
# Standard is to use plain R:

if(TRUE){
  calcG <- function(Z){
    p_i <- rowMeans(Z)
    Ztm <- Z - p_i
    p_i <- p_i/2
    A <- crossprod(Ztm) / (2 * sum(p_i*(1-p_i)))
    return(A)
  }
} else{
  calcG <- function(Z){
    A <- miraculix::relationshipMatrix(miraculix::genomicmatrix(Z), centered=TRUE, normalized=TRUE)
    return(A)
  }
}

G_chip <- calcG(Z_chip)
G_chip1 <- calcG(Z_chip1)
G_hbseq <- calcG(Z_hbseq)


nrep <- 1000

acc_test1 <- acc_train1 <- acc_test2 <- acc_train2 <- acc_test3 <-
  acc_train3  <-matrix(NA, nrow=nrep, ncol=9)

pheno <- as.matrix(pheno)
storage.mode(pheno) <- "numeric"

for(rep in 1:nrep){
  print("New rep:")
  print(rep)
  for(trait in 1:9){
    print(trait)
    training_set <- sample(1:321, 280)

    y <- y_real <- as.numeric(pheno[,trait+1])
    y[-training_set] <- NA

    test <- rrBLUP::mixed.solve(y, K = G_chip, method="REML", bounds = c(1e-9,1e9))
    acc_train1[rep, trait] <- cor(y_real[training_set], test$u[training_set], use = "pairwise.complete.obs")
    acc_test1[rep, trait] <- cor(y_real[-training_set], test$u[-training_set], use = "pairwise.complete.obs")


    test <- rrBLUP::mixed.solve(y, K = G_chip1, method="REML", bounds = c(1e-9,1e9))
    acc_train2[rep, trait] <- cor(y_real[training_set], test$u[training_set], use = "pairwise.complete.obs")
    acc_test2[rep, trait] <- cor(y_real[-training_set], test$u[-training_set], use = "pairwise.complete.obs")


    test <- rrBLUP::mixed.solve(y, K = G_hbseq, method="REML", bounds = c(1e-9,1e9))
    acc_train3[rep, trait] <- cor(y_real[training_set], test$u[training_set], use = "pairwise.complete.obs")
    acc_test3[rep, trait] <- cor(y_real[-training_set], test$u[-training_set], use = "pairwise.complete.obs")


  }

  save(file="HBimpute_predict_update4.RData",
       list=c("acc_test1", "acc_test2", "acc_test3", "acc_train1", "acc_train2", "acc_train3" ))

}
