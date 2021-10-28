## Error evaluation HBimpute:

# Data for chromosome 1 before imputation is available at
# https://github.com/tpook92/HBimpute/blob/master/MAZE_DH_chr1.vcf.gz

# All Data sets are available upon request
# Mail to torsten.pook@uni-goettingen.de

# Error rates overall
error_rate <- matrix(0, nrow=25, ncol=10)
# Error rates per marker (for allele frequency based analysis)
ps <- NULL

for(chromo in 10:1){

  # Load in assumed underlying truth
  load(paste0("Genetic_Datasets/Batch3_KEPE/PE_DH_chromo",chromo,".RData"))
  load(paste0("MAZE_depth/genodepth",chromo,"_rename.RData"))

  # make sure order of the genotypes is the same
  a <- lines
  keep <- which(duplicated(c(a,colnames(data)))[-(1:length(a))])
  dhm <- data[,keep]
  order <- numeric(length(keep))
  for(index in 1:length(keep)){
    order[index] <- which(colnames(dhm)==a[index])
  }
  dhm <- dhm[,order]

  map <- read.table(paste0("NRGene_MAZE_sequence/PEchromo",chromo,"map.map"))

  # Only consider three datasets here
  raw_imp <- vcfR::read.vcfR(paste0("hbseq_updatemaf_chr",chromo,"_tar095_temp.vcf.gz")) # BEAGLE 5.0 dataset
  hbarr_imp <- vcfR::read.vcfR(paste0("hbarray_updatemaf_chr",chromo,"_tar095.vcf.gz")) # HB-array
  hbseq_imp <- vcfR::read.vcfR(paste0("hbseq_updatemaf_chr",chromo,"_tar095.vcf.gz")) # HB-seq

  # Compute Genotypes from VCF-file

  g1 <- substr(raw_imp@gt[,-1], start=1, stop=1)
  g2 <- substr(raw_imp@gt[,-1], start=3, stop=3)
  storage.mode(g1) <- "integer"
  storage.mode(g2) <- "integer"
  geno_raw <- g1 + g2

  g1 <- substr(hbarr_imp@gt[,-1], start=1, stop=1)
  g2 <- substr(hbarr_imp@gt[,-1], start=3, stop=3)
  storage.mode(g1) <- "integer"
  storage.mode(g2) <- "integer"
  geno_hbarr <- g1 + g2

  g1 <- substr(hbseq_imp@gt[,-1], start=1, stop=1)
  g2 <- substr(hbseq_imp@gt[,-1], start=3, stop=3)
  storage.mode(g1) <- "integer"
  storage.mode(g2) <- "integer"
  geno_hbseq <- g1 + g2


  # Calculate overlapping markers between datasets
  joint <- intersect(intersect(intersect(as.numeric(raw_imp@fix[,2]), map[,4]), as.numeric(hbarr_imp@fix[,2])), as.numeric(hbseq_imp@fix[,2]))
  overlap1 <- which(duplicated(c(joint, map[,4]))[-(1:length(joint))])
  overlap2 <- which(duplicated(c( joint, raw_imp@fix[,2]))[-(1:length(joint))])
  overlap3 <- which(duplicated(c( joint, hbarr_imp@fix[,2]))[-(1:length(joint))])
  overlap4 <- which(duplicated(c( joint, hbseq_imp@fix[,2]))[-(1:length(joint))])

  overlap5 <- which(duplicated(c( joint, posi))[-(1:length(joint))])

  # Different classes of variant calls
  class3 <- is.na(geno_imputed)[overlap4,]
  class1 <- !is.na(geno)[overlap5,] & !class3
  class2 <- !class1 & !class3

  geno_reduced <- dhm[overlap1,] # this is assumed underlying truth
  geno_raw_reduced <- geno_raw[overlap2,] # BEAGLE 5.0 imputed
  geno_hbarr_reduced <- geno_hbarr[overlap3,] # HB-array imputed
  geno_hbseq_reduced <- geno_hbseq[overlap4,] # HB-seq imputed

  error_rate[1,chromo] <- mean(abs(geno_reduced-geno_raw_reduced))/2 # Error rates
  error_rate[2,chromo] <- mean(abs(geno_reduced-geno_hbarr_reduced))/2
  error_rate[3,chromo] <- mean(abs(geno_reduced-geno_hbseq_reduced))/2
  error_rate[4,chromo] <- nrow(geno_reduced) # Number of SNPs
  error_rate[5,chromo] <- mean(abs(geno_reduced-geno_raw_reduced)[class1])/2 # Errors in Class1
  error_rate[6,chromo] <- mean(abs(geno_reduced-geno_hbarr_reduced)[class1])/2
  error_rate[7,chromo] <- mean(abs(geno_reduced-geno_hbseq_reduced)[class1])/2
  error_rate[8,chromo] <- mean(abs(geno_reduced-geno_raw_reduced)[class2])/2 # Errors in Class2
  error_rate[9,chromo] <- mean(abs(geno_reduced-geno_hbarr_reduced)[class2])/2
  error_rate[10,chromo] <- mean(abs(geno_reduced-geno_hbseq_reduced)[class2])/2
  error_rate[11,chromo] <- mean(abs(geno_reduced-geno_raw_reduced)[class3])/2 # Errors in Class3
  error_rate[12,chromo] <- mean(abs(geno_reduced-geno_hbarr_reduced)[class3])/2
  error_rate[13,chromo] <- mean(abs(geno_reduced-geno_hbseq_reduced)[class3])/2
  error_rate[14,chromo] <- sum(class1) # Total number of cells looked at per class
  error_rate[15,chromo] <- sum(class2)
  error_rate[16,chromo] <- sum(class3)

  # Imputation accuracy evaluation

  corr <- function(A,B){
    am <- rowMeans(A)
    bm <- rowMeans(B)
    cov <- rowMeans((A-am) * (B-bm))
    cova <-  rowMeans((A-am) * (A-am))
    covb <-  rowMeans((B-bm) * (B-bm))
    cor <- cov/sqrt(cova*covb)
    cor[is.na(cor)] <- 0
    return(cor)
  }

  error_rate[17,chromo] <- mean( corr(geno_reduced, geno_raw_reduced))
  error_rate[18,chromo] <- mean( corr(geno_reduced, geno_hbarr_reduced))
  error_rate[19,chromo] <-  mean( corr(geno_reduced, geno_hbseq_reduced))

  # Evaluation of REF / ALT allele
  error_rate[20,chromo] <- mean(abs(geno_reduced-geno_raw_reduced)[geno_reduced==0])/2
  error_rate[21,chromo] <- mean(abs(geno_reduced-geno_hbarr_reduced)[geno_reduced==0])/2
  error_rate[22,chromo] <- mean(abs(geno_reduced-geno_hbseq_reduced)[geno_reduced==0])/2
  error_rate[23,chromo] <- mean(abs(geno_reduced-geno_raw_reduced)[geno_reduced==2])/2
  error_rate[24,chromo] <- mean(abs(geno_reduced-geno_hbarr_reduced)[geno_reduced==2])/2
  error_rate[25,chromo] <- mean(abs(geno_reduced-geno_hbseq_reduced)[geno_reduced==2])/2

  # Error rates per MAF
  geno_reduced_temp <- geno_reduced
  geno_raw_reduced_temp <- geno_raw_reduced
  geno_hbseq_reduced_temp <- geno_hbseq_reduced
  geno_hbarr_reduced_temp <- geno_hbarr_reduced

  geno_reduced_temp1 <- geno_reduced
  geno_raw_reduced_temp1 <- geno_raw_reduced
  geno_hbseq_reduced_temp1 <- geno_hbseq_reduced
  geno_hbarr_reduced_temp1 <- geno_hbarr_reduced

  geno_reduced_temp[geno_reduced==0] <- NA
  geno_raw_reduced_temp[geno_reduced==0] <- NA
  geno_hbseq_reduced_temp[geno_reduced==0] <- NA
  geno_hbarr_reduced_temp[geno_reduced==0] <- NA

  geno_reduced_temp1[geno_reduced==2] <- NA
  geno_raw_reduced_temp1[geno_reduced==2] <- NA
  geno_hbseq_reduced_temp1[geno_reduced==2] <- NA
  geno_hbarr_reduced_temp1[geno_reduced==2] <- NA

  ps <- rbind(ps, cbind(rowMeans(geno_reduced), rowMeans(geno_raw_reduced),
                        rowMeans(geno_hbseq_reduced), rowMeans(geno_hbarr_reduced),
                        rowMeans(abs(geno_reduced-geno_raw_reduced))/2,
                        rowMeans(abs(geno_reduced-geno_hbseq_reduced))/2,
                        rowMeans(abs(geno_reduced-geno_hbarr_reduced))/2,
                        rowMeans(abs(geno_reduced_temp - geno_raw_reduced_temp), na.rm=TRUE)/2,
                        rowMeans(abs(geno_reduced_temp - geno_hbseq_reduced_temp), na.rm=TRUE)/2,
                        rowMeans(abs(geno_reduced_temp - geno_hbarr_reduced_temp), na.rm=TRUE)/2,
                        rowMeans(abs(geno_reduced_temp1 - geno_raw_reduced_temp1), na.rm=TRUE)/2,
                        rowMeans(abs(geno_reduced_temp1 - geno_hbseq_reduced_temp1), na.rm=TRUE)/2,
                        rowMeans(abs(geno_reduced_temp1 - geno_hbarr_reduced_temp1), na.rm=TRUE)/2
  ))

  print(error_rate)
}

save(file="HBimpute_updatemaf.RData", list=c("ps", "error_rate"))
