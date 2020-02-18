'#
  Authors
Torsten Pook, torsten.pook@uni-goettingen.de

Copyright (C) 2019 -- 2020  Torsten Pook

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
'#


#' Breeding function
#'
#' Function to simulate a step in a breeding scheme
#' @param vcf Path of the file to impute (alternatively provide geno, depth, allele, lines, posi manually)
#' @param hetero Set TRUE when imputing a heterozygous organism (WORK IN PROGRESS)
#' @param hb_data Dataset to derive the HB-library on (if not provided use dataset itself)
#' @param hb_map Map-file for hb_data
#' @param max_hetero Only active if hetero=FALSE - remove marker with share of heterozygosity large max_hetero (default: 0.05)
#' @param maf Applied maf-filter (default: 0)
#' @param second_round Apply HB-imputation multiple times
#' @param hetero_is_missing Only active if hetero=FALSE - If TRUE set of heterozygous markers to NA (default: TRUE)
#' @param beagle_core Number of threats using in BEAGLE 5 (default: 10)
#' @param path_beaglejar Directory of the BEAGLE jar (default: "beagle5.jar")
#' @param estimate_del Set TRUE to estimate deletion
#' @param estimate_cnv Set TRUE to estimate copy number variation
#' @param use_del Set TRUE to include identified deletions in output vcf
#' @param use_cnv Set TRUE to include identified copy number variation in output vcf
#' @param del_freq minimum frequency of a deletion to be included in the data
#' @param cnv_min Minimum estimated depth to be counted as a CNV in the output vcf
#' @param cnv_freq Minimum frequency of a CNV to be included in the output vcf
#' @param max_depth Maximum read depth in the calling (Everything higher will be set to this value)
#' @param geno genotype dataset
#' @param depth read-depth
#' @param allele allele names (default: major "A", minor "C")
#' @param line line names
#' @param posi physical marker positions
#' @param window_size HB parameter - size of each window in the algorithm (default: 20)
#' @param target_coverage HB parameter - target Coverage in the blocklist
#' @param no_overwrite If TRUE genotypes will not be overwritten
#' @export
#'

impute <- function(vcf=NULL, chromo=NULL, out = "out", geno=NULL,
                   depth=NULL, allele=NULL, lines=NULL, posi=NULL,
                   hetero=FALSE, hb_data = NULL, hb_map = NULL, hb_maf=0,
                   activ_HB=TRUE, max_hetero=0.05, hetero_is_missing=TRUE,
                   remove_del=FALSE,
                   beagle_core = 10, path_beaglejar = "beagle5.jar",
                   estimate_del=TRUE, estimate_cnv=TRUE,
                   use_cnv=FALSE, use_del=FALSE,
                   del_freq=0.2, cnv_min=3, cnv_freq=0.1,
                   max_depth = 5, maf=0, second_round=FALSE,
                   cutoff=0.9999, beagle_ne=10000,
                   window_size=20, target_coverage=NULL,
                   min_majorblock = 5000,
                   no_overwrite=FALSE
                   ){

  {

    path_prebeagle <- paste0(out,"_temp1.vcf")
    out_temp <- paste0(out, "_temp")
    if(length(vcf)>0){

      data <- vcfR::read.vcfR(vcf)
      if(length(chromo)==0){
        chromo <- as.numeric(data@fix[2,1])
      }
      take <- which(data@fix[,1]==chromo)

      name <- strsplit((data@gt)[1,1], ":")[[1]]

      subdata <- strsplit(data@gt[take,-1], ":")

      subdata <- unlist(subdata)
      genotake <- 1:(length(subdata)/length(name)) * length(name) - length(name) + which(name=="GT")


      haplo1 <- matrix(substr(subdata[genotake], start=1, stop=1), ncol=ncol(data@gt)-1)
      haplo2 <- matrix(substr(subdata[genotake], start=3, stop=3), ncol=ncol(data@gt)-1)
      allele <- data@gt[take, 4:5]

      if(length(posi)==0){
        posi <- as.numeric(data@fix[take,2])
      }
      if(length(lines)==0){
        lines <- colnames(data@gt)[-1]
      }

      storage.mode(haplo1) <- "integer"
      storage.mode(haplo2) <- "integer"

      haplo1[which(haplo1>1)] <- 1
      haplo2[which(haplo2>1)] <- 1
      geno <- haplo1 + haplo2

      if(sum(name=="DP")>0){
        depthtake <- 1:(length(subdata)/length(name)) * length(name) - length(name) + which(name=="DP")
        depth <- matrix(subdata[depthtake], ncol=ncol(data@gt)-1)
        depth[depth=="."] <- 0
        storage.mode(depth) <- "integer"
      } else if(sum(name=="AD")>0){
        depthtake <- 1:(length(subdata)/length(name)) * length(name) - length(name) + which(name=="AD")
        calls <- subdata[genotake+1]
        calls <- strsplit(calls, ",")
        calls1 <- numeric(length(calls))
        for(index in 1:length(calls)){
          if(index%%100000==0){print(index)}
          calls1[index] <- max(as.numeric(calls[[index]]))
        }
        calls1[is.na(calls1)] <- 0
        depth <- matrix(calls1, ncol=ncol(data@gt)-1)

      } else{
        depth <- !is.na(geno)
        storage.mode(depth) <- "integer"
      }
    }


    ###########################################################
    ########### Quality control: ##############################
    ###########################################################

    if(max_hetero<1 && !hetero){
      depth_file <- rowMeans(geno==1, na.rm=TRUE)
      geno <- geno[depth_file<max_hetero,]
      allele <- allele[depth_file<max_hetero,]
      depth <- depth[depth_file<max_hetero,]
      posi <- posi[depth_file<max_hetero]
    }

    if(hetero_is_missing && !hetero){
      geno[geno==1] <- NA
      depth[is.na(geno)] <- 0
    }


    add <- rowMeans(is.na(geno))

    if(sum(add==1)>0){
      geno <- geno[add<1,]
      posi <- posi[add<1]
    }

    if(maf>=0){
      p_i <- rowMeans(geno, na.rm=TRUE) / 2
      p_i[p_i>0.5] <- 1 - p_i[p_i>0.5]
      geno <- geno[p_i>maf,]
      allele <- allele[p_i>maf,]
      depth <- depth[p_i>maf,]
      posi <- posi[p_i>maf]
    }

    depth[depth > max_depth] <- max_depth

    cat(paste0(nrow(geno), " markers survieved filtering.\n"))

    if(length(hb_data)==0){


      ref <- rep("A", nrow(geno))
      alt <- rep("C", nrow(geno))
      if(hetero){
        haplo <- geno[,sort(rep(1:ncol(geno),2))]
        haplo[,1:ncol(geno)*2][haplo[,1:ncol(geno)*2]==1] <- 0
      } else{
        haplo <- geno[,sort(rep(1:ncol(geno),2))]
      }

      haplo[haplo==2] <- 1
      haplo[is.na(haplo)] <- "."

      vcfgeno <- matrix(paste0(haplo[,(1:(ncol(haplo)/2))*2], "/", haplo[,(1:(ncol(haplo)/2))*2-1]), ncol=ncol(haplo)/2)


      map <- cbind(paste0("SNP", 1:nrow(geno)), as.numeric(posi))

      options(scipen=999)
      vcfgenofull <- cbind(1, map[,2], map[,1], ref, alt, ".", "PASS", ".", "GT", vcfgeno)
      vcfgenofull <- rbind(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", lines),vcfgenofull)

      headerfile <- rbind(
        "##fileformat=VCFv4.2",
        gsub("-", "", paste0("##filedate=",  Sys.Date())),
        paste0("##source='HBimpute_v0.0'"),
        "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>"
      )

      write.table(headerfile, file=path_prebeagle, quote=FALSE, col.names = FALSE, row.names = FALSE)
      write.table(vcfgenofull, file=path_prebeagle, quote=FALSE, col.names = FALSE, row.names = FALSE, append = TRUE, sep="\t")

      beagle_commandline <- paste0("java -Xss40g -jar ", path_beaglejar," ne=", beagle_ne, " gt=",path_prebeagle," out=",out_temp," nthreads=", beagle_core)
      system(beagle_commandline)

      data_temp <- vcfR::read.vcfR(paste0(out_temp,".vcf.gz"))

      g1 <- substr(data_temp@gt[,-1], start=1, stop=1)
      g2 <- substr(data_temp@gt[,-1], start=3, stop=3)
      storage.mode(g1) <- "integer"
      storage.mode(g2) <- "integer"
      if(hetero){
        hb_data <- cbind(g1,g2)
        hb_data <- hb_data[,(rep(c(0,ncol(g1)), ncol(g1))) + sort(rep(c(1:ncol(g1)), 2))] *2
        colnames(hb_data) <- (colnames(data_temp@gt)[-1])[sort(rep(1:ncol(g1),2))]
      } else{
        hb_data <- g1 + g2
        colnames(hb_data) <- colnames(data_temp@gt)[-1]
      }

      switch <- which(data_temp@fix[,4]=="C")
      hb_data[switch,] <- 2- hb_data[switch,]

      hb_map <- as.numeric(data_temp@fix[,2])
    } else{
      if(length(hb_map)==0){
        stop("No map file for hb_data file provided!")
      }
    }

    if(activ_HB){
      cat("Start derivation of HB library:\n")
      library(HaploBlocker)

      keep <- which(duplicated(c(lines,colnames(hb_data)))[-(1:length(lines))])
      dhm <- hb_data[,keep]
      order <- numeric(length(keep))
      if(hetero){
        for(index in 1:(length(keep)/2)){
          order[c(index*2-1, index*2)] <- which(colnames(dhm)==lines[index])
        }
      } else{
        for(index in 1:length(keep)){
          order[index] <- which(colnames(dhm)==lines[index])
        }
      }

      dhm <- cbind(dhm[,order], hb_data[,-keep])
      nmax <- length(order)

      if(hb_maf > 0){
        hb_pi <- rowMeans(dhm)/max(dhm)
        hb_pi[hb_pi>0.5] <- 1 - hb_pi[hb_pi>0.5]
        keep2 <- which(hb_pi>hb_maf)
        hb_map <- hb_map[keep2]
        dhm <- dhm[keep2,]
      }
      cat(paste0(length(hb_map), " SNPs used to derive haplotype library"))

      blocklist <- block_calculation(dhm, bp=hb_map,
                                     window_size=window_size,
                                     target_coverage=target_coverage,
                                     min_majorblock = min_majorblock)

      t <- coverage_test(blocklist)
      se <- blocklist_startend(blocklist, type="bp")
      t1 <- round(mean(t)*100, digits=2)
      le <- round(mean(se[,2]-se[,1])/1000000, digit=2)
      cat(paste0("Final haplotype library with ", t1, "% coverage.\n"))
      cat(paste0("Avg. block length is ", le, " MB\n"))
      if(t1<0.8){
        warning(paste0("Coverage in haplotype library is potentially too low at ", t1, "%"))
      }
      if(le < 0.1){
        warning(paste0("Avg.block length is only ", le, " MB. Potential problems with haplotype library!"))
      }

      if(hetero){
        geno_imputed <- matrix(NA, nrow=nrow(geno), ncol=ncol(geno)*2)
        new_depth <- estimated_cnv <- estimated_deletion <- matrix(0, nrow=nrow(geno), ncol=ncol(geno)*2)

      } else{
        geno_imputed <- matrix(NA, nrow=nrow(geno), ncol=ncol(geno))
        new_depth <- estimated_cnv <- estimated_deletion <- matrix(0, nrow=nrow(geno), ncol=ncol(geno))
      }


      loc <- as.numeric(posi)

      mean_depth <- mean(depth)

      cat("Start HB imputation:\n")
      pb <- utils::txtProgressBar(min = 0, max = ncol(geno) * (1 + hetero), style = 3)
      for(nr in 1:(ncol(geno)* (1 + hetero))){
        utils::setTxtProgressBar(pb, nr)
        se <- blocklist_startend(blocklist, type="bp")
        include <- which.block(blocklist, nr)
        se <- se[which(include>0),]
        if(length(se)==0){
          next
        }
        end_block  <- sort(unique(c(0,se[,1]-1, se[,2])))[-1]
        start_block <- c(1, end_block[1:(length(end_block)-1)]+1)

        for(index in 1:length(start_block)){
          take <- which(((loc>=start_block[index])+(loc<=end_block[index]))==2)
          indi <- which.indi(blocklist, nr, start=start_block[index], end=end_block[index])
          indi <- indi[indi<=nmax]
          if(length(indi)==0){
            indi <- nr
          }
          indi <- c(indi, nr, nr, nr, nr)
          if(length(indi)>0 && length(take)>0){
            if(hetero){
              ana <- hb_data[take,indi, drop=FALSE]
            } else{
              ana <- geno[take,indi, drop=FALSE]
            }

            if(hetero){
              anad <- depth[take,ceiling(indi/2), drop=FALSE]
            } else{
              anad <- depth[take,indi, drop=FALSE]
            }

            ana[is.na(ana)] <- -999

            zero <- rowSums((ana==0)*anad)
            two <- rowSums((ana==2)*anad)

            geno_imputed[take,nr][((zero>((two)*4)))*(zero>4)*(1:length(zero))] <- 0
            geno_imputed[take,nr][((4*(zero))<two) * (two>4)* (1:length(two))] <- 2



            if(estimate_del){
              mis <- rowSums((ana==(-999))) - 4 * (ana[,ncol(ana)]==(-999))
              p_zero <- exp(-mean_depth)
              total <- length(indi) - 4

              quali_na <- which(pbinom(size = total, prob = p_zero, q=mis) > cutoff & total > 5)
              estimated_deletion[take,nr][quali_na] <- 1
            }

            if(estimate_cnv){
              new_depth[take,nr] <- zero +  two
              estimated_cnv[take,nr] <- new_depth[take,nr] / length(indi) / mean_depth
            }
          }
        }
        close(pb)

      }

      if(remove_del && sum(estimated_deletion)>0){
        geno_imputed[estimated_deletion==1] <- NA
      }

      if(second_round){

        geno <- geno_imputed
        if(estimate_cnv){
          depth <- estimated_cnv
        }


        if(hetero){
          geno_imputed <- matrix(NA, nrow=nrow(geno), ncol=ncol(geno)*2)
          new_depth <- estimated_cnv <- estimated_deletion <- matrix(0, nrow=nrow(geno), ncol=ncol(geno)*2)

        } else{
          geno_imputed <- matrix(NA, nrow=nrow(geno), ncol=ncol(geno))
          new_depth <- estimated_cnv <- estimated_deletion <- matrix(0, nrow=nrow(geno), ncol=ncol(geno))
        }


        loc <- as.numeric(posi)

        mean_depth <- mean(depth)

        cat("Start HB imputation:\n")
        pb <- utils::txtProgressBar(min = 0, max = ncol(geno) * (1 + hetero), style = 3)
        for(nr in 1:(ncol(geno)* (1 + hetero))){
          utils::setTxtProgressBar(pb, nr)
          se <- blocklist_startend(blocklist, type="bp")
          include <- which.block(blocklist, nr)
          se <- se[which(include>0),]
          if(length(se)==0){
            next
          }
          end_block  <- sort(unique(c(0,se[,1]-1, se[,2])))[-1]
          start_block <- c(1, end_block[1:(length(end_block)-1)]+1)

          for(index in 1:length(start_block)){
            take <- which(((loc>=start_block[index])+(loc<=end_block[index]))==2)
            indi <- which.indi(blocklist, nr, start=start_block[index], end=end_block[index])
            indi <- indi[indi<=nmax]
            if(length(indi)==0){
              indi <- nr
            }
            indi <- c(indi, nr, nr, nr, nr)
            if(length(indi)>0 && length(take)>0){
              if(hetero){
                ana <- hb_data[take,indi, drop=FALSE]
              } else{
                ana <- geno[take,indi, drop=FALSE]
              }

              if(hetero){
                anad <- depth[take,ceiling(indi/2), drop=FALSE]
              } else{
                anad <- depth[take,indi, drop=FALSE]
              }

              ana[is.na(ana)] <- -999

              zero <- rowSums((ana==0)*anad)
              two <- rowSums((ana==2)*anad)

              geno_imputed[take,nr][((zero>((two)*4)))*(zero>4)*(1:length(zero))] <- 0
              geno_imputed[take,nr][((4*(zero))<two) * (two>4)* (1:length(two))] <- 2



              if(estimate_del){
                mis <- rowSums((ana==(-999))) - 4 * (ana[,ncol(ana)]==(-999))
                p_zero <- exp(-mean_depth)
                total <- length(indi) - 4

                quali_na <- which(pbinom(size = total, prob = p_zero, q=mis) > cutoff & total > 5)
                estimated_deletion[take,nr][quali_na] <- 1
              }

              if(estimate_cnv){
                new_depth[take,nr] <- zero +  two
                estimated_cnv[take,nr] <- new_depth[take,nr] / length(indi) / mean_depth
              }
            }
          }
          close(pb)

        }

        if(remove_del && sum(estimated_deletion)>0){
          geno_imputed[estimated_deletion==1] <- NA
        }
      }

      if(hetero){
        geno_imputed1 <- geno_imputed
        new_depth1 <- new_depth
        estimated_cnv1 <- estimated_cnv
        estimated_deletion1 <- estimated_deletion
        geno_imputed <- (geno_imputed[,1:ncol(geno)*2-1] + geno_imputed[,1:ncol(geno)*2])/2
        new_depth_temp <- new_depth[,1:ncol(geno)*2-1]
        new_depth <- new_depth[,1:ncol(geno)*2]
        new_depth[new_depth<new_depth_temp] <- new_depth_temp[new_depth<new_depth_temp]
        estimated_cnv <- (estimated_cnv[,1:ncol(geno)*2-1] + estimated_cnv[,1:ncol(geno)*2])/2
        estimated_deletion <- (estimated_deletion[,1:ncol(geno)*2-1] + estimated_deletion[,1:ncol(geno)*2])/2

        estimated_deletion[estimated_deletion>0] <- 1
      }

      if(no_overwrite){
        geno_imputed[!is.na(geno)] <- geno[!is.na(geno)]
      }

    } else{
      geno_imputed <- geno
      new_depth <- matrix(0, nrow=nrow(geno), ncol=ncol(geno))
      estimated_cnv <- matrix(0, nrow=nrow(geno), ncol=ncol(geno))
      estimated_deletion <- matrix(0, nrow=nrow(geno), ncol=ncol(geno))
    }


    cat(paste0("Genotype calls for ", round(mean(!is.na(geno_imputed))*100, digits=2), "% of the markers were obtained"))


    snpname <-  paste0("SNP", 1:nrow(geno))

    if(use_del){
      DEL <- rowSums(estimated_deletion)
      adddel <- which(DEL>(del_freq*ncol(geno)))
    } else{
      adddel <- NULL
    }
    if(use_cnv){
      CNV <- rowSums(estimated_cnv>cnv_min)
      addcnv <- which(CNV>(cnv_freq*ncol(geno)))
    } else{
      addcnv <- NULL
    }

    map <- rbind(cbind(paste0("SNP", 1:nrow(geno)), as.numeric(posi)) ,
                 cbind(paste0("SNP", 1:nrow(geno), "_DEL")[adddel], as.numeric(posi)[adddel] ),
                 cbind( paste0("SNP", 1:nrow(geno), "_CNV")[addcnv],  as.numeric(posi)[addcnv] ))

    geno_imputed <- rbind(geno_imputed, estimated_deletion[adddel,]*2, (estimated_cnv[addcnv,]>cnv_min)*2)

    order <- sort(as.numeric(map[,2]), index.return=TRUE)$ix
    map <- map[order,]
    geno_imputed <- geno_imputed[order,]

    # Just technical stuff to avoid to markers on the same bp
    while(sum(diff(as.numeric(map[,2]))==0)>0){
      up <- which(diff(as.numeric(map[,2]))==0)+1
      map[up,2] <- as.numeric(map[up,2])+1
    }

    #######################################
    ### Write input-vcf-file for BEAGLE ###
    #######################################

    ref <- rep("A", nrow(geno))
    alt <- rep("C", nrow(geno))

    if(hetero){
      haplo <- geno_imputed[,sort(rep(1:ncol(geno),2))]
      haplo[,1:ncol(geno)*2][haplo[,1:ncol(geno)*2]==1] <- 0
    } else{
      haplo <- geno_imputed[,sort(rep(1:ncol(geno),2))]

    }
    haplo[haplo==2] <- 1
    haplo[is.na(haplo)] <- "."

    vcfgeno <- matrix(paste0(haplo[,(1:(ncol(haplo)/2))*2], "/", haplo[,(1:(ncol(haplo)/2))*2-1]), ncol=ncol(haplo)/2)


    options(scipen=999)
    vcfgenofull <- cbind(1, map[,2], map[,1], ref, alt, ".", "PASS", ".", "GT", vcfgeno)
    vcfgenofull <- rbind(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", lines),vcfgenofull)

    headerfile <- rbind(
      "##fileformat=VCFv4.2",
      gsub("-", "", paste0("##filedate=",  Sys.Date())),
      paste0("##source='HBimpute_v0.0'"),
      "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>"
    )

    write.table(headerfile, file=path_prebeagle, quote=FALSE, col.names = FALSE, row.names = FALSE)
    write.table(vcfgenofull, file=path_prebeagle, quote=FALSE, col.names = FALSE, row.names = FALSE, append = TRUE, sep="\t")

    #####################################
    ##### BEAGLE Imputation #############
    #####################################
    beagle_commandline <- paste0("java -Xss40g -jar ", path_beaglejar," ne=", beagle_ne," gt=",path_prebeagle," out=",out," nthreads=", beagle_core)
    system(beagle_commandline)

  }
}

















