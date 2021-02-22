'#

  Authors
Torsten Pook, torsten.pook@uni-goettingen.de

Copyright (C) 2019 -- 2020  Torsten Pook

This software is free for non-commercial use. It may be copied,
modified, and redistributed provided that this copyright notice
is preserved on all copies. The intellectual property rights of
the algorithms used reside with the University of Goettingen
Animal Breeding and Genetics Group.

You may not use this software, in whole or in part, in support
of any commercial product without the express consent of the
author.

The HBimpute step implemented in this software is is patent pending
under application number EP20201121.9.
Patent applicants are KWS SAAT SE & Co. KGaA and University of Goettingen.
Inventors are Torsten Pook and Adnane Nemri.

There is no warranty or other guarantee of fitness of this
software for any purpose. It is provided solely "as is".
'#



#' Breeding function
#'
#' Function to simulate a step in a breeding scheme
#' @param vcf Path of the file to impute (alternatively provide geno, depth, allele, lines, posi manually)
#' @param vcf_ref Path of a reference panel
#' @param vcf_RData Path of an RData object to read instead of plain vcf-file
#' @param vcf_ref_RData Path of an RData object to read instead of plain vcf-file
#' @param chromo Manually provide a vector of chromosomes to perform imputation for (default: Process all available chromosomes)
#' @param out Path of the output file
#' @param hetero Set TRUE when imputing a heterozygous organism. The two haplotype of a individual are then processed separately
#' @param hb_data Dataset to derive the HB-library on (if not provided use dataset itself)
#' @param hb_map Map-file for hb_data
#' @param max_hetero Only active if hetero=FALSE - remove marker with share of heterozygosity large max_hetero (default: 0.05)
#' @param maf Applied maf-filter (default: 0)
#' @param hetero_is_missing Only active if hetero=FALSE - If TRUE set of heterozygous markers to NA (default: TRUE)
#' @param beagle_core Number of threats using in BEAGLE 5 (default: 1)
#' @param path_beaglejar Directory of the BEAGLE jar (default: "beagle5.jar")
#' @param beagle_ne Corresponds to ne parameter in BEAGLE
#' @param estimate_sv Set TRUE to activate structural variation detection
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
#' @param snpname SNP IDss
#' @param lines line names
#' @param posi physical marker positions
#' @param activ_HB Set to FALSE to not perform imputation via HBimpute (why do you use the package?!)
#' @param window_size HB parameter - size of each window in the algorithm (default: 20)
#' @param target_coverage HB parameter - target Coverage in the blocklist
#' @param overwrite_call Set to TRUE to always preserve original SNP calls
#' @param overwrite_call_min_depth Minimum read depth to overwrite a call (default: 1)
#' @param overwrite_na Set to FALSE to not replace NA calls in HB with the original call
#' @param overwrite_na_min_depth Minimum read depth to overwrite a call (default: 1)
#' @param extended_output Set to TRUE for all relevant outputs
#' @param zero_two_coding Set to TRUE in case you are providing an input in geno with 0/2 coded genotypes (and 1 heterozygous!)
#' @param log Path of a potential log file to be written (default: NULL will take vcf-file name)
#' @param hb_maf Minor allele frequency filter for the Haplotype library
#' @param min_majorblock HaploBlocker parameter
#' @param sv_cut1 sv_cut1
#' @param sv_cut2 sv_cut2
#' @param sv_window sv_window
#' @param ref_panel ref_panel
#' @param remove_del Set TRUE to set estimated deletions to NA calls
#' @param cutoff Old deletion algorithm binomial cutoff
#' @param min_confi Proportion of read the need to support one variant to be called (default: 4 - coding 4:1 ratio)
#' @param quali_filter Set TRUE to active quali filter after HBimpute step
#' @param max_na Cutoff for NA-filter after HBimpute step
#' @param min_depth Read-depth filter after HBimpute step
#' @param hb_max Maximum number of marker to consider in the haplotype library
#' @export
#'

impute <- function(vcf=NULL, vcf_ref=NULL, vcf_RData=NULL, vcf_ref_RData = NULL,
                   log = NULL, chromo=NULL, out = "out", geno=NULL,
                   depth=NULL, allele=NULL, lines=NULL, posi=NULL,
                   hetero=FALSE, hb_data = NULL, hb_map = NULL, hb_maf=0,
                   activ_HB=TRUE, max_hetero=0.01, hetero_is_missing=TRUE,
                   remove_del=FALSE,
                   beagle_core = 1, path_beaglejar = "beagle5.jar",
                   estimate_del=NULL, estimate_cnv=NULL,
                   use_cnv=FALSE, use_del=FALSE,
                   del_freq=0.1, cnv_min=2, cnv_freq=0.1,
                   max_depth = 10L, maf=0, zero_two_coding = FALSE,
                   cutoff=0.9999, beagle_ne=10000,
                   window_size=20, target_coverage=NULL,
                   min_majorblock = 5000,
                   overwrite_call=FALSE,
                   overwrite_call_min_depth = 1,
                   overwrite_na=TRUE,
                   overwrite_na_min_depth = 1,
                   min_confi = 4,
                   ref_panel=NULL,
                   extended_output= FALSE,
                   quali_filter = TRUE,
                   max_na = 0.5,
                   min_depth = NULL,
                   estimate_sv = FALSE,
                   sv_cut1 = 1.3,
                   sv_cut2 = 0.7,
                   sv_window = 250000,
                   hb_max = Inf,
                   early_remove = FALSE,
                   node_min_early = 2,
                   share_called = 0,
                   remove_lines = NULL
                   ){
  {


    options("scipen"=999)


    if(length(min_depth)==0){
      if(length(vcf_ref)>0){
        min_depth <- 0
      } else{
        min_depth <- 0.5
      }
    }
    if(length(estimate_cnv)==0 & use_cnv==FALSE){
      estimate_cnv <- FALSE
    } else if(length(estimate_cnv)==0){
      estimate_cnv <- TRUE
    }
    if(length(estimate_del)==0 & use_del==FALSE){
      estimate_del <- FALSE
    } else if(length(estimate_del)==0){
      estimate_del <- TRUE
    }

    out_temp <- paste0(out, "_temp")
    ref_temp <- paste0(out, "_reftemp")
    if(length(vcf)>0){
      if(length(vcf_RData)>0){
        load(vcf_RData)
      } else{
        data <- vcfR::read.vcfR(vcf)
      }

      if(length(remove_lines)>0){
        data@gt <- data@gt[,-(remove_lines+1)]
      }
      if(length(chromo)==0){
        chromo <- as.numeric(data@fix[2,1])
      }
      if(length(chromo)==1 && is.na(chromo)){
        chromo <- data@fix[2,1]
      }
    } else{
      if(length(chromo)==0){
        chromo <- 1
      }
    }

    if(length(vcf_ref)>0){
      if(length(vcf_ref_RData)>0){
        load(vcf_ref_RData)
      } else{
        data_ref <- vcfR::read.vcfR(vcf_ref)
      }
    }

    chromo_set <- chromo
    out_store <- out

    if(length(chromo_set)>1){
      hb_store <- hb_data
    }

    ##
    if(FALSE){
      data@fix <- data@fix[-(4:5),]
      data@gt <- data@gt[-(4:5),]
    }

    for(chromo in chromo_set){

      if(length(log)==0 && length(vcf)>0){
        log <- paste0(vcf, "_chromo_",chromo,".log")
      } else if(length(log)==0 && length(vcf)==0){
        log <- FALSE
      }

      if(log != FALSE){
        zz <- file(log, open="wt")
        sink(zz, append = TRUE, type = c("output"), split=TRUE)

        cat("#############################################################
        ############ HBimpute #############
        #############################################################")
        cat("\n")
        cat(paste0("Imputation started: ", Sys.time(), "\n"))
        cat(paste0("HBimpute version used: ", utils::sessionInfo()$otherPkgs$HBimpute$Version, "\n"))
        cat("Copyright (C) 2019-2021 Torsten Pook\n\n")
      }

      out <- out_store
      cat(paste0("Start processing chromosome ", chromo))
      if(length(chromo_set)>1){
        out <- paste0(out, "_", chromo)
        out_temp <- paste0(out, "_temp")
        ref_temp <- paste0(out, "_reftemp")
      }
      path_prebeagle1 <- paste0(out,"_temp1.vcf")
      path_prebeagle2 <- paste0(out,"_temp2.vcf")
      path_prebeagle3 <- paste0(out,"_temp3.vcf")
      path_prebeagle4 <- paste0(out,"_temp4.vcf")

      if(length(vcf)>0){

        ref_beagle <- FALSE

        take <- take_ref <- which(data@fix[,1]==chromo)
        if(length(vcf_ref)>0){
          take_ref <- which(data_ref@fix[,1]==chromo)
          cat(paste0("Reference panel with ", length(take_ref), " SNPs.\n"))
          cat(paste0("of which ", length(take), " SNPs are included in the data panel itself.\n"))

        }

        if(length(vcf_ref)==0){
          data_joined <- data
          is_ref <- rep(FALSE, ncol(data@gt)-1)
        } else{
          marker_panel <- data@fix[take,]
          gt_panel_set <- data@gt[take,]
          marker_ref <- data_ref@fix[take_ref,]
          data_ref_chromo <- data_ref
          data_ref_chromo@fix <- marker_ref
          data_ref_chromo@gt <- data_ref@gt[take_ref,]

          n_indi <- ncol(gt_panel_set)-1
          gt_panel <- matrix(NA, ncol= n_indi, nrow = nrow(marker_ref))
          replacement <- numeric(nrow(marker_panel))
          for(index in 1:nrow(marker_panel)){
            replacement[index] <- which(marker_ref[,2] == marker_panel[index,2])
          }
          gt_panel[replacement,] <- gt_panel_set[,-1]
          colnames(gt_panel) <- colnames(gt_panel_set)[-1]

          data_joined <- data_ref_chromo

          replace_na <- is.na(data_joined@fix[replacement,5])
          data_joined@fix[replacement,5][replace_na] <- marker_panel[,5][replace_na]
          if(sum(data_joined@fix[replacement,5]==data_joined@fix[replacement,4], na.rm=TRUE)>0){
            data_joined@fix[replacement,5][data_joined@fix[replacement,5]==data_joined@fix[replacement,4]] <-  marker_panel[data_joined@fix[replacement,5]==data_joined@fix[replacement,4],4]
          }


          stay <- which(rowSums(data_joined@fix[replacement,4:5]==marker_panel[,4:5], na.rm=TRUE)==2 | ((data_joined@fix[replacement,4]==marker_panel[,4]) + is.na(data_joined@fix[replacement,5]))==2)
          switch <- which(rowSums(data_joined@fix[replacement,4:5]==marker_panel[,5:4], na.rm=TRUE)==2)
          miss <- which(rowSums(is.na(marker_panel[,4:5]))>0)
          miss <- setdiff(miss,stay)
          remove <- (1:nrow(marker_panel))[-c(stay,switch,miss)]

          if(sum(!is.na(gt_panel[replacement,][miss,]))>0){
            cat("Check your inputs! Markers with entered genotypes were identified as missing.\n")
            gt_panel[replacement,][miss,] <- NA
          }
          if(length(remove)>0){
            cat(paste0("Remove ", length(remove), " SNP with missmatching alleles in the data panel and reference panel!\n"))
          }

          data_joined@gt <- cbind(data_joined@gt[,1], gt_panel, data_joined@gt[,-1])
          is_ref <- c(rep(FALSE, ncol(gt_panel)), rep(TRUE, ncol(data_ref_chromo@gt)-1))
        }


        name <- strsplit((data_joined@gt)[1,1], ":")[[1]]


        subdata_split <- list()

        cat("Start data import conversion:\n")
        pb <- utils::txtProgressBar(min = 0, max = ncol(data_joined@gt), style = 3)

        for(subindex in 2:ncol(data_joined@gt)){
          utils::setTxtProgressBar(pb, subindex)

          subdata <- strsplit(data_joined@gt[,subindex], ":", fixed=TRUE)
          lname <- length(name)
          for(index in 1:length(subdata)){
            if(length(subdata[[index]])<lname){
              subdata[[index]] <- c(subdata[[index]], rep(NA, length(name)-length(subdata[[index]])))
            }
          }
          subdata_split[[subindex]] <- subdata

        }
        close(pb)


        subdata <- unlist(subdata_split)
        rm(subdata_split)

        genotake <- 1:(length(subdata)/length(name)) * length(name) - length(name) + which(name=="GT")


        haplo1 <- matrix(substr(subdata[genotake], start=1, stop=1), ncol=ncol(data_joined@gt)-1)
        haplo2 <- matrix(substr(subdata[genotake], start=3, stop=3), ncol=ncol(data_joined@gt)-1)
        snpname <- data_joined@fix[,3]
        allele <- data_joined@fix[, 4:5]


        posi <- as.numeric(data_joined@fix[,2])
        if(length(lines)==0){
          lines <- colnames(data_joined@gt)[-1]
        }

        suppressWarnings(storage.mode(haplo1) <- "integer")
        suppressWarnings(storage.mode(haplo2) <- "integer")


        if(length(switch)>0 && sum(is_ref)>0){

          if(sum(haplo1[replacement[switch]]>1, na.rm=TRUE)>0){
            stop("Allele switch from reference panel called for non-binary variant?!")
          }
          haplo1[replacement[switch],] <- 1 - haplo1[replacement[switch],]
          haplo2[replacement[switch],] <- 1 - haplo2[replacement[switch],]
        }
        if(hetero){
          haplo1[which(haplo1>1)] <- 1
          haplo2[which(haplo2>1)] <- 1
          geno <- haplo1 + haplo2
        } else{
          geno <- haplo1
          geno[geno!=haplo2] <- NA

        }


        if(sum(name=="DP")>0){
          depthtake <- 1:(length(subdata)/length(name)) * length(name) - length(name) + which(name=="DP")
          depth <- matrix(subdata[depthtake], ncol=ncol(data@gt)-1)
          depth[depth=="."] <- 0L
          depth[is.na(depth)] <- 0L

          storage.mode(depth) <- "integer"
        } else if(sum(name=="AD")>0){
          depthtake <- 1:(length(subdata)/length(name)) * length(name) - length(name) + which(name=="AD")
          calls <- subdata[genotake+1]
          calls <- strsplit(calls, ",")
          calls1 <- numeric(length(calls))
          printer <- 250000
          for(index in 1:length(calls)){
            if(index==printer){
              print(index)
              printer <- printer + 250000
            }
            calls1[index] <- max(as.numeric(calls[[index]]))
          }
          calls1[is.na(calls1)] <- 0
          depth <- matrix(calls1, ncol=ncol(data@gt)-1)
          depth[is.na(depth)] <- 0
          storage.mode(depth) <- "integer"
        } else{
          depth <- !is.na(geno)
          storage.mode(depth) <- "integer"
        }


        if(length(remove)>0 && sum(is_ref)>0){
          remove <- replacement[remove]
          geno <- geno[-remove,]
          allele <- allele[-remove,]
          snpname <- snpname[-remove]
          depth <- depth[-remove,]
          posi <- posi[-remove]
        }


        if(sum(is_ref)>0){
          if(sum(is.na(geno[,is_ref]))){
            cat("Reference panel contains missing calls. Apply BEAGLE on Reference panel before imputation!\n")
            ref_beagle <- TRUE
          }
          is_phased <- substr(data_joined@gt[,which(is_ref)[1]+1], start=2, stop=2)
          if(sum(is_phased=="/", na.rm=TRUE) > sum(is_phased=="|", na.rm = TRUE)){
            cat("Reference panel is not phased (use of '/' instead of '|')! Apply BEAGLE on Reference panel before imputation!\n")
            ref_beagle <- TRUE
          }
        }

      }
      if(hetero){
        is_ref <- is_ref[sort(rep(1:length(is_ref),2))]
      }

    }

    if(length(lines)==0){
      lines <- paste0("ID_", 1:ncol(geno))
      colnames(geno) <- colnames(haplo1) <- colnames(haplo2) <- lines
    }

    ###########################################################
    ########### Quality control: ##############################
    ###########################################################

    if((max_hetero<1 && !hetero) && (length(vcf)>0 || zero_two_coding)){
      if(length(vcf)>0){
        check1 <- haplo1!=haplo2
      } else if(zero_two_coding){
        geno[geno>2] <- 1L
        check1 <- geno==1
        geno[check1] <- NA
        geno <- geno / 2L
      }

      depth_file <- rowMeans(check1, na.rm=TRUE)
      geno <- geno[depth_file<max_hetero,]
      allele <- allele[depth_file<max_hetero,]
      snpname <- snpname[depth_file<max_hetero]
      depth <- depth[depth_file<max_hetero,]
      posi <- posi[depth_file<max_hetero]
    }

    if(!hetero){
      depth[is.na(geno)] <- 0L
    }
    add <- rowMeans(is.na(geno))

    if(sum(add==1)>0){
      geno <- geno[add<1,]
      posi <- posi[add<1]
      allele <- allele[add<1,]
      snpname <- snpname[add<1]
      depth <- depth[add<1,]
    }

    if(maf>=0){
      if(hetero){
        p_i <- rowMeans(geno, na.rm=TRUE) / 2
      } else{
        p_i <- rowMeans(geno>0, na.rm=TRUE)
      }

      p_i[p_i>0.5] <- 1 - p_i[p_i>0.5]
      geno <- geno[p_i>maf,]
      allele <- allele[p_i>maf,]
      snpname <- snpname[p_i>maf]
      depth <- depth[p_i>maf,]
      posi <- posi[p_i>maf]
    }

    if(share_called>0){
      called <- rowMeans(!is.na(geno))

      geno <- geno[called>share_called,]
      posi <- posi[called>share_called]
      allele <- allele[called>share_called,]
      snpname <- snpname[called>share_called]
      depth <- depth[called>share_called,]

    }

    depth[depth > max_depth] <- max_depth


    if(length(allele)>0){
      ref <- allele[,1]
      alt <- allele[,2]
    } else{
      ref <- rep("A", nrow(geno))
      alt <- rep("C", nrow(geno))
    }

    if(length(snpname)!= nrow(allele)){
      snpname <-  paste0("SNP", 1:nrow(allele))
    }

    if(length(allele)>0 && (sum(is.na(allele))>0 || sum(allele[,1]==allele[,2])>0)){
      remove3 <- which(rowSums(is.na(allele))>0 | (allele[,1]==allele[,2]))
      geno <- geno[-remove3,]
      allele <- allele[-remove3,]
      snpname <- snpname[-remove3]
      depth <- depth[-remove3,]
      posi <- posi[-remove3]
      ref <- allele[,1]
      alt <- allele[,2]
    }

    cat(paste0(nrow(geno), " markers survieved filtering.\n"))

    ###########################################################
    ###### Generate dataset to derive Haplotype Library #######
    ############ Auxillary imputed SNP dataset ################
    ###########################################################

    if(length(chromo_set)>1){
      hb_data <- hb_store
    }
    if(length(hb_data)<=1){

      if(length(hb_data)==1){

        data_temp <- vcfR::read.vcfR(hb_data)

      } else{
        if(hetero){
          haplo <- geno[,sort(rep(1:ncol(geno),2))]
          haplo[,1:ncol(geno)*2][haplo[,1:ncol(geno)*2]==1] <- 0L
        } else{
          haplo <- geno[,sort(rep(1:ncol(geno),2))]
        }

        if(hetero){
          haplo[haplo==2] <- 1L
        }

        haplo[is.na(haplo)] <- "."


        separator <- rep("/", length(is_ref))
        if(!ref_beagle){
          separator[is_ref] <- "|"
        }

        vcfgeno <- matrix(paste0(haplo[,(1:(ncol(haplo)/2))*2], separator, haplo[,(1:(ncol(haplo)/2))*2-1]), ncol=ncol(haplo)/2)


        if(nrow(haplo)==length(snpname)){
          map <- cbind(snpname , as.numeric(posi))
        } else{
          map <- cbind(paste0("SNP", 1:nrow(haplo)), as.numeric(posi))
        }



        if(sum(is_ref)>0){
          options(scipen=999)
          ref_indi <- which(is_ref)
          vcfgenofull <- cbind(chromo, map[,2], map[,1], ref, alt, ".", "PASS", ".", "GT", vcfgeno[,ref_indi])
          vcfgenofull <- rbind(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", lines[ref_indi]),vcfgenofull)

          headerfile <- rbind(
            "##fileformat=VCFv4.2",
            gsub("-", "", paste0("##filedate=",  Sys.Date())),
            paste0("##source='HBimpute_v0.0'"),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>"
          )

          utils::write.table(headerfile, file=path_prebeagle4, quote=FALSE, col.names = FALSE, row.names = FALSE)
          utils::write.table(vcfgenofull, file=path_prebeagle4, quote=FALSE, col.names = FALSE, row.names = FALSE, append = TRUE, sep="\t")

          if(ref_beagle){
            beagle_commandline <- paste0("java -jar ", path_beaglejar," ne=", beagle_ne, " gt=",path_prebeagle4," out=",ref_temp," nthreads=", beagle_core)
            system(beagle_commandline)
          }
        }

        options(scipen=999)
        panel_indi <- which(!is_ref)
        vcfgenofull <- cbind(chromo, map[,2], map[,1], ref, alt, ".", "PASS", ".", "GT", vcfgeno[,panel_indi])
        vcfgenofull <- rbind(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", lines[panel_indi]),vcfgenofull)

        headerfile <- rbind(
          "##fileformat=VCFv4.2",
          gsub("-", "", paste0("##filedate=",  Sys.Date())),
          paste0("##source='HBimpute_v0.0'"),
          "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>"
        )

        utils::write.table(headerfile, file=path_prebeagle1, quote=FALSE, col.names = FALSE, row.names = FALSE)
        utils::write.table(vcfgenofull, file=path_prebeagle1, quote=FALSE, col.names = FALSE, row.names = FALSE, append = TRUE, sep="\t")

        if(sum(is_ref)>0){
          if(ref_beagle){
            ref_path <- paste0(ref_temp, ".vcf.gz")
          } else{
            ref_path <- path_prebeagle4
          }
          beagle_commandline <- paste0("java -jar ", path_beaglejar," ne=", beagle_ne, " gt=",path_prebeagle1," out=",out_temp," ref=", ref_path ," nthreads=", beagle_core)
        } else{
          beagle_commandline <- paste0("java -jar ", path_beaglejar," ne=", beagle_ne, " gt=",path_prebeagle1," out=",out_temp," nthreads=", beagle_core)
        }

        system(beagle_commandline)

        data_temp <- vcfR::read.vcfR(paste0(out_temp,".vcf.gz"))
        if(sum(is_ref)>0){
          data_temp_ref <-vcfR::read.vcfR(ref_path)

          data_temp@gt <- cbind(data_temp@gt, data_temp_ref@gt[,-1])
        }


      }



      g1 <- substr(data_temp@gt[,-1], start=1, stop=1)
      g2 <- substr(data_temp@gt[,-1], start=3, stop=3)
      storage.mode(g1) <- "integer"
      storage.mode(g2) <- "integer"
      if(hetero){
        hb_data <- cbind(g1,g2)
        g1m <- g1
        g2m <- g2
        g1m[depth==0L] <- NA
        g2m[depth==0L] <- NA
        hb_data_miss <- cbind(g1m,g2m)
        hb_data_miss <- hb_data_miss[,(rep(c(0,ncol(g1m)), ncol(g1m))) + sort(rep(c(1:ncol(g1m)), 2))] *2
        hb_data <- hb_data[,(rep(c(0,ncol(g1)), ncol(g1))) + sort(rep(c(1:ncol(g1)), 2))] *2
        colnames(hb_data) <- (colnames(data_temp@gt)[-1])[sort(rep(1:ncol(g1),2))]
      } else{
        hb_data <- g1 + g2
        colnames(hb_data) <- colnames(data_temp@gt)[-1]
      }


      hb_map <- as.numeric(data_temp@fix[,2])
    } else{
      if(length(hb_map)==0){
        stop("No map file for hb_data file provided!")
      }
    }

    ###########################################################
    ########################## HBimpute step ##################
    ###########################################################

    if(activ_HB){
      cat("Start derivation of HB library:\n")
      library(HaploBlocker)

      if(length(colnames(hb_data))==0){
        colnames(hb_data) <- lines
      }
      keep <- which(duplicated(c(lines,colnames(hb_data)))[-(1:length(lines))])
      dhm <- hb_data[,keep]
      order <- numeric(length(keep))
      to_remove <- NULL
      if(hetero){
        for(index in 1:(length(keep)/2)){
          order[c(index*2-1, index*2)] <- which(colnames(dhm)==lines[index])
        }
      } else{
        for(index in 1:length(keep)){
          keeper <- which(colnames(dhm)==lines[index])
          if(length(keeper)==1){
            order[index] <- keeper
          } else{
            to_remove <- c(to_remove, index)
            warning("Lines without haplotype library representative")
          }

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
      cat(paste0(length(hb_map), " SNPs used to derive haplotype library.\n"))

      if(FALSE){
        share_mis <- rowMeans(haplo==".")
        keep_hb <- which(share_mis < quantile(share_mis,0.5))

        dhm <- dhm[keep_hb,]
        hb_map <- hb_map[keep_hb]
      }
      if(hb_max < nrow(dhm)){
        remains <- sort(sample(1:nrow(dhm), hb_max))
        dhm <- dhm[remains,]
        hb_map <- hb_map[remains]
      }
      cat(paste0(length(hb_map), " SNPs used to derive haplotype library.\n"))

      blocklist <- HaploBlocker::block_calculation(dhm, bp=hb_map, node_min = 3,
                                                   edge_min = 3,
                                                   window_size=window_size,
                                                   target_coverage=target_coverage,
                                                   min_majorblock = min_majorblock,
                                                   weighting_length = 2,
                                                   early_remove = early_remove,
                                                   node_min_early = node_min_early)

      t <- HaploBlocker::coverage_test(blocklist)
      se <- HaploBlocker::blocklist_startend(blocklist, type="snp")
      print(mean(HaploBlocker::blocklist_size(blocklist)))
      cat(paste0("Avg. block length is ", round(mean(se[,2]-se[,1])), " SNPs\n"))
      se <- HaploBlocker::blocklist_startend(blocklist, type="bp")
      t1 <- round(mean(t)*100, digits=2)
      le <- round(mean(se[,2]-se[,1])/1000000, digits=2)
      cat(paste0("Final haplotype library with ", t1, "% coverage.\n"))
      cat(paste0("Avg. block length is ", le, " MB\n"))
      if(t1<80){
        warning(paste0("Coverage in haplotype library is potentially too low at ", t1, "%"))
      }
      if(le < 0.1){
        warning(paste0("Avg.block length is only ", le, " MB. Potential problems with haplotype library!"))
      }

      if(length(to_remove)>0){
        if(hetero){

        } else{
          geno <- geno[,-to_remove]
          depth <- depth[,-to_remove]
        }

      }

      if(hetero){
        geno_imputed <- matrix(NA, nrow=nrow(haplo), ncol=sum(!is_ref))
        new_depth <- hb_depth <-  estimated_cnv <- estimated_deletion <- matrix(0L, nrow=nrow(haplo), ncol=sum(!is_ref))
      } else{
        geno_imputed <- matrix(NA, nrow=nrow(geno), ncol=sum(!is_ref))
        new_depth <- hb_depth <-  estimated_cnv <- estimated_deletion <- matrix(0L, nrow=nrow(geno), ncol=sum(!is_ref))
      }


      loc <- as.numeric(posi)

      mean_depth <- mean(depth)


      if(hetero){

      } else{
        geno_temp <- geno
        geno_temp[is.na(geno_temp)] <- 0
        rmax <- RandomFieldsUtils::colMax(t(geno_temp))
      }

      cat("Start HB imputation:\n")
      pb <- utils::txtProgressBar(min = 0, max = sum(!is_ref), style = 3)
      for(nr in 1:sum(!is_ref)){
        utils::setTxtProgressBar(pb, nr)
        se <- se_temp <- HaploBlocker::blocklist_startend(blocklist, type="bp")
        include <- HaploBlocker::which.block(blocklist, nr)
        se <- se[which(include>0),,drop=FALSE]
        if(length(se)==0){
          next
        }
        end_block  <- sort(unique(c(0,se[,1]-1, se[,2], max(as.numeric(posi)))))[-1]
        start_block <- c(1, end_block[1:(length(end_block)-1)]+1)

        p2 <- numeric(length(start_block))
        activ <- 1
        for(index in 1:length(loc)){
          while(activ<=length(start_block) && loc[index]>=start_block[activ]){
            activ <- activ + 1
          }
          if(activ<=length(start_block)){
            p2[activ-1] <- index
          }
        }

        if(p2[1]==0){
          p2[1] <- -1
        }
        while(sum(p2==0)>0){
          p2[p2==0] <- p2[which(p2==0)-1]
        }
        if(p2[1]==(-1)){
          p2[1] <- 0
        }
        p2[length(p2)] <- sum(loc<=max(end_block))
        p1 <- c(1, p2[1:(length(p2)-1)]+1)

        for(index in (1:length(start_block))[p2>0]){
          take <- p1[index]:p2[index]
          indi <- HaploBlocker::which.indi(blocklist, nr, start=start_block[index], end=end_block[index], se = se_temp)
          indi <- indi[indi<=nmax]
          if(length(indi)==0){
            indi <- nr
          }
          indi <- c(indi, nr, nr, nr, nr)
          if(length(indi)>0 && length(take)>0){
            if(hetero){
              ana <- hb_data_miss[take,indi, drop=FALSE]
            } else{
              ana <- geno[take,indi, drop=FALSE]
            }

            if(hetero){
              anad <- depth[take,ceiling(indi/2), drop=FALSE]
            } else{
              anad <- depth[take,indi, drop=FALSE]
            }

            ana[is.na(ana)] <- -999L

            if(hetero){
              zero <- rowSums((ana==0L)*anad)
              two <- rowSums((ana==2L)*anad)
            } else{
              n_variant <- max(rmax[take]) + 1
              counter_list <- list()
              totalr <- numeric(length(take))
              for(index2 in 1:n_variant){
                counter_list[[index2]] <- rowSums((ana==as.integer(index2-1))*anad)
                totalr <- totalr +  counter_list[[index2]]
              }

            }

            if(hetero){
              geno_imputed[take,nr][((zero>((two)*min_confi)))*(zero>min_confi)*(1:length(zero))] <- 0L
              geno_imputed[take,nr][((min_confi*(zero))<two) * (two>min_confi)* (1:length(two))] <- 2L
            } else{
              for(index2 in 1:n_variant){
                geno_imputed[take,nr][((counter_list[[index2]] / totalr) >= (min_confi/ (min_confi+1)))*(counter_list[[index2]]>min_confi)*(1:length(totalr))] <- index2 - 1L
              }

            }





            hb_depth[take,nr] <- length(indi)

            if(estimate_del){
              mis <- rowSums((ana==(-999))) - 4 * (ana[,ncol(ana)]==(-999))
              p_zero <- exp(-mean_depth)
              total <- length(indi) - 4

              quali_na <- which(stats::pbinom(size = total, prob = p_zero, q=mis) > cutoff & total > 5)
              estimated_deletion[take,nr][quali_na] <- 1
            }

            if(hetero){
              new_depth[take,nr] <- zero +  two
            } else{
              new_depth[take,nr] <- totalr
            }


            if(estimate_cnv){
              estimated_cnv[take,nr] <- new_depth[take,nr] / length(indi) / mean_depth
            }
          }
        }
        close(pb)

      }

      if(remove_del && sum(estimated_deletion)>0){
        geno_imputed[estimated_deletion==1] <- NA
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

      if(overwrite_call){
        if(hetero){
          replaces <- which((geno_imputed[!is.na(geno)] != geno[!is.na(geno)]) |
                              (is.na(geno_imputed[!is.na(geno)]) ))
          reper <- geno[!is.na(geno)][replaces]
          reper[reper==1] <- 0
          geno_imputed1[,1:ncol(geno)*2-1][replaces] <- reper
          reper <- geno[!is.na(geno)][replaces]
          reper[reper==1] <- 2
          geno_imputed1[,1:ncol(geno)*2][replaces] <- reper
        }

        to_replace <- !is.na(geno) & depth >= overwrite_call_min_depth
        geno_imputed[to_replace] <- geno[to_replace]
      }
      if(overwrite_na & !hetero){
        to_replace <- is.na(geno_imputed) & !is.na(geno[,!is_ref]) & depth[,!is_ref] >= overwrite_na_min_depth
        geno_imputed[to_replace] <- geno[,!is_ref][to_replace]
      } else if(overwrite_na & hetero){
        replaces <- which(is.na(geno_imputed[!is.na(geno)]))
        reper <- geno[!is.na(geno)][replaces]
        reper[reper==1] <- 0
        geno_imputed1[,1:ncol(geno)*2-1][replaces] <- reper
        reper <- geno[!is.na(geno)][replaces]
        reper[reper==1] <- 2
        geno_imputed1[,1:ncol(geno)*2][replaces] <- reper


        geno_imputed[!is.na(geno)] <- geno[!is.na(geno)]
      }


    } else{
      geno_imputed <- geno
      new_depth <- estimated_cnv <- estimated_deletion <- hb_depth <- matrix(0, nrow=nrow(geno), ncol=ncol(geno))
    }

    if(hetero){
      cat(paste0("Genotype calls for ", round(mean(!is.na(geno_imputed1))*100, digits=2), "% of the markers were obtained"))
    } else{
      cat(paste0("Genotype calls for ", round(mean(!is.na(geno_imputed))*100, digits=2), "% of the markers were obtained"))
    }

    ### Quality filter for low imputing rates / low read depth

    if(quali_filter){
      quali_filter1 <- which(rowMeans(!is.na(geno_imputed))<max_na)
      quali_filter2 <- which( (rowMeans(new_depth) / rowMeans(hb_depth)) < (min_depth*mean_depth))
      quali_filter3 <- sort(unique(c(quali_filter1, quali_filter2)))
      if(length(quali_filter3)>0){
        geno_imputed <- geno_imputed[-quali_filter3,]
        estimated_deletion <- estimated_deletion[-quali_filter3,]
        estimated_cnv <- estimated_cnv[-quali_filter3,]
        posi <- posi[-quali_filter3]
        hb_depth <- hb_depth[-quali_filter3,]
        new_depth <- new_depth[-quali_filter3,]
        allele <- allele[-quali_filter3,]
        snpname <- snpname[-quali_filter3]
        ref <- ref[-quali_filter3]
        alt <- alt[-quali_filter3]
        depth <- depth[-quali_filter3]
      }

    }

    ####################################
    ######## CNV - Calling #############
    ####################################

    if(estimate_sv){
      cnv_window <- sv_window
      cnv_cutoff <- sv_cut1
      del_cutoff <- sv_cut2

      estimated_depth <- new_depth / hb_depth / mean_depth
      per_marker_depth <- rowMeans(estimated_depth)
      population_mean <- stats::ksmooth(as.numeric(posi), per_marker_depth, bandwidth = cnv_window, x.points = as.numeric(posi))
      for(row in 1:ncol(estimated_cnv)){
        activ_depth <- stats::ksmooth(as.numeric(posi), estimated_depth[,row], bandwidth = cnv_window, x.points = as.numeric(posi))
        estimated_cnv[,row] <- (activ_depth$y/population_mean$y)>cnv_cutoff
        estimated_deletion[,row] <- (activ_depth$y/population_mean$y)<del_cutoff
      }

    }


    if(use_del){
      DEL <- rowSums(estimated_deletion)
      adddel <- which(DEL>(del_freq*ncol(geno)))
    } else{
      adddel <- NULL
    }
    if(use_cnv){
      if(!estimate_cnv){
        CNV <- rowSums(estimated_cnv)
      } else{
        CNV <- rowSums(estimated_cnv>cnv_min)
      }

      addcnv <- which(CNV>(cnv_freq*ncol(geno)))
    } else{
      addcnv <- NULL
    }

    map <- rbind(cbind(snpname, as.numeric(posi)) ,
                 cbind(paste0(snpname, "_DEL")[adddel], as.numeric(posi)[adddel] ),
                 cbind( paste0(snpname, "_CNV")[addcnv],  as.numeric(posi)[addcnv] ))

    ref <- c(ref, rep("A", length(adddel)+ length(addcnv)))
    alt <- c(alt, rep("C", length(adddel)+ length(addcnv)))

    if(hetero){
      geno_imputed <- rbind(geno_imputed, estimated_deletion[adddel,]*2, (estimated_cnv[addcnv,])*2)
    } else{
      geno_imputed <- rbind(geno_imputed, estimated_deletion[adddel,], (estimated_cnv[addcnv,]))
    }


    order <- sort(as.numeric(map[,2]), index.return=TRUE)$ix
    map <- map[order,]
    ref <- ref[order]
    alt <- alt[order]
    geno_imputed <- geno_imputed[order,]

    # Just technical stuff to avoid to markers on the same bp
    while(sum(diff(as.numeric(map[,2]))==0)>0){
      up <- which(diff(as.numeric(map[,2]))==0)+1
      map[up,2] <- as.numeric(map[up,2])+1
    }

    #######################################
    ### Write input-vcf-file for BEAGLE ###
    #######################################


    if(hetero){
      haplo <- geno_imputed1
    } else{
      haplo <- geno_imputed[,sort(rep(1:sum(!is_ref),2))]

    }


    if(hetero){
      haplo[haplo==2] <- 1
    }

    haplo[is.na(haplo)] <- "."

    vcfgeno <- matrix(paste0(haplo[,(1:(ncol(haplo)/2))*2], "/", haplo[,(1:(ncol(haplo)/2))*2-1]), ncol=ncol(haplo)/2)


    options(scipen=999)
    vcfgenofull <- cbind(chromo, map[,2], map[,1], ref, alt, ".", "PASS", ".", "GT", vcfgeno)
    vcfgenofull <- rbind(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", lines[!is_ref]),vcfgenofull)

    headerfile <- rbind(
      "##fileformat=VCFv4.2",
      gsub("-", "", paste0("##filedate=",  Sys.Date())),
      paste0("##source='HBimpute_v0.0'"),
      "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>"
    )

    utils::write.table(headerfile, file=path_prebeagle2, quote=FALSE, col.names = FALSE, row.names = FALSE)
    utils::write.table(vcfgenofull, file=path_prebeagle2, quote=FALSE, col.names = FALSE, row.names = FALSE, append = TRUE, sep="\t")

    #####################################
    ##### BEAGLE Imputation #############
    #####################################


    if(sum(is_ref)>0){
      if(ref_beagle){
        ref_path <- paste0(ref_temp, ".vcf.gz")
      } else{
        ref_path <- path_prebeagle4
      }

      beagle_commandline <- paste0("java -jar ", path_beaglejar," ne=", beagle_ne," gt=",path_prebeagle2," out=",out," nthreads=", beagle_core," ref=", ref_path, " impute=false")
      system(beagle_commandline)
    } else{
      beagle_commandline <- paste0("java -jar ", path_beaglejar," ne=", beagle_ne," gt=",path_prebeagle2," out=",out," nthreads=", beagle_core)
      system(beagle_commandline)
    }


    #####################################################################
    ##### Reference panel (currently mostly experimental testing!) ######
    #####################################################################

    if(length(ref_panel)>0){

      ## Written only for heterozygous lines!



      data_panel <- vcfR::read.vcfR(paste0(out, ".vcf.gz"))

      if(FALSE){
        data_panel1 <- vcfR::read.vcfR("chromo7_3vcf.vcf.gz")
        k1 <- which(duplicated(c(data_panel@fix[,2], data_panel1@fix[,2]))[-(1:length(data_panel@fix[,2]))])

        data_panel1@gt <- data_panel1@gt[k1,]
        data_panel1@fix <- data_panel1@fix[k1,]

      }


      data_ref <- vcfR::read.vcfR(ref_panel)

      g1_ref <- substr(data_ref@gt[,-1], start=1, stop=1)
      g2_ref <- substr(data_ref@gt[,-1], start=3, stop=3)
      storage.mode(g1_ref) <- "integer"
      storage.mode(g2_ref) <- "integer"



      g1_ref  <- g1_ref*2
      g2_ref <- g2_ref*2

      g1_panel <- substr(data_panel@gt[,-1], start=1, stop=1)
      g2_panel <- substr(data_panel@gt[,-1], start=3, stop=3)
      storage.mode(g1_panel) <- "integer"
      storage.mode(g2_panel) <- "integer"




      loc_panel <- as.numeric(data_panel@fix[,2])

      ref_overlap <- duplicated(c( loc_panel, as.numeric(data_ref@fix[,2])))[-(1:length(loc_panel))]
      panel_overlap <- duplicated(c( as.numeric(data_ref@fix[,2]), loc_panel))[-(1:nrow(data_ref))]

      switch_coding <- which(data_panel@fix[panel_overlap,4] != data_ref@fix[ref_overlap,4])

      g1_panel[switch_coding,] <- 1- g1_panel[switch_coding,]
      g2_panel[switch_coding,] <- 1- g2_panel[switch_coding,]
      data_panel@fix[switch_coding,4:5] <- data_panel@fix[switch_coding,5:4]

      g1_panel  <- g1_panel*2
      g2_panel <- g2_panel*2
      geno <- (g1_panel + g2_panel)/2

      haplo_panel <- cbind(g1_panel, g2_panel)
      haplo_panel <- haplo_panel[,sort(rep(1:ncol(g1_panel),2))+ c(0, ncol(g1_panel))]


      dhm_ref <- cbind(haplo_panel[panel_overlap,], g1_ref[ref_overlap,], g2_ref[ref_overlap,])

      blocklist_ref <- HaploBlocker::block_calculation(dhm_ref, bp=loc_panel[panel_overlap],
                                                       window_size=window_size,
                                                       target_coverage=target_coverage,
                                                       min_majorblock = min_majorblock)

      full_panel <- matrix(NA, ncol=ncol(haplo_panel), nrow=nrow(g1_ref))

      full_panel[ref_overlap,] <- haplo_panel[panel_overlap,]
      full_panel <- cbind(full_panel, g1_ref, g2_ref)
      full_depth <- !is.na(full_panel)

      storage.mode(full_depth) <- "integer"

      t <- HaploBlocker::coverage_test(blocklist_ref)
      se <- HaploBlocker::blocklist_startend(blocklist_ref, type="bp")
      t1 <- round(mean(t)*100, digits=2)
      le <- round(mean(se[,2]-se[,1])/1000000, digits=2)
      cat(paste0("Final haplotype library Reference with ", t1, "% coverage.\n"))
      cat(paste0("Avg. block length is ", le, " MB\n"))
      if(t1<80){
        warning(paste0("Coverage in haplotype library is potentially too low at ", t1, "%"))
      }
      if(le < 0.1){
        warning(paste0("Avg.block length is only ", le, " MB. Potential problems with haplotype library!"))
      }

      geno_imputed <- matrix(NA, nrow=nrow(g1_ref), ncol=ncol(geno)*2)
      new_depth <- estimated_cnv <- estimated_deletion <- matrix(0, nrow=nrow(g1_ref), ncol=ncol(geno)*2)

      loc_ref <- as.numeric(data_ref@fix[,2])
      mean_depth_ref <- 1

      cat("Start HB imputation:\n")
      pb <- utils::txtProgressBar(min = 0, max = ncol(g1_panel) * (1 + hetero), style = 3)

      nmax <- ncol(full_panel)
      for(nr in 1:(ncol(g1_panel)* (1 + hetero))){
        utils::setTxtProgressBar(pb, nr)
        se <- HaploBlocker::blocklist_startend(blocklist_ref, type="bp")
        include <- HaploBlocker::which.block(blocklist_ref, nr)
        se <- se[which(include>0),]
        if(length(se)==0){
          next
        }
        end_block  <- sort(unique(c(0,se[,1]-1, se[,2])))[-1]
        start_block <- c(1, end_block[1:(length(end_block)-1)]+1)

        for(index in 1:length(start_block)){
          take <- which(((loc_ref>=start_block[index])+(loc_ref<=end_block[index]))==2)
          indi <- HaploBlocker::which.indi(blocklist_ref, nr, start=start_block[index], end=end_block[index])
          indi <- indi[indi<=nmax]
          if(length(indi)==0){
            indi <- nr
          }
          indi <- c(indi, nr, nr, nr, nr)
          if(length(indi)>0 && length(take)>0){
            if(hetero){
              ana <- full_panel[take,indi, drop=FALSE]
            } else{
              ana <- geno[take,indi, drop=FALSE]
            }

            if(hetero){
              anad <- full_depth[take,indi, drop=FALSE]
            } else{
              anad <- depth[take,indi, drop=FALSE]
            }

            ana[is.na(ana)] <- -999

            zero <- rowSums((ana==0)*anad)
            two <- rowSums((ana==2)*anad)

            geno_imputed[take,nr][((zero>((two)*min_confi)))*(zero>min_confi)*(1:length(zero))*(zero<10|two<10)] <- 0
            geno_imputed[take,nr][((min_confi*(zero))<two) * (two>min_confi)* (1:length(two))*(zero<10|two<10)] <- 2

            if(estimate_del){
              mis <- rowSums((ana==(-999))) - 4 * (ana[,ncol(ana)]==(-999))
              p_zero <- exp(-mean_depth)
              total <- length(indi) - 4

              quali_na <- which(stats::pbinom(size = total, prob = p_zero, q=mis) > cutoff & total > 5)
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

      if(hetero){
        geno_imputed2 <- geno_imputed
        new_depth2 <- new_depth
        estimated_cnv2 <- estimated_cnv
        estimated_deletion2 <- estimated_deletion
        geno_imputed <- (geno_imputed[,1:ncol(geno)*2-1] + geno_imputed[,1:ncol(geno)*2])/2
        new_depth_temp <- new_depth[,1:ncol(geno)*2-1]
        new_depth <- new_depth[,1:ncol(geno)*2]
        new_depth[new_depth<new_depth_temp] <- new_depth_temp[new_depth<new_depth_temp]
        estimated_cnv <- (estimated_cnv[,1:ncol(geno)*2-1] + estimated_cnv[,1:ncol(geno)*2])/2
        estimated_deletion <- (estimated_deletion[,1:ncol(geno)*2-1] + estimated_deletion[,1:ncol(geno)*2])/2

        estimated_deletion[estimated_deletion>0] <- 1
      }

      if(overwrite_na){
        if(hetero){
          replaces <- which(geno_imputed[ref_overlap,][!is.na(geno[panel_overlap,])] != geno[panel_overlap,][!is.na(geno[panel_overlap,])] |
                              (is.na(geno_imputed[ref_overlap,]) & !is.na(geno[panel_overlap,])))
          reper <- geno[panel_overlap,][!is.na(geno[panel_overlap,])][replaces]
          reper[reper==1] <- 0
          geno_imputed2[ref_overlap,1:ncol(geno)*2-1][replaces] <- reper
          reper <- geno[panel_overlap,][!is.na(geno[panel_overlap,])][replaces]
          reper[reper==1] <- 2
          geno_imputed2[ref_overlap,1:ncol(geno)*2][replaces] <- reper
        }

        geno_imputed[ref_overlap,][!is.na(geno[panel_overlap,])] <- geno[panel_overlap,][!is.na(geno[panel_overlap,])]

      }

      if(hetero){
        cat(paste0("Genotype calls for ", round(mean(!is.na(geno_imputed2))*100, digits=2), "% of the markers were obtained"))
      } else{
        cat(paste0("Genotype calls for ", round(mean(!is.na(geno_imputed))*100, digits=2), "% of the markers were obtained"))
      }




      if(use_del){
        DEL <- rowSums(estimated_deletion)
        adddel <- which(DEL>(del_freq*ncol(geno_imputed)))
      } else{
        adddel <- NULL
      }
      if(use_cnv){
        CNV <- rowSums(estimated_cnv>cnv_min)
        addcnv <- which(CNV>(cnv_freq*ncol(geno_imputed)))
      } else{
        addcnv <- NULL
      }

      map <- rbind(cbind(paste0("SNP", 1:nrow(geno_imputed)), as.numeric(loc_ref)) ,
                   cbind(paste0("SNP", 1:nrow(geno_imputed), "_DEL")[adddel], as.numeric(loc_ref)[adddel] ),
                   cbind( paste0("SNP", 1:nrow(geno_imputed), "_CNV")[addcnv],  as.numeric(loc_ref)[addcnv] ))

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

      ref <- data_ref@fix[,4]
      alt <- data_ref@fix[,5]

      if(hetero){
        haplo <- geno_imputed2
      } else{
        haplo <- geno_imputed[,sort(rep(1:ncol(geno),2))]

      }
      haplo[haplo==2] <- 1
      haplo[is.na(haplo)] <- "."


      vcfgeno <- matrix(paste0(haplo[,(1:(ncol(haplo)/2))*2], "/", haplo[,(1:(ncol(haplo)/2))*2-1]), ncol=ncol(haplo)/2)


      options(scipen=999)
      vcfgenofull <- cbind(chromo, map[,2], map[,1], ref, alt, ".", "PASS", ".", "GT", vcfgeno)
      vcfgenofull <- rbind(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", lines),vcfgenofull)

      headerfile <- rbind(
        "##fileformat=VCFv4.2",
        gsub("-", "", paste0("##filedate=",  Sys.Date())),
        paste0("##source='HBimpute_v0.0'"),
        "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>"
      )

      utils::write.table(headerfile, file=path_prebeagle3, quote=FALSE, col.names = FALSE, row.names = FALSE)
      utils::write.table(vcfgenofull, file=path_prebeagle3, quote=FALSE, col.names = FALSE, row.names = FALSE, append = TRUE, sep="\t")

      beagle_commandline <- paste0("java -jar ", path_beaglejar," ne=", beagle_ne," gt=",path_prebeagle3," out=",paste0(out,"_ref1")," nthreads=", beagle_core)
      system(beagle_commandline)


    }


    if(extended_output){
      save(file=paste0(out, ".RData"), list=c("geno_imputed", "new_depth", "estimated_cnv", "estimated_deletion", "hb_depth"))
    }

    if(log != FALSE){
      sink(zz, append = TRUE, type = c("message"))
      warnings()
      sink(NULL)
      sink(NULL, type=c("message"))
    }
  }


}

















