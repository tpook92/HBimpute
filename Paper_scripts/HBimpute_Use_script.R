# Data for chromosome 1 is available at
# https://github.com/tpook92/HBimpute/blob/master/MAZE_DH_chr1.vcf.gz

# All other data sets are available upon request
# Mail to torsten.pook@uni-goettingen.de



# HB-seq

library(HBimpute)
for(chromo in 1:10){

  impute(vcf = paste0("/mnt/ceph_space/shared_data/maze_bcf/MAZE_DH_chr",chromo,".vcf"),
         out = paste0("hbseq2_cnv_updatemaf_chr", chromo,"_tar095"),
         hetero=FALSE, # DH lines should be fully homozygous
         max_hetero = 0.01,
         chromo = chromo,
         min_majorblock=500, # This is a starting value for HaploBlocker that will automatically adapted to achieve the target coverage
         target_coverage = 0.95,
         overwrite_call = TRUE,
         overwrite_call_min_depth = 2,
         extended_output=TRUE,
         estimate_sv = TRUE, # SV detection mode (SV calls will affect final imputation in BEAGLE!)
         use_del = TRUE,
         use_cnv = TRUE,
         beagle_core = 10,
         maf=-1 # no MAF filtering
         )
}

# HB-array

for(chromo in 1:10){


  load(paste0("Genetic_Datasets/Batch3_KEPE/PE_DH_chromo",chromo,".RData")) # array data for the HaploBlocker haplotype library
  hb_data <- data

  hb_map <- read.table(paste0("Genetic_Datasets/Batch3_KEPE/PEchromo",chromo,"map.map"))[,4] ## PLINK-map file

  impute(vcf = paste0("/mnt/ceph_space/shared_data/maze_bcf/MAZE_DH_chr",chromo,".vcf"),
         out = paste0("hbarray2_cnv_updatemaf_chr", chromo,"_tar095"), hetero=FALSE,
         hb_data = hb_data,
         hb_map = hb_map,
         estimate_sv = TRUE,
         use_del = TRUE,
         use_cnv = TRUE,
         max_hetero = 0.01,
         chromo = chromo,
         min_majorblock=5000,
         target_coverage = 0.95,
         overwrite_call = TRUE,
         overwrite_call_min_depth = 2,
         extended_output=TRUE,
         beagle_core = 10,
         maf = -1)
}
