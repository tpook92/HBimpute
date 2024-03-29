% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute.R
\name{impute}
\alias{impute}
\title{Breeding function}
\usage{
impute(
  vcf = NULL,
  vcf_ref = NULL,
  vcf_RData = NULL,
  vcf_ref_RData = NULL,
  log = NULL,
  chromo = NULL,
  out = "out",
  geno = NULL,
  depth = NULL,
  allele = NULL,
  lines = NULL,
  posi = NULL,
  hetero = FALSE,
  hb_data = NULL,
  hb_map = NULL,
  hb_maf = 0,
  activ_HB = TRUE,
  max_hetero = 0.01,
  hetero_is_missing = TRUE,
  remove_del = FALSE,
  beagle_core = 1,
  path_beaglejar = "beagle5.jar",
  estimate_del = NULL,
  estimate_cnv = NULL,
  use_cnv = FALSE,
  use_del = FALSE,
  del_freq = 0.1,
  cnv_min = 2,
  cnv_freq = 0.1,
  max_depth = 10L,
  maf = 0,
  zero_two_coding = FALSE,
  cutoff = 0.9999,
  beagle_ne = 10000,
  window_size = 20,
  target_coverage = NULL,
  min_majorblock = 5000,
  overwrite_call = FALSE,
  overwrite_call_min_depth = 1,
  overwrite_na = TRUE,
  overwrite_na_min_depth = 1,
  min_confi = 4,
  ref_panel = NULL,
  extended_output = FALSE,
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
  remove_lines = NULL,
  hb_sets = NULL
)
}
\arguments{
\item{vcf}{Path of the file to impute (alternatively provide geno, depth, allele, lines, posi manually)}

\item{vcf_ref}{Path of a reference panel}

\item{vcf_RData}{Path of an RData object to read instead of plain vcf-file}

\item{vcf_ref_RData}{Path of an RData object to read instead of plain vcf-file}

\item{log}{Path of a potential log file to be written (default: NULL will take vcf-file name)}

\item{chromo}{Manually provide a vector of chromosomes to perform imputation for (default: Process all available chromosomes)}

\item{out}{Path of the output file}

\item{geno}{genotype dataset}

\item{depth}{read-depth}

\item{allele}{allele names (default: major "A", minor "C")}

\item{lines}{line names}

\item{posi}{physical marker positions}

\item{hetero}{Set TRUE when imputing a heterozygous organism. The two haplotype of a individual are then processed separately}

\item{hb_data}{Dataset to derive the HB-library on (if not provided use dataset itself)}

\item{hb_map}{Map-file for hb_data}

\item{hb_maf}{Minor allele frequency filter for the Haplotype library}

\item{activ_HB}{Set to FALSE to not perform imputation via HBimpute (why do you use the package?!)}

\item{max_hetero}{Only active if hetero=FALSE - remove marker with share of heterozygosity large max_hetero (default: 0.05)}

\item{hetero_is_missing}{Only active if hetero=FALSE - If TRUE set of heterozygous markers to NA (default: TRUE)}

\item{remove_del}{Set TRUE to set estimated deletions to NA calls}

\item{beagle_core}{Number of threats using in BEAGLE 5 (default: 1)}

\item{path_beaglejar}{Directory of the BEAGLE jar (default: "beagle5.jar")}

\item{estimate_del}{Set TRUE to estimate deletion}

\item{estimate_cnv}{Set TRUE to estimate copy number variation}

\item{use_cnv}{Set TRUE to include identified copy number variation in output vcf}

\item{use_del}{Set TRUE to include identified deletions in output vcf}

\item{del_freq}{minimum frequency of a deletion to be included in the data}

\item{cnv_min}{Minimum estimated depth to be counted as a CNV in the output vcf}

\item{cnv_freq}{Minimum frequency of a CNV to be included in the output vcf}

\item{max_depth}{Maximum read depth in the calling (Everything higher will be set to this value)}

\item{maf}{Applied maf-filter (default: 0)}

\item{zero_two_coding}{Set to TRUE in case you are providing an input in geno with 0/2 coded genotypes (and 1 heterozygous!)}

\item{cutoff}{Old deletion algorithm binomial cutoff}

\item{beagle_ne}{Corresponds to ne parameter in BEAGLE}

\item{window_size}{HB parameter - size of each window in the algorithm (default: 20)}

\item{target_coverage}{HB parameter - target Coverage in the blocklist}

\item{min_majorblock}{HaploBlocker parameter}

\item{overwrite_call}{Set to TRUE to always preserve original SNP calls}

\item{overwrite_call_min_depth}{Minimum read depth to overwrite a call (default: 1)}

\item{overwrite_na}{Set to FALSE to not replace NA calls in HB with the original call}

\item{overwrite_na_min_depth}{Minimum read depth to overwrite a call (default: 1)}

\item{min_confi}{Proportion of read the need to support one variant to be called (default: 4 - coding 4:1 ratio)}

\item{ref_panel}{ref_panel}

\item{extended_output}{Set to TRUE for all relevant outputs}

\item{quali_filter}{Set TRUE to active quali filter after HBimpute step}

\item{max_na}{Cutoff for NA-filter after HBimpute step}

\item{min_depth}{Read-depth filter after HBimpute step}

\item{estimate_sv}{Set TRUE to activate structural variation detection}

\item{sv_cut1}{sv_cut1}

\item{sv_cut2}{sv_cut2}

\item{sv_window}{sv_window}

\item{hb_max}{Maximum number of marker to consider in the haplotype library}

\item{snpname}{SNP IDss}
}
\description{
Function to simulate a step in a breeding scheme
}
