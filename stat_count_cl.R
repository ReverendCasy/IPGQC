#!/usr/bin/env Rscript

pacman::p_load(argparser,
               dplyr,
               ggtree,
               magrittr,
               shipunov,
               stringr)

####### Argument parsing #######################################
ap <- arg_parser("StatCount: a standalone script for deducing species attribution by Identical Protein Group composition")
ap <- add_argument(ap,
                   "--ipg_mapping",
                   help = "A file with IPG-to-species mapping",
                   type = "character",
                   short = "-m")
ap <- add_argument(ap,
                   "--all_ipgs",
                   help = "A single-column list of all IPGs stored in the dataset",
                   type = "character",
                   short = "-a")
ap <- add_argument(ap,
                   "--input",
                   help = "Raw IPGQC output file",
                   type = "character",
                   short = "-i")
ap <- add_argument(ap,
                   "--bootstrap",
                   help = "Launch boostrapped t-test (default: off)",
                   flag = TRUE,
                   short = "-b")
ap <- add_argument(ap,
                   "--bootstrap_rounds",
                   help = "Rounds of bootstrap to perform (default: 100)",
                   type = "numeric",
                   default = 100,
                   short = "-r")
argv <- parse_args(ap)
################################################################

####### Data preprocessing #####################################
pa_table <-  read.delim(argv$ipg_mapping,
                        sep = '\t',
                        header = FALSE,
                        stringsAsFactors = FALSE)
all_ipgs <- scan(argv$all_ipgs,
                 character())
species_ipgs <- readLines(argv$input)

total_proteins <- species_ipgs[1] %>% 
  str_split(' ') %>% 
  unlist() %>% 
  last() %>% 
  regexPipes::gsub('\\.', '') %>%
  as.numeric()

detected_ipgs <- species_ipgs[2] %>% 
  str_split(' ') %>% 
  unlist() %>% 
  last() %>% 
  regexPipes::gsub(':', '') %>% 
  as.numeric()

## Terminate the script... softly
if (detected_ipgs == 0) {
  stop("No IPGs were found. Exiting.")
}

ipgs_ids <-  species_ipgs %>% 
  .[3] %>% 
  str_split(',') %>% 
  unlist()

constructSpeciesTable <- function(specs) {
  species_table <- specs %>% 
    str_split('\t') %>% 
    list(Species = sapply(., function(x) x[1]),
         Count = sapply(., function(x) gsub('[)(]', '', x[2]) %>% 
                          as.numeric())) %>% 
    .[c('Species', 'Count')] %>% 
    as.data.frame()
  return(species_table)
}

species_table <- constructSpeciesTable(species_ipgs[5:length(species_ipgs)])

################################################################



##### Fisher test evalutation ##################################

## a function to count two-field one-sided Fisher test
ftest <- function(k, m, exp_spec, exp_non) {
  contingency <- matrix(c(k, m, exp_spec, exp_non), 2)
  fres <- fisher.test(contingency,
                      alternative = 'greater')
  return(fres$p.value)
}

## a master function for Fisher test
countFisher <- function(patable = pa_table,
                        spec_table = species_table,
                        all.ipgs = all_ipgs) {
  adjacency_list <- patable %>% 
    `colnames<-`(c('Species', 'IPGs')) %>% 
    mutate(IPGs = sapply(IPGs, 
                         function(x) str_count(x, 
                                               ',') + 1))
  spec_table$m <- sapply(spec_table$Count,
                         function(x) total_proteins - x)
  
  spec_table$n <- sapply(spec_table$Species,
                         function(x) patable[patable$V1 == x,2] %>% 
                           strsplit(',') %>% 
                           unlist() %>% 
                           length())
  spec_table$l <- sapply(spec_table$Species,
                         function(x) setdiff(all.ipgs,
                                             patable[patable$V1 == x, 2] %>% 
                                               strsplit(',') %>% 
                                               unlist()) %>% 
                           length())
  spec_table$exp_spec <- round(spec_table$n * total_proteins / length(all.ipgs))
  spec_table$exp_non <- round(spec_table$l * total_proteins / length(all.ipgs))
  p <- sapply(1:nrow(spec_table),
                         function(x) ftest(spec_table$Count[x], 
                                           spec_table$m[x], 
                                           spec_table$exp_spec[x], 
                                           spec_table$exp_non[x]))
  
  padj <- p.adjust(p, 
                   method = 'BH')
  return(list(Expected = spec_table$exp_spec,
              fisherp = p,
              fisherpadj = padj))
}
#################################################################


######### Bootstrap evaluation ##################################
## a function to subsample initial IPG set
subsample <- function(all_ipgs = all_ipgs,
                      spec_ipgs = species_ipgs,
                      num = argv$bootstrap_rounds) {
  subs <- sample(all_ipgs, num)
  return(length(intersect(subs, spec_ipgs)))
}

## a function to emulate bootstrapping
bootstrapEmulation <- function(n, func, ...) {
  opt <- list(...)
  return(replicate(n , do.call(func, opt)))
}

## a function to calculate one-sided t-test for the species 
tFunc <- function(all_ipgs = all_ipgs,
                  spec_ipgs,
                  num,
                  n_boot,
                  mu) {
  boot_vec <- bootstrapEmulation(n_boot, 
                                 subsample, 
                                 all_ipgs, 
                                 spec_ipgs, 
                                 num)
  return(t.test(boot_vec, 
                mu = mu, 
                alternative = 'less')$p.value)
}

################################################################

fisher_res <- countFisher()
species_table %<>% cbind(as.data.frame(fisher_res))
if (isTRUE(argv$bootstrap)) {
  boot_p <- sapply(1:nrow(species_table),
                     function(x) tFunc(all_ipgs = all_ipgs,
                                       spec_ipgs = pa_table[pa_table$V1 == species_table$Species[x],2] %>% 
                                         strsplit(',') %>% 
                                         unlist(),
                                       num = detected_ipgs,
                                       n_boot = argv$bootstrap_rounds,
                                       mu = species_table$Count[x]))
  boot_padj <- p.adjust(boot_p, 
                        method = 'BH')
  species_table %<>% cbind(list(bootp = boot_p,
                               bootpadj = boot_padj))
}

print(species_table)

