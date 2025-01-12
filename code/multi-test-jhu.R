# Test on HDL GWAS data from GLGC
library(optparse)
library(parallel)
library(readr)
library(dplyr)
library(data.table)
library(scales)
library(caret) # for findCorrelation to train ensemble PRS
library(stringr)
library(inline)
library(bigreadr)
library(bigsnpr)

options(stringsAsFactors=F)
option_list = list(
  make_option("--PennPRS_path", action = "store", default = NA, type = "character",
              help="Path to the PennPRS folder [Required]"),
  make_option("--userID", action = "store", default = NA, type = "character",
              help="User account ID [Required]"),
  make_option("--submissionID", action = "store", default = NA, type = "character",
              help="Job ID [Required]"),
  make_option("--methods", action = "store", default = 'PROSPER', type = "character",
              help="Options: a subset of methods from PROSPER, PRS-CSx, and MUSSEL, divided by comma"),
  make_option("--trait", action = "store", default = NA, type = "character",
              help="trait name [Optional]"),
  make_option("--races", action = "store", default = NA, type = "character",
              help="Races of the training GWAS data. Options: a subset (need to have at least two) of: EUR (European), AFR (African),
              AMR (Mixed American, Hispanic/Latio), EAS (East Asian), or SAS (South Asian), divided by comma [Required]"),
  make_option("--LDrefpanel", action = "store", default = '1kg', type = "character",
              help="LD reference panel. Options: '1kg' (1000 Genomes Project Phase 3) or 'ukbb' (UK Biobank) [Optional]"),
  make_option("--k", action = "store", default = 2, type = "numeric",
              help = "k-fold Monte Carlo Cross Validation (MCCV) for PUMAS. Options: any integer greater than or equal to 2 [Optional]"),
  
  make_option("--partitions", action = "store", default = '0.8,0.2', type = "character",
              help="Partitions for PUMAS subsampling. 
              Format: '% training, % testing' (equal to 1 - % training), divided by comma [Optional]"),
  
  make_option("--ndelta", action = "store", default = 5, type = "integer",
              help="Number of candidate values of the shrinkage parameter in L2 regularization. Options: 
              any positive integer number [Optional]"),
  make_option("--nlambda", action = "store", default = 5, type = "integer",
              help="Number of different candidate values for lambda (shrinkage parameter in the L1 regularization.
              Options: any positive integer [Optional]"),
  make_option("--lambda.min.ratio", action = "store", default = 0.01, type = "numeric",
              help="Ratio between the lowest and highest candidate values of lambda.
              Options: any value in (0,1) [Optional]"),
  make_option("--Ll", action="store", default=5, type='integer',
              help="Length of path for the tuning parameter lambda in the PROSPER step [default: %default]"),
  make_option("--Lc", action="store", default=5, type='integer',
              help="Length of path for the tuning parameter c in the PROSPER step [default: %default]"),
  
  make_option("--phi", action = "store", default = '1e+00,1e-02,1e-04,1e-6', type = "character",
              help="Global shrinkage parameter phi. For GWAS with limited sample sizes (including most of the current disease GWAS), 
              fixing phi to 1e-2 (for highly polygenic traits) or 1e-4 (for less polygenic traits), or doing a small-scale grid search 
              (e.g., phi=1e-6, 1e-4, 1e-2, 1) to find the optimal phi value in the validation dataset often improves perdictive performance.
              Alternatively, phi can be learnt from the data using a fully Bayesian approach. This works well for polygenic traits 
              with very large GWAS sample sizes (hundreds of thousands of subjects), but is computationally too intensive and thus is not included in our pipeline.
              Default (default values in the PRS-CS algorithm, Nov 21, 2024 version): grid search with candidate values 1e+00,1e-02,1e-04,1e-6.
              Options: candidate values in (0,1], divided by comma [Optional]"),
  
  make_option("--verbose", action="store", default=1, type="integer",
              help="Print logfile? 0 = no; 1 = yes [default: %default]")
)
opt = parse_args(OptionParser(option_list=option_list))


# Input: ---------
userID = opt$userID
submissionID = opt$submissionID
methods = str_split(opt$methods,",")[[1]]
trait = opt$trait
races = str_split(opt$races,",")[[1]]
LDrefpanel = opt$LDrefpanel
# Parameters for subsampling
k = opt$k
# ----------------

# Example Manual Input (for testing the code): ---------
userID = 'jin'
submissionID = 'multi'
methods = c('PRS-CSx', 'MUSSEL') # We might add other multi-ancestry methods in the future, but we will test PROSPER for now
trait = 'dMRI1'
races = c('EUR','EAS') # ----------------

# Optional input parameters (for PUMAS subsampling):
partitions <- opt$partitions
PennPRS_path = opt$PennPRS_path
# ----------------
homedir = '/depot/feixue/data/PennPRS/Files/'
package = '/depot/feixue/data/PennPRS/software/PROSPER'
plink_path = '/depot/feixue/data/PennPRS/software/'
PRScs_path = paste0(PennPRS_path, 'software/PRScs/')
PRScsx_path = paste0(PennPRS_path, 'software/PRScsx/')


type = 'multi-ancestry' # 'multi-ancestry' or 'single-ancestry'
ld_path <- paste0(homedir)
PUMAS_path = paste0(homedir,'code/')
threads = 1
ld_path0 <- paste0(ld_path, 'LD_1kg/') # set to the /LD_1kg folder under /LD/
if (LDrefpanel == '1kg'){
  eval_ld_ref_path <- paste0(ld_path, '/1KGref_plinkfile/') # set to the /1KGref_plinkfile folder under /LD/
  path_precalLD <- paste0(ld_path, '/LDpred2_lassosum2_corr_1kg/') # set to the /LDpred2_lassosum2_corr_1kg folder under /LD/
}
NCORES = 5 # don't provide this option to users, o.w. it will become too complicated, especially if we add other multi-ancestry methods.
# ----------------

K = length(races)

if ('PROSPER' %in% methods){
  ndelta = opt$ndelta # number of candidate values of the shrinkage parameter in L2 regularization
  nlambda = opt$nlambda # number of different candidate values for lambda (shrinkage parameter in the L1 regularization). Default in lassosum2 pipeline: 30, which may lead to issues when using PUMAS subsampling to tune parameters
  lambda.min.ratio = opt$lambda.min.ratio # Ratio between the lowest and highest candidate values of lambda. Dandidate values in (0,Inf), divided by comma
  Ll = opt$Ll
  Lc = opt$Lc
}
if ('PRS-CSx' %in% methods){
  PHI = opt$phi; phi.vals = as.numeric(str_split(PHI,",")[[1]])
  N_THREADS = opt$N_THREADS
}
# Job name/ID: e.g., trait_race_method_userID_submissionID
jobID = paste(c(trait,paste0(races,collapse = '.'),paste0(methods,collapse = '.'),userID,submissionID), collapse = '_')
# Create a job-specific (trait, race, method, userID, jobID) directory to save all the outputs, set the working directory to this directory
workdir = paste0('/scratch/bell/bingxin/PennPRS/',jobID,'/')
suppressWarnings(dir.create(workdir))
setwd(workdir)

source(paste0(PUMAS_path, 'PennPRS_functions.R')) # please save the PennPRS_functions.R file to the /PUMAS/code/ directory


# please change this directory to the directory specified by jobID
gwas_path <- paste0(workdir, 'sumdata/')
output_path <- paste0(workdir, 'output/')
input_path <- paste0(workdir, 'input_for_eval/')
PennPRS_finalresults_path <- paste0(workdir, 'PennPRS_results/')
dir.create(gwas_path, showWarnings = F)
dir.create(ld_path, showWarnings = F)
dir.create(input_path, showWarnings = F)
dir.create(output_path, showWarnings = F)
dir.create(PennPRS_finalresults_path, showWarnings = F)
# Create a separate directory 'PRS_model_training/' to store input for training PRS models
prsdir0 = paste0(workdir, 'PRS_model_training/')
dir.create(prsdir0, showWarnings = F)

for (method in methods) prsdir = paste0(prsdir0, method,'/')
dir.create(prsdir, showWarnings = F)

summdata = paste0(prsdir, 'summdata/')
dir.create(summdata, showWarnings = F)
rscriptsdir = paste0(prsdir, 'rscripts/')
dir.create(rscriptsdir, showWarnings = F)
logfiledir = paste0(rscriptsdir, 'logfile/')
dir.create(logfiledir, showWarnings = F)
path_out = paste0(prsdir, 'output/')
dir.create(path_out, showWarnings = F)
path_out_lassosum2 = paste0(prsdir, 'output/lassosum2')
dir.create(path_out_lassosum2, showWarnings = F)
for (race in races) dir.create(paste0(path_out_lassosum2, '/', race), showWarnings = F)
path_out_PROSPER = paste0(prsdir, 'output/PROSPER/')
dir.create(path_out_PROSPER, showWarnings = F)
for (race in races) dir.create(paste0(path_out_PROSPER, '/', race), showWarnings = F)
output_path_eval = paste0(workdir, 'output_for_eval/')
dir.create(output_path_eval, showWarnings = F)


# copy the input GWAS summary data, {Ancestry}_{Trait}.txt, to the /sumdata/folder
for (race in races){
  trait_name = paste0(race,'_',trait)
  system(paste0('cp -r /depot/feixue/data/PennPRS/inputfiles/',trait_name,'.txt ', workdir, 'sumdata/'))
}

