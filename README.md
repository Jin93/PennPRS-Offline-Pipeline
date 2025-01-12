# PennPRS Offline Pipeline 

[PennPRS](https://pennprs.org/) is a cloud-based platform dedicated to online PRS model training. On this Github page, we provide an offline version of PennPRS.

To use the tool, please follow the instructions in [Getting Started](#gettingstarted) to download the required files, then try our example code in [Example](#example). Please refer to the [paper](https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00095-8) or contact our team at pennprs@googlegroups.com for details.
</br>



## Version History
- [ ] __January 5, 2025:__  The PennPRS offline pipeline was made available on Github.
</br>



## Getting Started

1. Clone the Github repository by `git clone https://github.com/Jin93/PennPRS.git`. From now on, we will refer to the folder as `/PennPRS/`.

2. Create a folder `/LD/` under `/PennPRS/`. Download LD reference data files for different populations and save the uncompressed folder(s) in `/PennPRS/LD/`.

Each compressed LD data file contains the following folders:
  1. ~ 1.2 million [HapMap3](https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3) and [MEGA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5156387/) SNPs.
  2. For five ancestry groups: European (EUR), African/African America (AFR), Hispanic/Latino/Admixed American (AMR), East Asian (EAS), and South Asian (SAS).
  3. Based on the 1000 Genomes samples (recommended when GWAS training sample sizes are relatively small, e.g., N<sub>GWAS</sub> < 50K for all ancestry groups).





### 2. LD information generated based on the 1000 Genomes Project phase 3 samples 

- 498 EUR, 659 AFR, 347 AMR, 503 EAS, 487 SAS.

[EUR LD information](https://drive.google.com/file/d/1448w1cMdFuBmhnKEpamUBkGWX93WcaI-/view?usp=drive_link) (~26.33G): Google Drive File ID: `1448w1cMdFuBmhnKEpamUBkGWX93WcaI-`; decompress by `tar -zxvf EUR.tar.gz`

[AFR LD information](https://drive.google.com/file/d/12j9y0brD004HNakC7LoEot8RHFhToVAq/view?usp=drive_link) (~36.81G): Google Drive File ID: `12j9y0brD004HNakC7LoEot8RHFhToVAq`; decompress by `tar -zxvf AFR.tar.gz`

[AMR LD information](https://drive.google.com/file/d/1ItDxm8mllEAoaUyddd_wz8UTp0clNmLW/view?usp=drive_link) (~32.82G): Google Drive File ID: `1ItDxm8mllEAoaUyddd_wz8UTp0clNmLW`; decompress by `tar -zxvf AMR.tar.gz`

[EAS LD information](https://drive.google.com/file/d/1z4bztGLssmXmOiA4AHq4jzRQ5u9A-w7-/view?usp=drive_link) (~19.90G): Google Drive File ID: `1z4bztGLssmXmOiA4AHq4jzRQ5u9A-w7-`; decompress by `tar -zxvf EAS.tar.gz`

[SAS LD information](https://drive.google.com/file/d/1PCiy_rMDGxPpBZpQzRSVNocc9ZwoHt0-/view?usp=drive_link) (~21.61G): Google Drive File ID: `1PCiy_rMDGxPpBZpQzRSVNocc9ZwoHt0-`; decompress by `tar -zxvf SAS.tar.gz`

Each unzipped folder contains the following LD information that is needed for analysis of various supported PRS methods:

`./1KGref_plinkfile/`: raw LD reference genotype data.

`./LD_1kg/`: LD information by LD block required for subsampling.

`./LDpred2_lassosum2_corr_1kg/`: Pre-computed LD matrices and SNP information by LD block for implementing [lassosum2](https://privefl.github.io/bigsnpr/articles/LDpred2.html), [LDpred2](https://privefl.github.io/bigsnpr/articles/LDpred2.html), [PROSPER](https://github.com/Jingning-Zhang/PROSPER), and [MUSSEL](https://github.com/Jin93/MUSSEL).

`./LD/`: Pre-computed LD matrices and SNP information by LD block, which are input files for [MUSSEL](https://github.com/Jin93/MUSSEL).

`./map/`: SNP information (SNP ID, alleles) for mapping alleles between LD reference data and GWAS summary data for implementing [lassosum2](https://privefl.github.io/bigsnpr/articles/LDpred2.html), [LDpred2](https://privefl.github.io/bigsnpr/articles/LDpred2.html), [PROSPER](https://github.com/Jingning-Zhang/PROSPER), and [MUSSEL](https://github.com/Jin93/MUSSEL).



If PLINK or PLINK2 is not working, please follow the instructions for [PLINK1.9](https://www.cog-genomics.org/plink/) and [PLINK2](https://www.cog-genomics.org/plink/2.0/) to re-install them under `/PennPRS/software/`.

- Launch R and install required libraries:

``` r
install.packages(c('RISCA','optparse','bigreadr','bigsnpr','bigparallelr', 'bigmemory','stringr','caret','Rcpp', 'RcppArmadillo','RcppTN','inline','doMC','foreach','doParallel','data.table','readr','MASS','reshape','parallel',
'devtools','genio','dplyr','pryr','Matrix','lavaan','xtable','SuperLearner'))
```

3. Prepare and save the input GWAS summary data file(s) in `\${input_GWAS_path}` with file name `${ancestry}_${trait}.txt`.

Please refer to the `example_data` in [Examples](#examples) to prepare the input data files.

An example of the GWAS summary data format:
```
CHR	SNP	A1	A2	MAF	BETA	SE	P	N
1	rs3131969	A	G	0.129396	-0.00478692	0.0105404	0.64972	32586
1	rs2286139	C	T	0.136484	0.0018321	0.0105153	0.861684	31359
1	rs12562034	A	G	0.103985	-0.0016805	0.0114914	0.883732	32827
```
The following columns (and column names) are required for the GWAS summary data files:

 1. CHR: chromosome, 1, 2, ..., or 22.
 2. SNP: SNP RSID (in the format of rsXXXX). If only the position information is available, please impute RSID using reference genotype data of the same genome build.
 3. A1: effect allele (the allele which BETA corresponds to).
 4. A2: alternative/Other allele.
 5. BETA: SNP effect. For binary traits, beta is the log of odds ratio (logOR) from logistic regressions.
 6. SE: standard error of beta.
 7. P: p-value.
 8. N: GWAS sample size by SNP. For binary traits, it is the effective sample size: 4 / (1 / N_control + 1 / N_case); and for continuous traits, it is the total sample size.


Before running the MUSSEL pipeline, please consider applying the following quality control (QC) steps to the GWAS summary data:

 1. Only keep the biallelic HapMap3 + MEGA SNPs (SNP IDs can be found in the second column of `ref_bim.txt`) to avoid troubles caused by reading huge files (e.g., > 8 million SNPs) in R.
 2. Remove SNPs with minor allele frequencies (MAF) lower than 1% in all populations.
 3. Remove SNPs with very small GWAS sample sizes (e.g., < 90% of the total GWAS sample size). This step can be omitted if too many SNPs are removed.
 4. In the case where each training GWAS samples consist of multiple ancestry groups. Ideally, a customized LD reference dataset should be created with matched ancestral composition. Please contact for conducting suck a task.

</br>






# Use the pipeline 

The PennPRS offline pipeline provides three modules: 

[Single Ancestry Analysis with Pseudo Training](#option-1:single-ancestry-analysis-with-pseudo-training)

([LDpred2](#step-0:-run-ldpred2-by-ancestry))

([LDpred2](#step-0:-run-ldpred2-by-ancestry))

## Option 1: Single-Ancestry Analysis with Pseudo Training

    Rscript <path_to_single-ancestry.R> \
      --userID <userID> \
      --submissionID <submissionID> \
      --method <method> \
      --trait <trait> \
      --race <race> \
      --LDrefpanel <LDrefpanel> \
      --k <k> \
      --partitions <partitions> \
      --ndelta <ndelta> \
      --nlambda <nlambda> \
      --lambda.min.ratio <lambda.min.ratio> \
      --alpha <alpha> \
      --p_seq <p_seq> \
      --sparse <sparse> \
      --kb <kb> \
      --Pvalthr <Pvalthr> \
      --R2 <R2> \
      --ensemble <ensemble> \
      --verbose <verbose> \
      --temp_path <temp_path>


__Note:__ 

there are several command lines that may need to be customized by users because of discrepancies in server:

- The linux commands for submitting jobs on server (lines 88 - 94 and line 110: "sbatch --mem=30G" in `LDpred2_jobs.R`, and lines 120 - 126 and lines 145 - 146 (e.g., "sbatch --mem=35G"), in `MUSS_jobs.R`), may need to be modified according to the server used (e.g., "module load conda_R", "module load R/4.3", etc.).
- For the commands in line 110 in `LDpred2_jobs.R`: "sbatch --mem=25G", the memory is required for running LDpred2 sequentially on two ancestry groups with NCORES=11 cores. If more ancestries are included, a larger memory (e.g., 40G for $K=5$) may need to requested (mainly for loading the LD information). Similarly, for the commands in lines 145 - 146 in `MUSS_jobs.R`: the required memory should be customized according to the number of training ancestry groups ($K$). The default memory requirement in MUSS_jobs.R is for jointly modeling 2 ancestry groups. For modeling more ancestry groups, the requested memory can be estimated as a linear function of $K$.



  **Parameters:** <br>
--userID: User account ID. This is a unique identifier for the user running the job.<br>
--submissionID: Job ID. This is a unique identifier for the job submission.<br>
--methods: Specifies the methods to be used for analysis. Multiple methods can be provided as a comma-separated list.<br>
--trait: The name of the trait being analyzed.<br>
--race: Specifies the race of the individuals in the training GWAS. Available options are EUR (European), AFR (African), AMR (Mixed American, Hispanic/Latino), EAS (East Asian), or SAS (South Asian).<br>
--LDrefpanel: Specifies the linkage disequilibrium (LD) reference panel to be used. Options include '1kg' (1000 Genomes Project Phase 3) or 'ukbb' (UK Biobank).<br>
--k: The number of folds for k-fold Monte Carlo Cross Validation (MCCV) used in PUMAS. Must be an integer greater than or equal to 2.<br>
--partitions: Specifies the partitioning of the data for subsampling in PUMAS. The format should be '% training, % testing', where % testing equals 1 - % training.<br>
--delta: Candidate values for the shrinkage parameter in L2 regularization. Multiple values can be provided as a comma-separated list.<br>
--nlambda: The number of different candidate values for the lambda parameter (shrinkage parameter in L1 regularization).<br>
--lambda.min.ratio: Ratio between the lowest and highest candidate values of lambda.<br>
--alpha: A scaling factor for heritability, where H_2 = alpha * H_20. The default values recommended for the LDpred2 algorithm can be used.<br>
--p_seq: A sequence of candidate values for the proportion of causal SNPs. Values should be provided as a comma-separated list.<br>
--sparse: Whether to consider a sparse effect size distribution, where the majority of SNPs have effects shrunk to zero.<br>
--kb: SNPs within this range of the index SNP are considered for the Clumping (C) step in the analysis.<br>
--Pvalthr: Specifies the p-value thresholds for the Thresholding (T) step. Multiple values can be provided as a comma-separated list.<br>
--R2: SNPs with squared correlation higher than the specified R2 value with the index SNPs will be removed.<br>
--ensemble: Whether to train a weighted combination of the single PRS models.<br>
--verbose: Controls the verbosity of the output. Set to 0 for no output or 1 to print a log file.<br>
--temp_path: Specifies the path to the GWAS data after quality control.<br>



## Multi Ancestry

    Rscript /home/ubuntu/pennprs/backend/code/multi_ancestry.R \
      --userID <userID> \
      --submissionID <submissionID> \
      --method <method> \
      --trait <trait> \
      --race <race> \
      --LDrefpanel <LDrefpanel> \
      --k <k> \
      --partitions <partitions> \
      --ndelta <ndelta> \
      --nlambda <nlambda> \
      --lambda.min.ratio <lambda.min.ratio> \
      --Ll <Ll> \
      --Lc <Lc> \
      --verbose <verbose> \
      --temp_path <temp_path>

**Parameters:**
--userID: User account ID. This is a unique identifier for the user running the job.<br>
--submissionID: Job ID. This is a unique identifier for the job submission.<br>
--methods: Specifies the methods to be used for analysis. The default method is 'PROSPER'.<br>
--trait: The name of the trait being analyzed.<br>
--races: Specifies the races of the individuals in the training GWAS data. At least two races must be specified, separated by commas. Available options are EUR (European), AFR (African), AMR (Mixed American, Hispanic/Latino), EAS (East Asian), or SAS (South Asian).<br>
--LDrefpanel: Specifies the linkage disequilibrium (LD) reference panel to be used. Options include '1kg' (1000 Genomes Project Phase 3) or 'ukbb' (UK Biobank).<br>
--k: The number of folds for k-fold Monte Carlo Cross Validation (MCCV) used in PUMAS. Must be an integer greater than or equal to 2.<br>
--partitions: Specifies the partitioning of the data for subsampling in PUMAS. The format should be '% training, % testing', where % testing equals 1 - % training.<br>
--ndelta: The number of candidate values for the shrinkage parameter in L2 regularization. Must be a positive integer.<br>
--nlambda: The number of different candidate values for the lambda parameter (shrinkage parameter in L1 regularization). Must be a positive integer.<br>
--lambda.min.ratio: Ratio between the lowest and highest candidate values of lambda.<br>
--Ll: Specifies the length of the path for the tuning parameter lambda in the PROSPER step.<br>
--Lc: Specifies the length of the path for the tuning parameter c in the PROSPER step.<br>
--verbose: Controls the verbosity of the output. Set to 0 for no output or 1 to print a log file.<br>
--temp_path: Specifies the path to the GWAS data after quality control.<br>




**Data requirements：**

| CHR | A1 | A2 | SNP        | MAF       | BETA                  | SE        | P       | N     |
|-----|----|----|------------|-----------|-----------------------|-----------|---------|-------|
| 1   | G  | A  | rs12562034 | 0.085736  | 0.012104              | 0.051204  | 0.8131  | 16162 |
| 1   | C  | A  | rs4970383  | 0.367982  | 0.0032909             | 0.027646  | 0.9052  | 17079 |
| 1   | C  | T  | rs4475691  | 0.402124  | -0.0021203000000000003| 0.026989  | 0.9374  | 17079 |


**Required columns in the files:**

SNP: SNP ID in the format of rsXXXX.<br>
CHR: chromosome number in the format of 1,2,...,22.<br>
BETA: SNP effect. Note that for binary traits, beta is the coefficient in logistic regression, i.e. log(OR).<br>
SE: Standard error of beta.<br>
A1: effective allele (counted allele in regression).<br>
A2: alternative allele (non-A1 allele).<br>
MAF: Minor Allele Frequency.<br>
P: P value.<br>
N: Sample size per variant. Note that for binary traits, it is effective sample sizes = 4 / (1 / N_control + 1 / N_case); and for continuous traits, it is simply the sample size.<br>



# Output

The scripts single_ancestry.R will create directory 
$PATH_out/trait_race_method_userID_submissionID/ and writes trait.method.PRS.txt and PRS_model_training_info.txt.

The scripts multi_ancestry.R will create directory 
$PATH_out/trait_race_method_userID_submissionID/ and writes trait.method.PRS.txt and PRS_model_training_info.txt.








## PennPRS Manual



### Step 0: run LDpred2 by ancestry

This step is to obtain the estimated causal SNP proportion ($p_k, k=1,2,\ldots,K$) and heritability ($h^2_k, k=1,2,\ldots,K$) parameters in LDpred2 for each of $K$ training ancestry groups. These parameters will be used to specify the prior causal SNP proportions and heritability parameters in [MUSS](#step-1:-muss).


LDpred2_jobs.R: submit LDpred2 jobs by chromosome.
```r
LDpred2_jobs.R --PATH_package --PATH_data --PATH_LDref --PATH_out --FILE_sst --pop --bfile_tuning --NCORES
```

LDpred2_tuning.R: obtain estimated LDpred2 parameters.
```r
LDpred2_tuning.R --PATH_package --PATH_out --PATH_plink --FILE_sst --pop --chrom 1-22 --bfile_tuning --pheno_tuning --bfile_testing --pheno_testing --testing --NCORES
```

### Step 1: MUSS

MUSS: a Bayesian model that jointly models the GWAS summary data across all training populations to obtain a total of $L \times K$ PRS models under $L$ different tuning parameter settings for $Pr⁡(δ_{1j},…,δ_{Kj})$ (functions of $p_k$s) and $\rho_{k_1,k_2}$s across all $K$ training populations.

```r
MUSS_jobs.R --PATH_package --PATH_data --PATH_LDref --PATH_out --FILE_sst --pop --LDpred2_params --chrom --bfile_tuning --NCORES
```

### Step 2: MUSSEL

For each target population, apply the Super Learning (SL) algorithm (default base learners: elastic net regression, ridge regression, and linear regression) to train an “optimal” linear combination of the ($L \times K$) PRS models, which we call the MUSSEL PRS model, based on the tuning set of the target population. 

Optional: the prediction performance of the final MUSSEL PRS model can be reported on an independent testing set, if the testing set is provided as an input.

```r
MUSSEL.R --PATH_package --PATH_out --PATH_plink --FILE_sst --pop --chrom --bfile_tuning --pheno_tuning --bfile_testing --pheno_testing --testing --NCORES
```

- PATH_package (required): path to the directory where the downloaded files (decompressed) are saved.

- PATH_data (required): path to the directory where the training data by ancestry group are saved.

- PATH_LDref (required): path to the directory where the LD reference data by ancestry group and chromosome are saved.

- PATH_out (required): path to the output directory where the results are saved.

- PATH_plink (required): path to plink2.

- FILE_sst (required): paths followed by file names of the population-specific GWAS summary statistics, separated by comma. Required columns: chr, rsid, pos, a0, a1, beta, beta_se, n_eff.

- pop (required: populations of the GWAS samples, separated by comma.

- chrom (required): the chromosome on which the model is fitted, input in the format of 1-22 or 1,2,3. Default: 1-22

- p: candidate values for tuning parameter p (causal SNP proportion). Default:
1.0e-05, 1.8e-05, 3.2e-05, 5.6e-05, 1.0e-04, 1.8e-04, 3.2e-04, 5.6e-04, 1.0e-03, 1.8e-03, 3.2e-03, 5.6e-03, 1.0e-02, 1.8e-02, 3.2e-02, 5.6e-02, 1.0e-01, 1.8e-01, 3.2e-01, 5.6e-01, 1.0e+00.

- H2: candidate values for tuning parameter H2 (heritability = H2 * h2_est from LDSC). Default: 0.3, 0.7, 1, 1.4.

- sparse: whether to consider a sparse model: 0, 1, or 0,1. Default: 0.

- bfile_tuning (required): path to PLINK binary input file prefix (excluding ".bed"/".bim"/".fam"") for tuning, save by chromosome.

- pheno_tuning (optional): path to phenotype file (PLINK format) for tuning.

- covar_tuning (optional): path to quantitative covariates (PLINK format) for tuning.

- testing (required): whether to perform testing in seperate dataset. Default: F.

- trait_type (required): Type of phenotype, continuous or binary. Default: 'continuous'.

- bfile_testing (optional): path to PLINK binary input file prefix (.bed/.bim/.fam) for testing, save by chromosome.

- pheno_testing (optional): path to phenotype file (PLINK format) for testing.

- covar_testing (optional): path to quantitative covariates (PLINK format) for testing.

- verbose: how much chatter to print: 0=nothing; 1=minimal; 2=all. Default: 1.

- cleanup: cleanup temporary files or not. Default: T.

- NCORES: how many cores to use. (Default: 13 for LDpred2_jobs.R, 5 for MUSS_jobs.R, and 1 for LDpred2_tuning.R and MUSSEL.R)

- LDpred2_params (required): path to the directory where the tuned LDpred2 parameters (population-specific causal SNP proportions, heritability and whether or not a sparse model is used) are saved, separated by comma.

- cors_additional (optional): additional candidate values for tuning parameter: genetic correlation across ancestry groups, example: 3 groups with label 1,2,3, want to add two additional settings: cor_setting1(1,2),cor_setting1(1,3),cor_setting1(2,3);cor_setting2(1,2),cor_setting2(1,3),cor_setting2(2,3).

- ps_additional (optional): typically not necessary. Additional candidate values for tuning parameter: ancestry-specific causal SNP proportions, example: 3 groups with label 1,2,3, want to add two additional settings: p1_setting1,p2_setting1,p3_setting1;p1_setting2,p2_setting2,p3_setting2.

- SL_library (optional): the base learners implemented in SuperLearner, separated by comma. Default: SL.glmnet,SL.ridge,SL.lm.

- linear_score (optional): whether the trained linear models will be saved. If not, only the Super Learner model will be saved. Note: some models in SL_library are non-linear. In this case, linear score file cannot be generated.

- target_pop (required): Target population (used to save output).

</br>


## Example
Download [example data](https://www.dropbox.com/scl/fi/bne781g2qsq67p0r9y9cl/example.zip?rlkey=6aw5tnfpbnjc3ieq1ee4do2ay&dl=0), decompress it by `tar -zxvf example.tar.gz` and save the files under the directory ${path_example}. Download the 1000 Genomes reference data and save the decompressed files in ${path_LDref}. Create a new folder `path_out` (e.g., in this example, `/dcs04/nilanjan/data/jjin/MUSSEL/test`) to save the output. Run the example code below with your own data directories and check if the results/outputs (saved in ${path_out}) are consistent with the results/outputs here: [example results](https://www.dropbox.com/scl/fi/44z9rkku5nsc45yvdm6hk/MEBayesSL_example_data_results.zip?rlkey=rrmgjwu4at0nbr4spjj6836l5&dl=0).

```r 
module load R

package='/dcs04/nilanjan/data/jjin/MUSSEL'
path_data='/dcs04/nilanjan/data/jjin/example'
path_LDref='/dcs04/nilanjan/data/jjin/LD_1kg'
path_out='/dcs04/nilanjan/data/jjin/MUSSEL/test'
path_plink='/dcl01/chatterj/data/jin/software/plink2'
target_pop='EUR,AFR'
trait_type='continuous'
```
Note: load the R version for which the required R packages were installed, in this example, R Version 4.2.2 Patched (2023-03-01 r83924).


### Step 1: Run LDpred2 by chromosome (22 jobs, each for one chromosome). 

In each job, the algorithm will run under different tuning parameter settings in parallel.
Note: as side products, $K$ LDpred2 PRS models trained based on GWAS data for each ancestry group will be generated by this step.

``` r
Rscript ${package}/R/LDpred2_jobs.R \
--PATH_package ${package} \
--PATH_data ${path_data} \
--PATH_LDref ${path_LDref} \
--PATH_out ${path_out} \
--FILE_sst ${path_data}/summdata/EUR.txt,${path_data}/summdata/AFR.txt \
--pop EUR,AFR \
--bfile_tuning ${path_data}/sample_data/EUR/tuning_geno,${path_data}/sample_data/AFR/tuning_geno \
--NCORES 11

```




### Step 2: Obtain tuned parameters from LDpred2

Wait until all LDpred2 jobs are completed to run this step. Tuned LDpred2 effect size estimates and the optimal tuning parameters are saved in ${path_out}. This step also generates the LDpred2 PRS models for each ancestry group as by products, which will be saved in `${path_out}/LDpred2/{race}_LDpred2_beta.txt`.

``` r
Rscript ${package}/R/LDpred2_tuning.R \
--PATH_package ${package} \
--PATH_out ${path_out} \
--PATH_plink ${path_plink} \
--FILE_sst ${path_data}/summdata/EUR.txt,${path_data}/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 1-22 \
--bfile_tuning ${path_data}/sample_data/EUR/tuning_geno,${path_data}/sample_data/AFR/tuning_geno \
--pheno_tuning ${path_data}/sample_data/EUR/pheno.txt,${path_data}/sample_data/AFR/pheno.txt \
--covar_tuning ${path_data}/sample_data/EUR/covar.txt,${path_data}/sample_data/AFR/covar.txt \
--bfile_testing ${path_data}/sample_data/EUR/testing_geno,${path_data}/sample_data/AFR/testing_geno \
--pheno_testing ${path_data}/sample_data/EUR/pheno.txt,${path_data}/sample_data/AFR/pheno.txt \
--covar_testing ${path_data}/sample_data/EUR/covar.txt,${path_data}/sample_data/AFR/covar.txt \
--trait_type ${trait_type} \
--testing TRUE \
--NCORES 1

```
Note: 
 1. the `pheno.txt` dataset contains phenotype information for both tuning and testing individuals. 
 2. For continuous traits, instead of using the `--covar_tuning` and `--covar_testing` option to adjust for covariates, an alternative approach is to first regress the phenotype on covariates, then save the residual (adjusted phenotype data) in `pheno.txt`.
 3. To try customized tuning parameter settings, please specify the corresponding inputs in the .sh files.

### Step 3: Run MUSS by chromosome (submit 22 jobs simultaneously, each for one chromosome). 

``` r
Rscript ${package}/R/MUSS_jobs.R \
--PATH_package ${package} \
--PATH_data ${path_data} \
--PATH_LDref ${path_LDref} \
--PATH_out ${path_out} \
--FILE_sst ${path_data}/summdata/EUR.txt,${path_data}/summdata/AFR.txt \
--pop EUR,AFR \
--LDpred2_params ${path_out}/LDpred2/EUR_optim_params.txt,${path_out}/LDpred2/AFR_optim_params.txt \
--chrom 1-22 \
--bfile_tuning ${path_data}/sample_data/EUR/tuning_geno,${path_data}/sample_data/AFR/tuning_geno \
--NCORES 5

```

### Step 4: Combine PRS models generated under different parameter settings with a Super Learner (SL) algorithm to obtain the final ensembled MUSSEL PRS model. 

If a testing dataset is provided, the prediction $R^2$ or $AUC$ of the final MUSSEL PRS model will be reported on the testing set.

``` r
Rscript ${package}/R/MUSSEL.R \
--PATH_package ${package} \
--PATH_out ${path_out} \
--PATH_plink ${path_plink} \
--pop EUR,AFR \
--target_pop ${target_pop} \
--chrom 1-22 \
--bfile_tuning ${path_data}/sample_data/EUR/tuning_geno,${path_data}/sample_data/AFR/tuning_geno \
--pheno_tuning ${path_data}/sample_data/EUR/pheno.txt,${path_data}/sample_data/AFR/pheno.txt \
--covar_tuning ${path_data}/sample_data/EUR/covar.txt,${path_data}/sample_data/AFR/covar.txt \
--bfile_testing ${path_data}/sample_data/EUR/testing_geno,${path_data}/sample_data/AFR/testing_geno \
--pheno_testing ${path_data}/sample_data/EUR/pheno.txt,${path_data}/sample_data/AFR/pheno.txt \
--covar_testing ${path_data}/sample_data/EUR/covar.txt,${path_data}/sample_data/AFR/covar.txt \
--trait_type ${trait_type} \
--testing TRUE \--NCORES 1

```

<<<<<<< HEAD
   install.packages(c('optparse','bigreadr','readr','stringr', 'caret', 'SuperLearner', 'glmnet', 'MASS', 'Rcpp', 'RcppArmadillo', 'inline', 'doMC', ‘foreach'))



## Questions

Please report any issues on the Wiki page and contact pennprs@googlegroups.com if you have any questions. Our team will check the reported issues and questions weekly.


## Citation


=======
# PennPRS
Code used in PennPRS manuscript.
>>>>>>> 819840b12b93b6df89bcc3f5ef61d6e9e59d80eb
