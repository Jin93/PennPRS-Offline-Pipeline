# Example Manual Input (for testing the code): ---------
PennPRS_path = '/depot/feixue/data/PennPRS/'
homedir = '/scratch/bell/bingxin/PennPRS/'
input_GWAS_path = '/depot/feixue/data/PennPRS/inputfiles/'
submissionID = 'single_ans'
methods = c('C+T', 'PRS-CS')
trait = 'dMRI1'
race = 'EUR'
ensemble = TRUE


# ----------------

##### Single-Ancestry:
opt$PennPRS_path='/scratch/bell/bingxin/PennPRS/PennPRS_dir/'
opt$homedir='/scratch/bell/bingxin/PennPRS/'
opt$input_GWAS_path='/depot/feixue/data/PennPRS/inputfiles/'
opt$submissionID='single_ancestry_test'
opt$methods='C+T,lassosum2,LDpred2'
opt$trait='dMRI1'
opt$race='EUR'
opt$ensemble=TRUE
opt$NCORES=17



PennPRS_path='/scratch/bell/bingxin/PennPRS/PennPRS_dir/'
homedir='/scratch/bell/bingxin/PennPRS/'
input_GWAS_path='/depot/feixue/data/PennPRS/inputfiles/'
submissionID='single_ancestry_test'
methods='C+T,lassosum2,LDpred2'
trait='dMRI1'
race='EUR'
ensemble=TRUE
NCORES=17
nlambda=2
p_seq='1.0e-05,1.0e-04'

sbatch /depot/feixue/data/PennPRS/job_submission/Single-Ancestry.sh ${PennPRS_path} ${homedir} ${input_GWAS_path} ${submissionID} ${methods} ${trait} ${race} ${ensemble} ${NCORES} ${nlambda} ${p_seq}

sbatch -t 24:00:00 --nodes=1 --ntasks=1 --cpus-per-task=17 -A bingxin --mem=100G /depot/feixue/data/PennPRS/job_submission/Single-Ancestry.sh ${PennPRS_path} ${homedir} ${input_GWAS_path} ${submissionID} ${methods} ${trait} ${race} ${ensemble} ${NCORES}

#!/bin/bash

#SBATCH -A bingxin 
#SBATCH --time=24:00:00
#SBATCH --job-name=pumas-buxin
#SBATCH --mail-user=jin.jin@pennmedicine.upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=100g                # 50GB memory per job
#SBATCH --nodes=1                # One node per job
#SBATCH --cpus-per-task=17

# Load modules
module load r/4.3.1

PennPRS_path=$1
homedir=$2
input_GWAS_path=$3
submissionID=$4
methods=$5
trait=$6
race=$7
ensemble=$8
NCORES=$9
nlambda=${10}
p_seq=${11}

Rscript ${PennPRS_path}/code/single-ancestry-step1.R \
--PennPRS_path ${PennPRS_path} \
--homedir ${homedir} \
--input_GWAS_path ${input_GWAS_path} \
--submissionID ${submissionID} \
--methods ${methods} \
--trait ${trait} \
--race ${race} \
--ensemble ${ensemble} \
--NCORES ${NCORES} \
--nlambda "${nlambda}" \
--p_seq "${p_seq}" \

if [ $? -eq 0 ]; then
# If the first script succeeded, run the second script
Rscript ${PennPRS_path}code/single-ancestry-step2.R \
--PennPRS_path ${PennPRS_path} \
--homedir ${homedir} \
--submissionID ${submissionID} \
--methods ${methods} \
--trait ${trait} \
--race ${race} \
--ensemble ${ensemble}
--nlambda "${nlambda}" \
--p_seq "${p_seq}" \
else
  # If the first script failed, print an error message and do not proceed
  echo "Error: First script failed. Second script will not run."
fi










##### PRS-CS:
opt$PennPRS_path='/scratch/bell/bingxin/PennPRS/PennPRS_dir/'
opt$homedir='/scratch/bell/bingxin/PennPRS/'
opt$input_GWAS_path='/depot/feixue/data/PennPRS/inputfiles/'
opt$submissionID='prscstest'
opt$trait='dMRI1'
opt$race='EUR'
opt$ensemble=TRUE
opt$N_THREADS=11
opt$type='grid'
opt$phi='1e-04,1e-6'


PennPRS_path='/scratch/bell/bingxin/PennPRS/PennPRS_dir/'
homedir='/scratch/bell/bingxin/PennPRS/'
input_GWAS_path='/depot/feixue/data/PennPRS/inputfiles/'
submissionID='prscstest'
trait='dMRI1'
race='EUR'
N_THREADS=11
type='grid'
phi='1e-04,1e-6'

sbatch /depot/feixue/data/PennPRS/job_submission/PRS-CS.sh ${PennPRS_path} ${homedir} ${input_GWAS_path} ${submissionID} ${trait} ${race} ${N_THREADS} ${type} ${phi}


#!/bin/bash

#SBATCH -A highmem 
#SBATCH --time=24:00:00
#SBATCH --job-name=pumas-buxin
#SBATCH --mail-user=jin.jin@pennmedicine.upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=100g                # memory per job
#SBATCH --nodes=1                # One node per job
#SBATCH --cpus-per-task=11

# Load modules
module load r/4.1.2
module load anaconda

PennPRS_path=$1
homedir=$2
input_GWAS_path=$3
submissionID=$4
trait=$5
race=$6
N_THREADS=$7
type=$8
phi=$9

export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

Rscript ${PennPRS_path}/code/PRS-CS.R \
--PennPRS_path ${PennPRS_path} \
--homedir ${homedir} \
--input_GWAS_path ${input_GWAS_path} \
--submissionID ${submissionID} \
--trait ${trait} \
--race ${race} \
--N_THREADS ${N_THREADS} \
--type ${type} \
--phi ${phi} \








##### Tuning Parameter Free Methods:
opt$PennPRS_path='/scratch/bell/bingxin/PennPRS/PennPRS_dir/'
opt$homedir='/scratch/bell/bingxin/PennPRS/'
opt$input_GWAS_path='/depot/feixue/data/PennPRS/inputfiles/'
opt$submissionID='tuning_free'
opt$methods='LDpred2-auto,DBSLMM'
opt$trait='dMRI1'
opt$race='EUR'
opt$NCORES=17


PennPRS_path='/scratch/bell/bingxin/PennPRS/PennPRS_dir/'
homedir='/scratch/bell/bingxin/PennPRS/'
input_GWAS_path='/depot/feixue/data/PennPRS/inputfiles/'
submissionID='tuning_free'
methods='LDpred2-auto,DBSLMM'
trait='dMRI1'
race='EUR'
NCORES=17

sbatch -t 4:00:00 --nodes=1 --ntasks=1 --cpus-per-task=17 -A standby --mem=100G /depot/feixue/data/PennPRS/job_submission/Tuning-Parameter-Free.sh ${PennPRS_path} ${homedir} ${input_GWAS_path} ${userID} ${submissionID} ${methods} ${trait} ${race} ${NCORES}

sbatch /depot/feixue/data/PennPRS/job_submission/Tuning-Parameter-Free.sh ${PennPRS_path} ${homedir} ${input_GWAS_path} ${submissionID} ${methods} ${trait} ${race} ${NCORES}


#!/bin/bash

#SBATCH -A standby 
#SBATCH --time=4:00:00
#SBATCH --job-name=pumas-buxin
#SBATCH --mail-user=jin.jin@pennmedicine.upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=100g                # memory per job
#SBATCH --nodes=1                # One node per job
#SBATCH --cpus-per-task=17

# Load modules
module load r/4.1.2

PennPRS_path=$1
homedir=$2
input_GWAS_path=$3
submissionID=$4
methods=$5
trait=$6
race=$7
NCORES=$8

Rscript ${PennPRS_path}/code/Tuning-Parameter-Free.R \
--PennPRS_path ${PennPRS_path} \
--homedir ${homedir} \
--input_GWAS_path ${input_GWAS_path} \
--submissionID ${submissionID} \
--methods ${methods} \
--trait ${trait} \
--race ${race} \
--NCORES ${NCORES} \













##### Multi-ancestry:

#!/bin/bash

#SBATCH -A bingxin 
#SBATCH --time=24:00:00
#SBATCH --job-name=pumas-buxin
#SBATCH --mail-user=jin.jin@pennmedicine.upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=100g                # memory per job
#SBATCH --nodes=1                # One node per job
#SBATCH --cpus-per-task=11

export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

# Load modules
module load r/4.3.1

PennPRS_path=$1
homedir=$2
input_GWAS_path=$3
submissionID=$4
methods=$5
trait=$6
races=$7
NCORES=$8
phi=$9

Rscript ${PennPRS_path}/code/multi-ancestry-step1.R \
--PennPRS_path ${PennPRS_path} \
--homedir ${homedir} \
--input_GWAS_path ${input_GWAS_path} \
--submissionID ${submissionID} \
--methods ${methods} \
--trait ${trait} \
--races ${races} \
--N_THREADS ${N_THREADS} \
--phi "${phi}" \

if [ $? -eq 0 ]; then
# If the first script succeeded, run the second script
Rscript ${PennPRS_path}code/multi-ancestry-step2.R 
--PennPRS_path ${PennPRS_path} \
--homedir ${homedir} \
--submissionID ${submissionID} \
--methods ${methods} \
--trait ${trait} \
--races ${races} \
--N_THREADS ${N_THREADS} \
--phi "${phi}" \
else
  # If the first script failed, print an error message and do not proceed
  echo "Error: First script failed. Second script will not run."
fi


# PROSPER
opt$PennPRS_path='/scratch/bell/bingxin/PennPRS/PennPRS_dir/'
opt$homedir='/scratch/bell/bingxin/PennPRS/'
opt$input_GWAS_path='/depot/feixue/data/PennPRS/inputfiles/'
opt$submissionID='jin_multi'
opt$methods='PROSPER'
opt$trait='dMRI1'
opt$races='EUR,EAS'
opt$NCORES=15


PennPRS_path='/scratch/bell/bingxin/PennPRS/PennPRS_dir/'
homedir='/scratch/bell/bingxin/PennPRS/'
input_GWAS_path='/depot/feixue/data/PennPRS/inputfiles/'
submissionID='jin_multi'
methods='PROSPER'
trait='dMRI1'
races='EUR,EAS'
NCORES=11

sbatch /depot/feixue/data/PennPRS/job_submission/PROSPER.sh ${PennPRS_path} ${homedir} ${input_GWAS_path} ${submissionID} ${methods} ${trait} ${races} ${NCORES}


# PRS-CSx
opt$PennPRS_path='/scratch/bell/bingxin/PennPRS/PennPRS_dir/'
opt$homedir='/scratch/bell/bingxin/PennPRS/'
opt$input_GWAS_path='/depot/feixue/data/PennPRS/inputfiles/'
opt$submissionID='jin_multi'
opt$methods='PRS-CSx'
opt$trait='dMRI1'
opt$races='EUR,EAS'
opt$N_THREADS=10
opt$phi='1e-4,1e-6'


PennPRS_path='/scratch/bell/bingxin/PennPRS/PennPRS_dir/'
homedir='/scratch/bell/bingxin/PennPRS/'
input_GWAS_path='/depot/feixue/data/PennPRS/inputfiles/'
submissionID='jin_multi'
methods='PRS-CSx'
trait='dMRI1'
races='EUR,EAS'
N_THREADS=10
phi='1e-4,1e-6'

sbatch /depot/feixue/data/PennPRS/job_submission/PRS-CSx.sh ${PennPRS_path} ${homedir} ${input_GWAS_path} ${submissionID} ${methods} ${trait} ${race} ${N_THREADS} ${phi}



# MUSSEL:
opt$PennPRS_path='/scratch/bell/bingxin/PennPRS/PennPRS_dir/'
opt$homedir='/scratch/bell/bingxin/PennPRS/'
opt$input_GWAS_path='/depot/feixue/data/PennPRS/inputfiles/'
opt$submissionID='multijob1'
opt$methods='MUSSEL'
opt$trait='dMRI1'
opt$races='EUR,EAS'
opt$NCORES=17



PennPRS_path='/scratch/bell/bingxin/PennPRS/PennPRS_dir/'
homedir='/scratch/bell/bingxin/PennPRS/'
input_GWAS_path='/depot/feixue/data/PennPRS/inputfiles/'
submissionID='jin_multi'
methods='MUSSEL'
trait='dMRI1'
races='EUR,EAS'
NCORES=10


