#!/bin/bash -l
#SBATCH -J Gen_seed3000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=camilo.suarez@uni.lu
#SBATCH --nodes=1
###SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH -c 128
#SBATCH -p batch
#SBATCH -t 0-23:59:59
#SBATCH --qos=normal
#SBATCH -o %x-%j.log

### Submission info
echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"

### Loading modules for jobs 
module purge
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
module load toolchain/foss
module load tools/Singularity

###echo SRUN_CPUS_PER_TASK
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

job_name=${SLURM_JOB_NAME}
input_file=1000_07152024.py
###input_file=1000_fr_stl.py
scratch_dir=$PWD

### execution
srun singularity exec --bind ${scratch_dir} ContGeoGen_V4.sif python3 ${input_file}

echo "== Ending run at $(date)"
