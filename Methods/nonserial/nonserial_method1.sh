#!/bin/bash

#SBATCH --account=def-cdlin
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mem=1gb

module load StdEnv/2020
module load r/4.2.1

Rscript nonserial_method1.R 
