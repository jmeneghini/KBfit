#!/bin/bash
#SBATCH --job-name=isosinglet_doublystrange_D200      # Job name
#SBATCH --output=output_%j.txt         # Output file (%j will be replaced with the job ID)
#SBATCH --error=error_%j.txt           # Error file
#SBATCH --time=02:00:00                # Time limit (hh:mm:ss)
#SBATCH --ntasks=1                     # Number of tasks (usually 1 for a single Python script)
#SBATCH --cpus-per-task=32              # Number of CPU cores per task

source /home/jmeneghini/.bashrc
shopt -s expand_aliases

cd /home/jmeneghini/KBfit/new_testing/test
KBfit ../input_print_Kpi.xml
