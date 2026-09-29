#!/bin/bash
#SBATCH --job-name=microplan
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=microplan-%j.log
# Illustrative requests, not benchmarked requirements; adjust for your site.
set -euo pipefail
source /path/to/environment/bin/activate
microplan plan --config /path/to/study.json --out "/path/to/results-${SLURM_JOB_ID}"
