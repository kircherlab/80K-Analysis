#!/bin/bash
#SBATCH --job-name=ngn2_epi_marks_sampling_1000
#SBATCH --output=ngn2_epi_marks_sampling_1000_%j.out
#SBATCH --error=ngn2_epi_marks_sampling_1000_%j.err
#SBATCH --time=2-00:00:00            # 2 days
#SBATCH --cpus-per-task=65
#SBATCH --mem=100G                      # use all available memory on node
#SBATCH --mail-type=BEGIN,END,FAIL         # optional: email on end/failure
#SBATCH --mail-user=kilian.salomon@bih-charite.de  # optional: your email

source /home/kisa11/miniforge3/etc/profile.d/conda.sh
conda activate updated_r_kernel

# Run R script
Rscript /sc-projects/sc-proj-bih-reg-seqs/users/kisa11/projects/80k_region_enrichment_regioneReloaded/regioneReloaded_NGN2_with_and_without_epigentic_marks_sampling_results_enhancer_silencer_destinguished.R
