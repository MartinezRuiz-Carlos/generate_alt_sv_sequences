#!/bin/bash
#SBATCH --job-name=test_sv_sequence   # Job name
#SBATCH --ntasks=1                    # Run six tasks
#SBATCH --mem=32G                     # Job Memory
#SBATCH --time=36:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs/sv_sequence_call_%j%t.log   # Standard output and error log

SEQUENCE_SIZE=90000
REFERENCE=/nemo/project/proj-brcax/data/assets/reference/t2t/hs1.fa
python generate_sv_sequence.py -i input/03-005-0130_03-5994_nl_sniffles.vcf.gz -c sniffles -r ${REFERENCE} -s ${SEQUENCE_SIZE} -o results/test/sniffles
# python generate_sv_sequence.py -i input/03-005-0130_03-5994_nl_svim.vcf.gz -c svim -r ${REFERENCE} -s ${SEQUENCE_SIZE} -o results/test/svim
python generate_sv_sequence.py -i input/03-005-0130_03-5994_nl_savana.vcf.gz -c savana -r ${REFERENCE} -s ${SEQUENCE_SIZE} -o results/test/savana