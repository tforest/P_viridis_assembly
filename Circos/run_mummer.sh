#!/bin/sh

#SBATCH -J mummPic
#SBATCH -o mummPic.out
#SBATCH -e mummPic.err
#SBATCH --mem=84G
#SBATCH --mail-type=FAIL
#SBATCH -p workq
#SBATCH -t 04-00:00:00
#SBATCH -c 32

module load bioinfo/MUMmer/4.0.0rc1

nucmer -c 100 -p piCar_100 GCA_015227895.2_Caur_TTU_1.0_genomic.fna CHR.fa 
### map_coords.mum file DESCRIPTION
#"S1" and "E1" are "start" and "end" positions of the match on the first sequence. 
#"S2" and "E2" are "start" and "end" positions of the match on the 2nd sequence. 
#"LEN" is the length of the match, "% IDY" percent identity of the alignment. 
###

# Filter of length > 10000bp
delta-filter -m -i 80 -l 10000 piCar_100.delta > piCar_filtered.delta
show-coords -l -q -T piCar_filtered.delta > map_coords_L10000_80i.mum

