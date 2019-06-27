#!/bin/bash

for i in {0..15}
do
 cp ~/scratch/tabula_muris/lane_1/10X_P7_${i}/outs/molecule_info.h5 ~/projects/rrg-hsn/rfarouni/index_hopping/data/novaseq_l1/input/10X_P7_${i}.h5
 cp ~/scratch/tabula_muris/lane_2/10X_P7_${i}/outs/molecule_info.h5 ~/projects/rrg-hsn/rfarouni/index_hopping/data/novaseq_l2/input/10X_P7_${i}.h5
done
