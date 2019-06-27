#!/bin/bash

for i in Normal Tumor;
do
 cp /home/rfarouni/scratch/hiseq4000_data/${i}/outs/molecule_info.h5 ~/projects/rrg-hsn/rfarouni/index_hopping/data/hiseq4000_nonplexed/input/${i}_nonplexed.h5
 cp /home/rfarouni/scratch/hiseq4000_data/${i}_Multiplexed/outs/molecule_info.h5 ~/projects/rrg-hsn/rfarouni/index_hopping/data/hiseq4000_plexed/input/${i}_plexed.h5
done
