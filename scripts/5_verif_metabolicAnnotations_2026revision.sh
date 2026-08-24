#!/bin/bash

module load muscle/
muscle -align /shared/projects/feob_iron_oxidizing_bacteria_/META/ANALYSIS_MG/5_MAG_treatment/metabolic_verifications/LSHF-pmoAcompleteCluster_and_amoA_pmoA_References.faa -output pmoAcompleteClust_w_amoA_pmoA_Refs.afa

module load raxml-ng/
raxml-ng --msa pmoAcompleteClust_w_amoA_pmoA_Refs.afa --model LG+G --prefix pmoAcpltClust_w_amoA_pmoA_Refs --seed 12345 --bs-trees 1000 --threads 6 --all

#module load iqtree/
#iqtree3 -s pmoAcompleteClust_w_amoA_pmoA_Refs.afa -m MFP -B 1000 -T 8 --alrt 1000 --seed 535238
