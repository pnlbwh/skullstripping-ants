#!/bin/sh

input=input.nrrd
output=output_mask.nrrd

/projects/schiz/software/skullstripping-ants/build/multiAtlasANTS/mainANTSAtlasWeightedOutputProbability $input $output t1s.txt masks.txt
