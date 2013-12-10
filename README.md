How to Build
------------


First, download and compile ANTS:

    git clone https://github.com/stnava/ANTs
    mkdir ANTS-build && cd ANTS-build
    ccmake ../ANTS
    make

If this doesn't work on your system then read the README.txt in the
`ANTS/` folder for directions.

Next, we can build the skullstripping-ants code:

    git clone https://github.com/pnlbwh/skullstripping-ants.git
    mkdir skullstripping-ants-build && cd skullstripping-ants-build
    ccmake .. -DITK_DIR=/path/to/ANTS-build/ITKv4-build/ -DANTS_BUILD=/path/to/ANTS-build -DANTS_SRC=/path/to/ANTS/trunk/
    make


How to Run
----------

You need to create 2 text files, which contain the paths to your training set: 

1) text file with list of structural images, e.g. `t1s.txt`
2) text file with list of structural image masks, e.g. `masks.txt`.  

Requirements:

- Each labelmap must have the same label number

Then you can compute the new mask:

    skullstripping-ants-build/multiAtlasANTS/mainANTSAtlasWeightedOutputProbability /path/to/structural_image /path/to/output_atlas_mask t1s.txt masks.txt

Finally, you will need to threshold the mask.  E.g. threshold a nrrd at 50:

    unu 2op gt /path/to/output_atlas_mask 50 | unu save -e gzip -f nrrd -o /path/to/output_mask


TODO
----

Add how to threshold a nifti in README
