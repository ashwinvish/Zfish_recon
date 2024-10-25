# Zfish_recon
Contains scripts to analyze the reconstructions from the hindbrain dataset. 
-   These scripts were used for the pre-print https://www.biorxiv.org/content/10.1101/2020.10.28.359620v3. 
-   Requires the installation of the treestoolbox v1.15 (https://www.mathworks.com/matlabcentral/fileexchange/68886-trees-toolbox/) developed by Hermann Cuntz. 
-  Required Statistics and Machine Learning provided by mathworks
- Some hardcoded paths in the scrips might need to be updated before running the scripts.
- `\scripts` - contains all the scripts needed to generate the figures in the pre-print.
		- Third party scripts provided are colorcet.m, adjustedRandIndex.m, compute_confusion.m, VOI.m, locallog,m
- `\matFiles` - Important .matFiles, are used in analysis and figures. 
		- -Run multiScaleModularity.m before plotMultiScaleModularity.m or modularityScaling_500_OM_0_01_cleaned.mat 
- `\skeletons` - contains most if not all skeletons used for analysis in the paper. Any missing skeletons will be preserved in the private google drive folder used by the lab.
- `fig_mFiles` - contains some of the m files, that generate the figures from the manuscript.

## Terminology
-   `scw` is the common extension for representing neuronal skeletons.
    -   `Root` - root node of the tree, typically is the soma. For trees without somas in the reconstructed volume, the first plane of the 3D volume is considered as the root node. 
-   `Zbrain` - refers to the Zbrain atals for larval zebrafish developed by the Engert lab.

# Large files
-	Some large files to run the funcuton `isRhombomere.m` are located here https://drive.google.com/drive/folders/1eYR6s-Lb6w0D4RkY8tmQwHO5b5oLUsO8?usp=drive_link
