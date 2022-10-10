# Zfish_recon
Contains scripts to analyze the reconstructions from the hindbrain dataset. 
-   These scripts were used for the pre-print https://www.biorxiv.org/content/10.1101/2020.10.28.359620v3. 
-   Requires the installation of the treestoolbox (https://www.mathworks.com/matlabcentral/fileexchange/68886-trees-toolbox/) developed by Hermann Cuntz. 
- Some hardcoded paths in the scrips might need to be updated before running the scripts.
- `\scripts` - contains all the scritps needed to generate the figures in the pre-print.
- `\matFiles` - Important .matFiles, are used in analysis and figures. 
- `\skeletons` - contains most if not all skeltons used for analysis in the paper. Any missing skeletons will be preserved in the private google drive folder used by the lab.
- `fig_mFiles` - contains some of the m files, that generate the figures from the manuscript.

## Terminology
-   `scw` is the commin extension for representing neurnal skeletons.
    -   `Root` - root node of the tree, typicall is the soma. For trees without somas in the reconstructed volume, the first plane of the 3D volume is considered as the root node. 
-   `Zbrain` - refers to the Zbrain atals for larval zebrafish developed by the Engert lab.