## list of matfiles needed to generate figures

1. `slopes_CO_top500_unscaled_08062020.mat` --> K values from model. Model is at Values.ipynb (written by Alex.S). Left side of all figures with model in black

2. `slopes_many_CO_top500_unscaled_08062020.mat` --> 1000 variants of K values with randomized weights/ supplementary figure (synapse detection error) generated by k_values.ipynb (Alex S)

3. Slopes for potential connections. 
`slopes_CO_top500_unscaled_Pot_2.mat` --> 2um jitter
`slopes_CO_top500_unscaled_Pot_5.mat` --> 5um jitter
`slopes_CO_top500_unscaled_Pot_10.mat` --> 10um jitter

4. `slopesThresh_09102020.csv` --> experimental K values are loaded from, where column 4 - ABD; 5 - ABDi, 6 - Do, 7 - threshold; 8 - slopes

5. `fitCellsBasedOnLocationDO.mat` --> experimental results for DO neurons
	
6.`MatOrder_CO_top500_2blocks_gamma038_08062020.mat` -->  Neurons ordered according to Figure2 (center+periphery) Post Louvan Clustering

7. `ConnMatrix_CO_top500_2blocks_gamma038_08062020.mat` --> Connectivity matrix Figure2 (center+periphery) Post Louvan Clustering

8. `cellIDType_CO_top500_2blocks_gamma038_08062020.mat` --> Identity of cells from the above ordering (Axl - axial (modA), Int - Integrator (modO), DO - vest, vSPNs, ABD)   

9. ` df_cleaned.mat` --> Synapse edge list at data frame
		 `df_cleaned.mat` was not included due to its size>100MB. Please check the google drive for access.
10. `soma_loc_zbrainSpace_MatOrder_CO_top500_2blocks_gamma038_08062020.mat` --> soma locations of all neurons in the center, identified in `MatOrder_CO_top500_2blocks_gamma038_08062020.mat`. All coordinates are in zbrain atlas space.

11. `tformRough-LowEMtoHighEM-set2-Elavl3-Mnx-SB-set4.mat` --> transformation matrix to convert from neuroglancer space to z-brain atlas space.

12. `AllCells.mat` --> List of call cells reconstructed

13. `block1_cleaned_gamma038_08072020.mat` --> list of cells in block1 (Axial module) from the center matrix

14. `block2_cleaned_gamma038_08072020.mat` --> list of cells in block2 (ocular module) from the center matrix



