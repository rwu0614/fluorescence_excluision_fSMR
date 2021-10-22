# fluorescence_excluision_fSMR
Analysis pipeline for raw fSMR data from fluorescence exclusion experiments

## Update on 2021-10-22
1) PMT_readout_combo.m file is compatible with and without fluorescence exclusion signals. User can use this to analyze samples pured labeled with positive fluorescence markers. For details please see openning description on thresholding in the .m file
2) Compensations from fluorescence exclusion channel to downstream channel are automated in the PMT_readou_combo.m so users do not need to add additional compensation due to fxm. ***User still need to use compensation_matrix.m and apply_compensation.m to compute and apply compensation factors from other positive fluorescence markers. 
