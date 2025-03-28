# fluorescence_excluision_fxSMR
Analysis pipeline for raw fxSMR data from fluorescence exclusion experiments

## Update on 2021-10-22
1) PMT_readout_combo.m file is compatible with or without fluorescence exclusion signals. User can use this to analyze samples purely labeled with positive fluorescence markers. For details please see then openning description on thresholding in the .m file
2) Compensations from fluorescence exclusion channel to downstream channels are automated in  PMT_readout_combo.m so users do not need to add additional compensation cuased by fxm. **User still need to use compensation_matrix.m and apply_compensation.m to compute and apply compensation factors from other positive fluorescence markers.** 
