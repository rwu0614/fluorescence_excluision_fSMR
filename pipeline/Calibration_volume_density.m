clear all
close all 
clc
currentFolder = pwd;
addpath('plotting_functions\');
addpath('analysis_functions\');
%%
% Input UI to grab path to a single readout paired sample txt file
fprintf('\nGetting paired sample...\n')
[input_info.sample_filename, input_info.sample_dir, ~] = uigetfile('../*.*','Select Readout_paired_[sample name].txt',' ');
sample_path = [input_info.sample_dir,'\',input_info.sample_filename];
opts = detectImportOptions(sample_path,'ReadVariableNames',true,'VariableNamingRule','preserve','Delimiter','\t');
sample = readtable(sample_path,opts);
name_split = strsplit(input_info.sample_dir,'\');   
sample_name = name_split{end-1};   
sample_name= strrep(sample_name,'_',' ');

% Input UI to grab path to volume calibration factor
fprintf('\nGetting volume calibration factor...\n')
[input_info.vol_filename, input_info.vol_dir, ~] = uigetfile('../*.*','Select Calibration_factor .txt',' ');
vol_path = [input_info.vol_dir,'\',input_info.vol_filename];
opts = detectImportOptions(vol_path,'ReadVariableNames',true,'VariableNamingRule','preserve','Delimiter','\t');
vol_calibration = readtable(vol_path,opts);
refine_calibration_factor=vol_calibration.calibration_factor_fLoverAU;
%% Compute volume and density
fluid_density = input('\nInput fluid density (g/mL):');
sample.volume_fL = sample.vol_au.*refine_calibration_factor;
sample.density_gcm3 = sample.buoyant_mass_pg./sample.volume_fL+fluid_density;

%%
cd(input_info.sample_dir)
out_file_name = ['Calibrated_readout_paired_' sample_name '.txt'];
writetable(sample,out_file_name, 'delimiter', '\t');
cd(currentFolder)













