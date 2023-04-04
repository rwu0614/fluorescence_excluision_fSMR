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


%% Compute volume and density
fluid_density = input('\nInput fluid density (g/mL):'); %%1.01189056 for 10mg/ml dex-RPMIT
bead_density = input('\nInput bead density (g/mL):'); %%1.025185752 g/mL

calibration_factor = median(sample.buoyant_mass_pg./sample.vol_au)/(bead_density-fluid_density);


volume_calibration_result = table();
volume_calibration_result.calibration_factor_fLoverAU = calibration_factor;

%         cd(input_info.coulter_dir)
%         out_file_name = ['Calibration_factor_' coulter_sample_name '.txt'];
%         writetable(volume_calibration_result,out_file_name, 'delimiter', '\t');
%         cd(currentFolder)

cd(input_info.sample_dir)
out_file_name = ['Calibration_factor_from_34hydrogelbead' sample_name '.txt'];
writetable(volume_calibration_result,out_file_name, 'delimiter', '\t');
cd(currentFolder)
disp(strcat(out_file_name," ",string(calibration_factor)))


