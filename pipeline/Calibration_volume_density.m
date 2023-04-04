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


%%
% figure('Position',[744,630,213.800000000000,250.200000000000],'Color','w')
% tiledlayout(2,1,'Padding','tight')
% 
% dot_color = [0.43921568627451         0.435294117647059         0.435294117647059];
% 
% nexttile
% scatter(sample.buoyant_mass_pg,sample.volume_fL,3,'filled','MarkerFaceColor',dot_color,'MarkerFaceAlpha',0.25)
% 
% set(gca,'XColor','none');
% ylim([400,2200])
% xlim([20,80])
% ylabel('Volume (fL)')
% nexttile
% scatter(sample.buoyant_mass_pg,sample.density_gcm3,3,'filled','MarkerFaceColor',dot_color,'MarkerFaceAlpha',0.25)
% 
% ylim([1.01,1.1])
% xlim([20,80])
% ylabel('Density (g/mL)')
% xlabel('Buoyant mass (pg)')
% legend("n="+string(height(sample)),'Location','southeast')
% legend box off
