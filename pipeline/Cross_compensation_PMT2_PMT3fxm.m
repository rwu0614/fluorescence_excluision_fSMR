clear 
close all
format long;
currentFolder = pwd;
addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');
%%
% Input UI to grab path to fullflex pmt readout fullflex txt file
fprintf('\nGetting fullflex pmt fullflex...\n')
[input_info.fullflex_filename, input_info.fullflex_dir, ~] = uigetfile('../*.*','Select Readout_pmt_fxm_fullflex_[sample name].txt',' ');
fullflex_path = [input_info.fullflex_dir,'\',input_info.fullflex_filename];
opts = detectImportOptions(fullflex_path,'ReadVariableNames',true,'VariableNamingRule','preserve','Delimiter','\t');
fullflex = readtable(fullflex_path,opts);
name_split = strsplit(input_info.fullflex_dir,'\');   
sample = name_split{end-1};   
sample_name= strrep(sample,'_',' ');
%%

% Set pmt2 to pmt3 spillover factor
s = 0.038;

amp_detect_ch2 = fullflex.pmt2_amp_mV;
amp_detect_ch3 = fullflex.pmt3_amp_mV;
base_detect_ch2 = fullflex.pmt2_baseline_mV;
base_detect_ch3 = fullflex.pmt3_baseline_mV;

amp_true_ch2 = (base_detect_ch3.*amp_detect_ch2-base_detect_ch2.*amp_detect_ch3)./(base_detect_ch3-s*base_detect_ch2);
amp_true_ch3 = amp_detect_ch3-s*amp_true_ch2;


pmt_out = fullflex;

pmt_out.vol_au = abs(amp_true_ch3-base_detect_ch3)./base_detect_ch3;
pmt_out.pmt1_mV = abs(pmt_out.pmt1_amp_mV - pmt_out.pmt1_baseline_mV);
pmt_out.pmt2_mV = amp_true_ch2;
pmt_out.pmt4_mV = abs(fullflex.pmt4_amp_mV-fullflex.pmt4_baseline_mV.*...
             (amp_true_ch3./base_detect_ch3));
pmt_out.pmt5_mV = abs(fullflex.pmt5_amp_mV-fullflex.pmt5_baseline_mV.*...
             (amp_true_ch3./base_detect_ch3));



figure(1)
scatter(pmt_out.pmt2_mV,(pmt_out.buoyant_mass_pg./(abs(amp_detect_ch3-base_detect_ch3)./base_detect_ch3)),5,'filled')
hold on
scatter(pmt_out.pmt2_mV,(pmt_out.buoyant_mass_pg./pmt_out.vol_au),5,'filled')
symlog()

figure(2)
scatter(pmt_out.pmt2_mV,(pmt_out.buoyant_mass_pg./(abs(pmt_out.pmt5_amp_mV-pmt_out.pmt5_baseline_mV)./pmt_out.pmt5_baseline_mV)),5,'filled')
symlog()
%%
figure(2)
scatter(abs(pmt_out.pmt2_amp_mV-pmt_out.pmt2_baseline_mV),pmt_out.vol_au,5,'filled')
hold on
scatter(pmt_out.pmt2_mV,pmt_out.vol_au,5,'filled')
symlog()
%%
figure(3)
scatter(abs(pmt_out.pmt2_amp_mV-pmt_out.pmt2_baseline_mV),pmt_out.vol_au,5,'filled')
hold on
scatter(pmt_out.pmt2_mV,pmt_out.vol_au,5,'filled')
symlog()
















