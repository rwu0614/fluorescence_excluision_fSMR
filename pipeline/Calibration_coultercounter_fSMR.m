clear all
close all 
clc
currentFolder = pwd;
addpath('plotting_functions\');
addpath('analysis_functions\');


%% Input .#m4 file from Coulter counter, and a sample fSMR Readout_paired file
[input_info.coulter_filename, input_info.coulter_dir,~] = uigetfile('../*.*','Select Coulter Counter .#m4 File', ' ');
file = fullfile( input_info.coulter_dir, input_info.coulter_filename);
[tempDir, tempFile] = fileparts(file); 
status = copyfile(file, fullfile(tempDir, [tempFile, '.txt']));
Coulter_data = readtable(strcat(tempDir,"\", [tempFile, '.txt']),'Delimiter',' ');
name_split = strsplit(input_info.coulter_filename,'.');   
coulter_sample_name = name_split{end-1};   

% Input UI to grab path to a single readout paired sample txt file
fprintf('\nGetting paired sample...\n')
[input_info.sample_filename, input_info.sample_dir, ~] = uigetfile('../*.*','Select Readout_paired_[sample name].txt',' ');
sample_path = [input_info.sample_dir,'\',input_info.sample_filename];
opts = detectImportOptions(sample_path,'ReadVariableNames',true,'VariableNamingRule','preserve','Delimiter','\t');
sample = readtable(sample_path,opts);

%% Grabbing coulter volume data and creating mock single cell volume array
ind_size_bin_start  = find(Coulter_data.Var1 == "[#Bindiam]");
ind_size_bin_end  = find(Coulter_data.Var1 == "[Binunits]");
ind_count_bin_start  = find(Coulter_data.Var1 == "[#Binheight]");
ind_count_bin_end  = find(Coulter_data.Var1 == "[SizeStats]");

Vol_data = table();
Vol_data.diameter= str2double(string(Coulter_data.Var1(ind_size_bin_start+1:ind_size_bin_end-1)));
Vol_data.count= str2double(string(Coulter_data.Var1(ind_count_bin_start+1:ind_count_bin_end-1)));
Vol_data.volume_fL = (4*pi*(Vol_data.diameter./2).^3)/3;

mock_sc_vol = zeros(sum(Vol_data.count),1);
sc_count_start_ind = 1;
for i = 1:height(Vol_data)
    bin_count = Vol_data.count(i);
    mock_sc_vol(sc_count_start_ind:sc_count_start_ind+bin_count) = Vol_data.volume_fL(i);
    sc_count_start_ind = sc_count_start_ind+bin_count+1;
end
%% Set single cell volume range for calibration
vol_ceiling =(4*pi*(40./2).^3)/3; % for 40um filter
figure(1)
set(gcf, 'Position',[236,478,1500,420])
subplot(1,2,1)
histogram('BinEdges',[0;Vol_data.volume_fL]','BinCounts',Vol_data.count)
% hold on
% histogram(mock_sc_vol,'BinWidth',20)
xlim([0,vol_ceiling])
set(gca,'Xscale','log')
xlabel('Volume (fL)')
ylabel('Count')
title('Coulter counter all size bins')
subplot(1,2,2)
vol_range = input('\nCoulter volume range cutoff [ , ]:');
vol_cut_low = vol_range(1);
vol_cut_high = vol_range(2);
ind_vol_range = find(Vol_data.volume_fL>vol_cut_low & Vol_data.volume_fL<vol_cut_high);
target_vol_data = Vol_data(ind_vol_range,:);
subplot(1,2,2)
histogram('BinEdges',[0;Vol_data.volume_fL]','BinCounts',Vol_data.count)
hold on
histogram('BinEdges',[vol_cut_low;target_vol_data.volume_fL]','BinCounts',target_vol_data.count)
xlim([0,vol_ceiling])
set(gca,'Xscale','log')
xlabel('Volume (fL)')
ylabel('Count')
title('Coulter volume range for fSMR calibration')
legend(["All bins","Selected bins"])
%% Perform refined calibration

% Find a starting calibration factor by calibration fSMR to median cell
% volume from coulter counter data
target_vol_data.cum_count(1) = target_vol_data.count(1);
for i = 2:height(target_vol_data)
    target_vol_data.cum_count(i) = target_vol_data.cum_count(i-1)+target_vol_data.count(i);
end
target_vol_data.normalized_count = round(target_vol_data.count.*(height(sample)/target_vol_data.cum_count(end)))*1.2;

median_ind = median(1:1:target_vol_data.cum_count(end));
[~,vol_ind] = min(abs(target_vol_data.cum_count-median_ind));
median_volume_fL = target_vol_data.volume_fL(vol_ind);

sample_calibration = sample;

sample_median_vol_au = median(sample_calibration.vol_au);
pre_pass_calibration_factor = median_volume_fL/sample_median_vol_au;

% Refine calibration factor by sweeping and finding most optimal
% calibration factor defined by a scoring test
% Score is defined by the sum of difference between fSMR and coulter volume
% at each percentile from 5 to 95%
ind_mock_range = find(mock_sc_vol>vol_cut_low & mock_sc_vol<vol_cut_high);
mock_p_test_sc = mock_sc_vol(ind_mock_range);

mock_prt = prctile(mock_p_test_sc,5:1:95);

pre_pass_calibration_factor_test = pre_pass_calibration_factor-3:0.01:pre_pass_calibration_factor+3;
score = zeros(size(pre_pass_calibration_factor_test));
for i = 1:length(pre_pass_calibration_factor_test)
    pre_pass_sample_vol_fL = sample_calibration.vol_au.*pre_pass_calibration_factor_test(i);
    pre_pass_sample_prt = prctile(pre_pass_sample_vol_fL,5:1:95);
    score(i) = sum(abs(pre_pass_sample_prt-mock_prt));
end
[~,score_ind] = min(abs(pre_pass_calibration_factor_test-pre_pass_calibration_factor));
pre_pass_calibration_factor_score = score(score_ind);
[min_score,refine_ind] = min(score);
refine_calibration_factor = pre_pass_calibration_factor_test(refine_ind);

figure(2)
set(gcf, 'Position',[236,478,1500,420])
subplot(1,2,1)
    plot(pre_pass_calibration_factor_test,score)
    hold on
    plot(pre_pass_calibration_factor,pre_pass_calibration_factor_score,'.','MarkerSize',15)
    hold on
    plot(refine_calibration_factor,min_score,'.','MarkerSize',15)

    title("Calibration refinement")
    ylabel('Sum of difference by percentile')
    xlabel('Calibration factor (fL/vol_{au})')
    legend(["Score","Calibration by median","Refined calibration"])
subplot(1,2,2)
    pre_pass_sample_vol_fL = sample_calibration.vol_au.*refine_calibration_factor;
    histogram(pre_pass_sample_vol_fL,'BinEdges',[0;target_vol_data.volume_fL]')
    hold on
    histogram('BinEdges',[0;target_vol_data.volume_fL]','BinCounts',target_vol_data.normalized_count)
    xlim([vol_cut_low,vol_cut_high])
    xlabel('Volume (fL)')
    set(gca,'Xscale','log')
    legend(["fSMR volume with refined calibration","Coulter Counter"])

refinement_pass = input("\nDoes volume calibration looks good? Input 1 if yes:");

%%
if refinement_pass ==1
    volume_calibration_result = table();
    volume_calibration_result.calibration_factor_fLoverAU = refine_calibration_factor;
    
    cd(input_info.coulter_dir)
    out_file_name = ['Calibration_factor_' coulter_sample_name '.txt'];
    writetable(volume_calibration_result,out_file_name, 'delimiter', '\t');
    cd(currentFolder)
end




