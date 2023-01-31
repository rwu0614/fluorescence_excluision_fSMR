clear all
close all 
% clc
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

color_uninfect = [0.8500 0.3250 0.0980];

    vol_ceiling =(4*pi*(40./2).^3)/3; % for 40um filter
    figure('Position',[236,754.600000000000,990.600000000000,143.400000000000],'color','w')
    tiledlayout(1,2,'Padding','compact')
    nexttile
        histogram('BinEdges',[0;Vol_data.volume_fL]','BinCounts',Vol_data.count)
        % hold on
        % histogram(mock_sc_vol,'BinWidth',20)
        xlim([0,vol_ceiling])
        set(gca,'Xscale','log')
        xlabel('Volume (fL)')
        ylabel('Count')
        title('Coulter counter all size bins')
    nexttile
%         vol_range = input('\nCoulter volume range cutoff [ , ]:');
         vol_range = [400,40000];
        vol_cut_low = vol_range(1);
        vol_cut_high = vol_range(2);
        ind_vol_range = find(Vol_data.volume_fL>vol_cut_low & Vol_data.volume_fL<vol_cut_high);
        target_vol_data = Vol_data(ind_vol_range,:);
        histogram('BinEdges',[vol_cut_low;target_vol_data.volume_fL]','BinCounts',target_vol_data.count,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor',color_uninfect)
        xlim([0,vol_ceiling])
        set(gca,'Xscale','log')
        xlabel('Volume (fL)')
        ylabel('Count')
 

%
