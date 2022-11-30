clear 
close all
clc
format long;
currentFolder = pwd;
addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');
%%
input_dir = uigetdir('input directory');
input_dir = string([input_dir,'\']);
%%
 %fprintf('\nGetting SMR data...\n')
    input_info.smr_dir = input_dir;
    smr_sample_path = strsplit(input_info.smr_dir,'\');
    smr_sample_name = smr_sample_path(end-1);
    S = dir(fullfile(input_info.smr_dir ,smr_sample_name));
    input_info.smr_filename = S.name;
    smr_file_ID = fopen(strcat(input_info.smr_dir, input_info.smr_filename), 'r');
    %fprintf('\n%s selected for analysis\n', smr_filename)

    %fprintf('\nGetting Time data...\n')
    input_info.smr_time_dir = input_dir;
    S = dir(fullfile(input_info.smr_time_dir ,sprintf('*SMR_time*')));
    input_info.smr_time_filename = S.name;
    smr_time_file_ID  = fopen(strcat(input_info.smr_time_dir, input_info.smr_time_filename), 'r', 'b');
    %fprintf('\n%s selected for analysis\n', smr_time_filename)
    %grab sample name from input file name
    name_split = strsplit(input_info.smr_dir,'\');
    sample_name = name_split{end-1};
    sample_name= strrep(sample_name,'_',' '); % this variable is for annotating output filename

%THIS IS FOR FEEDBACK
rawdata_smr = fread(smr_file_ID, 'int32','b')/2^32;
rawdata_smr(1:129:end)=[];        
rawdata_smr = rawdata_smr';
rawdata_smr = diff(atan(rawdata_smr(1:2:end)./rawdata_smr(2:2:end)));                
rawdata_smr(rawdata_smr > pi/2) = rawdata_smr(rawdata_smr > pi/2) - pi;
rawdata_smr(rawdata_smr <= -pi/2) = rawdata_smr(rawdata_smr <= -pi/2) + pi;
rawdata_smr = rawdata_smr*10000/(2*pi);
fclose(smr_file_ID);

%Read time file. Removing the first element in the time file to initialize the 
%rawdata_smr_time array to be the same size of rawdata_smr
rawdata_smr_time = fread(smr_time_file_ID,'float64=>double');
rawdata_smr_time(1)=[];

%%
n_pmt_channel=5;
input_info.pmt_filename = strings(n_pmt_channel,1); %initialize PMT filenames
input_info.pmt_dir = strings(n_pmt_channel,1);%initialize PMT directory paths
pmt_file_ID = zeros(n_pmt_channel,1);%initialize PMT local file IDs
    %loop through to get each PMT local file ID
    for i = 1:n_pmt_channel
        %fprintf('\nGetting PMT Channel %d readout...', i)
        input_info.pmt_dir(i) = input_dir;
        target_pmt_channel = sprintf('*PMT_ch%i*.bin',i);
        S = dir(fullfile(input_info.pmt_dir(i),target_pmt_channel));
        input_info.pmt_filename(i) = S.name;
        %fprintf('\n%s selected for analysis\n', input_info.pmt_filename(i))
        pmt_file_ID(i) = fopen(strcat(input_info.pmt_dir(i), input_info.pmt_filename(i)), 'r', 'b'); %open PMT file and create a local file ID for extracting data downstream

    end

    %grabing time local file ID
    %fprintf('\nGetting time data...\n')
    input_info.time_dir = input_dir;
    S = dir(fullfile(input_info.time_dir,sprintf('*PMT_time*.bin')));
    input_info.time_filename = S.name;
    time_file_ID = fopen(strcat(input_info.time_dir, input_info.time_filename), 'r', 'b');
    %fprintf('\n%s selected for analysis\n', input_info.time_filename)
    
%%
% Input UI to grab path to a single readout paired sample txt file
fprintf('\nGetting paired sample...\n')
[input_info.sample_filename, input_info.sample_dir, ~] = uigetfile('../*.*','Select Readout_paired_[sample name].txt',' ');
sample_path = [input_info.sample_dir,'\',input_info.sample_filename];
opts = detectImportOptions(sample_path,'ReadVariableNames',true,'VariableNamingRule','preserve','Delimiter','\t');
sample = readtable(sample_path,opts);
%%
figure(1)
scatter(sample.buoyant_mass_pg,sample.vol_au./sample.buoyant_mass_pg,5,'filled')

dcm_obj = datacursormode();
set(dcm_obj,'UpdateFcn',@myupdatefcn);


%%
time_ind = 460;
time_pick = sample.real_time_sec(time_ind);

[~,smr_time_pick_ind] = min(abs(rawdata_smr_time-time_pick));
%%
n = 0;
pmt_seg_datasize = 2e4;
rawdata_pmt=[];
rawdata_time_pmt=[];
while(fseek(pmt_file_ID(2), n*8*pmt_seg_datasize, 'bof') == 0)
    % flip forward 8*analysis_params.datasize bytes repeatedly until file ends
    fseek(pmt_file_ID(2),n*8*pmt_seg_datasize, 'bof');
    fseek(time_file_ID,n*8*pmt_seg_datasize, 'bof');
    
    % read raw pmt and time file
    rawdata_pmt = fread(pmt_file_ID(2),pmt_seg_datasize,'float64=>double');
    rawdata_time_pmt = fread(time_file_ID,pmt_seg_datasize,'float64=>double');
    
    if rawdata_time_pmt(end)>time_pick
        break
    end
    n = n + 1;
end

[~,pmt_time_pick_ind] = min(abs(rawdata_time_pmt-time_pick));
%%
smr_ind_series = smr_time_pick_ind-300:smr_time_pick_ind+500;
pmt_ind_series = pmt_time_pick_ind-5:pmt_time_pick_ind+3000;
figure(2)
tiledlayout(2,1,'Padding','compact')
nexttile
plot(rawdata_smr_time(smr_ind_series),rawdata_smr(smr_ind_series))
nexttile
plot(pmt_ind_series,rawdata_pmt(pmt_ind_series))

%%
function txt = myupdatefcn(obj,event_obj)
pos = get(event_obj,'Position');
x = pos(1); y = pos(2); I = get(event_obj, 'DataIndex');
txt = {['X: ',num2str(x)],...
       ['Y: ',num2str(y)],...
       ['I: ',num2str(I)]};
%print(txt);
end










