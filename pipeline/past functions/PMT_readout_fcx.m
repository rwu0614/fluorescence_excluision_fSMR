clear all
close all
clc
currentFolder = pwd;
addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimize using following parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Peak_length = 100; % estimated peak length
    datasize = 2e5;   % establish a segment size (~32Mbytes)
    
    analysis_params.Baseline_rough_cutoff = -20; % default is -20 for populational fSMR experiments where light source is always on
    analysis_params.med_filt_length = 5; %full PMT data median filter window size, default 50
    analysis_params.moving_average_window_size = 5; %full PMT data moving-average filter window size, default 5    
    analysis_params.med_filt_window_size = 3*Peak_length ; % baseline median filter window size, sampling distance for extrapolating flat baseline   
    analysis_params.min_distance_btw_peaks = 50; % minimum distance between peaks, for identifying unique peaks
    analysis_params.uni_peak_range_ext = 5; % number of data points from each side of detection cutoff to be considered as part of the peak
    analysis_params.uni_peak_baseline_window_size = 100; % length of data points from each side of detection cutoff to compute the local baseline
    
    % Below is for choosing which side of the baseline to use when analyzing
    % fluorescence exclusion signal. Fxm baseline is flow rate dependent so
    % there might be systemetic differences from on side to the other.
    % Recommendation: when fluorescence-detection region is close to SMR
    % cantiliver, choose the side where cell are travelling the fastest
    % (i.e. steady state flow)
    % left baseline -> 1
    % right baseline -> 2
    % average baseline from both side -> 3
    analysis_params.fxm_baseline_choice = 2;
    
    
    % This is for setting detection threshold for each PMT channel,
    % threshold is in the unit of standard deviation of baseline amplitude
    % i.e. noise level
    % *************** IMPORTANT******************
    % For fluorescence exclusion threshold, always use a negative value,
    % and still postitive threshold for downstream channels
    analysis_params.detect_thresh_pmt(1) = 10; 
    analysis_params.detect_thresh_pmt(2) = 10; 
    analysis_params.detect_thresh_pmt(3) = -10;
    analysis_params.detect_thresh_pmt(4) = 10;
    analysis_params.detect_thresh_pmt(5) = 10;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---------------- Run mode and upstream compensation determination ----------- %%
% To detemine if the user wants to analyze data including fluorescence
% exclusion or purly for positive labeling

fxm_channel = find(analysis_params.detect_thresh_pmt<0, 1);

if isempty(find(analysis_params.detect_thresh_pmt<0, 1))
    analysis_params.fxm_mode = 0;
else
    analysis_params.fxm_mode = 1;
end

% To detemine if the user wants to input compensation factor(s) from upstream
% channels to compensate for fxm channel

if fxm_mode == 1
    analysis_params.upstream_compen = input('Apply upstream compensation for fxm channel? Yes-1 No-0\n');
    if analysis_params. upstream_compen == 1
        compen_flag = 0;
        while compen_flag == 0
            compen_factor = input('Input compensation factor(s) as an array (array must have entries for every upstream channel of fxm channel):\n');
            if length(compen_factor) ~= find(analysis_params.detect_thresh_pmt<0, 1)-1
                warning('Length of the compensation array does not match with expected input length, please re-enter.')
            else
                compen_flag = 1;
            end
        end
    end
end

%% ------------------initialization-------------------%%

%writing global variables to connect with other .m files in the pipeline
global elapsed_time
global elapsed_index
global elapsed_peak_count;
elapsed_time=0;
elapsed_index=0;
elapsed_peak_count=0;


%Grabbing file path for PMT data
n_pmt_channel = 5; %number of pmt channels
input_info.pmt_filename = strings(n_pmt_channel,1); %initialize PMT filenames
input_info.pmt_dir = strings(n_pmt_channel,1);%initialize PMT directory paths
pmt_file_ID = zeros(n_pmt_channel,1);%initialize PMT local file IDs

%loop through to get each PMT local file ID
for i = 1:n_pmt_channel
    fprintf('\nGetting PMT Channel %d readout...', i)
    [input_info.pmt_filename(i), input_info.pmt_dir(i), exist] = uigetfile('../*.*','Select PMT Channel Data File', ' '); % exist variable =! 0 if a file is selected

    if(exist == 0)
        fprintf('Quitting analysis program now...')
        return
    else
        fprintf('\n%s selected for analysis\n', input_info.pmt_filename(i))
        pmt_file_ID(i) = fopen(strcat(input_info.pmt_dir(i), input_info.pmt_filename(i)), 'r', 'b'); %open PMT file and create a local file ID for extracting data downstream
    end
end

%grabing time local file ID
fprintf('\nGetting time data...\n')
[input_info.time_filename, input_info.time_dir, exist_time] = uigetfile('../*.*','Select time File',' ');
if(exist_time ~= 0)
    time_file_ID = fopen(strcat(input_info.time_dir, input_info.time_filename), 'r', 'b');
    fprintf('\n%s selected for analysis\n', input_info.time_filename)
else
    fprintf('Quitting analysis program now...\n')
    return
end
name_split = strsplit(input_info.time_dir,'\');
sample_name = name_split{end-1};
sample_name= strrep(sample_name,'_',' ');

% select analysis mode, initialize readout and display
analysis_mode = input('\nRapid analysis mode? (1 = Yes, 0 = No):    ');
if analysis_mode == 1
    disp_progress = input('\nDisplay progress? (1 = Yes, 0 = No):      ');
else
    disp_progress = 1;
end

if disp_progress==1
    %set analysis display screen size
    scrsize = get(0, 'Screensize');
    figure('OuterPosition',[0 0.05*scrsize(4) scrsize(3) 0.95*scrsize(4)])
end
disp_params=[analysis_mode, disp_progress];

%get how many segments to conduct analysis 
n = 1;
while(fseek(pmt_file_ID(1), n*8*datasize, 'bof') == 0)
    % flip forward 8*datasize bytes repeatedly until file ends
    n = n + 1;
end
num_segments = n - 1; % total number of segments = length of file in segments


%% Main analysis on looping data segments
n_out_features = 6; % total number of output features per detected signal, ie time, pmt, length, shape factors

rawdata_pmt = cell(1,n_pmt_channel);

time_of_detection  = [];
voltage_pmt1 = [];
voltage_pmt2 = [];
voltage_pmt3 = [];
voltage_pmt4 = [];
voltage_pmt5 = [];

segment_loop=0;
flag = 0;
while(flag==0)
    % seek data for current segement, datatype int is 8bytes
    for channel = 1:n_pmt_channel
        fseek(pmt_file_ID(channel),segment_loop*8*datasize, 'bof');
    end
    fseek(time_file_ID,segment_loop*8*datasize, 'bof');
    
    % read raw pmt and time file
    for channel = 1:n_pmt_channel
        rawdata_pmt{1,channel} = fread(pmt_file_ID(channel),datasize,'float64=>double');
    end
    rawdata_time_pmt = fread(time_file_ID,datasize,'float64=>double');
    
    [seg_readout_pmt] = P1_peakanalysis_pmt(segment_loop,num_segments,rawdata_pmt, rawdata_time_pmt, disp_params, analysis_params,input_info);

    if ~isempty(seg_readout_pmt.time_of_detection)
        time_of_detection = [time_of_detection, seg_readout_pmt.time_of_detection];
        voltage_pmt1= [voltage_pmt1, seg_readout_pmt.amplitude{1}];
        voltage_pmt2= [voltage_pmt2, seg_readout_pmt.amplitude{2}];
        voltage_pmt3= [voltage_pmt3, seg_readout_pmt.amplitude{3}];
        voltage_pmt4= [voltage_pmt4, seg_readout_pmt.amplitude{4}];
        voltage_pmt5= [voltage_pmt5, seg_readout_pmt.amplitude{5}];
    end
    
    segment_loop=segment_loop+1;
    
    if length(rawdata_pmt{1,1}) < datasize
        flag = 1;
    end
end

% remove first zero element, conversion from unit of V to mV
readout_pmt.time_of_detection= time_of_detection;
readout_pmt.voltage_pmt1 = voltage_pmt1*1000;
readout_pmt.voltage_pmt2 = voltage_pmt2*1000;
readout_pmt.voltage_pmt3 = voltage_pmt3*1000;
readout_pmt.voltage_pmt4 = voltage_pmt4*1000;
readout_pmt.voltage_pmt5 = voltage_pmt5*1000;

%% Quality check to remove low-quality signals


%% Apply compensation to fxm channel if needed by user


%% Apply compensation to fxm - downstream channel to remove effect from fxm spillover


%% Generate PMT readout output file
%format follows: [time of detection(computer real time), PacificBlue(mV),FITC(mV), PE(mV), APC(mV), Cy7(mV)]
output_pmt = [readout_pmt.time_of_detection',readout_pmt.voltage_pmt1',...
    readout_pmt.voltage_pmt2',readout_pmt.voltage_pmt3',readout_pmt.voltage_pmt4'...
    readout_pmt.voltage_pmt5'];

cd(input_info.pmt_dir(1))
out_file_name = ['readout_pmt_uncompensated_' sample_name '.csv'];
dlmwrite(out_file_name, output_pmt, 'delimiter', ',', 'precision', 25);
cd(currentFolder)

% generate analysis report
report_dir = [input_info.pmt_dir{1} '\' sample_name '_report\PMT_report\'];
mkdir(report_dir)
PMT_readout_report_v1(report_dir,input_info,input_info.pmt_dir(1),sample_name, analysis_params, readout_pmt);

