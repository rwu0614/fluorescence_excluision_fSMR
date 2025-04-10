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
    datasize = 20e4;   % establish a segment size (~32Mbytes)
    
    analysis_params.Baseline_rough_cutoff = -20; % default is -20 for populational fSMR experiments where light source is always on
    analysis_params.med_filt_length = 5; %full PMT data median filter window size, default 50
    analysis_params.moving_average_window_size = 5; %full PMT data moving-average filter window size, default 5    
    analysis_params.med_filt_window_size = 3*Peak_length ; % baseline median filter window size, sampling distance for extrapolating flat baseline   
    analysis_params.min_distance_btw_peaks = 50; % minimum distance between peaks, for identifying unique peaks
    analysis_params.uni_peak_range_ext = 30; % number of data points from each side of detection cutoff to be considered as part of the peak
    analysis_params.uni_peak_baseline_window_size = 200; % length of data points from each side of detection cutoff to compute the local baseline
    
    %calcien+annexin detection stragety: prioritize calcien
    analysis_params.detect_thresh_pmt(1) = 20; 
    analysis_params.detect_thresh_pmt(2) = 5; 
    analysis_params.detect_thresh_pmt(3) = 20;
    analysis_params.detect_thresh_pmt(4) = 20;
    analysis_params.detect_thresh_pmt(5) = 20;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
%%
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
baseline_left_slope = [];
baseline_right_slope = [];
baseline_left_height = [];
baseline_right_height = [];

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
        baseline_left_slope = [baseline_left_slope, seg_readout_pmt.baseline_left_slope{2}];
        baseline_right_slope = [baseline_right_slope, seg_readout_pmt.baseline_right_slope{2}];
        baseline_left_height = [baseline_left_height, seg_readout_pmt.baseline_left_height{2}];
        baseline_right_height = [baseline_right_height, seg_readout_pmt.baseline_right_height{2}];
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
readout_pmt.baseline_left_slope = baseline_left_slope;
readout_pmt.baseline_right_slope = baseline_right_slope;
readout_pmt.baseline_left_height = baseline_left_height;
readout_pmt.baseline_right_height =baseline_right_height;

%%
thresh_baselineDiff_over_sig = 0.08; % cutoff for left-right baseline height difference normalized by the signal amplitude
thresh_base_slope = 0.00004; % cutoff for left-right baseline slopes

all_cell_baselineDiff_over_sig = abs(baseline_left_height-baseline_right_height)./voltage_pmt2;
%all_cell_baselineDiff_over_sig = all_cell_baselineDiff_over_sig(all_cell_baselineDiff_over_sig<1);
base_diff_pass_ind =  find(all_cell_baselineDiff_over_sig < thresh_baselineDiff_over_sig);
base_leftslope_pass_ind = find(abs(baseline_left_slope) < thresh_base_slope);
base_rightslope_pass_ind = find(abs(baseline_right_slope) < thresh_base_slope);

cell_pass_ind = intersect(base_diff_pass_ind, intersect(base_leftslope_pass_ind,base_rightslope_pass_ind));

%%
histogram(all_cell_baselineDiff_over_sig,50)
scatter(all_cell_baselineDiff_over_sig,baseline_left_slope')

%%
ex_x = readout_pmt.time_of_detection(cell_pass_ind);
ex_x = ex_x-ex_x(1);

ex_y = readout_pmt.voltage_pmt2(cell_pass_ind);
% ex_y_base = medfilt1(ex_y,10);
% ex_y = ex_y-ex_y_base;
%cell_1_ind = find(ex_x>0&ex_x<200&ex_y>100);
%cell_1_ind = find(ex_x>200&ex_x<325&ex_y>0);
cell_1_ind = find(ex_x>0&ex_x<250&ex_y>60&ex_y<150);
%cell_1_ind = find(ex_x>2000&ex_x<2200&ex_y>100&ex_y<400);

cell_1_ind = cell_1_ind(1:1:end);
cell_1_x = ex_x(cell_1_ind);
cell_1_y = ex_y(cell_1_ind);

cell_1_y_base = medfilt1(cell_1_y,100);
cell_1_y_base(1) = cell_1_y_base(2);
cell_1_y_base = cell_1_y_base-cell_1_y_base(1);

cell_1_y_mod = cell_1_y;

scatter(ex_x,ex_y,'.');
hold on
scatter(cell_1_x,cell_1_y,'.');



cell_1_std = std(cell_1_y);
cell_1_average = mean(cell_1_y);
cell_1_cv = cell_1_std/cell_1_average;
disp(cell_1_cv*100)








%% Generate PMT readout output file
%format follows: [time of detection(computer real time), PacificBlue(mV),FITC(mV), PE(mV), APC(mV), Cy7(mV)]
output_pmt = [readout_pmt.time_of_detection(cell_pass_ind)',readout_pmt.voltage_pmt1(cell_pass_ind)',...
    readout_pmt.voltage_pmt2(cell_pass_ind)',readout_pmt.voltage_pmt3(cell_pass_ind)',readout_pmt.voltage_pmt4(cell_pass_ind)'...
    readout_pmt.voltage_pmt5(cell_pass_ind)'];

cd(input_info.pmt_dir(1))
out_file_name = ['readout_pmt_uncompensated_' sample_name '.csv'];
dlmwrite(out_file_name, output_pmt, 'delimiter', ',', 'precision', 25);
cd(currentFolder)

% generate analysis report
report_dir = [input_info.pmt_dir{1} '\' sample_name '_report\PMT_report\'];
mkdir(report_dir)
PMT_readout_report_v1(report_dir,input_info,input_info.pmt_dir(1),sample_name, analysis_params, readout_pmt);

