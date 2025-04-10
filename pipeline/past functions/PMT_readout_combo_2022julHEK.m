function pass_rate = PMT_readout_combo_2022julHEK(varargin)
% input should be empty for one-time execution on a single sample. For
% batch analysis, use PMT_readout_combo(input_directory,batch_analysis_params)
% batch_analysis_params should have the same field as the analysis_params
% structure below

%% ------------------initialization-------------------%%
if nargin ==0
    clear 
    close all
    clc
end
currentFolder = pwd;
addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');

%writing global variables to connect with other .m files in the pipeline
global elapsed_time
global elapsed_index
global elapsed_peak_count;
elapsed_time=0;
elapsed_index=0;
elapsed_peak_count=0;


%%
if nargin==0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimize using following parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    analysis_params.Peak_length = 100; % estimated peak length
    analysis_params.datasize = 2e4;   % establish a segment size (~32Mbytes)
    
    analysis_params.Baseline_rough_cutoff = -20; % default is -20 for populational fSMR experiments where light source is always on
    analysis_params.med_filt_length = 5; %full PMT data median filter window size, default 50
    analysis_params.moving_average_window_size = 5; %full PMT data moving-average filter window size, default 5    
    analysis_params.med_filt_window_size = 3*analysis_params.Peak_length ; % baseline median filter window size, sampling distance for extrapolating flat baseline   
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
    analysis_params.detect_thresh_pmt(2) = 3; 
    analysis_params.detect_thresh_pmt(3) = -3;
    analysis_params.detect_thresh_pmt(4) = 10;
    analysis_params.detect_thresh_pmt(5) = 10;
    
    % for upstream compensation
    analysis_params. upstream_compen = 0; % 0- no compensation from upstream channel of fxm channel to initialize
    % For signal QC filtering
    analysis_params.thresh_baselineDiff_over_sig = 0.05; % cutoff for left-right baseline height difference normalized by the signal amplitude
    analysis_params.thresh_base_slope = 2*10^-3; % cutoff for left-right baseline slopes
    analysis_params.thresh_base_height_range = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    input_dir = varargin{1};
    batch_params = varargin{2};
    analysis_params = batch_params;
end
%% ---------------- Run mode and upstream compensation determination ----------- %%
% To detemine if the user wants to analyze data including fluorescence
% exclusion or purly for positive labeling

fxm_channel = find(analysis_params.detect_thresh_pmt<0, 5);

if isempty(find(analysis_params.detect_thresh_pmt<0, 5))
    analysis_params.fxm_mode = 0;
else
    analysis_params.fxm_mode = 1;
    analysis_params.fxm_channel = fxm_channel;
end
fxm_mode = analysis_params.fxm_mode;

if nargin ==0
% To detemine if the user wants to input compensation factor(s) from upstream
% channels to compensate for fxm channel
    if fxm_mode == 1
        fprintf('Fluorescence exclusion (fxm) analysis mode entered \n')
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
else
    if fxm_mode == 1
        if analysis_params. upstream_compen == 1
            compen_flag = 0;
            while compen_flag == 0
                compen_factor = batch_params.compen_factor;
                if length(compen_factor) ~= find(analysis_params.detect_thresh_pmt<0, 1)-1
                    warning('Length of the compensation array does not match with expected input length, please re-enter.')
                else
                    compen_flag = 1;
                end
            end
        end
    end
end

%% ---------------- Load PMT raw data ----------------- %%
n_pmt_channel = 5; %number of pmt channels
input_info.pmt_filename = strings(n_pmt_channel,1); %initialize PMT filenames
input_info.pmt_dir = strings(n_pmt_channel,1);%initialize PMT directory paths
pmt_file_ID = zeros(n_pmt_channel,1);%initialize PMT local file IDs
if nargin ==0
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
    
else
    analysis_mode = 1;
    disp_progress = 0;
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

    name_split = strsplit(input_info.time_dir,'\');
    sample_name = name_split{end-1};
    sample_name= strrep(sample_name,'_',' ');
end


%% ---------------- Runtime Display settings ----------------- %%
if disp_progress==1
    %set analysis display screen size
    scrsize = get(0, 'Screensize');
    figure('OuterPosition',[0 0.05*scrsize(4) scrsize(3) 0.95*scrsize(4)])
end
disp_params=[analysis_mode, disp_progress];

%get how many segments to conduct analysis 
n = 1;
while(fseek(pmt_file_ID(1), n*8*analysis_params.datasize, 'bof') == 0)
    % flip forward 8*analysis_params.datasize bytes repeatedly until file ends
    n = n + 1;
end
num_segments = n - 1; % total number of segments = length of file in segments

%% -------------  Main analysis on looping data segments  ------------ %%
% display progress bar
progress_bar = waitbar(0,'Starting analysis...');
progress_bar.Position=[600,400,290,170];
pause(0.5)
n_out_features = 6; % total number of output features per detected signal, ie time, pmt, length, shape factors

rawdata_pmt = cell(1,n_pmt_channel);

segment_loop=0;
flag = 0;
full_readout_initialized = 0;

while(flag==0)

    % seek data for current segement, datatype int is 8bytes
    for channel = 1:n_pmt_channel
        fseek(pmt_file_ID(channel),segment_loop*8*analysis_params.datasize, 'bof');
    end
    fseek(time_file_ID,segment_loop*8*analysis_params.datasize, 'bof');
    
    % read raw pmt and time file
    for channel = 1:n_pmt_channel
        rawdata_pmt{1,channel} = fread(pmt_file_ID(channel),analysis_params.datasize,'float64=>double');
    end
    rawdata_time_pmt = fread(time_file_ID,analysis_params.datasize,'float64=>double');
  
    [seg_readout_pmt,progress_msg] = P1_peakanalysis_pmt(segment_loop,num_segments,rawdata_pmt, rawdata_time_pmt, disp_params, analysis_params,input_info);
  
    if ~isempty(seg_readout_pmt)
        if full_readout_initialized == 0
            full_readout_pmt = seg_readout_pmt;
            full_readout_initialized = 1;
        else        
            full_readout_pmt = vertcat(full_readout_pmt,seg_readout_pmt);
        end
        waitbar(segment_loop/num_segments,progress_bar,{progress_msg.line0,progress_msg.line1,progress_msg.line2,progress_msg.line3});
        pause(0.01)
    end
    
%     if segment_loop == 1498
%         waitbar(1,progress_bar,{'Finishing'});
%         pause(0.5)
%         flag = 1;
%     end
%     
    segment_loop=segment_loop+1;
    
    if length(rawdata_pmt{1,1}) < analysis_params.datasize
        waitbar(1,progress_bar,{'Finishing'});
        pause(0.5)
        flag = 1;
    end
end


%% Quality check to remove low-quality signals
if fxm_mode == 1
    non_nan_ind = find(~isnan(full_readout_pmt.baseline(:,fxm_channel)));
    base_med = median(full_readout_pmt.baseline(non_nan_ind(1:round(length(non_nan_ind)*0.3)),fxm_channel));
    all_base_med_norm = abs(full_readout_pmt.baseline(:,fxm_channel)-base_med);
    base_amp_pass_ind = find(all_base_med_norm < base_med*analysis_params.thresh_base_height_range);
    all_cell_baselineDiff_over_sig = abs(full_readout_pmt.baseline_left_height(:,fxm_channel)-full_readout_pmt.baseline_right_height(:,fxm_channel))...
        ./abs(full_readout_pmt.amplitude(:,fxm_channel)-full_readout_pmt.baseline(:,fxm_channel));
    base_diff_pass_ind = find(all_cell_baselineDiff_over_sig<analysis_params.thresh_baselineDiff_over_sig);
    base_leftslope_pass_ind = find(abs(full_readout_pmt.baseline_left_slope(:,fxm_channel)) < analysis_params.thresh_base_slope);
    base_rightslope_pass_ind = find(abs(full_readout_pmt.baseline_right_slope(:,fxm_channel)) < analysis_params.thresh_base_slope);

    cell_pass_ind =intersect(base_amp_pass_ind,intersect(base_diff_pass_ind, intersect(base_leftslope_pass_ind,base_rightslope_pass_ind)));
    pass_rate = 100*length(cell_pass_ind)/height(full_readout_pmt.amplitude);
    QC_msg = sprintf('%% %0.2f of detected signals passed QC check',pass_rate);
else
    pass_rate = 100;
    cell_pass_ind = 1:1:height(full_readout_pmt);
    QC_msg = sprintf('QC check are passed');
end
waitbar(1,progress_bar,QC_msg);
pause(0.5)
% hist(all_cell_baselineDiff_over_sig,1000)
%  scatter(abs(full_readout_pmt.baseline_right_slope(:,fxm_channel)),abs(full_readout_pmt.baseline_left_slope(:,fxm_channel)))
%% Apply compensation to fxm channel if needed by user
if fxm_mode == 1
    fxm_compen_amp = full_readout_pmt.amplitude(:,fxm_channel);
    if analysis_params. upstream_compen == 1
        for i = 1:length(compen_factor)
            fxm_compen_amp = (abs(full_readout_pmt.amplitude(:,i)-full_readout_pmt.baseline(:,i)))*compen_factor(i)+fxm_compen_amp;
        end
        full_readout_pmt.amplitude(:,fxm_channel) = fxm_compen_amp;
        upComp_msg = sprintf('Fxm upstream compensation applied');
    else
        upComp_msg = sprintf('Fxm upstream compensation skipped');
    end
else
    upComp_msg = sprintf('Fxm upstream compensation skipped');
end

waitbar(1,progress_bar,{QC_msg,upComp_msg});
pause(0.5)
%% Apply compensation to fxm - downstream channel to remove effect from fxm spillover
full_readout_pmt.amplitude_mV = full_readout_pmt.amplitude*1000; %Convert to mV
full_readout_pmt.baseline_mV = full_readout_pmt.baseline*1000; %Convert to mV
% Set pmt2 to pmt3 spillover factor
s = 0.042;

amp_detect_ch2 = full_readout_pmt.amplitude_mV(:,2);
amp_detect_ch3 = full_readout_pmt.amplitude_mV(:,3);
base_detect_ch2 = full_readout_pmt.baseline_mV(:,2);
base_detect_ch3 = full_readout_pmt.baseline_mV(:,3);

amp_true_ch2 = (base_detect_ch3.*amp_detect_ch2-base_detect_ch2.*amp_detect_ch3)./(base_detect_ch3-s*base_detect_ch2);
amp_true_ch3 = amp_detect_ch3-s*amp_true_ch2;

full_readout_pmt.signal(:,3) = 1000*abs(amp_true_ch3-base_detect_ch3)./base_detect_ch3;
full_readout_pmt.signal(:,1) = abs(full_readout_pmt.amplitude_mV(:,1) - full_readout_pmt.baseline_mV(:,1));
full_readout_pmt.signal(:,2) = amp_true_ch2;
full_readout_pmt.signal(:,4) = abs(full_readout_pmt.amplitude_mV(:,4) - full_readout_pmt.baseline_mV(:,4).*...
             (amp_true_ch3./base_detect_ch3));
full_readout_pmt.signal(:,5) = abs(full_readout_pmt.amplitude_mV(:,5) - full_readout_pmt.baseline_mV(:,5).*...
             (amp_true_ch3./base_detect_ch3));

downComp_msg = sprintf('Fxm downstream compensation applied');
waitbar(1,progress_bar,{QC_msg,upComp_msg,downComp_msg});
pause(0.5)
%% Generate PMT readout output file
%format follows: 
%for only positive labeling: [time of detection(computer real time), PacificBlue(mV),FITC(mV), PE(mV), APC(mV), Cy7(mV)]
%for Fxm, the fxm channel have raw volume signal as the unitless number, it
%is the dip in baseline normlized by the baseline height
output_msg = 'Generating output';
waitbar(1,progress_bar,{QC_msg,upComp_msg,downComp_msg,output_msg});
pause(0.5)

output_pmt = [full_readout_pmt.time_of_detection(cell_pass_ind),full_readout_pmt.signal((cell_pass_ind),:)];
output_pmt_table = array2table(output_pmt);
output_pmt_table.Properties.VariableNames = {'real_time_sec','pmt1_mV','pmt2_mV','vol_au','pmt4_mV','pmt5_mV'};
% if fxm_mode == 1
%     output_pmt_table.Properties.VariableNames{fxm_channel+1} = 'vol_au'; %% rename fxm channel with 1 offset colume of timestamp
% end

cd(input_info.pmt_dir(1))
if fxm_mode == 0
    out_file_name = ['readout_pmt_crosscompensated_' sample_name '.txt'];
else
    out_file_name = ['readout_pmt_crosscompensated_' sample_name '.txt'];
end
writetable(output_pmt_table, out_file_name, 'delimiter', '\t');
cd(currentFolder)

output_msg = 'Generating output... done';
waitbar(1,progress_bar,{QC_msg,upComp_msg,downComp_msg,output_msg});
pause(0.5)
report_msg = 'Generating report...';
waitbar(1,progress_bar,{QC_msg,upComp_msg,downComp_msg,output_msg,report_msg});
pause(0.5)
% generate analysis report
report_dir = [input_info.pmt_dir{1} '\' sample_name '_report\PMT_report\'];
mkdir(report_dir)
PMT_readout_report_v1(report_dir,input_info,input_info.pmt_dir(1),sample_name, analysis_params, full_readout_pmt(cell_pass_ind,:),QC_msg);
report_msg = 'Generating report... done';
waitbar(1,progress_bar,{QC_msg,upComp_msg,downComp_msg,output_msg,report_msg});

%% generate analysis params txt files 
cd(input_info.pmt_dir{1})
    writetable(struct2table(analysis_params),'readout_pmt_analysis_params.txt','Delimiter',' ')
cd(currentFolder)

if nargin ~=0
    delete(progress_bar)
end
%%
figure(1)
scatter(full_readout_pmt.time_of_detection(:),full_readout_pmt.baseline(:,fxm_channel)...
    ,5,'filled',"MarkerFaceAlpha",0.5)
hold on
scatter(full_readout_pmt.time_of_detection(cell_pass_ind),full_readout_pmt.baseline(cell_pass_ind,fxm_channel)...
    ,5,'filled',"MarkerFaceAlpha",0.5)

%%
figure(2)
scatter(full_readout_pmt.time_of_detection(cell_pass_ind),full_readout_pmt.signal(cell_pass_ind,fxm_channel)...
    ,5,'filled',"MarkerFaceAlpha",0.5)
%%
% scatter(abs(full_readout_pmt.signal(cell_pass_ind,fxm_channel)),abs(full_readout_pmt.amplitude(cell_pass_ind,4)-full_readout_pmt.baseline(cell_pass_ind,4))*1000)
% hold on
% tt = (full_readout_pmt.amplitude(cell_pass_ind,4)-(full_readout_pmt.baseline(cell_pass_ind,4)).*...
%             (full_readout_pmt.amplitude(cell_pass_ind,fxm_channel)./full_readout_pmt.baseline(cell_pass_ind,fxm_channel)));
% scatter(abs(full_readout_pmt.signal(cell_pass_ind,fxm_channel)),tt*1000)
% %%
% co_check_ind = find(full_readout_pmt.signal(cell_pass_ind,fxm_channel)>10&full_readout_pmt.signal(cell_pass_ind,fxm_channel)<300);
% co_check_ind = cell_pass_ind(co_check_ind);
% x = abs(full_readout_pmt.signal(co_check_ind,fxm_channel));
% y =abs((full_readout_pmt.amplitude(co_check_ind,3)-full_readout_pmt.baseline(co_check_ind,3))./full_readout_pmt.baseline(co_check_ind,3))*1000;
% scatter(x,y)
% temp_stats = regstats(y,x,'linear');
% disp(temp_stats.beta(2))
% 
% scatter(full_readout_pmt.signal(cell_pass_ind,fxm_channel),full_readout_pmt.signal(cell_pass_ind,4))


end

