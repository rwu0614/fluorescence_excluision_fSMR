function [readout_pmt,progress_msg] = P1_peakanalysis_pmt(n_segment,num_segments,...
    rawdata_pmt,rawdata_time_pmt, disp_params, analysis_params,input_info)

global elapsed_time;
global elapsed_index;
global elapsed_peak_count;

%Allocating peak detection parameters
n_pmt_channel =5;
analysis_mode = disp_params(1);
disp_progress = disp_params(2);

Baseline_rough_cutoff = analysis_params.Baseline_rough_cutoff;
med_filt_length = analysis_params.med_filt_length;
moving_average_window_size = analysis_params.moving_average_window_size;
med_filt_window_size = analysis_params.med_filt_window_size;
min_distance_btw_peaks = analysis_params.min_distance_btw_peaks;
detect_thresh_pmt = analysis_params.detect_thresh_pmt;
uni_peak_range_ext = analysis_params.uni_peak_range_ext;
uni_peak_baseline_window_size = analysis_params.uni_peak_baseline_window_size;

fxm_mode = analysis_params.fxm_mode;
fxm_channel = find(analysis_params.detect_thresh_pmt<0);
fxm_baseline_choice = analysis_params.fxm_baseline_choice;

%initialize cells structures for analysis and output
Data.read = cell(1, n_pmt_channel);
Data.normalized = cell(1, n_pmt_channel);
Data.filtered_med = cell(1, n_pmt_channel);
Data.filtered_med_ave = cell(1, n_pmt_channel);

Process.baseline = cell(1, n_pmt_channel);
Process.peak_threshold = cell(1, n_pmt_channel);
Process.peak_indices = cell(1, n_pmt_channel);

Peak.start = cell(1, n_pmt_channel);
Peak.end = cell(1, n_pmt_channel);
Peak.count = cell(1, n_pmt_channel);
Peak.amplitude = cell(1, n_pmt_channel);
Peak.amp_over_base = cell(1, n_pmt_channel);
Peak.location = cell(1, n_pmt_channel);
Peak.time = cell(1, n_pmt_channel);
Peak.width = cell(1, n_pmt_channel);
Peak.baseline = cell(1, n_pmt_channel);
Peak.time_of_detection = [];


% Rough quality check on PMT data being above minimum expected voltage
idx_rangeofinterest = find(rawdata_pmt{1,1} > Baseline_rough_cutoff);

% Existing segment when too few signals are above baseline cutoff (for long
% experiment when light source is turned off when not used
if length(idx_rangeofinterest) < med_filt_window_size
    readout_pmt = [];
    progress_msg = [];
    %fprintf('Too short for peak selection...exiting this segment\n')
    elapsed_time = elapsed_time + (rawdata_time_pmt(end) - rawdata_time_pmt(1));
    elapsed_index = elapsed_index + length(rawdata_pmt{1});
    return;
end

% Allocate raw pmt and time data
for i = 1:n_pmt_channel
    Data.read{1,i} = rawdata_pmt{1,i}(idx_rangeofinterest);
end
t = rawdata_time_pmt(idx_rangeofinterest); 

% Nomarlize each channel by its median to remove channel-dependent
% baseline voltage difference

%The median normalization step is skipped to make sure fluorescence exclusion baseline does not get erased
% Data.normalized = cellfun(@(x) x-median(x),Data.read,...
%     'UniformOutput',false);
Data.normalized = Data.read; 

%===== Apply median filter then moving-average filter =============%

Data.filtered_med = cellfun(@(x) medfilt1(x,med_filt_length),Data.normalized,...
    'UniformOutput',false);

%setting up parameter for moving average filter
b = (1/moving_average_window_size)*ones(1,moving_average_window_size); a=1; 

Data.filtered_med_ave = cellfun(@(x) filter(b,a,x), Data.filtered_med,...
    'UniformOutput',false);



%================Median filtered baseline creation================%
%select x-axis index for pmt baseline fitting using defined sampling distance
x_axis_Ind = [1:length(Data.normalized{1})]';
med_filt_x_axis_Ind = 1:med_filt_window_size:x_axis_Ind(end); 

%initialize fitted pmt baselines
med_filt_pmt_base = zeros(n_pmt_channel, length(med_filt_x_axis_Ind));

%applying median filter between each baseline sampling distance
for channel = 1:n_pmt_channel
    for i= 1:length(med_filt_x_axis_Ind)-1
        base_fit_range = med_filt_x_axis_Ind(i):med_filt_x_axis_Ind(i+1);
        med_filt_pmt_base(channel,i) = median(Data.filtered_med{channel}(base_fit_range));
    end
end
% tying up ends of fitted baselines
if med_filt_x_axis_Ind(end) == x_axis_Ind(end)
    med_filt_pmt_base(:,end) = med_filt_pmt_base(:,end-1);
else % if ends are different
    med_filt_pmt_base(:,end) = cellfun(@(x) median(x(med_filt_x_axis_Ind(end):x_axis_Ind(end))), Data.filtered_med)';
end

% linear extrapolation to get full time domain baselines from median fitted
% baselines
for i = 1:n_pmt_channel
    Process.baseline{i}=interp1(med_filt_x_axis_Ind, med_filt_pmt_base(i,:), x_axis_Ind, 'linear', 'extrap');
end



%======================Find rough peak indices=====================%   
%creating baseline-matched thresholds for peak detection from user-defined
%delta-voltage thresholds

for i = 1:n_pmt_channel
% using absolute threshold cutoff values
%   Process.peak_threshold{i}= Process.baseline{i} + detect_thresh_pmt(i);
% using multiplier of noise amplitude as cutoff values
    % compute noise of baseline by taking the differential of baseline and
    % remove 2.5% of lower and higher outliers, and take standard deviation
    noise_data_keep = rmoutliers(diff(Data.normalized{1,i}),'percentiles',[1,99]);
    Process.peak_threshold{i}= Process.baseline{i} + detect_thresh_pmt(i)*std(noise_data_keep);
end


% Find all indices where median filtered data passes threshold. Note that
% for fxm it's the points that are below the threshold that are consided a
% peak
Process.peak_indices = cellfun(@(x,y) find(x > y), Data.filtered_med, Process.peak_threshold,"UniformOutput",false);
if fxm_mode == 1
    Process.peak_indices{fxm_channel} = find(Data.filtered_med{fxm_channel} < Process.peak_threshold{fxm_channel});
end


%%=================== FINDING PEAK RANGE===================%%
for i=1:n_pmt_channel
    if ~isempty(Process.peak_indices{i})
        temp_ind = Process.peak_indices{i}';
        %find indecies for the starting points of all unique peaks
        Peak.start{i} = [temp_ind(1) temp_ind([0 diff(temp_ind)]>min_distance_btw_peaks)];
        %find indecies for the end points of all unique peaks
        Peak.end{i} = fliplr([temp_ind(end) fliplr(temp_ind(fliplr([0 diff(fliplr(temp_ind))]<-min_distance_btw_peaks)))]);
        Peak.count{i} = length(Peak.end{i});
    else
        Peak.start{i}=[];
        Peak.end{i}=[];
        Peak.count{i}=length(Peak.end{i});
    end
end
%exiting the script if no peaks were found
if sum([Peak.count{:}])==0
    %fprintf('\nNo peak found in this section...exiting this segment\n');
    readout_pmt = [];
    progress_msg = [];
    elapsed_time = elapsed_time + (rawdata_time_pmt(end) - rawdata_time_pmt(1));
    elapsed_index = elapsed_index + length(rawdata_pmt{1,i});
    return;
end

%======================Finding individual peaks=====================%

if disp_progress==1
    %allocating space to plot real-time analysis
    plotspace.n_row = 5; 
    plotspace.n_col = n_pmt_channel;

    plotspace.segment_start = 1+(0:plotspace.n_col-1)*plotspace.n_row;
    plotspace.segment_end = 3+(0:plotspace.n_col-1)*plotspace.n_row;
    plotspace.peak = 4+(0:plotspace.n_col-1)*plotspace.n_row;
    plotspace.peakshape = 5+(0:plotspace.n_col-1)*plotspace.n_row;
    channel_name = ["Pacific Blue","FITC","PE","APC","CY7"];
    for i = 1:plotspace.n_col
        subplot(plotspace.n_col,plotspace.n_row,[plotspace.segment_start(i) plotspace.segment_end(i)]);
            plot(x_axis_Ind, Data.filtered_med_ave{1,i}, 'r-');hold on
            plot(x_axis_Ind, Process.peak_threshold{1,i}, 'k'); 
            plot(x_axis_Ind, Process.baseline{1,i}, 'LineStyle', ':');
            ylabel(append(channel_name(i),' (V)'));
            title(append('Segment view [segment #',num2str(n_segment),']'));
        subplot(plotspace.n_col,plotspace.n_row,plotspace.peak(channel));
            ylabel(append(channel_name(channel),' (V)'));
            title(append('Peak view'));
    end
    

end
% initializing loop parameters
seg_num_peaks=0; 
exitflag = 0;  
while (exitflag ~=1)
    seg_num_peaks=seg_num_peaks+1;

    [Peak, tempstart, tempend] = P2_find_pmt_peaks(Peak,uni_peak_range_ext);

    maxrange = max(tempstart,1):min(tempend, x_axis_Ind(end));
    
    plotrange = max(tempstart-uni_peak_baseline_window_size,1):min(tempend+uni_peak_baseline_window_size, x_axis_Ind(end));
    %plotrange = max(tempend,1):min(tempend+uni_peak_baseline_window_size, x_axis_Ind(end)); % to correct for scattering on the left side
    left_base_range = max(tempstart-uni_peak_baseline_window_size,1):tempstart;
    right_base_range = tempend:min(tempend+uni_peak_baseline_window_size, x_axis_Ind(end));
%     local_baseline = cellfun(@(x) median(x(setxor(plotrange, maxrange))),Data.filtered_med);
    local_baseline = cellfun(@(x) median(x(plotrange)),Data.filtered_med);
    
    if fxm_baseline_choice==1 && fxm_mode==1
        plotrange_fxm = max(tempstart-uni_peak_baseline_window_size,1):min(tempend, x_axis_Ind(end)); % grab only left side baseline
        local_baseline(fxm_channel) = median(Data.filtered_med{fxm_channel}(plotrange_fxm));
    elseif fxm_baseline_choice ==2 && fxm_mode==1
        plotrange_fxm = max(tempend,1):min(tempend+uni_peak_baseline_window_size, x_axis_Ind(end)); % grab only right side baseline
        local_baseline(fxm_channel) = median(Data.filtered_med{fxm_channel}(plotrange_fxm));
    end

    i=seg_num_peaks;

    for channel = 1:n_pmt_channel
         [~, Peak.location{channel}(i)] = cellfun(@(x) max(abs(x(maxrange)-local_baseline(channel))), Data.filtered_med_ave(channel));
        Peak.location{channel}(i) = Peak.location{channel}(i) + maxrange(1) - 1;
        Peak.amplitude{channel}(i) = Data.filtered_med_ave{channel}(Peak.location{channel}(i)) ; % absolute pmt voltage at the peak
        Peak.baseline{channel}(i) = local_baseline(channel); % absolute pmt voltage at the baseline near the peak
        Peak.time{channel}(i)=t(Peak.location{channel}(i));
        Peak.width{channel}(i)= length(maxrange);
        
        left_base_fit = polyfit(1:1:length(left_base_range),Data.normalized{channel}(left_base_range),1);
       
        right_base_fit = polyfit(1:1:length(right_base_range),Data.normalized{channel}(right_base_range),1);
        
        Peak.baseline_left_slope{channel}(i) = left_base_fit(1);
        Peak.baseline_right_slope{channel}(i) = right_base_fit(1);
        Peak.baseline_left_height{channel}(i) = median(Data.normalized{channel}(left_base_range));
        Peak.baseline_right_height{channel}(i) = median(Data.normalized{channel}(right_base_range));

        
    end
    if analysis_mode==0 && fxm_mode==1
       disp("Fxm left baseline fitted slope") 
       disp(Peak.baseline_left_slope{fxm_channel}(i))
       disp("Fxm right_base_fitted slope ") 
       disp(Peak.baseline_right_slope{fxm_channel}(i))
       disp("Fxm baseline height difference (% peak height)") 
       disp(abs(Peak.baseline_left_height{fxm_channel}(i)-Peak.baseline_right_height{fxm_channel}(i))/...
           abs(Peak.amplitude{fxm_channel}(i)-Peak.baseline{fxm_channel}(i)));
    end
    
    % Grab time of detection from highist intensity signal
    temp_all_peak_time = [Peak.time{1}(i),Peak.time{2}(i),Peak.time{3}(i),...
        Peak.time{4}(i),Peak.time{5}(i)];
    [~,channelID_max_amp] = max([Peak.amplitude{1}(i),Peak.amplitude{2}(i),...
        Peak.amplitude{3}(i),Peak.amplitude{4}(i),Peak.amplitude{5}(i)]);
    Peak.time_of_detection(i) = temp_all_peak_time(channelID_max_amp);

    % Display analysis process (peak recognition)
    if disp_progress==1
        % Enter peak shape analysis
    %Peak_shape = P2_peak_shape(t,Data,Peak,channelID_max_amp,plotrange,maxrange,disp_progress,plotspace,channel_name,input_info); 
        for channel = 1:plotspace.n_col
            subplot(plotspace.n_col,plotspace.n_row,[plotspace.segment_start(channel) plotspace.segment_end(channel)]);
                plot(Peak.location{channel}, Data.filtered_med_ave{channel}(Peak.location{channel}), '*g');hold on;
                if channel ==plotspace.n_col(end)
                    axP = get(gca,'Position');
                    legend(["Median + moving average filter",...
                        "Peak detection threshold","Baseline","Peak location (running)"]...
                    ,'Location','SouthOutside');
                    set(gca, 'Position', axP);
                end
            subplot(plotspace.n_col,plotspace.n_row,plotspace.peak(channel));
%                 plot(plotrange, Data.filtered_med{channel}(plotrange)-local_baseline(channel), '.', 'color', 'blue'); hold on;
                plot(plotrange, Data.normalized{channel}(plotrange)-local_baseline(channel), '.', 'color', 'blue'); hold on;
                plot(left_base_range, Data.normalized{channel}(left_base_range)-local_baseline(channel), '.', 'color', 'yellow'); hold on;
                plot(right_base_range, Data.normalized{channel}(right_base_range)-local_baseline(channel), '.', 'color', 'cyan'); hold on;
                plot(maxrange, Data.filtered_med_ave{channel}(maxrange)-local_baseline(channel), '-r');
                plot(Peak.location{channel}(i), Peak.amplitude{channel}(i)-local_baseline(channel), 'or');
                plot(plotrange, zeros(1,length(plotrange)), 'k-');
                hold off;
                if channel ==plotspace.n_col(end)
                    axP = get(gca,'Position');
                    legend(["Raw data","Left baseline","Left baseline","Median + moving average filter",...
                        "Peak","Baseline"]...
                    ,'Location','SouthOutside');
                    set(gca, 'Position', axP);
                end
            drawnow;
        end

    end

    
    
    
    for i=1:n_pmt_channel
        Peak.count{i} = length(Peak.end{i});
    end

    if sum([Peak.count{:}])==0
        exitflag = 1;
    end
    if analysis_mode==0
        pass = input('pass? yes-1,no-0');
        if pass==0
            Peak.time_of_detection(i)=-1;
        end
    end
end

%% Generate output
readout_pmt = table();
readout_pmt.time_of_detection = Peak.time_of_detection';
readout_pmt.amplitude=cell2mat(Peak.amplitude')';
readout_pmt.baseline=cell2mat(Peak.baseline')';
readout_pmt.width=cell2mat(Peak.width(1)')';
readout_pmt.location = cell2mat(cellfun(@(x) x + elapsed_index,Peak.location,"UniformOutput",false)')';
readout_pmt.baseline_left_slope = cell2mat(Peak.baseline_left_slope')';
readout_pmt.baseline_right_slope = cell2mat(Peak.baseline_right_slope')';
readout_pmt.baseline_left_height = cell2mat(Peak.baseline_left_height')';
readout_pmt.baseline_right_height =cell2mat(Peak.baseline_right_height')';


if disp_progress==1
    for i = 1:plotspace.n_col
            subplot(plotspace.n_col,plotspace.n_row,[plotspace.segment_start(i) plotspace.segment_end(i)]);
            hold off;
    end
end

elapsed_time = elapsed_time + (rawdata_time_pmt(end) - rawdata_time_pmt(1));
elapsed_index = elapsed_index + length(rawdata_pmt{1});
elapsed_peak_count = elapsed_peak_count + length(Peak.time_of_detection);

progress_msg.line0 = sprintf('   Segment#%1.0f of %1.0f\n',n_segment,num_segments);
progress_msg.line1 = sprintf('  Segment peak count = %1.0f, Total peak count = %1.0f \n',seg_num_peaks,elapsed_peak_count);
progress_msg.line2 = sprintf('  %%-- elapsed time   : %1.2f minutes --%%\n', elapsed_time/60);
progress_msg.line3 = sprintf('  %%-- Average detection throughput: %1.2f events/minutes --%%\n\n', elapsed_peak_count/(elapsed_time/60));

end