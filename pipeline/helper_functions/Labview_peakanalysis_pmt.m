
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimize using following parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Peak_length = 100; % estimated peak length
    
    analysis_params.Baseline_rough_cutoff = -20; % default is -20 for populational fSMR experiments where light source is always on
    analysis_params.med_filt_length = 5; %full PMT data median filter window size, default 50
    analysis_params.moving_average_window_size = 5; %full PMT data moving-average filter window size, default 5    
    analysis_params.med_filt_window_size = 8*Peak_length ; % baseline median filter window size, sampling distance for extrapolating flat baseline   
    analysis_params.min_distance_btw_peaks = 50; % minimum distance between peaks, for identifying unique peaks
    analysis_params.uni_peak_range_ext = 5; % number of data points from each side of detection cutoff to be considered as part of the peak
    analysis_params.uni_peak_baseline_window_size = 100; % length of data points from each side of detection cutoff to compute the local baseline
    
    %calcien+annexin detection stragety: prioritize calcien
    analysis_params.detect_thresh_pmt(1) = thresh_multi*noise(1); 
    analysis_params.detect_thresh_pmt(2) = thresh_multi*noise(2);
    analysis_params.detect_thresh_pmt(3) = thresh_multi*noise(3);
    analysis_params.detect_thresh_pmt(4) = thresh_multi*noise(4);
    analysis_params.detect_thresh_pmt(5) = thresh_multi*noise(5);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Allocating peak detection parameters
n_pmt_channel =5;

Baseline_rough_cutoff = analysis_params.Baseline_rough_cutoff;
med_filt_length = analysis_params.med_filt_length;
moving_average_window_size = analysis_params.moving_average_window_size;
med_filt_window_size = analysis_params.med_filt_window_size;
min_distance_btw_peaks = analysis_params.min_distance_btw_peaks;
detect_thresh_pmt = analysis_params.detect_thresh_pmt;
uni_peak_range_ext = analysis_params.uni_peak_range_ext;
uni_peak_baseline_window_size = analysis_params.uni_peak_baseline_window_size;

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
Peak.location = cell(1, n_pmt_channel);
%Peak.time = cell(1, n_pmt_channel);
Peak.width = cell(1, n_pmt_channel);
Peak.baseline = cell(1, n_pmt_channel);
%Peak.time_of_detection = [];

%if labview enqueued pmt matrix is a 5xn numerical matrix, then
for i = 1:n_pmt_channel
    rawdata_pmt{1,i} = pmt_que(i,:);
end

% Rough quality check on PMT data being above minimum expected voltage
idx_rangeofinterest = find(rawdata_pmt{1,1} > Baseline_rough_cutoff);

% Allocate raw pmt and time data
for i = 1:n_pmt_channel
    Data.read{1,i} = rawdata_pmt{1,i}(idx_rangeofinterest);
end
%t = rawdata_time_pmt(idx_rangeofinterest); 

% Nomarlize each channel by its median to remove channel-dependent
% baseline voltage difference
Data.normalized = cellfun(@(x) x-median(x),Data.read,...
    'UniformOutput',false);

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
Process.peak_threshold{i}= Process.baseline{i} + detect_thresh_pmt(i);
end

Process.peak_indices = cellfun(@(x,y) find(x > y), Data.filtered_med, Process.peak_threshold,"UniformOutput",false);



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
    %readout_pmt.time_of_detection = [];
    readout_pmt.amplitude=[];
    readout_pmt.location = [];
    readout_pmt.baseline=[];
    
else

    %======================Finding individual peaks=====================%

    % initializing loop parameters
    seg_num_peaks=0; 
    exitflag = 0;  
    while (exitflag ~=1)
        seg_num_peaks=seg_num_peaks+1;

        [Peak, tempstart, tempend] = P2_find_pmt_peaks(Peak,uni_peak_range_ext);

        maxrange = max(tempstart,1):min(tempend, x_axis_Ind(end));
        plotrange = max(tempstart-uni_peak_baseline_window_size,1):min(tempend+uni_peak_baseline_window_size, x_axis_Ind(end));
        local_baseline = cellfun(@(x) median(x(setxor(plotrange, maxrange))),Data.filtered_med);

        i=seg_num_peaks;

        for channel = 1:n_pmt_channel
            [Peak.amplitude{channel}(i), Peak.location{channel}(i)] = cellfun(@(x) max(x(maxrange)), Data.filtered_med_ave(channel));
            Peak.location{channel}(i) = Peak.location{channel}(i) + maxrange(1) - 1;
            Peak.amplitude{channel}(i) = Peak.amplitude{channel}(i) - local_baseline(channel); %correct for local baseline
            %Peak.time{channel}(i)=t(Peak.location{channel}(i));
            Peak.width{channel}(i)= length(maxrange);
            Peak.baseline{channel}(i) = local_baseline(channel) + median(Data.read{channel}); 
        end


        % Grab time of detection from highist intensity signal
%         temp_all_peak_time = [Peak.time{1}(i),Peak.time{2}(i),Peak.time{3}(i),...
%             Peak.time{4}(i),Peak.time{5}(i)];
%         [~,channelID_max_amp] = max([Peak.amplitude{1}(i),Peak.amplitude{2}(i),...
%             Peak.amplitude{3}(i),Peak.amplitude{4}(i),Peak.amplitude{5}(i)]);
%         Peak.time_of_detection(i) = temp_all_peak_time(channelID_max_amp);

        for i=1:n_pmt_channel
            Peak.count{i} = length(Peak.end{i});
        end

        if sum([Peak.count{:}])==0
            exitflag = 1;
        end

    end

    % Generate output

    %readout_pmt.time_of_detection = Peak.time_of_detection;
    readout_pmt.amplitude=Peak.amplitude;
    readout_pmt.location = cellfun(@(x) x + elapsed_index,Peak.location,"UniformOutput",false);
    readout_pmt.baseline=Peak.baseline;
end

    %time_of_detection = readout_pmt.s;
    voltage_pmt1= readout_pmt.amplitude{1};
    voltage_pmt2= readout_pmt.amplitude{2};
    voltage_pmt3= readout_pmt.amplitude{3};
    voltage_pmt4= readout_pmt.amplitude{4};
    voltage_pmt5= readout_pmt.amplitude{5};


