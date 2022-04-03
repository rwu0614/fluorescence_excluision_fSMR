function data = S1_PeakAnalysis_time(x, t, datalast, sectionnumber, analysisparam)

global xdata
global ydata
global left_base
global right_base
global segmentbound
global elapsed_time;
global samplepeak
global sampletime 
global sidelength;

analysismode = analysisparam(1);
dispprogress = analysisparam(2);
vv = zeros(length(x),1);
% clear all S1 variables
pkidx_poly = [];
apkidx_poly = [];
pk_t = [];
apk_t = [];
pkht_poly = [];
apkht_poly = [];
peakwidth = [];
pk_leftbase = [];
pk_rightbase = [];
baselinedist = [];
baselineslope = [];
htdiff_poly = [];
ahtdiff_poly = [];                                                                                          
pknum = [];
sectnum = [];
pkorder = [];



%t=t-t(1);

if length(x)<10
    %disp('No GOOD peak found in this section: data Length< 10');
    data=zeros(13,1); 
    %disp(' '); 
elapsed_time = elapsed_time + t(end);  return;
end



                                                             % these are just indices


%% ================== estimated data points per frequency peak & noise in system ==============
%added 03/22/2019
estimated_datapoints = 600; %for full transit; deafult
estimated_noise = 0.6; %Hz


sgolay_length_idx = 1; %this will make default length of 5 for 50 idx transits
sgolay_length = 2*round(sgolay_length_idx*(estimated_datapoints/50))+1;
fprintf('sgolay length for datapoints per transit of %4.0f will be %2.0f', estimated_datapoints, sgolay_length); disp(' ');
if elapsed_time==0
    xtest = x([1:min([1e5, length(x)])]) - x(1);
figure(2); plot(xtest); hold on;
plot(sgolayfilt(xtest, 3,5));
if sgolay_length>3
plot(sgolayfilt(xtest, 3,sgolay_length));
legend('raw', 'sgolay 3-5', strcat('sgolay 3-', num2str(sgolay_length)));
else
legend('raw', 'sgolay 3-5');
end

%noise_estimate_frame = input('What is the x-axis frame for estimating noise? [_,_]');
%measured_noise = std(xtest(noise_estimate_frame));
measured_noise = std(xtest([1:10000]));
fprintf('estimated noise level is %2.3f', measured_noise); disp(' ');
input('check if the filtering is okay and your noise level');
estimated_noise = measured_noise;
end    



%%
if sgolay_length>3
ydata = sgolayfilt(x,3,sgolay_length); 
else
ydata = sgolayfilt(x,3,5);
end
xdata = [1:length(ydata)]';    
%disp('Finding peak indices...')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline selection
% Optimize using below parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_threshold = 0.005;      % Find extremely flat part of curve (sys1)
med_filt_wd = 200;           % window of median filter, which removes the flat part in the anti-node
bs_dev_thres = 0.5;         % baseline dev_threshold; threshold used to remove the flat part in the anti-node
unqPeakDist=150;            % distance over which is a unique 2nd mode peaks,default 200
offset_input = 3;             % baseline offset to select for peaks

%% added 03222019 - compensate for number of data points & noise leve;
diff_threshold = diff_threshold*((estimated_noise/0.1)^(1/2))/(estimated_datapoints/400);
med_filt_wd = round(med_filt_wd * estimated_datapoints/400);
bs_dev_thres = bs_dev_thres*((estimated_noise/0.1)^(1/2));
unqPeakDist = round(unqPeakDist*estimated_datapoints/400);


idx = find(abs(diff(ydata))<diff_threshold); %remove fast varying points to get baseline. i.e., remove cell frequency                                                  


% ignore the flat points found over the anti-node
mf_ydata_thres=medfilt1(ydata(idx), med_filt_wd);


idx_f=find(abs(ydata(idx)-mf_ydata_thres)<bs_dev_thres);

idx=idx(idx_f);

if length(idx)<=1                   %% add on 01242022 by Ye
    if isempty(idx)
        disp('No GOOD peak found in this section: idx Length == 0');
    else
        disp('No GOOD peak found in this section: idx Length == 1');
    end
    data=zeros(13,1); 
    disp(' '); 
    elapsed_time = elapsed_time + t(end);  return;
end

ydata_thres=interp1(xdata(idx), ydata(idx), xdata);

size(ydata_thres)

ydata_thres=ydata_thres-offset_input;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = find(ydata < ydata_thres);                                                  % find all y_diff values that pass this threshold - these are the main peaks

repeatflag = 1;
while(repeatflag == 1)
    
    if ~isempty(idx)
        idx_ends = find(abs(diff(idx)) > 1);                                                    % find all breaks in y_diff indices; each segment is a peak of interest, with end on idx_end and start on (idx_end + 1)
        idx_ends = [0 idx_ends' length(idx)];                                                   % cap the beginning and end of idx with markers

        peak_idx = zeros(1,length(idx_ends) - 1);                                               % set empty matrix for loop efficiency

        for i = 1:length(idx_ends) - 1
            ydata_segment = ydata(idx(idx_ends(i)+1:idx_ends(i + 1)));                        % focus on the piece of y_diff between two consecutive idx_end markers
            segment_idx_max = min(find(ydata_segment == min(ydata_segment)));                      % find index of the maximum deviation, i.e. the apex within the piece
            global_idx_max = idx(idx_ends(i) + 1) + segment_idx_max - 1;
            peak_idx(i) = global_idx_max(1);                                                       % convert the apex index to the global index within y_diff      
        end
        
        unique_peaks = diff([peak_idx length(xdata)]) > unqPeakDist;
        peak_idx = peak_idx(unique_peaks);                                                      % which for single-cell 2nd-mode peaks are 3n+1, and thus we are saving the each third peak
                                                                                                %%% now we have the global indices of the main peaks in this segment of
                                                                                                %%% the entire frequency data 
        repeatflag = 0;
        %disp('...Done.')
    else
        %disp('No peak found in this section');
        repeatflag = 0;
        data=zeros(13,1); 
        elapsed_time = elapsed_time + t(end);  
        return;
    end
end



%disp(' ')
%disp('Segmentizing the dataset...')

segmentbound = zeros(1,length(peak_idx));

segmentbound(1) = peak_idx(1) + round((peak_idx(2) - peak_idx(1))/2);
for i = 2:length(peak_idx)-1
    segmentbound(i) = peak_idx(i) + round((peak_idx(i+1) - peak_idx(i))/2);
end
segmentbound(end) = xdata(end);
segmentbound = [1 segmentbound];
segment_med_wd = round(median(diff(segmentbound))/10);

set(gca,'XLim',[0 length(xdata)]);
set(gca,'YLim',[min(ydata) - 10 max(ydata) + 10]);
%disp('Done.')

%%=======================================================================
%% ==================== INDIVIDUAL PEAK ANALYSIS================ %%
%% ==============================================================%%
% disp(' ')
% disp('Beginning individual peak analysis.')
i = 0;

%-------------- Segmentizing the dataset -------------------------%%
segment_threshold=200; %number of data points from peaks. 
segment_threshold = segment_threshold*estimated_datapoints/400;
sidelength = 3000; %number of more indices on each side of baseline 
 
for i=1:length(peak_idx)
    
    if min(abs(peak_idx(i)-segmentbound(i)), abs(peak_idx(i)-segmentbound(i+1))) > segment_threshold
    
        local_xdata = xdata(segmentbound(i):segmentbound(i + 1)) - xdata(segmentbound(i)) + 1;
        local_ydata = ydata(segmentbound(i):segmentbound(i + 1));
    
    
        
    
    %%%%%%%% USER ADJUSTABLE (Below is assuming 400 data points per transit) %%%%%%%%
    edgethres = 0.12;               % choose the first point left/right of the secondary peaks 40% percent of the average baseline freqvalue
    stdevmultiplier = 3;         % allow 102% of the minimum standard deviation
    diffmultiplier = 0.9;           % allow 90% of the deviation from mean frequency closest to 2ndary peaks
    winsize = 150;                  % number of points searching for baseline collection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    baseparams = [stdevmultiplier diffmultiplier edgethres winsize diff_threshold med_filt_wd bs_dev_thres];
    
    %%% Identify primary and secondary peaks within segment:
%     disp(' ')
%     disp('Locating segment peaks...')
    peaks = S2_PeaksetFinder(local_xdata, local_ydata, offset_input, baseparams, analysismode);  
%     disp('...Done.')
    
    if numel(peaks)==3
        
        
    peakdist_temp = peaks(end)-peaks(1);
    %%% Identify left- and right-hand baselines within segment:
%     disp(' ')
%     disp('Identifying baseline...')
    
    [left_base, right_base, edgeidx] = S2_BaselineFinder(local_xdata, local_ydata, peaks, baseparams, peakdist_temp, analysismode);  
    if(numel(left_base) == 0 || numel(right_base) == 0 || numel(edgeidx) == 0 || numel(peaks) == 0)
    i=i+1;   
    else
        
        
        
        local_peakwidth = diff(edgeidx);  
        
        pk_xdata = local_xdata(left_base(1):right_base(end)) - local_xdata(left_base(1)) + 1;
        pk_ydata = local_ydata(left_base(1):right_base(end));
        
        
       
        
        local_peaks = peaks - left_base(1) + 1;
        local_baseline = [left_base - left_base(1) + 1, right_base - left_base(1) + 1];
        %disp('...Done.')
        
%         input('continue?');
        %%% Perform polynomial fits to refine peak location and to measure peak height:
        %disp(' ')
        %disp('Performing polynomial fit on segment peaks...')
        [local_pkidx_poly, local_pkht_poly, local_apkidx_poly, local_apkht_poly, local_baselineslope, local_htdiff_poly, local_ahtdiff_poly] ...
            = S2_PeakFitter(pk_xdata, pk_ydata, local_baseline, local_peaks, local_peakwidth, dispprogress); 
        %disp('...Done.')
        %disp(' ')
        
         %% ---------------- added by JK 09/18/14 ------
        %here we are trying to save the mode shape for each peak detected
     
       
        %% Input Experiment parameters
            Experiment.R = 32768;
            Experiment.decimation = 1;
            Experiment.Fs = 100e6 / Experiment.R / Experiment.decimation;

       
        temppeak = local_ydata(left_base(1): right_base(end));
        
        temptime = [t(segmentbound(i) + left_base(1))+1/Experiment.Fs:1/Experiment.Fs:t(segmentbound(i) + left_base(1))+1/Experiment.Fs*length(temppeak)];
        samplepeak=[samplepeak temppeak' 1e3 i sectionnumber]; 
        sampletime = [sampletime temptime 0 i sectionnumber];
        % 0 is to distinguish different peaks

        format shortg
        fprintf('------- Section %1.0f ------ \n', sectionnumber);
        fprintf('-- Data for Segment %1.0f -- \n', i);
        fprintf('Baseline slope: %1.5f \n', local_baselineslope);
        fprintf('2nd mode %%diff: %1.5f \n', local_htdiff_poly);
        fprintf('-------------------------- \n')
        disp(' ')

        pkidx_poly = [pkidx_poly local_pkidx_poly + segmentbound(i) + left_base(1)];
        apkidx_poly = [apkidx_poly local_apkidx_poly + segmentbound(i) + left_base(1)];
        pk_t = [pk_t t(local_pkidx_poly + segmentbound(i) + left_base(1))'];
        apk_t = [apk_t t(local_apkidx_poly + segmentbound(i) + left_base(1))' 0];
        pkht_poly = [pkht_poly local_pkht_poly];
        apkht_poly = [apkht_poly local_apkht_poly 0];
        peakwidth = [peakwidth local_peakwidth*ones(1,length(peaks))];
        pk_leftbase = [pk_leftbase mean(local_ydata(left_base))*ones(1,length(peaks))];
        pk_rightbase = [pk_rightbase mean(local_ydata(right_base))*ones(1,length(peaks))];
        baselinedist = [baselinedist (right_base(end) - left_base(1) + 1)*ones(1,length(peaks))];
        baselineslope = [baselineslope local_baselineslope*ones(1,length(peaks))];
        htdiff_poly = [htdiff_poly local_htdiff_poly*ones(1,length(peaks))];
        ahtdiff_poly = [ahtdiff_poly local_ahtdiff_poly*ones(1,length(peaks))];
        pknum = [pknum 1:length(peaks)];
        pkorder = [pkorder i*ones(1,length(peaks))];
        sectnum = [sectnum sectionnumber*ones(1,length(peaks))];

       
         if dispprogress == 1
          hold off; drawnow;
            figure(1);
       
        subplot(2,2,[1 2]); plot(xdata, ydata, '-'); 
        hold on;
        subplot(2,2,[1 2]); plot(peak_idx, ydata(peak_idx), '.r');
        subplot(2,2,[1 2]); plot(left_base + segmentbound(i), ydata(left_base + segmentbound(i)), '.g', right_base + segmentbound(i), ydata(right_base + segmentbound(i)), '.g')
%         if (left_base(1)-segment_med_wd<1)
%             subplot(2,2,[1 2]); plot(left_base-left_base(1)+1, ydata(left_base-left_base(1)+1), '.g', ...
%                 right_base + peak_idx(i)-segment_med_wd, ydata(right_base + peak_idx(i)-segment_med_wd), '.g')
%         else
%             subplot(2,2,[1 2]); plot(left_base + peak_idx(i)-segment_med_wd, ydata(left_base + peak_idx(i)-segment_med_wd), '.g', ...
%                 right_base + peak_idx(i)-segment_med_wd, ydata(right_base + peak_idx(i)-segment_med_wd), '.g')
%         end
%         subplot(2,2,[1 2]); plot(left_base + peak_idx(i)-segment_med_wd, ydata(left_base + peak_idx(i)-segment_med_wd), '.g', ...
%             right_base + peak_idx(i)-segment_med_wd, ydata(right_base + peak_idx(i)-segment_med_wd), '.g')
        subplot(2,2,[1 2]); plot(xdata(peak_idx(i)), ydata(peak_idx(i)), '.c');                 % mark last analyzed peakset
        subplot(2,2,[1 2]); line([segmentbound' segmentbound'], [min(ydata) - 100 max(ydata) + 100], 'Color', 'k', 'LineStyle', ':')
        subplot(2,2,[1 2]); plot(xdata, ydata_thres, 'k')
        set(gca,'XLim',[0 length(xdata)]);
        set(gca,'YLim',[min(ydata) - 10 max(ydata) + 10]); hold off;
   
                                                                                            % Plot all peak heights over time in bottom right graph
        if sum(datalast)==0
            subplot(2,2,4); plot(elapsed_time + pk_t, pkht_poly, '.')                       % add each new analyzed point to this existing plot
            hold on
        else
            subplot(2,2,4); plot(datalast(1,2:end), datalast(2,2:end), '.')                     % plot all previous data points
            hold on
            subplot(2,2,4); plot(elapsed_time + pk_t, pkht_poly, '.')                       % add each new analyzed point to this existing plot
        end
        
        subplot(2,2,4); plot(pk_t((end - length(peaks) + 1):end) + elapsed_time, pkht_poly((end - length(peaks) + 1):end), 'or') % highlight last point
        hold off
        
         else
            clc
         end
        
        if(analysismode ~= 1)
            input('Hit ENTER to continue with analysis, CTRL+C to stop analysis......  ');
        end

       
        
        i=i+1;
    end
    
    else
        %disp('no good peaks in this segment or number of peaks is not 3');
     
    end
    
    else
       %disp('JK: skipping this segment. Too short for baseline selection');
    end
    
end
%%% assign peak number in data matrix (1,2,3)

if (~isempty(pk_t))
    data(1,:) = pk_t ; %LabVIEW record computer real-time
    data(2,:) = pkht_poly; %peak height
    data(3,:) = peakwidth; %peak width
    data(4,:) = pk_leftbase;
    data(5,:) = pk_rightbase;
    data(6,:) = baselinedist; %baseline left to right distance
    data(7,:) = baselineslope; %1st order baseline fit slope before flattening
    data(8,:) = apkht_poly;    %node deviation
    data(9,:) = ahtdiff_poly;   %FWHM
    data(10,:) = sectnum;
    data(11,:) = pknum;
    data(12,:) = pkorder;
    data(13,:) = vv(pkidx_poly);
    %save('data.csv','data','-ascii','-double','-append');                                   % save the current processing
else
    %disp('No GOOD peak found in this section');
    data=zeros(13,1);
end

elapsed_time = elapsed_time + t(end);                                            % update initial time for the next iteration

%disp(' ')

return

end