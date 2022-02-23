% SMR_readout
% This file is for processing raw SMR .bin files and outputs a csv file that
% contains peak heights of every detected smr signal as well as the time stamp
% of the peak. User input required for both the smr frequency .bin file and 
% smr time .bin file that should have been simulaneously written by the fSMR
% LabVIEW DataReciever VI. This pipeline can not be used for only analyzing
% SMR raw data without matched time file.
% 
% Output file will be saved in the same directory as the raw smr frequency
% file. In the output matrix, each row contains features of the same 
% detected peak, while each column contains a defined feature.
% Column features are as following:
%     1st column = time stamp from labview recorded computer real-time
%     2nd column = peak height in Hz of the first peak from the smr signal

%%
clear all
close all
format long;
currentFolder = pwd;
addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');


global datafull
global x
global elapsed_time;
global samplepeak;
global sampletime;
global datasmr;
elapsed_time = 0;
samplepeak=[];
sampletime=[];


%% Initialization

%=========================== Read data============================%
% user input for smr frequency file location
fprintf('\nGetting SMR data...\n')
[input_info.smr_filename, input_info.smr_dir, exist_smr] = uigetfile('../*.*','Select SMR File',' ');
if(exist_smr ~= 0)
    smr_file_ID = fopen(strcat(input_info.smr_dir, input_info.smr_filename), 'r');
    fprintf('\n%s selected for analysis\n', input_info.smr_filename)
else
    fprintf('Quitting analysis program now...')
    return
end

%grab sample name from input file name
name_split = strsplit(input_info.smr_dir,'\');
sample_name = name_split{end-1};
sample_name= strrep(sample_name,'_',' '); % this variable is for annotating output filename

% user input for smr time file location
fprintf('\nGetting time data...\n')
[input_info.smr_time_filename, input_info.smr_time_dir, exist_time] = uigetfile('../*.*','Select time File',' ');
if(exist_time ~= 0)
    smr_time_file_ID = fopen(strcat(input_info.smr_time_dir, input_info.smr_time_filename), 'r', 'b');
    fprintf('\n%s selected for analysis\n', input_info.smr_time_filename)
else
    fprintf('\nContinuing analysis without time data...\n')
end
%=================================================================%


%%=======================Raw data pre-processing=======================%
%%This part converts the raw smr binary file into 1D time-series 
%%freqeuncy readout array. Depending on the LabVIEW SMRDriverClient VI
%%setting that created the raw data (ie PLL or Feedback mode), only part of this
%%code block is used while the other part is commented out

%THIS IS FOR PLL
%rawdata_smr = fread(smr_file_ID, 'int32','b')/2^32*1.25e7;
%rawdata_smr(1:129:end)=[];
%fclose(smr_file_ID);

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

%set the size of each data block to be analyzed one at a time       
datasize = 2e6;
num_blocks = ceil(length(rawdata_smr)/datasize);
%=================================================================%


%=========================Analysis and display mode =======================%
%%User inputs to choose analysis and display mode, rapid analysis (if yes) means the 
%%pipeline will go through the entire dataset without seeking user approval 
%%for passing an identified smr signal. Displaying progress (if yes) will
%%display individual detected peak shapes as well as all detected peaks from
%%the whole data block under analysis(data block size defined by "datasize" variable)

analysismode = input('Rapid analysis mode? (1 = Yes, 0 = No):    ');
if analysismode == 1
    dispprogress = input('Display progress? (1 = Yes, 0 = No):       ');
else
    dispprogress = 1;
end

% This is for passing parameters to downstream functions
analysisparams.analysismode = analysismode;       
analysisparams.dispprogress = dispprogress;

% Set real-time reporting screensize
if dispprogress == 1
    scrsize = get(0, 'Screensize');
    figure('OuterPosition',[0 0.05*scrsize(4) scrsize(3) 0.95*scrsize(4)])
end
%======================================================================%


%% Main analysis

%=======================Main analysis loop=======================%
% Initialize global variable to store all features of all detected smr
% signal from the entire dataset
datatest = zeros(13,1);
datafull = zeros(13,1); %total of 13 features

%Initialize loop to go through each data block (size defined by "datasize"
%variable)
loop=0;

%Initialize estimated number of datapoints want to test
estimated_datapoints = 50:50:800;

datasizetest = 5e5;
num_peaks_compiled = [];

for i = 1:length(estimated_datapoints)
    % Creat a column vector of the frequency data in the current block being analyzed  
    rawdata_smr_current_test = rawdata_smr([loop*datasizetest+1:min(length(rawdata_smr), (loop+1)*datasizetest)]);
    % Creat a column vector of the time data in the current block being analyzed 
    rawdata_smr_time_current_test = rawdata_smr_time([loop*datasizetest+1:min(length(rawdata_smr), (loop+1)*datasizetest)]);
    % Detect peaks from current block
    datalasttest = S1_PeakAnalysis_time(-rawdata_smr_current_test', rawdata_smr_time_current_test, datatest, loop, analysisparams, estimated_datapoints(i));
    num_peaks_compiled(i) = length(datalasttest)./3;
end

[~,optimized_idx] = max(num_peaks_compiled);
estimated_datapoints_best = estimated_datapoints(optimized_idx);

% estimated_datapoints_best = 150;

while(1)
%     if loop ==2
%         loop = loop+2;
%     end
    
    % Creat a column vector of the frequency data in the current block being analyzed  
    rawdata_smr_current = rawdata_smr([loop*datasize+1:min(length(rawdata_smr), (loop+1)*datasize)]);
    % Creat a column vector of the time data in the current block being analyzed 
    rawdata_smr_time_current = rawdata_smr_time([loop*datasize+1:min(length(rawdata_smr), (loop+1)*datasize)]);
    % Detect peaks from current block
    [datalast, analysis_params] = S1_PeakAnalysis_time(-rawdata_smr_current', rawdata_smr_time_current, datafull, loop, analysisparams, estimated_datapoints_best);
    analysis_params.estimated_datapoints_optimized = estimated_datapoints_best;
    
    
    % Compile all detected indiviaul peaks from each segement through the
    % loop. NOTE: one smr single from 2nd mode vibration will have 3 peaks
    % and thus 3 consecutive rows in the datafull matrix
    datafull = [datafull datalast]; 
     
    loop= loop+1; %iterate to next data block
    
    % STOP if loop reaches end of main file
    if length(rawdata_smr_current) < datasize
        break
    end
end
%=================================================================%
%%
%=======================Merging and filtering=======================%
%Run the merge function to get cell by peak-feature matrix and save
S3_Merge;
%save('data.mat', 'datasmr');    
%Run select_smrpeaks_fast to filter out bad peaks and save filtered matrix
select_smrpeaks_fast;
%save('data.mat', 'datasmr_good');    
%=================================================================%


%% Output readout csv file
% datasmr_good format is [tm' tm'/60 mm' bm' bs' m1' m2' m3' nd1' nd2' ndm' w' bd' vs' sectnum' tm'/3600 mm'/2 pkorder' ndm'./mm'];

output_smr = [datasmr_good(:,1),datasmr_good(:,6)]; %time, first peakhight

cd(input_info.smr_dir)
out_file_name = ['readout_smr_' sample_name '.csv'];
%Write csv file, NOTE: a high precision number is required to output
%full length time-stamps of each signal. This is crucial for accurate
%pairing to the pmt signals
dlmwrite(out_file_name, output_smr, 'delimiter', ',', 'precision', 25);
cd(currentFolder)

%% Generate report file

report_dir = [input_info.smr_dir '\' sample_name '_report\SMR_report\'];
mkdir(report_dir)
SMR_readout_report_v1(report_dir,input_info,sample_name, datasmr_good, output_smr, analysis_params)


%% Create report figures
scrsize = get(0, 'Screensize');
out_fig = figure('OuterPosition',[0 0.05*scrsize(4) scrsize(3) 0.95*scrsize(4)]);
%plotting time-series mean buoyant mass in Hz as well as the baseline amplitude
sub_1 = subplot(3,1,1);
    [ax, h1, h2]=plotyy(datasmr_good(:,1), datasmr_good(:,3), datasmr_good(:,1), datasmr_good(:,4));
    set(h1, 'LineStyle', '-.');
    set(h2, 'LineStyle', '-');
    xlabel('Time (computor real-time)')
    ylabel('Frequency shifts (Hz)')
    hold on
    yyaxis right
    ylabel('Baseline amplitude (Hz)')
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',-90,'VerticalAlignment','bottom')

%plotting 
Rough_datarate = 20000; %For Feedback mode CIC rate = 10000, datarate is around 20000
dt_data = 1000/Rough_datarate;% unit = ms, delta t between two consecutive raw data point
keep_ind = find(datafull(3,:)<10000); %filter on transit_time
transit_fil_datafull = datafull(:,keep_ind);
keep_ind_2 = find(transit_fil_datafull(2,:)<8000); %filter on buoyant mass, 400pg cap
mass_transit_fil_datafull = transit_fil_datafull(:,keep_ind_2);
mass_transit_fil_datafull(:,1) =[];
peak_time = mass_transit_fil_datafull(1,:)-mass_transit_fil_datafull(1,1);
peak_time_s = peak_time/60;

sub_2 = subplot(3,1,2);
    scatter(peak_time_s,mass_transit_fil_datafull(3,:)*dt_data,5,mass_transit_fil_datafull(2,:),'fill');c=colorbar;
    set(get(c,'title'),'string','Buoyant mass(pg)','Rotation',0);
    xlabel('Time (min)')
    ylabel('Transit time (ms)')
    hold on
    yyaxis right
    est_timedomain = 0:round(peak_time_s(end)/6):peak_time_s(end);
    lower_bound =  est_timedomain(1:end-1);
    upper_bound = est_timedomain(2:end);
    for i=1:length(lower_bound)
        lower_t = lower_bound(i);
        upper_t = upper_bound(i);
        cell_count = (1/3)*length(find(peak_time_s >lower_t&peak_time_s <upper_t));
        est_30min = 30*cell_count/(upper_t-lower_t);
    %     fprintf('\n#cells per min (%1.2f to %1.2f min): %1.0f, estimation for 30min= %1.0f\n ',lower_t,upper_t,...
    %         cell_count(i),est_30min(i));
        plot(lower_t:0.01:upper_t,est_30min*ones(length(lower_t:0.01:upper_t),1),'-','color','[0.8500, 0.3250, 0.0980]');
        hold on
    end
    ylabel('Throughput - #cell per 30min')
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',-90,'VerticalAlignment','bottom')

sub_3 = subplot(3,1,3);
hist(mass_transit_fil_datafull(2,:),150)
xlabel('Frequency shift (Hz)')
ylabel('Counts')

%% Hz to pg calibration
beads_radi = 4*10^-4; % 4um to cm
beads_dens = 1.049; % beads density
fluid_dens = 1.003; % PBS density
BM_beads = (10^12)*(beads_dens-fluid_dens)*((4/3)*pi*(beads_radi^3)); % in pg
scaling_factor = BM_beads/median(mass_transit_fil_datafull(2,:));

% chip from 04/12/21 scaling factor is 0.60779
 
disp(estimated_datapoints_best)
 