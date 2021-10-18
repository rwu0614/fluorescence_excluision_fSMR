clear all
close all
clc
currentFolder = pwd;
addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');


%%

%Grabbing file path for PMT data
n_pmt_channel = 1; %number of pmt channels
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
fprintf('\nGetting PMT time data...\n')
[input_info.pmt_time_filename, input_info.pmt_time_dir, exist_pmt_time] = uigetfile('../*.*','Select time File',' ');
if(exist_pmt_time ~= 0)
    time_file_ID = fopen(strcat(input_info.pmt_time_dir, input_info.pmt_time_filename), 'r', 'b');
    fprintf('\n%s selected for analysis\n', input_info.pmt_time_filename)
else
    fprintf('Quitting analysis program now...\n')
    return
end

fprintf('\nGetting SMR data...\n')
[smr_filename, smr_dir, exist_smr] = uigetfile('../*.*','Select SMR File',' ');
if(exist_smr ~= 0)
    smr_file_ID = fopen(strcat(smr_dir, smr_filename), 'r');
    fprintf('\n%s selected for analysis\n', smr_filename)
else
    fprintf('Quitting analysis program now...')
    return
end

%grab sample name from input file name
name_split = strsplit(smr_dir,'\');
sample_name = name_split{end-1};
sample_name= strrep(sample_name,'_',' '); % this variable is for annotating output filename

% user input for smr time file location
fprintf('\nGetting time data...\n')
[smr_time_filename, smr_time_dir, exist_time] = uigetfile('../*.*','Select time File',' ');
if(exist_time ~= 0)
    smr_time_file_ID = fopen(strcat(smr_time_dir, smr_time_filename), 'r', 'b');
    fprintf('\n%s selected for analysis\n', smr_time_filename)
else
    fprintf('\nContinuing analysis without time data...\n')
end


%% Main analysis

% seek data for current segement, datatype int is 8bytes

fseek(pmt_file_ID(1),200000, 'bof');    
rawdata_pmt= fread(pmt_file_ID(1),200000,'float64=>double');

fseek(time_file_ID,200000 ,'bof');
rawdata_time_pmt = fread(time_file_ID,200000,'float64=>double');


rawdata_smr = fread(smr_file_ID, 'int32','b')/2^32;
rawdata_smr(1:129:end)=[];        
rawdata_smr = rawdata_smr';
rawdata_smr = diff(atan(rawdata_smr(1:2:end)./rawdata_smr(2:2:end)));                
rawdata_smr(rawdata_smr > pi/2) = rawdata_smr(rawdata_smr > pi/2) - pi;
rawdata_smr(rawdata_smr <= -pi/2) = rawdata_smr(rawdata_smr <= -pi/2) + pi;
rawdata_smr = rawdata_smr*10000/(2*pi);
fclose(smr_file_ID);

rawdata_smr_time = fread(smr_time_file_ID,'float64=>double');
rawdata_smr_time(1)=[];
%%
pmt_time_diff = diff(rawdata_time_pmt);
pmt_time_stamp_ind = find(pmt_time_diff~=0);
time_x_ind = 1:1:length(rawdata_time_pmt);
pmt_time_interp = interp1(time_x_ind(pmt_time_stamp_ind),rawdata_time_pmt(pmt_time_stamp_ind),time_x_ind);
pmt_time_interp(1) = pmt_time_interp(2);
pmt_time_interp(end-5:end) = pmt_time_interp(end-6) ;
%%
figure(1)
initial_ind = 1:1:200000;
plot(rawdata_smr_time(200000:1:300000),rawdata_smr(200000:1:300000))
ylim([-500,400])
hold on
yyaxis right
plot(pmt_time_interp(initial_ind),rawdata_pmt(initial_ind))
ylim([0,5])
%%
fig = figure;
left_color = [0 0.4470 0.7410];
right_color = [0.4660 0.6740 0.1880];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
rawdata_smr_norm = rawdata_smr-median(rawdata_smr(1:1:200000));
t_entrance = 3.713265608673*10^12;
syn_smr_time = rawdata_smr_time*1000-t_entrance;
[~,t_smr_enter_ind] = min(abs(syn_smr_time));
syn_pmt_time = pmt_time_interp*1000-t_entrance;
[~,t_pmt_enter_ind] = min(abs(syn_pmt_time));
signal_frame_smr = 500;
signal_frame_pmt = signal_frame_smr*2;
signal_offset_smr =100;
signal_offset_pmt = signal_offset_smr*2;

yyaxis left
plot(syn_smr_time(t_smr_enter_ind-signal_offset_smr:(t_smr_enter_ind+signal_frame_smr)),...
    rawdata_smr_norm(t_smr_enter_ind-signal_offset_smr:(t_smr_enter_ind+signal_frame_smr)),'LineWidth',2)
ylim([-200,180])
ylabel('SMR resonance frequency (Hz)')
hold on
yyaxis right
plot(syn_pmt_time(t_pmt_enter_ind-signal_offset_pmt:(t_pmt_enter_ind+signal_frame_pmt)),...
    rawdata_pmt(t_pmt_enter_ind-signal_offset_pmt:(t_pmt_enter_ind+signal_frame_pmt)),'LineWidth',2)
ylim([0.90,1.65])
ylabel('FITC PMT detection voltage (V)')
hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',-90,'VerticalAlignment','bottom')
xlim([-signal_offset_smr/10-1,signal_frame_smr/10])
xlabel('Time from entrance to SMR cantilever (ms)')
set(gca,'FontSize',15)
title("Example SMR-PMT signal of a RPE cell")
%% Noise characterization
segment_size = 100000;
rawdata_pmt = rawdata_pmt(end-segment_size:end);
baseline = median(rawdata_pmt);
x_norm = rawdata_pmt - median(rawdata_pmt);
noise_level = rms(x_norm);
disp(baseline)
disp(noise_level)
% disp(noise_level/baseline)
%% filtering
figure(1)
med_filt_length = 300;
x_medfilt = medfilt1(x_norm,med_filt_length);

%Ploting
plot(rawdata_pmt)
hold on
plot(x_medfilt)

%% FFT on noise signal
figure(2)
Fs = 10000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(rawdata_pmt);             % Length of signal
t = (0:L-1)*T;        % Time vector

Y = fft(x_norm);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
















 