clear all;
clc;
close all;
currentFolder = pwd;
scrsize = get(0, 'Screensize');
addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');

    
%% Loading SMR and PMT readout
fprintf('\nGetting SMR data...\n')
[input_info.smr_filename, input_info.smr_dir, exist_smr] = uigetfile('../*.*','Select SMR File',' ');
cd(input_info.smr_dir)
smr_input = readtable(input_info.smr_filename);
cd(currentFolder)

smr_data.time=smr_input{:,1}*1000; % convert to ms

name_split = strsplit(input_info.smr_dir,'\');   
sample_name = name_split{end-1};   
sample_name= strrep(sample_name,'_',' ');   

fprintf('\nGetting PMT data...\n')
[input_info.pmt_filename, input_info.pmt_dir, exist_pmt] = uigetfile('../*.*','Select PMT File',' ');
cd(input_info.pmt_dir)
pmt_input = readtable(input_info.pmt_filename);
cd(currentFolder)

pmt_data.time=pmt_input{:,1}*1000; % ms
pmt_data.pmt{1}=pmt_input{:,2};
pmt_data.pmt{2}=pmt_input{:,3};
pmt_data.pmt{3}=pmt_input{:,4};
pmt_data.pmt{4}=pmt_input{:,5};
pmt_data.pmt{5}=pmt_input{:,6};

%% SMR data conversion from hz to pg
smr_data.chipID = 'G5W2 505-B15-L';
smr_data.Hz2pg_conversion_factor = 0.59553;
smr_data.smr=smr_input{:,2}*smr_data.Hz2pg_conversion_factor; % convert from Hz to pg

% fil_smr_ind = find(smr_data.smr>50);
% smr_data.smr = smr_data.smr(fil_smr_ind);
% smr_data.time = smr_data.time(fil_smr_ind);

%% Window of PMT detection with SMR time-cordinates as origins
format long g
trace=[];

window = 300;

n= 1;
m=[];
m_i=1;
c =1;
hit = 0;
n = 1;
no_match_smr =[];
for j = 1:length(smr_data.time)
    for i= 1:length(pmt_data.time)
        run_trace_diff = smr_data.time(j) - pmt_data.time(i);
            if run_trace_diff<window && run_trace_diff>=-window
                trace(c) = smr_data.time(j)-pmt_data.time(i);
                c =c+1;
                hit = 1;
            end
    end
    if hit==0
        no_match_smr(n) = j;
        n =n+1;
    end
    hit = 0;
end

match_smr = smr_data.smr(~no_match_smr);
%%
figure('OuterPosition',[0.2*scrsize(3) 0.5*scrsize(4) 0.6*scrsize(3) 0.3*scrsize(4)]);
%[counts,edges] = histcounts(trace,2*window);
[counts,edges] = histcounts(trace,500);
x = edges(1:end-1);
y = zeros(1,length(x));
col = counts;  % This is the color, vary with x in this case.

surface([x;x],[y;y],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',100);
xlim([-80,80])
xlabel('Delta T from PMT to SMR signals (ms)')
set(gca, 'YTick', [])
c=colorbar;
set(get(c,'title'),'string','PMT occurrence','Rotation',0);


%% Pairing

%set window limit for PMT-SMR pairing
min_time_threshold = input('\nLower threshold for PMT to SMR transit time:');
max_time_threshold = input('\nUpper threshold for PMT to SMR transit time:');


paired_pmt_ind = zeros(length(pmt_data.time),1);
paired_smr_ind = zeros(length(pmt_data.time),1);
paired_delta_t = zeros(length(pmt_data.time),1);
multi_count =0;
n= 1;
for i= 1:length(smr_data.time)
    for j = 1:length(pmt_data.time)
        run_diff = pmt_data.time(j)-smr_data.time(i);
            if run_diff> -max_time_threshold && run_diff<-min_time_threshold   
                paired_smr_ind(n) = i;
                paired_pmt_ind(n) = j;     
                paired_delta_t(n) = smr_data.time(i)- pmt_data.time(j);
                n=n+1;   
            end
       
    end
end

%%================ HIGHLY CRITICAL STEP========================%% 
%%======doublets and multiplets removal from paired array======%%

% logging raw unfiltered paired array
pre_filt_paired_pmt_ind = paired_pmt_ind(1:n-1);
pre_filt_paired_smr_ind = paired_smr_ind(1:n-1);
pre_filt_paired_delta_t = paired_delta_t(1:n-1);

% Step 1: Grabing SMR indices that only appeared once in the raw paired array
[~,ia_smr,ic_smr] = unique(paired_smr_ind(1:n-1));
 a_counts_smr = accumarray(ic_smr,1);
 value_counts_smr = [ia_smr, a_counts_smr];
 uni_paried_ind_smr = ia_smr(a_counts_smr==1);
 
% Step 2: updateing paired array with unique-smr indices
uni_paired_pmt_ind = paired_pmt_ind(uni_paried_ind_smr);
uni_paired_smr_ind = paired_smr_ind(uni_paried_ind_smr);
uni_paired_delta_t = paired_delta_t(uni_paried_ind_smr);

% Step 3: Grabing pmt indices that only appeared once in unique-smr paired array
[~,ia_smr_pmt,ic_smr_pmt] = unique(uni_paired_pmt_ind);
 a_counts_smr_pmt = accumarray(ic_smr_pmt,1);
 value_counts_smr_pmt = [ia_smr_pmt, a_counts_smr_pmt];
 uni_paried_ind_smr_pmt = ia_smr_pmt(a_counts_smr_pmt==1);

% Step 4: updateing paired array with unique-smr-pmt indices
paired_pmt_ind = uni_paired_pmt_ind(uni_paried_ind_smr_pmt);
paired_smr_ind = uni_paired_smr_ind(uni_paried_ind_smr_pmt);
paired_delta_t = uni_paired_delta_t(uni_paried_ind_smr_pmt);
n = length(paired_pmt_ind);
multiplet_count = length(pre_filt_paired_pmt_ind)-n;
%% Generate output
%format follows: [time of detection(computer real time), smr(pg), PacificBlue(mV),FITC(mV), PE(mV), APC(mV), Cy7(mV),PMTtoSMR transit time(ms)]
readout_paired = [smr_data.time(paired_smr_ind),smr_data.smr(paired_smr_ind),...
    pmt_input{:,2}(paired_pmt_ind),pmt_input{:,3}(paired_pmt_ind),...
    pmt_input{:,4}(paired_pmt_ind),pmt_input{:,5}(paired_pmt_ind),...
    pmt_input{:,6}(paired_pmt_ind),paired_delta_t];
cd(input_info.pmt_dir)
out_file_name = ['Readout_paired_' sample_name '.csv'];
dlmwrite(out_file_name, readout_paired, 'delimiter', ',', 'precision', 25);
cd(currentFolder)

% generate analysis report
report_dir = [input_info.pmt_dir '\' sample_name '_report\Pairing_report\'];
mkdir(report_dir)
Readout_pairing_report_v1(report_dir,input_info,sample_name,smr_data,pmt_data,trace, window,min_time_threshold,max_time_threshold,readout_paired,multiplet_count);


% %%
% figure(1)
% voltage_plot_lim_lower = 1e-2;
% voltage_plot_lim_higher = 1e+4;
% smr_high_exclude_pct = 99.5;
% smr_plot_lim_lower = 0;
% smr_plot_lim_higher = prctile(smr_data.smr,smr_high_exclude_pct);
% 
% obj2plot =smr_data.smr;
% bin_n = round(length(obj2plot)/80);
% obj2plot = obj2plot(obj2plot<smr_plot_lim_higher & obj2plot>smr_plot_lim_lower);
% hh=histogram(obj2plot,bin_n );
% hold on
% obj2plot =readout_paired(:,2);
% obj2plot = obj2plot(obj2plot<smr_plot_lim_higher & obj2plot>smr_plot_lim_lower);
% hh=histogram(obj2plot,bin_n);
% xlabel("Buoyant mass (pg)")
% title("Paired SMR signal distribution")
% xlim([smr_plot_lim_lower,smr_plot_lim_higher])
% hold off
