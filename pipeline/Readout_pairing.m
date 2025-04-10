function pairing_stats=Readout_pairing(varargin)

if nargin ==0
    clear;
    clc;
    close all;
else
    input_dir = varargin{1};
    analysis_params_pair = varargin{2};
    output_dir = varargin{3};
end

currentFolder = pwd;
scrsize = get(0, 'Screensize');
addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');

    
%% Loading SMR and PMT readout
if nargin ==0
    fprintf('\nGetting SMR data...\n')
    [input_info.smr_filename, input_info.smr_dir, exist_smr] = uigetfile('../*.*','Select SMR File',' ');
    cd(input_info.smr_dir)
    smr_input = readtable(input_info.smr_filename);
    cd(currentFolder)

    name_split = strsplit(input_info.smr_dir,'\');   
    sample_name = name_split{end-1};   
    sample_name= strrep(sample_name,'_',' ');   

    fprintf('\nGetting PMT data...\n')
    [input_info.pmt_filename, input_info.pmt_dir, exist_pmt] = uigetfile('../*.*','Select PMT File',' ');
    cd(input_info.pmt_dir)
    pmt_input = readtable(input_info.pmt_filename);
    cd(currentFolder)

    % Ask user for PMT analysis params txt file
    fprintf('\nGetting PMT analysis params txt file...')
    [input_info.params_filename, input_info.params_dir,~] = uigetfile('../*.*','Select PMT analysis params .txt File', ' ');
    fprintf('\n%s selected for analysis\n', input_info.params_filename)
    analysis_params_file_fullname = strcat(input_info.params_dir, input_info.params_filename); 
    analysis_params = readtable(analysis_params_file_fullname,'ReadRowNames',true,'ReadVariableNames',true,'Delimiter',' ');
else
    %fprintf('\nGetting SMR data...\n')
    input_info.smr_dir = input_dir;
    S = dir(fullfile(input_info.smr_dir ,sprintf('*readout_smr*.txt')));
    input_info.smr_filename = S.name;
    cd(input_info.smr_dir)
    smr_input = readtable(input_info.smr_filename);
    cd(currentFolder)

    name_split = strsplit(input_info.smr_dir,'\');   
    sample_name = name_split{end-1};   
    sample_name= strrep(sample_name,'_',' ');   

    %fprintf('\nGetting PMT data...\n')
    input_info.pmt_dir = input_dir;
    S = dir(fullfile(input_info.pmt_dir  ,sprintf('*compensate*.txt')));
    input_info.pmt_filename = S.name;
    cd(input_info.pmt_dir)
    pmt_input = readtable(input_info.pmt_filename);
    cd(currentFolder)
    
    %fprintf('\nGetting PMT analysis params data...\n')
    input_info.params_dir = input_dir;
    S = dir(fullfile(input_info.params_dir ,sprintf('*readout_pmt_analysis_params*.txt')));
    input_info.params_filename = S.name;
    analysis_params_file_fullname = strcat(input_info.params_dir, input_info.params_filename); 
    analysis_params = readtable(analysis_params_file_fullname,'ReadRowNames',true,'ReadVariableNames',true,'Delimiter',' ');
end

pmt_timestamp=pmt_input.real_time_sec*1000; % ms
smr_timestamp=smr_input.real_time_sec*1000; % ms
%% SMR data conversion from hz to pg
if nargin ==0
    smr_data.chipID = input("Chip ID?\n",'s');
    smr_data.Hz2pg_conversion_factor = input("Hz to pg conversion factor?\n");
else
    smr_data.chipID = analysis_params_pair.chip_id;
    smr_data.Hz2pg_conversion_factor = analysis_params_pair.hz2pg_factor;
end

smr_input.buoyant_mass_pg = smr_input.buoyant_mass_hz*smr_data.Hz2pg_conversion_factor; % convert from Hz to pg


%% Window of PMT detection with SMR time-cordinates as origins
format long g
trace=[];

window = 3000;

c =1;
hit = 0;
n = 1;
for j = 1:length(smr_timestamp)
    for i= 1:length(pmt_timestamp)
        run_trace_diff = smr_timestamp(j) - pmt_timestamp(i);
            if run_trace_diff<window && run_trace_diff>=-window
                trace(c) = smr_timestamp(j)-pmt_timestamp(i);
                c =c+1;
                hit = 1;
            end
    end
    if hit==0
        n =n+1;
    end
    hit = 0;
end

%%

if nargin ==0
    figure('OuterPosition',[0.2*scrsize(3) 0.3*scrsize(4) 0.6*scrsize(3) 0.5*scrsize(4)]);
    subplot(2,1,1)
    %[counts,edges] = histcounts(trace,2*window);
    [counts,edges] = histcounts(trace,500);
    x = edges(1:end-1);
    y = zeros(1,length(x));
    col = counts;  % This is the color, vary with x in this case.

    surface([x;x],[y;y],[col;col],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',100);
    xlim([-800,800])
    xlabel('Delta T from PMT to SMR signals (ms)')
    set(gca, 'YTick', [])
    c=colorbar;
    set(get(c,'title'),'string','PMT occurrence','Rotation',0);

    subplot(2,1,2)
    plot_trace_mode_offset = 15;
    histogram(trace,'BinWidth',0.5)
    if abs(mode(trace))<100
        xlim([mode(trace)-plot_trace_mode_offset,mode(trace)+plot_trace_mode_offset])
    else
        xlim([-50,50])
    end
    xlabel('Delta T from PMT to SMR signals (ms)')
    ylabel('Counts')
end
%% Pairing

%set window limit for PMT-SMR pairing
if nargin ==0
    min_time_threshold = input('\nLower threshold for PMT to SMR transit time:');
    max_time_threshold = input('\nUpper threshold for PMT to SMR transit time:');
else
    min_time_threshold = analysis_params_pair.min_time_threshold;
    max_time_threshold = analysis_params_pair.max_time_threshold;   
end

paired_pmt_ind = zeros(length(pmt_timestamp),1);
paired_smr_ind = zeros(length(pmt_timestamp),1);
paired_delta_t = zeros(length(pmt_timestamp),1);
multi_count =0;
n= 1;
for i= 1:length(smr_timestamp)
    for j = 1:length(pmt_timestamp)
        run_diff = pmt_timestamp(j)-smr_timestamp(i);
            if run_diff> -max_time_threshold && run_diff<-min_time_threshold   
                paired_smr_ind(n) = i;
                paired_pmt_ind(n) = j;     
                paired_delta_t(n) = smr_timestamp(i)- pmt_timestamp(j);
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

%%
% paired_smr_ind = paired_smr_ind(end-1500:end);
% paired_pmt_ind = paired_pmt_ind(end-1500:end);
%% Generate output
prt_paired_smr = 100*length(paired_smr_ind)/height(smr_input);
prt_paired_pmt = 100*length(paired_pmt_ind)/height(pmt_input);
dropout_rate = 100*multiplet_count/length(pre_filt_paired_pmt_ind);

pairing_stats = [prt_paired_smr,prt_paired_pmt,dropout_rate,length(paired_smr_ind)];

%format follows: SMR input table merged with PMT input table, each row is a
%paired particle
readout_paired = [smr_input(paired_smr_ind,:),pmt_input(paired_pmt_ind,2:end)];
readout_paired.elapsed_time_min = (smr_input.real_time_sec(paired_smr_ind)-smr_input.real_time_sec(paired_smr_ind(1)))/60;
readout_paired.pmt2smr_transit_time_ms = paired_delta_t;
% datetime(3.715875052615539e+09, 'ConvertFrom','epochTime','Epoch','1904-01-01')
cd(input_info.pmt_dir)
out_file_name = ['Readout_paired_' sample_name '.txt'];
writetable(readout_paired,out_file_name, 'delimiter', '\t');
cd(currentFolder)

% generate pairing window plot for batch analysis
if nargin ~=0
    cd(output_dir)
    
    % create pairing window plot
    figure('color','w','position',[669,593,647,285])
    plot_trace_mode_offset = 10;
    histogram(trace,'BinWidth',0.5)
    hold on
    xline(min_time_threshold,'Color','r')
    hold on
    xline(max_time_threshold,'Color','r')
    if abs(mode(trace))<100
        xlim([min_time_threshold-plot_trace_mode_offset,max_time_threshold+plot_trace_mode_offset])
    else
        xlim([-80,80])
    end
    xlabel('Delta T from PMT to SMR signals (ms)')
    ylabel('Counts')

    out_file_name = [sample_name ' pairing window' '.png'];
    print(gcf,out_file_name,'-dpng','-r1200')
    cd(currentFolder)
end

% generate analysis report
if nargin ~=0
    input_info.pmt_dir = convertStringsToChars(input_info.pmt_dir);
end
report_dir = [input_info.pmt_dir '\' sample_name '_report\Pairing_report\'];
mkdir(report_dir)

Readout_pairing_report_v1(report_dir,input_info,sample_name,smr_data,smr_input,pmt_input,analysis_params,trace, window,min_time_threshold,max_time_threshold,readout_paired,multiplet_count);

end
