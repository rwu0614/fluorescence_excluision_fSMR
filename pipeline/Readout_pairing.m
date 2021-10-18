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

window = 3000;

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
xlim([-800,800])
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



%%
ex_x = smr_data.time(paired_smr_ind)/1000;
ex_x = ex_x-ex_x(1);

ex_y = pmt_input{:,3}(paired_pmt_ind);

ex_smr = smr_data.smr(paired_smr_ind);

%cell_1_ind = find(ex_x>0&ex_x<200&ex_y>100);
cell_1_ind = find(ex_x>6200&ex_x<6800&ex_y<80);
%cell_1_ind = find(ex_x>325&ex_x<600&ex_y>0&ex_y<300);
%cell_1_ind = find(ex_x>2000&ex_x<2200&ex_y>100&ex_y<400);

cell_1_ind = cell_1_ind(1:1:end);
cell_1_x = ex_x(cell_1_ind);
cell_1_y = ex_y(cell_1_ind);

% cell_1_y_base = medfilt1(cell_1_y,100);
% cell_1_y_base(1) = cell_1_y_base(2);
% cell_1_y_base = cell_1_y_base-cell_1_y_base(1);

cell_1_y_mod = cell_1_y;

figure(1)
scatter(ex_x,ex_y,'.');
hold on
scatter(cell_1_x,cell_1_y,'.');
hold off

figure(2)
scatter(cell_1_x,ex_smr((cell_1_ind)),'.');

figure(3)
scatter(cell_1_x,ex_smr((cell_1_ind))./cell_1_y,'.');

cell_1_std = std(cell_1_y);
cell_1_average = mean(cell_1_y);
cell_1_cv = cell_1_std/cell_1_average;
disp(cell_1_cv*100)



%%
%2.4485634
%8.8485634
median_vol_real = 1100; %fL for L1210
median_vol_au = median(pmt_input{:,3}(paired_pmt_ind));
PMT_to_pL_conversion_factor = median_vol_real/median_vol_au; %um3
real_vol = pmt_input{:,3}(paired_pmt_ind)*PMT_to_pL_conversion_factor;
real_dia = (real_vol*6/pi).^(1/3);
real_density = smr_data.smr(paired_smr_ind)./real_vol+1.005584+0.018;
figure (1)
scatter(real_vol,real_density)
figure (2)
scatter(real_vol,smr_data.smr(paired_smr_ind))
figure (3)
scatter(smr_data.smr(paired_smr_ind),real_density)

%%
 figure(10) %density histo
        obj2plot =real_density;
        bin_n = round(length(obj2plot)/40);
        h1=histogram(obj2plot);
        h1.Normalization = 'probability';
        h1.BinWidth = 0.001;
        h1.FaceColor = [0.8500 0.3250 0.0980];
        hold on
        obj2plot =fc;
        h2=histogram(obj2plot);
        h2.Normalization = 'probability';
        h2.BinWidth = 0.001;
        h2.FaceColor = [0.3010 0.7450 0.9330];
        hold off
        legend('Volume exclusion','Fluid exchange')
        xlabel('Density (g/cm^{3})')
        ylabel('Probability density estimation')
xlim([1.04,1.09])
title('L1210 density measurement')

%%
 figure(11) %density histo
        obj2plot =smr_data.smr(paired_smr_ind)*1.2;
        bin_n = round(length(obj2plot)/40);
        h1=histogram(obj2plot);
        h1.Normalization = 'probability';
        h1.BinWidth = 5;
        h1.FaceColor = [0.8500 0.3250 0.0980];
        hold on
        obj2plot =fbm;
        h2=histogram(obj2plot);
        h2.Normalization = 'probability';
        h2.BinWidth = 5;
        h2.FaceColor = [0.3010 0.7450 0.9330];
        hold off
        legend('Volume exclusion','Fluid exchange')
        xlabel('Density (g/cm^{3})')
        ylabel('Probability density estimation')

%%
 figure(12) %density histo
        obj2plot =real_vol;
        bin_n = round(length(obj2plot)/40);
        h1=histogram(obj2plot);
        h1.Normalization = 'probability';
        h1.BinWidth = 50;
        h1.FaceColor = [0.8500 0.3250 0.0980];
        hold on
        obj2plot =fv;
        h2=histogram(obj2plot);
        h2.Normalization = 'probability';
        h2.BinWidth = 50;
        h2.FaceColor = [0.3010 0.7450 0.9330];
        hold off
        legend('Volume exclusion','Fluid exchange')
        xlabel('Density (g/cm^{3})')
        ylabel('Probability density estimation')
%%
figure(13)
scatter(real_vol,real_density,10,"filled", 'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',.5) 
hold on
scatter(fv,fc,10,"filled",'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerFaceAlpha',0.8) 
ylabel('Density (g/cm^{3})')
xlabel('Volume (fL)')
legend('Volume exclusion','Fluid exchange')
title('L1210 density measurement')

%%
CV_volexclusion = std(real_density)/mean(real_density);
disp(CV_volexclusion)
CV_fluidexchange = std(fc)/mean(fc);
disp(CV_fluidexchange)


