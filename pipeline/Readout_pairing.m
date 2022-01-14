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

% Ask user for PMT analysis params txt file
fprintf('\nGetting PMT analysis params txt file...')
[input_info.params_filename, input_info.params_dir,~] = uigetfile('../*.*','Select PMT analysis params .txt File', ' ');
fprintf('\n%s selected for analysis\n', input_info.params_filename)
analysis_params_file_fullname = strcat(input_info.params_dir, input_info.params_filename); 
analysis_params = readtable(analysis_params_file_fullname,'ReadRowNames',true,'ReadVariableNames',true,'Delimiter',' ');

pmt_data.time=pmt_input{:,1}*1000; % ms
pmt_data.pmt{1}=pmt_input{:,2};
pmt_data.pmt{2}=pmt_input{:,3};
pmt_data.pmt{3}=pmt_input{:,4};
pmt_data.pmt{4}=pmt_input{:,5};
pmt_data.pmt{5}=pmt_input{:,6};

%% SMR data conversion from hz to pg
smr_data.chipID = '505_20211119';
smr_data.Hz2pg_conversion_factor = 0.64982;
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

%%
% paired_smr_ind = paired_smr_ind(end-1500:end);
% paired_pmt_ind = paired_pmt_ind(end-1500:end);
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
Readout_pairing_report_v1(report_dir,input_info,sample_name,smr_data,pmt_data,analysis_params,trace, window,min_time_threshold,max_time_threshold,readout_paired,multiplet_count);



%% Below are for a quick look at the paired data from fxm experiment with FITC dextran

%% For trapping
% ex_x = smr_data.time(paired_smr_ind)/1000;
% ex_x = ex_x-ex_x(1);
% 
% ex_y = pmt_input{:,3}(paired_pmt_ind);
% 
% ex_smr = smr_data.smr(paired_smr_ind);
% 
% %cell_1_ind = find(ex_x>0&ex_x<200&ex_y>100);
% cell_1_ind = find(ex_x>6200&ex_x<6800&ex_y<80);
% %cell_1_ind = find(ex_x>325&ex_x<600&ex_y>0&ex_y<300);
% %cell_1_ind = find(ex_x>2000&ex_x<2200&ex_y>100&ex_y<400);
% 
% cell_1_ind = cell_1_ind(1:1:end);
% cell_1_x = ex_x(cell_1_ind);
% cell_1_y = ex_y(cell_1_ind);
% 
% % cell_1_y_base = medfilt1(cell_1_y,100);
% % cell_1_y_base(1) = cell_1_y_base(2);
% % cell_1_y_base = cell_1_y_base-cell_1_y_base(1);
% 
% cell_1_y_mod = cell_1_y;
% 
% figure(1)
% scatter(ex_x,ex_y,'.');
% hold on
% scatter(cell_1_x,cell_1_y,'.');
% hold off
% 
% figure(2)
% scatter(cell_1_x,ex_smr((cell_1_ind)),'.');
% 
% figure(3)
% scatter(cell_1_x,ex_smr((cell_1_ind))./cell_1_y,'.');
% 
% cell_1_std = std(cell_1_y);
% cell_1_average = mean(cell_1_y);
% cell_1_cv = cell_1_std/cell_1_average;
% disp(cell_1_cv*100)
% 


%%
%2.4485634
%8.8485634
median_vol_real = 1100; %fL for L1210
median_vol_au = median(pmt_input{:,3}(paired_pmt_ind));
PMT_to_pL_conversion_factor = median_vol_real/median_vol_au; %um3
%PMT_to_pL_conversion_factor=22.5483499895012;
real_vol = pmt_input{:,3}(paired_pmt_ind)*PMT_to_pL_conversion_factor;
real_dia = (real_vol*6/pi).^(1/3);
real_density = smr_data.smr(paired_smr_ind)./real_vol+1.005584+0.018;
figure (1)
scatter(real_vol,real_density)
xlabel('Volume (fL)')
ylabel('Density (g/cm3)')
figure (2)
scatter(real_vol,smr_data.smr(paired_smr_ind))
xlabel('Volume (fL)')
ylabel('Buoyant mass (pg)')
figure (3)
scatter(smr_data.smr(paired_smr_ind),real_density)
xlabel('Buoyant mass (pg)')
ylabel('Density (g/cm3)')


%%
%  figure(10) %density histo
%         obj2plot =real_density;
%         bin_n = round(length(obj2plot)/40);
%         h1=histogram(obj2plot);
%         h1.Normalization = 'probability';
%         h1.BinWidth = 0.001;
%         h1.FaceColor = [0.5 0.5 0.5];
%         hold on 
%          obj2plot =FC.density_g_cm_3_;
%         bin_n = round(length(obj2plot)/40);
%         h2=histogram(obj2plot);
%         h2.Normalization = 'probability';
%         h2.BinWidth = 0.001;
%         h2.FaceColor = [0.9290 0.6940 0.1250];
%         xlabel('Density (g/cm^{3})')
%         ylabel('Probability density estimation')
% xlim([1.04,1.09])
% set(gca,'FontSize',15)
% legend('Volume exclusion','Fluid exhange')
% %title('Density measurement')
% 
%  figure(11) %density histo
%         obj2plot =smr_data.smr(paired_smr_ind)*1.2;
%         bin_n = round(length(obj2plot)/40);
%         h1=histogram(obj2plot);
%         h1.Normalization = 'probability';
%         h1.BinWidth = 0.5;
%         h1.FaceColor = [0.8500 0.3250 0.0980];
%         xlabel('Buoyant mass (pg)')
%         ylabel('Probability density estimation')
%         xlim([0,40])
%         set(gca,'FontSize',15)
%  figure(12) %density histo
%         obj2plot =real_vol;
%         bin_n = round(length(obj2plot)/40);
%         h1=histogram(obj2plot);
%         h1.Normalization = 'probability';
%         h1.BinWidth = 30;
%         h1.FaceColor = [0.8500 0.3250 0.0980];
%         
%         xlabel('Volume (fL)')
%         ylabel('Probability density estimation')
%         xlim([0,1500])
%         set(gca,'FontSize',15)
%%
% figure(13)
% scatter(smr_data.smr(paired_smr_ind)*1.2,real_vol,20,"filled", 'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',.4) 
%  hold on
%  scatter(FC.BuoyantMassInNormalRPMI_pg_,FC.volume_fL_,20,"filled", 'MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerFaceAlpha',1) 
% %dscatter(real_vol,real_density,'SMOOTHING',10,'BINS',[3000,2000],'PLOTTYPE','scatter','msize',8) 
% 
% ylabel('Volume (fL)')
% xlabel('Buoyant mass (pg)')
%  xlim([20,110])
%   ylim([400,2000])
% legend('Volume exclusion','Fluid exhange')
% set(gca,'FontSize',15)
% pbaspect([1 1 1])
%%
% figure(15)
% 
% scatter(real_vol,real_density,6,smr_data.smr(paired_smr_ind)*1.2,'filled','MarkerFaceAlpha',0.8)
% 
% cmap = crameri('roma');
% cmap = colormap(parula(5));
% %cmap =  tab20(3);
% %set(gca,'ColorScale','log')
% colormap(hsv(100));
% c=colorbar;
% caxis([6 35]);
% set(get(c,'title'),'string','Buoyant mass (pg)','Rotation',0);
% ylabel('Density (g/cm^{3})')
% xlabel('Volume (fL)')
% 
% xlim([0,1500])
% ylim([1.02,1.1])
% set(gca,'FontSize',15)
% legend('n = 5547')
%%
figure(14)
scatter(pmt_input{:,1}(paired_pmt_ind),real_density,10,"filled", 'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',.5) 
apc = pmt_input{:,4}(paired_pmt_ind);
%dscatter(smr_data.smr(paired_smr_ind),log2(apc),'SMOOTHING',15,'BINS',[3000,2000],'PLOTTYPE','scatter') 

ylabel('Density (g/cm^{3})')
xlabel('Volume (fL)')
legend('Volume exclusion')
title('L1210 density measurement')
ylim([1.0,1.2])
% %%
% CV_volexclusion = std(real_density)/mean(real_density);
% disp(CV_volexclusion)
% CV_fluidexchange = std(fc)/mean(fc);
% disp(CV_fluidexchange)


