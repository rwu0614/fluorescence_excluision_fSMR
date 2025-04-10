clear all
close all 
clc
currentFolder = pwd;
addpath('plotting_functions\');
addpath('analysis_functions\');
%%
% Input UI to grab path to a single readout paired sample txt file
fprintf('\nGetting paired sample...\n')
[input_info.sample_filename, input_info.sample_dir, ~] = uigetfile('../*.*','Select Readout_paired_[sample name].txt',' ');
sample_path = [input_info.sample_dir,'\',input_info.sample_filename];
opts = detectImportOptions(sample_path,'ReadVariableNames',true,'VariableNamingRule','preserve','Delimiter','\t');
sample = readtable(sample_path,opts);
name_split = strsplit(input_info.sample_dir,'\');   
sample_name = name_split{end-1};   
sample_name= strrep(sample_name,'_',' ');

%%
figure('Color','w')
scatter(sample.buoyant_mass_pg,sample.pmt2_mV,5,'filled','MarkerFaceAlpha',0.5)
set(gca,'yscale','log')

%%
figure('Color','w')
scatter(sample.buoyant_mass_pg,sample.pmt5_mV,5,'filled','MarkerFaceAlpha',0.5)
set(gca,'yscale','log')
live_draq7_thresh = 15;
live_in = find(sample.pmt5_mV<live_draq7_thresh);

%%
figure('Color','w')
histogram(log10(sample.pmt2_mV(live_in)),'NumBins',50)

%%
gfp_thresh = 0.9;

livsc_gfp_pos_ind = find(sample.pmt5_mV<live_draq7_thresh&sample.pmt2_mV>gfp_thresh);
livsc_gfp_neg_ind = find(sample.pmt5_mV<live_draq7_thresh&sample.pmt2_mV<=gfp_thresh);
figure('Color','w')
scatter(sample.buoyant_mass_pg(livsc_gfp_pos_ind),sample.pmt2_mV(livsc_gfp_pos_ind),5,'filled','MarkerFaceAlpha',0.5)
hold on
scatter(sample.buoyant_mass_pg(livsc_gfp_neg_ind),sample.pmt2_mV(livsc_gfp_neg_ind),5,'filled','MarkerFaceAlpha',0.5)

set(gca,'yscale','log')

%%
figure('Color','w')
scatter(sample.buoyant_mass_pg(livsc_gfp_pos_ind),sample.pmt5_mV(livsc_gfp_pos_ind),5,'filled','MarkerFaceAlpha',0.5)
hold on
scatter(sample.buoyant_mass_pg(livsc_gfp_neg_ind),sample.pmt5_mV(livsc_gfp_neg_ind),5,'filled','MarkerFaceAlpha',0.5)

set(gca,'yscale','log')

%% break to bm bins
bin_num=30;
log_bm_bin_edge = linspace(prctile(log10(sample.buoyant_mass_pg(live_in)),1.5),prctile(log10(sample.buoyant_mass_pg(live_in)),99.8),bin_num+1);
bm_bin_edge = 10.^log_bm_bin_edge;


freq_matrix = zeros(1,bin_num);
count_matrix = zeros(1,bin_num);
count_matrix_gfppos = zeros(1,bin_num);
count_matrix_gfpneg = zeros(1,bin_num);

    [count_matrix(:),~] = histcounts(log10(sample.buoyant_mass_pg(live_in)),log_bm_bin_edge);
    [freq_matrix(:),~] = histcounts(log10(sample.buoyant_mass_pg(live_in)),log_bm_bin_edge,'Normalization','pdf');

    [count_matrix_gfppos(:),~] = histcounts(log10(sample.buoyant_mass_pg(livsc_gfp_pos_ind)),log_bm_bin_edge);

    [count_matrix_gfpneg(:),~] = histcounts(log10(sample.buoyant_mass_pg(livsc_gfp_neg_ind)),log_bm_bin_edge);

figure('Color','w')

    histogram(log10(sample.buoyant_mass_pg(live_in)),'binedges',log_bm_bin_edge,'Normalization','pdf');


%% plot normalized gfp + vs - bm histo
figure('Color','w')

    
    histogram(sample.buoyant_mass_pg(live_in),'binedges',bm_bin_edge,'Normalization','pdf');
    hold on
    histogram(sample.buoyant_mass_pg(livsc_gfp_pos_ind),'binedges',bm_bin_edge,'Normalization','pdf');
    set(gca,'XScale','log')

    xlabel('Buoyant mass (pg)')
    ylabel('Probability density')

%% stacked histogram of gfp+ and - counts
figure('Color','w')

    bar(log_bm_bin_edge(1:end-1)',[count_matrix_gfppos(:),count_matrix_gfpneg(:)],1,'stacked');
    xlabel('Log10 (Buoyant mass)')
    ylabel('Count')


%% ratio of gfp+ cells per bin

ratio_matrix_gfppos = 100*count_matrix_gfppos./count_matrix;

figure('Color','w')

    bar(log_bm_bin_edge(1:end-1)',ratio_matrix_gfppos(:),1);
    xlabel('Log10 (Buoyant mass)')
    ylabel('%GFP+ cells')
    ylim([0,80])
    
%%
figure('Color','w')
scatter(sample.buoyant_mass_pg,sample.vol_au./sample.buoyant_mass_pg,5,log10(sample.pmt5_mV),'filled','MarkerFaceAlpha',0.5)
colorbar
colormap jet
caxis([1.5, 2]);
    
%%
figure('Color','w')
scatter(sample.buoyant_mass_pg(livsc_gfp_pos_ind),sample.vol_au(livsc_gfp_pos_ind)./sample.buoyant_mass_pg(livsc_gfp_pos_ind),5,log10(sample.pmt5_mV(livsc_gfp_pos_ind)),'filled','MarkerFaceAlpha',0.5)
colorbar
colormap jet
caxis([1.5, 2]);
%%
figure('Color','w')
histogram(sample.vol_au(livsc_gfp_pos_ind)./sample.buoyant_mass_pg(livsc_gfp_pos_ind),'BinWidth',0.02,'Normalization','pdf');
hold on
histogram(sample.vol_au(livsc_gfp_neg_ind)./sample.buoyant_mass_pg(livsc_gfp_neg_ind),'BinWidth',0.02,'Normalization','pdf');

legend(["GFP+","GFP-"])
xlim([0.25,1.5])

%%
figure('Color','w')
histogram(sample.node_deviation_hz(livsc_gfp_pos_ind),'BinWidth',0.1,'Normalization','pdf');
hold on
histogram(sample.node_deviation_hz(livsc_gfp_neg_ind),'BinWidth',0.1,'Normalization','pdf');

legend(["GFP+","GFP-"])
xlim([-2,2])
%%
figure('Color','w')
scatter(sample.pmt2_mV(livsc_gfp_pos_ind),sample.vol_au(livsc_gfp_pos_ind)./sample.buoyant_mass_pg(livsc_gfp_pos_ind),5,'filled','MarkerFaceAlpha',0.5)

set(gca,'XScale','log')

%%
figure('Color','w')
scatter(sample.pmt2_mV(livsc_gfp_pos_ind),sample.buoyant_mass_pg(livsc_gfp_pos_ind),5,'filled','MarkerFaceAlpha',0.5)

set(gca,'XScale','log')
