function rpt_log_sec3 = Readout_pairing_report_v1(report_dir,input_info,sample_name,smr_data,pmt_data,trace,window,min_time_threshold,max_time_threshold,readout_paired,multiplet_count)

% Readout_pairing_report.m creates a pdf report of the pairing results
% between smr_readout and pmt_readou. The goal for the report is to
% bookmark key parameters from the original analysis for future replication
% and reflection of the analysis.
close all;
import mlreportgen.report.*
import mlreportgen.dom.*

currentFolder = pwd;

addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');

cd (report_dir)
time_stamp = datestr(datetime("now"),"yyyymmddTHHMMSS");
rpt = Report(append(report_dir,"Readout_pairing_report_",time_stamp),'pdf');
open(rpt)
%% Initialize report structure
append(rpt,TableOfContents);
chrpt_title = Heading1('Part ');
chrpt_title.Style = {CounterInc('sect1'),CounterReset('sect2'),...
     WhiteSpace('preserve')...
     Color('black'),...
     Bold, FontSize('18pt'),OuterMargin("0pt", "0pt","0pt","0pt")};
append(chrpt_title,AutoNumber('sect1'));
append(chrpt_title,'. ');

sectrpt_title = Heading2();
sectrpt_title.Style = {CounterInc('sect2'),Bold,...
     WhiteSpace('preserve'),OuterMargin("0pt", "0pt","30pt","0pt")};
append(sectrpt_title,AutoNumber('sect1'));
append(sectrpt_title,'.');
append(sectrpt_title,AutoNumber('sect2'));
append(sectrpt_title,'. ');

paramrpt_title = Heading3();
paramrpt_title.Style = {CounterInc('sect3'),...
WhiteSpace('preserve')};
% append(paramrpt_title,AutoNumber('sect1'));
% append(paramrpt_title,'.');
% append(paramrpt_title,AutoNumber('sect2'));
% append(paramrpt_title,'. ');
% append(paramrpt_title,AutoNumber('sect3'));
% append(paramrpt_title,'. ');

paragraph_style = {Color('black'),FontFamily('Arial'),FontSize('10pt')...
    ,OuterMargin("0pt", "0pt","0pt","0pt")};
caption_style = {Color('black'),FontFamily('Arial'),FontSize('10pt')...
    ,OuterMargin("0pt", "0pt","10pt","20pt")};


%% Part 1 Analysis log
rpt_title = clone(chrpt_title);
append(rpt_title,'Analysis log');
ch1 = Chapter('title',rpt_title);

% Section 1 Input and output
rpt_log_sec1=[];
rpt_log_sec1.Time = sprintf('This report is genereated on %s\n', datetime(now,'ConvertFrom','datenum'));
rpt_log_sec1.Sample = sprintf('Sample directory: "%s" from path "%s" \n',sample_name,input_info.smr_dir);
rpt_log_sec1.SMR_input = sprintf('Input SMR data: "%s" from path "%s" \n',input_info.smr_filename,input_info.smr_dir);
rpt_log_sec1.PMT_input = sprintf('Input PMT data: "%s" from path "%s" \n',input_info.pmt_filename,input_info.pmt_dir);
fn = fieldnames(rpt_log_sec1);

rpt_title = clone(sectrpt_title());
append(rpt_title,'Input and output');
append(ch1,Section('title',rpt_title));
for i = 1:numel(fn)
    rpt_title = clone(paramrpt_title());
    append(rpt_title,fn{i});
    p = Paragraph(rpt_log_sec1.(fn{i}));
    p.Style = paragraph_style;
    append(ch1,Section('title',rpt_title,'Content',p));
    %append(rpt,p);
end

% Section 2 Key analysis parameters
rpt_log_sec2=[];
rpt_log_sec2.SMR_chip_ID = sprintf('SMR chip ID: %s\n', smr_data.chipID);
rpt_log_sec2.SMR_HztoPg_conversion_factor = sprintf('SMR hz-to-pg conversion factor: %1.5f pg/hz\n', smr_data.Hz2pg_conversion_factor);
rpt_log_sec2.Pairing_lower_cutoff = sprintf('Pairing window - lower PMT to SMR transit time cutoff = %1.5f ms\n', min_time_threshold);
rpt_log_sec2.Pairing_higher_cutoff = sprintf('Pairing window - higher PMT to SMR transit time cutoff = %1.5f ms\n', max_time_threshold);
fn = fieldnames(rpt_log_sec2);

rpt_title = clone(sectrpt_title());
append(rpt_title,'Key analysis parameters');
append(ch1,Section('title',rpt_title));
for i = 1:numel(fn)
    rpt_title = clone(paramrpt_title());
    append(rpt_title,fn{i});
    p = Paragraph(rpt_log_sec2.(fn{i}));
    p.Style = paragraph_style;
    append(ch1,Section('title',rpt_title,'Content',p));
    %append(rpt,p);
end

% Section 3 Pairing results
rpt_log_sec3 =[];
rpt_log_sec3.Paired_events = sprintf('There are %1.0f paired events\n',length(readout_paired(:,1)));
rpt_log_sec3.Percent_SMR_paired = sprintf('SMR:  %%%1.2f paired. [%1.0f out of %1.0f signals]\n', ...
    100*length(readout_paired(:,1))/length(smr_data.smr),length(readout_paired(:,1)),length(smr_data.smr));
rpt_log_sec3.Percent_PMT_paired = sprintf('PMT:  %%%1.2f paired. [%1.0f out of %1.0f signals]\n', ...
    100*length(readout_paired(:,1))/length(pmt_data.time),length(readout_paired(:,1)),length(pmt_data.time));
rpt_log_sec3.Dropout_rate = sprintf('Dropout rate (due to non-unique pairing): %%%1.2f \n', 100*multiplet_count/(length(readout_paired(:,1))+multiplet_count));
fn = fieldnames(rpt_log_sec3);

rpt_title = clone(sectrpt_title());
append(rpt_title,'Pairing results');
append(ch1,Section('title',rpt_title));
for i = 1:numel(fn)
    rpt_title = clone(paramrpt_title());
    append(rpt_title,fn{i});
    p = Paragraph(rpt_log_sec3.(fn{i}));
    p.Style = paragraph_style;
    append(ch1,Section('title',rpt_title,'Content',p));
end
disp(rpt_log_sec3.Percent_SMR_paired);
disp(rpt_log_sec3.Percent_PMT_paired);
disp(rpt_log_sec3.Dropout_rate);

% Section 4 Key plotting parameters
voltage_plot_lim_lower = 1e-2;
voltage_plot_lim_higher = 1e+4;

smr_high_exclude_pct = 99.5;
smr_plot_lim_lower = 0;
smr_plot_lim_higher = prctile(smr_data.smr,smr_high_exclude_pct);

rpt_log_sec4 =[];
rpt_log_sec4.Fluorescence_plotting_lower_cutoff = sprintf('Lower fluorescence plotting limit: %1.5f mV\n',voltage_plot_lim_lower);
rpt_log_sec4.Fluorescence_plotting_higher_cutoff = sprintf('Higher fluorescence plotting limit: %1.5f mV\n',voltage_plot_lim_higher);
rpt_log_sec4.Buoyant_mass_plotting_percetile = sprintf('Lower %1.2f pecentile of buoyant mass signals were plotted to exclude extreme outliers.\n',smr_high_exclude_pct);
rpt_log_sec4.Buoyant_mass_plotting_lower_cutoff = sprintf('Lower buoyant mass plotting limit: %1.5f pg\n',smr_plot_lim_lower);
rpt_log_sec4.Buoyant_mass_plotting_higher_cutoff = sprintf('Higher buoyant mass plotting limit: %1.5f pg\n',smr_plot_lim_higher);
fn = fieldnames(rpt_log_sec4);
 
rpt_title = clone(sectrpt_title());
append(rpt_title,'Key plotting parameters');
append(ch1,Section('title',rpt_title));
for i = 1:numel(fn)
    rpt_title = clone(paramrpt_title());
    append(rpt_title,fn{i});
    p = Paragraph(rpt_log_sec4.(fn{i}));
    p.Style = paragraph_style;
    append(ch1,Section('title',rpt_title,'Content',p));
end
append(rpt,ch1);

%% Part 2 Report figures
rpt_title = clone(chrpt_title());
append(rpt_title,'Report figures');
ch2 = Chapter('title',rpt_title);
n_pmt_channel = 5;
scrsize = get(0, 'Screensize');

%% Figure 1 Plotting timestamps from all detected smr and pmt events before pairing
rpt_title = clone(sectrpt_title());
append(rpt_title,'Plotting timestamps from all detected smr and pmt events before pairing');
append(ch2,Section('title',rpt_title));

p = Paragraph(['Expect to have largely overlapping timestamps ',...
    'between pmt and smr. Long gap(s) in one domain while the other one has ',...
    'continueous timestamps could result in lower pairing efficiency in the contineous domain']);
    p.Style = caption_style; 
append(ch2,p);

dot_size = 8;
figure('OuterPosition',[0 0.5*scrsize(4) 0.95*scrsize(3) 0.3*scrsize(4)]);
    scatter((10^-3/60)*(smr_data.time-smr_data.time(1)),ones(size(smr_data.time)),dot_size,'filled','r')
    hold on
    scatter((10^-3/60)*(pmt_data.time-smr_data.time(1)),-1*ones(size(pmt_data.time)),dot_size,'filled','b')
    ylim([-2,2])
    title('Real time SMR vs PMT stamps')
    legend(['SMR signal';'PMT signal'],'Location','eastout')
    xlabel('Elasped time (min)')
    set(gca, 'YTick', [], 'FontSize',20)
    

print(gcf,'-dpng','-r400','temp_high_res_image_1')
figReporter1 =Image(which('temp_high_res_image_1.png'));
figReporter1.Style = {ScaleToFit};
add(ch2,figReporter1);

%% Figure 2 Plotting collapesed PMT occurances on normalized smr detection-time cordinate
rpt_title = clone(sectrpt_title());
append(rpt_title,'Plotting collapesed PMT occurances on normalized smr detection-time cordinate');
append(ch2,Section('title',rpt_title));
p = Paragraph(['This plot is used to help user decide on the window for unique pairing between a smr event '...
        ,'and a pmt event. For each smr event, all pmt events that occured within +-'...
        ,'[see below]ms of the smr detection are recorded. The differences in detection time '...
        ,'between the identified pmt event and the smr event of interest are hence '...
        ,'calculated. After this process has been iterated for every smr event, a '...
        ,'1D heatmap of all the pmt-to-smr transit time is plotted with the color '...
        ,'indicating the number of pmt event occurances at a particular delta_t to '...
        ,'the smr signal.'...
        ,'Cutoffs for the window of unique pairing that the user has selected to generate this '...
        ,'analysis are indicated under the Key parameters section']);
p.Style = caption_style; 
append(ch2,p);

figure('OuterPosition',[0.2*scrsize(3) 0.2*scrsize(4) 0.6*scrsize(3) 0.5*scrsize(4)]);
    subplot(2,1,1);
        [counts,edges] = histcounts(trace,2*window);
        x = edges(1:end-1);
        y = zeros(1,length(x));
        col = counts;  % This is the color, vary with x in this case.
        surface([x;x],[y;y],[col;col],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',100);
        xlim([-window,window])
        xlabel('Delta T from PMT to SMR signals (ms)')
        set(gca, 'YTick', [],'FontSize',13)
        c=colorbar;
        set(get(c,'title'),'string','PMT occurrence','Rotation',0);

    subplot(2,1,2);
        [counts,edges] = histcounts(trace,2*window);
        x = edges(1:end-1);
        y = zeros(1,length(x));
        col = counts;  % This is the color, vary with x in this case.
        surface([x;x],[y;y],[col;col],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',100);
        xlim([-80,80])
        xlabel('Delta T from PMT to SMR signals (ms)')
        title('Zoomed in')
        set(gca, 'YTick', [],'FontSize',13)
        c=colorbar;
        set(get(c,'title'),'string','PMT occurrence','Rotation',0);

print(gcf,'-dpng','-r400','temp_high_res_image_2')
figReporter2 =Image(which('temp_high_res_image_2.png'));
figReporter2.Style = {ScaleToFit};
add(ch2,figReporter2);


%% Figure 3 Plotting all paire-wise comparisions between channels from paired data
rpt_title = clone(sectrpt_title());
append(rpt_title,'Plotting all paire-wise comparisions between channels from paired data');
append(ch2,Section('title',rpt_title));

figure('OuterPosition',[0*scrsize(3) 0.2*scrsize(4) 0.92*scrsize(3) 0.7*scrsize(4)]);

x_axis_pmt_channel = [1 1 1 1 2 2 2 3 3 4];
y_axis_pmt_channel = [2 3 4 5 3 4 5 4 5 5];
rpt_title_color_lab = ["Pacific Blue","FITC","PE","APC","Cy7"];

ha = tight_subplot(2,n_pmt_channel,[.15 .05],[.1 .05]);

for i = 1:length(x_axis_pmt_channel)
    axes(ha(i));
        x_raw = readout_paired(:,x_axis_pmt_channel(i)+2); % +2 to account for time and smr data column
        y_raw = readout_paired(:,y_axis_pmt_channel(i)+2); 
        x_filter_ind = find(x_raw<voltage_plot_lim_higher & x_raw>voltage_plot_lim_lower);
        y_filter_ind = find(y_raw<voltage_plot_lim_higher & y_raw>voltage_plot_lim_lower);
        [filter_ind,~] = intersect(x_filter_ind,y_filter_ind);
        x_filtered = x_raw(filter_ind);
        y_filtered = y_raw(filter_ind);
        scatter(x_filtered,y_filtered, 3,'filled')
        symlog()
        title(append(rpt_title_color_lab(x_axis_pmt_channel(i)), ' vs ' ,rpt_title_color_lab(y_axis_pmt_channel(i))))
        xlabel(append(rpt_title_color_lab(x_axis_pmt_channel(i)), ' (mV)'))
        ylabel(append(rpt_title_color_lab(y_axis_pmt_channel(i)), ' (mV)'))
        set(gca,'FontSize',12)
        legend(['n' '=' int2str(length(y_filtered))],'location',"northeast",'FontSize',9)
end
print(gcf,'-dpng','-r400','temp_high_res_image_3')
figReporter3 =Image(which('temp_high_res_image_3.png'));
figReporter3.Style = {ScaleToFit};
add(ch2,figReporter3);


%% Figure 4 Plotting fluorescence vs buoyant mass for all channels
rpt_title = clone(sectrpt_title());
append(rpt_title,'Plotting fluorescence vs buoyant mass for all channels');
append(ch2,Section('title',rpt_title));

figure('OuterPosition',[0*scrsize(3) 0.2*scrsize(4) 0.92*scrsize(3) 0.7*scrsize(4)]);

color_rbg = [0.3010 0.7450 0.9330;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.8500 0.3250 0.0980;
    0.6350 0.0780 0.1840];

ha = tight_subplot(2,n_pmt_channel,[.15 .05],[.1 .05]);

for i = 1:n_pmt_channel
    axes(ha(i));
%     subplot(2,n_pmt_channel,i);
        x_raw = readout_paired(:,2); 
        y_raw = readout_paired(:,i+2); % +2 to account for time and smr data column
        x_filter_ind = find(x_raw<smr_plot_lim_higher & x_raw>smr_plot_lim_lower);
        y_filter_ind = find(y_raw<voltage_plot_lim_higher & y_raw>voltage_plot_lim_lower);
        [filter_ind,~] = intersect(x_filter_ind,y_filter_ind);
        x_filtered = x_raw(filter_ind);
        y_filtered = y_raw(filter_ind);
        scatter(x_filtered,y_filtered, 3,color_rbg(i,:),'filled')
        symlog('y')
        ylabel(append(rpt_title_color_lab(i), ' (mV)'))
        xlabel("Buoyant mass (pg)")
        title(append(sample_name, ' ', rpt_title_color_lab(i), ' vs BM'))
        legend(['n' '=' int2str(length(y_filtered))],'location',"northwest",'FontSize',9)
        set(gca,'FontSize',12)
end
for i = 1:n_pmt_channel
%     subplot(2,n_pmt_channel,i+n_pmt_channel);
    axes(ha(i+n_pmt_channel));
        x_raw = readout_paired(:,2); 
        y_raw = readout_paired(:,i+2); % +2 to account for time and smr data column
        x_filter_ind = find(x_raw<smr_plot_lim_higher & x_raw>smr_plot_lim_lower);
        y_filter_ind = find(y_raw<voltage_plot_lim_higher & y_raw>voltage_plot_lim_lower);
        [filter_ind,~] = intersect(x_filter_ind,y_filter_ind);
        x_filtered = x_raw(filter_ind);
        y_filtered = y_raw(filter_ind);
        scatter(x_filtered,y_filtered, 3,color_rbg(i,:),'filled')
        symlog('xy')
        ylabel(append(rpt_title_color_lab(i), ' (mV)'))
        xlabel("Buoyant mass (pg)")
        title(append(sample_name, ' ', rpt_title_color_lab(i), ' vs BM'))
        set(gca,'FontSize',12)
        legend(['n' '=' int2str(length(y_filtered))],'location',"northwest",'FontSize',9)
        
end

print(gcf,'-dpng','-r400','temp_high_res_image_4')
figReporter4 =Image(which('temp_high_res_image_4.png'));
figReporter4.Style = {ScaleToFit};
add(ch2,figReporter4);

%% Figure 6 Plotting distribution of smr and pmt signals pre and post-pairing
rpt_title = clone(sectrpt_title());
append(rpt_title,'Plotting distribution of smr and pmt signals pre and post-pairing');
append(ch2,Section('title',rpt_title));
p = Paragraph(['This plot is intended for a rough check on whether certain signals are ',...
'biased towards being lost during the pairing process. Extremly low ',...
'buoyant mass events are prone to not be paired likely because ',...
'they might be debris that are not labeled with any fluoresence markers. ',...
'Low intensity fluorescence signals are also prone to be drop off during ',...
'the pairing process, likely because they might be higher noise signals ',...
'that passed the pmt detection threshold and got recorded as pmt events.']);
p.Style = caption_style; 
append(ch2,p);

n_bin_bm = zeros(1,2);
bm_size = [length(smr_data.smr), length(readout_paired(:,1))];
for i = 1:length(n_bin_bm)
    if bm_size(i) < 1000
       n_bin_bm(i) = 20;
    elseif bm_size(i) < 3000
       n_bin_bm(i) = 40;
    elseif bm_size(i) < 6000
       n_bin_bm(i) = 50;
    elseif bm_size(i) < 9000
       n_bin_bm(i) = 60;
    else
       n_bin_bm(i) = 80;
    end
end

n_bin_pmt = zeros(1,2);
pmt_size = [length(pmt_data.time), length(readout_paired(:,1))];
for i = 1:length(n_bin_pmt)
    if pmt_size(i) < 1000
       n_bin_pmt(i) = 10;
    elseif pmt_size(i) < 3000
       n_bin_pmt(i) = 20;
    elseif pmt_size(i) < 6000
       n_bin_pmt(i) = 40;
    elseif pmt_size(i) < 9000
       n_bin_pmt(i) = 45;
    else
       n_bin_pmt(i) = 50;
    end
end

figure('OuterPosition',[0.1*scrsize(3) 0.05*scrsize(4) 0.4*scrsize(3) 0.95*scrsize(4)]);
color_rbg = [0.3010 0.7450 0.9330;
        0.4660 0.6740 0.1880;
        0.9290 0.6940 0.1250;
        0.8500 0.3250 0.0980;
        0.6350 0.0780 0.1840];
ha = tight_subplot(n_pmt_channel+1,2,[.05 .08],[.05 .03],[.1 .05]);
    
for i = 1:n_pmt_channel+1
    axes(ha(i*2));
    obj2plot =readout_paired(:,i+1); %+1 to account for time column
    if i ==1
        obj2plot = obj2plot(obj2plot<smr_plot_lim_higher & obj2plot>smr_plot_lim_lower);
        hh=histogram(obj2plot,n_bin_bm(2));
        xlabel("Buoyant mass (pg)")
        title("Paired SMR signal distribution")
        xlim([smr_plot_lim_lower,smr_plot_lim_higher])
    else
        obj2plot = obj2plot(obj2plot<voltage_plot_lim_higher & obj2plot>voltage_plot_lim_lower);
        [~,edges] = histcounts(log10(obj2plot),n_bin_pmt(2));
        hh = histogram(obj2plot,10.^edges);
        xlim([voltage_plot_lim_lower,voltage_plot_lim_higher])
        set(gca, 'xscale','log')
        xlabel(append(rpt_title_color_lab(i-1), ' (mV)'))
        title(append(' Paired ', rpt_title_color_lab(i-1), ' signal distribution'))
        hh.FaceColor = color_rbg(i-1,:);
    end
        
        ylabel("Counts")
        legend(['n' '=' int2str(length(obj2plot))],'location',"northeast")
        hh.EdgeColor = 'w';
end

for i = 1:n_pmt_channel+1
    axes(ha(i*2-1));
    if i ==1
        obj2plot =smr_data.smr;
        obj2plot = obj2plot(obj2plot<smr_plot_lim_higher & obj2plot>smr_plot_lim_lower);
        hh=histogram(obj2plot,n_bin_bm(1));
        xlabel("Buoyant mass (pg)")
        title("Raw SMR signal distribution")
        xlim([smr_plot_lim_lower,smr_plot_lim_higher])
    else
        obj2plot =pmt_data.pmt{i-1}; %-1 to account for smr plotting
        obj2plot = obj2plot(obj2plot<voltage_plot_lim_higher & obj2plot>voltage_plot_lim_lower);
        [~,edges] = histcounts(log10(obj2plot),n_bin_pmt(1));
        hh = histogram(obj2plot,10.^edges);
        set(gca, 'xscale','log')
        xlabel(append(rpt_title_color_lab(i-1), ' (mV)'))
        title(append(' Raw ', rpt_title_color_lab(i-1), ' signal distribution'))
        hh.FaceColor = color_rbg(i-1,:);
        xlim([voltage_plot_lim_lower,voltage_plot_lim_higher])
    end
        
        ylabel("Counts")
        legend(['n' '=' int2str(length(obj2plot))],'location',"northeast")
        hh.EdgeColor = 'w';
end

print(gcf,'-dpng','-r400','temp_high_res_image_6')
figReporter6 =Image(which('temp_high_res_image_6.png'));
figReporter6.Style = {ScaleToFit};
add(ch2,figReporter6);

%% Figure 7 Plotting PMT to SMR transit time vs bm and fluorescence
rpt_title = clone(sectrpt_title());
append(rpt_title,'Plotting PMT to SMR transit time vs bm and fluorescence');
append(ch2,Section('title',rpt_title));
% p = Paragraph(['This plot is intended for a rough check on whether certain signals are ',...
% 'biased towards being lost during the pairing process. Extremly low ',...
% 'buoyant mass events are prone to not be paired likely because ',...
% 'they might be debris that are not labeled with any fluoresence markers. ',...
% 'Low intensity fluorescence signals are also prone to be drop off during ',...
% 'the pairing process, likely because they might be higher noise signals ',...
% 'that passed the pmt detection threshold and got recorded as pmt events.']);
% p.Style = caption_style; 
% append(ch2,p);

figure('OuterPosition',[0*scrsize(3) 0.2*scrsize(4) 0.98*scrsize(3) 0.7*scrsize(4)]);

ha = tight_subplot(2,n_pmt_channel,[.15 .04],[.1 .05],[.03 .005]);

for i = 1:n_pmt_channel
    axes(ha(i));
        x_raw = readout_paired(:,2); % +1 to account for time and smr data column
        y_raw = readout_paired(:,8); 
        x_filter_ind = find(x_raw<smr_plot_lim_higher & x_raw>smr_plot_lim_lower);
        x_filtered = x_raw(x_filter_ind);
        y_filtered = y_raw(x_filter_ind);
        color_overlay = readout_paired(:,i+2);
        color_overlay = color_overlay(x_filter_ind);
        scatter(x_filtered,y_filtered,5,color_overlay,'filled')
        ylabel('PMT to SMR transit time (ms)')
        xlabel("Buoyant mass (pg)")
        cmap = tab20(5);
        set(gca,'ColorScale','log','FontSize',10)
        colormap(cmap);
        c=colorbar;
%         caxis([prctile(color_overlay,10) prctile(color_overlay,90)]);
        title(append(sample_name, ' ', 'Transit time', ' vs BM'))
        legend(['n' '=' int2str(length(y_filtered))],'location',"southeast")
        set(get(c,'title'),'string',append(rpt_title_color_lab(i), ' (mV)'),'Rotation',0);
end
for i = 1:n_pmt_channel
    axes(ha(i+n_pmt_channel));
        x_raw = readout_paired(:,2); % +1 to account for time and smr data column
        y_raw = readout_paired(:,8); 
        x_filter_ind = find(x_raw<smr_plot_lim_higher & x_raw>smr_plot_lim_lower);
        x_filtered = x_raw(x_filter_ind);
        y_filtered = y_raw(x_filter_ind);
        color_overlay = readout_paired(:,i+2);
        color_overlay = color_overlay(x_filter_ind);
        scatter(x_filtered,y_filtered,5,color_overlay,'filled')
        ylabel('PMT to SMR transit time (ms)')
        xlabel("Buoyant mass (pg)")
        cmap = tab20(5);
        symlog('x')
        set(gca,'ColorScale','log','FontSize',10)
        colormap(cmap);
        c=colorbar;
%        caxis([prctile(color_overlay,10) prctile(color_overlay,90)]);
        title(append(sample_name, ' ', 'Transit time', ' vs BM'))
        legend(['n' '=' int2str(length(y_filtered))],'location',"southeast")
        set(get(c,'title'),'string',append(rpt_title_color_lab(i), ' (mV)'),'Rotation',0);
end

print(gcf,'-dpng','-r400','temp_high_res_image_7')
figReporter7 =Image(which('temp_high_res_image_7.png'));
figReporter7.Style = {ScaleToFit};
add(ch2,figReporter7);
%%
append(rpt,ch2);
close(rpt)
close all
delete *.png
%rptview(rpt)
cd (currentFolder)
end