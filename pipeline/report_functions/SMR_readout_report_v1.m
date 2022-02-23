function  SMR_readout_report_v1(report_dir,input_info,sample_name, datasmr_good, output_smr, analysis_params)
% PMT_report.m creates a html report of a newly completed analysis from PMT_readout.m
% where fluorescent events were recognized from raw PMT data. The goal for the report is to
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
rpt = Report(append(report_dir,"SMR_readout_report_",time_stamp),'pdf');
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
rpt_log_sec1.Sample = sprintf('This report is genereated on sample "%s" from path "%s" \n',sample_name,input_info.smr_dir);
rpt_log_sec1.SMR_input = sprintf('SMR input is from "%s" from path "%s" \n',input_info.smr_filename,input_info.smr_dir);
rpt_log_sec1.SMR_time_input = sprintf('SMR time input is from "%s" from path "%s" \n',input_info.smr_time_filename,input_info.smr_time_dir);
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
rpt_log_sec2 = [];
rpt_log_sec2.Analysis_parameter_1= sprintf('analysis_params.diff_threshold = %1.5f \n', analysis_params.diff_threshold);
rpt_log_sec2.Analysis_parameter_2 = sprintf('analysis_params.med_filt_wd = %1.5f \n', analysis_params.med_filt_wd);
rpt_log_sec2.Analysis_parameter_3 = sprintf('analysis_params.bs_dev_thres = %1.5f \n', analysis_params.bs_dev_thres);
rpt_log_sec2.Analysis_parameter_4 = sprintf('analysis_params.unqPeakDist = %1.5f \n', analysis_params.unqPeakDist);
rpt_log_sec2.Analysis_parameter_5 = sprintf('analysis_params.offset_input = %1.5f \n', analysis_params.offset_input);
rpt_log_sec2.Analysis_parameter_6 = sprintf('analysis_params.edgethres = %1.5f \n', analysis_params.edgethres);
rpt_log_sec2.Analysis_parameter_7 = sprintf('analysis_params.stdevmultiplier = %1.5f \n', analysis_params.stdevmultiplier);
rpt_log_sec2.Analysis_parameter_8 = sprintf('analysis_params.diffmultiplier = %1.5f \n', analysis_params.diffmultiplier);
rpt_log_sec2.Analysis_parameter_9 = sprintf('analysis_params.winsize = 150 = %1.5f \n', analysis_params.winsize);
rpt_log_sec2.Analysis_parameter_10 = sprintf('analysis_params.estimated_datapoints_optimized = %1.5f \n', analysis_params.estimated_datapoints_optimized);
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

% Section 3 QC parameters
% rpt_log_sec3 =[];
% rpt_log_sec3.QC_thresh_baselineDiff_over_sig = sprintf('analysis_params.thresh_baselineDiff_over_sig = %1.5f \n', analysis_params.thresh_baselineDiff_over_sig);
% rpt_log_sec3.QC_thresh_base_slope = sprintf('analysis_params.thresh_base_slope = %1.5f \n', analysis_params.thresh_base_slope);
% rpt_log_sec3.QC_thresh_base_height_range = sprintf('analysis_params.thresh_base_height_range = %1.5f \n',analysis_params.thresh_base_height_range);
% rpt_log_sec3.QC_pass_rate = QC_msg;
% fn = fieldnames(rpt_log_sec3);
% 
% rpt_title = clone(sectrpt_title());
% append(rpt_title,'QC parameters and pass rate');
% append(ch1,Section('title',rpt_title));
% for i = 1:numel(fn)
%     rpt_title = clone(paramrpt_title());
%     append(rpt_title,fn{i});
%     p = Paragraph(rpt_log_sec3.(fn{i}));
%     p.Style = paragraph_style;
%     append(ch1,Section('title',rpt_title,'Content',p));
%     %append(rpt,p);
% end

% Section 4 Key plotting parameters
% voltage_plot_lim_lower = 1e-2;
% voltage_plot_lim_higher = 1e+4;
% 
% rpt_log_sec4 =[];
% rpt_log_sec4.Fluorescence_plotting_lower_cutoff = sprintf('Lower fluorescence plotting limit: %1.5f mV\n',voltage_plot_lim_lower);
% rpt_log_sec4.Fluorescence_plotting_higher_cutoff = sprintf('Higher fluorescence plotting limit: %1.5f mV\n',voltage_plot_lim_higher);
% fn = fieldnames(rpt_log_sec4);
%  
% rpt_title = clone(sectrpt_title());
% append(rpt_title,'Key plotting parameters');
% append(ch1,Section('title',rpt_title));
% for i = 1:numel(fn)
%     rpt_title = clone(paramrpt_title());
%     append(rpt_title,fn{i});
%     p = Paragraph(rpt_log_sec4.(fn{i}));
%     p.Style = paragraph_style;
%     append(ch1,Section('title',rpt_title,'Content',p));
% end
append(rpt,ch1);

%% Part 2 Report figures
rpt_title = clone(chrpt_title());
append(rpt_title,'Report figures');
ch2 = Chapter('title',rpt_title);
scrsize = get(0, 'Screensize');

global datafull

% making summary matrix that includes all pmt channels for systemic
% plotting
% n_pmt_channel =5;
% peak_amp_summary = readout_pmt.signal';
% scrsize = get(0, 'Screensize');

%% Figure 1 plotting baseline amplitude and time-series mean buoyant mass in Hz
rpt_title = clone(sectrpt_title());
append(rpt_title,'Plotting baseline amplitude and time-series mean buoyant mass in Hz');
append(ch2,Section('title',rpt_title));
figure('OuterPosition',[0 0.5*scrsize(4) 0.95*scrsize(3) 0.3*scrsize(4)]);

[ax, h1, h2]=plotyy(datasmr_good(:,1), datasmr_good(:,3), datasmr_good(:,1), datasmr_good(:,4));
set(h1, 'LineStyle', '-.');
set(h2, 'LineStyle', '-');
xlabel('Time (computer real-time)')
ylabel('Frequency shifts (Hz)')
hold on
yyaxis right
ylabel('Baseline amplitude (Hz)')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',-90,'VerticalAlignment','bottom')
    
print(gcf,'-dpng','-r400','temp_high_res_image_1')
figReporter1 =Image(which('temp_high_res_image_1.png'));
figReporter1.Style = {ScaleToFit};
add(ch2,figReporter1);
%% Figure 2 plotting individual transit time and throughput across measurement duration
rpt_title = clone(sectrpt_title());
append(rpt_title,'Plotting individual transit time and throughput across measurement duration');
append(ch2,Section('title',rpt_title));
p = Paragraph(['This plot is intended to check on the consistency of transit time ',...
'over the course of the whole measurement time. The expected ',...
'If step wise increases or decreases are seen in the transity time ',...
'plot, it could mean issues with pressure that could result from the system itself ',...
'or a clog. The throughput plot gives some indication of how many cells were measured ',...
'throighout the measurement.']);
p.Style = caption_style; 
append(ch2,p);

figure('OuterPosition',[0 0.5*scrsize(4) 0.95*scrsize(3) 0.3*scrsize(4)]);

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

print(gcf,'-dpng','-r400','temp_high_res_image_2')
figReporter2 =Image(which('temp_high_res_image_2.png'));
figReporter2.Style = {ScaleToFit};
add(ch2,figReporter2);
%% Figure 3 Plotting histogram of detected frequency shifts
rpt_title = clone(sectrpt_title());
append(rpt_title,'Plotting histogram of detected frequency shifts');
append(ch2,Section('title',rpt_title));

figure('OuterPosition',[0 0.5*scrsize(4) 0.95*scrsize(3) 0.3*scrsize(4)]);

histogram(mass_transit_fil_datafull(2,:),150)
xlabel('Frequency shift (Hz)')
ylabel('Counts')

print(gcf,'-dpng','-r400','temp_high_res_image_3')
figReporter3 =Image(which('temp_high_res_image_3.png'));
figReporter3.Style = {ScaleToFit};
add(ch2,figReporter3);

%% Figure 4 Plotting detected frequency shifts (peaks) over time
rpt_title = clone(sectrpt_title());
append(rpt_title,'Plotting detected frequency shifts (peaks) over time');
append(ch2,Section('title',rpt_title));

figure('OuterPosition',[0 0.5*scrsize(4) 0.95*scrsize(3) 0.3*scrsize(4)]);

scatter(output_smr(:,1),output_smr(:,2),5,'filled')
xlabel('Time (computer real-time)')
ylabel('Frequency shift (Hz)')

print(gcf,'-dpng','-r400','temp_high_res_image_4')
figReporter4 =Image(which('temp_high_res_image_4.png'));
figReporter4.Style = {ScaleToFit};
add(ch2,figReporter4);

%% Figure 5 Plotting histogram of detected peak width
rpt_title = clone(sectrpt_title());
append(rpt_title,'Plotting histogram of detected peak width');
append(ch2,Section('title',rpt_title));

figure('OuterPosition',[0 0.5*scrsize(4) 0.95*scrsize(3) 0.3*scrsize(4)]);

histogram(datasmr_good(:,20),150)
xlabel('Peak width')
ylabel('Counts')

print(gcf,'-dpng','-r400','temp_high_res_image_5')
figReporter5 =Image(which('temp_high_res_image_5.png'));
figReporter5.Style = {ScaleToFit};
add(ch2,figReporter5);
%%
append(rpt,ch2);
close(rpt)
close all
delete *.png
%rptview(rpt)
cd (currentFolder)
end