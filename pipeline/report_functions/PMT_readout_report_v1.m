function  PMT_readout_report_v1(report_dir,input_info,sample_dir,sample_name, analysis_params, readout_pmt)
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
rpt = Report(append(report_dir,"PMT_readout_report_",time_stamp),'pdf');
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
rpt_log_sec1.Sample = sprintf('This report is genereated on sample "%s" from path "%s" \n',sample_name,sample_dir);
rpt_log_sec1.PMT1_input = sprintf('PMT channel 1 (Pacific Blue) input is from "%s" from path "%s" \n',input_info.pmt_filename{1},input_info.pmt_dir{1});
rpt_log_sec1.PMT2_input = sprintf('PMT channel 2 (FITC) input is from "%s" from path "%s" \n',input_info.pmt_filename{2},input_info.pmt_dir{2});
rpt_log_sec1.PMT3_input = sprintf('PMT channel 3 (PE) input is from "%s" from path "%s" \n',input_info.pmt_filename{3},input_info.pmt_dir{3});
rpt_log_sec1.PMT4_input = sprintf('PMT channel 4 (APC) input is from "%s" from path "%s" \n',input_info.pmt_filename{4},input_info.pmt_dir{4});
rpt_log_sec1.PMT5_input = sprintf('PMT channel 5 (Cy7) input is from "%s" from path "%s" \n',input_info.pmt_filename{5},input_info.pmt_dir{5});
rpt_log_sec1.PMT_time_input = sprintf('PMT time input is from "%s" from path "%s" \n',input_info.time_filename,input_info.time_dir);
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
rpt_log_sec2.Analysis_parameter_1= sprintf('analysis_params.Baseline_rough_cutoff = %1.5f \n', analysis_params.Baseline_rough_cutoff);
rpt_log_sec2.Analysis_parameter_2 = sprintf('analysis_params.med_filt_length = %1.5f \n', analysis_params.med_filt_length);
rpt_log_sec2.Analysis_parameter_3 = sprintf('analysis_params.moving_average_window_size = %1.5f \n', analysis_params.moving_average_window_size);
rpt_log_sec2.Analysis_parameter_4 = sprintf('analysis_params.med_filt_window_size = %1.5f \n', analysis_params.med_filt_window_size);
rpt_log_sec2.Analysis_parameter_5 = sprintf('analysis_params.min_distance_btw_peaks = %1.5f \n', analysis_params.min_distance_btw_peaks);
rpt_log_sec2.Analysis_parameter_6 = sprintf('analysis_params.uni_peak_range_ext = %1.5f \n', analysis_params.uni_peak_range_ext);
rpt_log_sec2.Analysis_parameter_7 = sprintf('analysis_params.uni_peak_baseline_window_size = %1.5f \n', analysis_params.uni_peak_baseline_window_size);
rpt_log_sec2.Analysis_parameter_8 = sprintf('analysis_params.detect_thresh_pmt(1) = %1.5f (multiplier of noise amplitude)\n', analysis_params.detect_thresh_pmt(1));
rpt_log_sec2.Analysis_parameter_9 = sprintf('analysis_params.detect_thresh_pmt(2) = %1.5f (multiplier of noise amplitude)\n', analysis_params.detect_thresh_pmt(2));
rpt_log_sec2.Analysis_parameter_10 = sprintf('analysis_params.detect_thresh_pmt(3) = %1.5f (multiplier of noise amplitude)\n', analysis_params.detect_thresh_pmt(3));
rpt_log_sec2.Analysis_parameter_11 = sprintf('analysis_params.detect_thresh_pmt(4) = %1.5f (multiplier of noise amplitude)\n', analysis_params.detect_thresh_pmt(4));
rpt_log_sec2.Analysis_parameter_12 = sprintf('analysis_params.detect_thresh_pmt(5) = %1.5f (multiplier of noise amplitude)\n', analysis_params.detect_thresh_pmt(5));
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
% Section 4 Key plotting parameters
voltage_plot_lim_lower = 1e-2;
voltage_plot_lim_higher = 1e+4;

rpt_log_sec4 =[];
rpt_log_sec4.Fluorescence_plotting_lower_cutoff = sprintf('Lower fluorescence plotting limit: %1.5f mV\n',voltage_plot_lim_lower);
rpt_log_sec4.Fluorescence_plotting_higher_cutoff = sprintf('Higher fluorescence plotting limit: %1.5f mV\n',voltage_plot_lim_higher);
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

% making summary matrix that includes all pmt channels for systemic
% plotting
n_pmt_channel =5;
peak_amp_summary = readout_pmt.signal';
scrsize = get(0, 'Screensize');

%% Figure 1 plotting distributions of all fluorescence signal in each channel 
rpt_title = clone(sectrpt_title());
append(rpt_title,'plotting distributions of all fluorescence signal in each channel');
append(ch2,Section('title',rpt_title));
figure('OuterPosition',[0.3*scrsize(3) 0.05*scrsize(4) 0.4*scrsize(3) 0.95*scrsize(4)]);

n_histcount_bin = round(length(peak_amp_summary(1,:))/50);
title_color_lab = ["Pacific Blue","FITC","PE","APC","Cy7"];
color_rbg = [0.3010 0.7450 0.9330;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.8500 0.3250 0.0980;
    0.6350 0.0780 0.1840];

for i = 1:n_pmt_channel
    subplot(n_pmt_channel,1,i);
        y_raw= peak_amp_summary(i,:);
        filter_ind = y_raw<voltage_plot_lim_higher & y_raw>voltage_plot_lim_lower;
        y_filtered = y_raw(filter_ind);
        [~,edges] = histcounts(log10(y_filtered),n_histcount_bin);
        hh = histogram(y_filtered,10.^edges);
        set(gca, 'xscale','log')
        xlim([voltage_plot_lim_lower,voltage_plot_lim_higher])
        xlabel(append(title_color_lab(i), ' (mV)'))
        ylabel("Counts")
        title(append(sample_name, ' ', title_color_lab(i), ' PMT voltage distribution'))
        legend(['n' '=' int2str(length(y_filtered))],'location',"northeast")
        hh.FaceColor = color_rbg(i,:);
        hh.EdgeColor = 'w';
end
print(gcf,'-dpng','-r400','temp_high_res_image_1')
figReporter1 =Image(which('temp_high_res_image_1.png'));
figReporter1.Style = {ScaleToFit};
add(ch2,figReporter1);
%% Figure 2 plotting fluorecence level over the time course of the whole measurement
rpt_title = clone(sectrpt_title());
append(rpt_title,'plotting fluorecence level over the time course of the whole measurement');
append(ch2,Section('title',rpt_title));
p = Paragraph(['This plot is intended to check on the consistency of fluorescence ',...
'intensities over the course of the whole measurement time. The expected ',...
'pattern of a successful measurement would be a relatively flat line with no ',...
'slope, meaning that the sample fluorescene level are consistant throught ',...
'out the measurement. If step wise increases or decreases are seen in this ',...
'plot, it could mean that the pmt reference voltage were changed during an ',...
'active measurement, which is a huge red flag. If this occurs, the dataset ',...
'must be truncated using time domain cutoffs to ensure all pmt events used ',...
'for downstream analysis were recorded from the same pmt reference voltage settings ',...
'ie. from the flatten part of this plot.']);
p.Style = caption_style; 
append(ch2,p);

figure('OuterPosition',[0.3*scrsize(3) 0.05*scrsize(4) 0.4*scrsize(3) 0.95*scrsize(4)]);

for i = 1:n_pmt_channel
    subplot(n_pmt_channel,1,i);
        y_raw= peak_amp_summary(i,:);
        filter_ind = y_raw<voltage_plot_lim_higher & y_raw>voltage_plot_lim_lower;
        y_filtered = y_raw(filter_ind);
        scatter((readout_pmt.time_of_detection(filter_ind)-...
            readout_pmt.time_of_detection(1))/60,y_filtered,5,color_rbg(i,:),'filled');
        set(gca, 'yscale','log')
        ylabel(append(title_color_lab(i), ' (mV)'))
        xlabel("Time (min)")
        title(append(sample_name, ' ', title_color_lab(i), ' PMT voltage vs time'))
        legend(['n' '=' int2str(length(y_filtered))],'location',"northeast")
end

print(gcf,'-dpng','-r400','temp_high_res_image_2')
figReporter2 =Image(which('temp_high_res_image_2.png'));
figReporter2.Style = {ScaleToFit};
add(ch2,figReporter2);
%% Figure 3 Plotting all paire-wise comparisions between channels from all detected pmt events (extreme outliers excluded)
rpt_title = clone(sectrpt_title());
append(rpt_title,'Plotting all paire-wise comparisions between channels from all detected pmt events (extreme outliers excluded)');
append(ch2,Section('title',rpt_title));

figure('OuterPosition',[0.02*scrsize(3) 0.1*scrsize(4) 0.97*scrsize(3) 0.7*scrsize(4)]);

x_axis_pmt_channel = [1 1 1 1 2 2 2 3 3 4];
y_axis_pmt_channel = [2 3 4 5 3 4 5 4 5 5];
title_color_lab = ["Pacific Blue","FITC","PE","APC","Cy7"];
for i = 1:length(x_axis_pmt_channel)
    subplot(2,5,i);
        x_raw = peak_amp_summary(x_axis_pmt_channel(i),:);
        y_raw = peak_amp_summary(y_axis_pmt_channel(i),:);
        x_filter_ind = find(x_raw<voltage_plot_lim_higher & x_raw>voltage_plot_lim_lower);
        y_filter_ind = find(y_raw<voltage_plot_lim_higher & y_raw>voltage_plot_lim_lower);
        [filter_ind,~] = intersect(x_filter_ind,y_filter_ind);
        x_filtered = x_raw(filter_ind);
        y_filtered = y_raw(filter_ind);
        scatter(x_filtered,y_filtered, 1.5,'filled')
        symlog()
        title(append(title_color_lab(x_axis_pmt_channel(i)), ' vs ' ,title_color_lab(y_axis_pmt_channel(i))))
        xlabel(append(title_color_lab(x_axis_pmt_channel(i)), ' (mV)'))
        ylabel(append(title_color_lab(y_axis_pmt_channel(i)), ' (mV)'))
        legend(['n' '=' int2str(length(filter_ind))],'location',"southeast")
end
print(gcf,'-dpng','-r400','temp_high_res_image_3')
figReporter3 =Image(which('temp_high_res_image_3.png'));
figReporter3.Style = {ScaleToFit};
add(ch2,figReporter3);
%%
append(rpt,ch2);
close(rpt)
close all
delete *.png
%rptview(rpt)
cd (currentFolder)
end