function  apply_compensation_report(sample_name,input_info,compen_mat_input,fs_data_uncomp,fs_data_comp)
% PMT_report.m creates a html report of a newly completed analysis from PMT_readout.m
% where fluorescent events were recognized from raw PMT data. The goal for the report is to
% bookmark key parameters from the original analysis for future replication
% and reflection of the analysis.
close all;
addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');

%% Report log
fprintf('This report is genereated on %s\n', datetime(now,'ConvertFrom','datenum'))
fprintf('The compensation matrix is from "%s" from path "%s" \n',...
    input_info.compen_mat_filename,input_info.compen_mat_filedir)
fprintf('The full stain uncompensated input is from "%s" from path "%s" \n',...
    input_info.fullstain_filename,input_info.fullstain_filedir)
%% Key analysis parameters
fprintf('\nCompensation matrix applied:\n')
disp(compen_mat_input)


%% Key plotting parameters
voltage_plot_lim_lower = -1e+2;
voltage_plot_lim_higher = 1e+4;

fprintf('Fluorescence intensity limits for plotting pmt voltage readouts to exclude extreme outliers:\n')
fprintf('Low intensity limit: %1.5f mV\n',voltage_plot_lim_lower)
fprintf('High intensity limit: %1.5f mV\n',voltage_plot_lim_higher)

% making summary matrix that includes all pmt channels for systemic
% plotting
n_pmt_channel =5;
peak_amp_summary = fs_data_comp';
title_color_lab = ["Pacific Blue","FITC","PE","APC","Cy7"];
scrsize = get(0, 'Screensize');

%% Figure 1 all pair-wise comparision before and after compensation
% Blue - before compensation, Red-after compensation
dot_trans = 0.06;
plot_ct =1;
figure('OuterPosition',[0.2*scrsize(3) 0.05*scrsize(4) 0.55*scrsize(3) 0.95*scrsize(4)]);
legend(["Raw","Compensated"],'location',"northwest")
for input_channel = 1:n_pmt_channel   
        for i = 1:n_pmt_channel
                subplot(n_pmt_channel,n_pmt_channel,plot_ct);
                    % Excluding outliers for uncompensated data
                    x_raw = fs_data_uncomp(:,input_channel); 
                    y_raw = fs_data_uncomp(:,i); % +2 to account for time and smr data column
                    x_filter_ind = find(x_raw<voltage_plot_lim_higher & x_raw>voltage_plot_lim_lower);
                    y_filter_ind = find(y_raw<voltage_plot_lim_higher & y_raw>voltage_plot_lim_lower);
                    [filter_ind,~] = intersect(x_filter_ind,y_filter_ind);
                    x_filtered = x_raw(filter_ind);
                    y_filtered = y_raw(filter_ind);
                    % Plot
                    scatter(x_filtered,y_filtered, 3,'b','filled',...
                        'MarkerFaceAlpha',dot_trans ,'MarkerEdgeAlpha',dot_trans ); hold on;
                    
                    % Excluding outliers for uncompensated data
                    x_raw = fs_data_comp(:,input_channel); 
                    y_raw = fs_data_comp(:,i); % +2 to account for time and smr data column
                    x_filter_ind = find(x_raw<voltage_plot_lim_higher & x_raw>voltage_plot_lim_lower);
                    y_filter_ind = find(y_raw<voltage_plot_lim_higher & y_raw>voltage_plot_lim_lower);
                    [filter_ind,~] = intersect(x_filter_ind,y_filter_ind);
                    x_filtered = x_raw(filter_ind);
                    y_filtered = y_raw(filter_ind);
                    scatter(x_filtered,y_filtered, 3,'r','filled'...
                        ,'MarkerFaceAlpha',dot_trans ,'MarkerEdgeAlpha',dot_trans ); hold off;
                    symlog('xy',1.3)
                    ylabel(append(title_color_lab(i), ' (mV)'))
                    xlabel(append(title_color_lab(input_channel), ' (mV)'))
                    title(sample_name)
                    plot_ct = plot_ct+1;
        end
end

%% Figure 2 plotting distributions of all compensated fluorescence signal in each channel 

figure('OuterPosition',[0.3*scrsize(3) 0.05*scrsize(4) 0.4*scrsize(3) 0.95*scrsize(4)]);
% n_histcount_bin = round(length(peak_amp_summary(1,:))/50);
n_histcount_bin = 100;
color_rbg = [0.3010 0.7450 0.9330;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.8500 0.3250 0.0980;
    0.6350 0.0780 0.1840];

for i = 1:n_pmt_channel
    subplot(n_pmt_channel,1,i);
        y_raw= peak_amp_summary(i,:);
        filter_ind = y_raw<voltage_plot_lim_higher & y_raw>0.01;
        y_filtered = y_raw(filter_ind);
        [~,edges] = histcounts(log10(y_filtered),n_histcount_bin);
        hh = histogram(y_filtered,10.^edges);
        set(gca, 'xscale','log')
        xlim([0.01,voltage_plot_lim_higher])
        xlabel(append(title_color_lab(i), ' (mV)'))
        ylabel("Counts")
        title(append(sample_name, ' ', title_color_lab(i), ' Compensated PMT distribution'))
        legend(['n' '=' int2str(length(y_filtered))],'location',"northeast")
        hh.FaceColor = color_rbg(i,:);
        hh.EdgeColor = 'w';
end

%% Figure 3 Plotting all paire-wise comparisions between channels from all compensated pmt events (extreme outliers excluded)
figure('OuterPosition',[0.02*scrsize(3) 0.1*scrsize(4) 0.97*scrsize(3) 0.7*scrsize(4)]);
dot_trans = 0.1;
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
        scatter(x_filtered,y_filtered, 5,'filled'...
            ,'MarkerFaceAlpha',dot_trans ,'MarkerEdgeAlpha',dot_trans)
        symlog('xy',1.3)
        title(append(title_color_lab(x_axis_pmt_channel(i)), ' vs ' ,title_color_lab(y_axis_pmt_channel(i))))
        xlabel(append(title_color_lab(x_axis_pmt_channel(i)), ' (mV)'))
        ylabel(append(title_color_lab(y_axis_pmt_channel(i)), ' (mV)'))
        legend(['n' '=' int2str(length(filter_ind))],'location',"southeast")
end


end