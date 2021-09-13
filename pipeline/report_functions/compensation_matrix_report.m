function compensation_matrix_report(input_info,output_mat)

% compensation_matrix_report.m creates a report on the newly generated
% compensation matrix from single color control samples. The goal for the report is to
% bookmark key parameters from the original analysis for future replication
% and reflection of the analysis.
close all;
currentFolder = pwd;
addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');

%% Report log
fprintf('This report is genereated on %s\n', datetime(now,'ConvertFrom','datenum'))
fprintf('This report is genereated on experiment from path "%s" \n',input_info.experiment_dir)
for i = 1:length(input_info.channel)
    fprintf('%s single color control pmt data is from "%s" from path "%s" \n',...
        input_info.channel{i},input_info.single_color_ctrl_filename{i},input_info.single_color_ctrl_filedir{i})
end

%% Key analysis parameters
channel_w_comp = find(~cellfun(@isempty,input_info.single_color_compen_array));

fprintf('\nAll user-defined compensation factors as following:\n')
for i = 1:length(channel_w_comp)
    factor_name = input_info.single_color_compen_array_indicator{channel_w_comp(i)};
    factor_array = input_info.single_color_compen_array{channel_w_comp(i)};
    for j = 1:length(input_info.channel)
        fprintf('%s:%1.5f\t',factor_name(j),factor_array(j))
    end
    fprintf('\n')
end

fprintf('\nFinal output compensation matrix:\n')
disp(output_mat)

%% Key plotting parameters
voltage_plot_lim_lower = 1e-2;
voltage_plot_lim_higher = 1e+4;

fprintf('Fluorescence intensity limits for plotting pmt voltage readouts to exclude extreme outliers:\n')
fprintf('Low intensity limit: %1.5f mV\n',voltage_plot_lim_lower)
fprintf('High intensity limit: %1.5f mV\n',voltage_plot_lim_higher)

scrsize = get(0, 'Screensize');
dot_trans = 0.4;
n_pmt_channel=5;
channels = ["Pacific Blue","FITC","PE","APC","Cy7"];
%% Compensation result(s) from single color control(s)
for m = 1:length(channel_w_comp)
    input_channel = channel_w_comp(m);
    cd(input_info.single_color_ctrl_filedir{input_channel})
    sc_input = readtable(input_info.single_color_ctrl_filename{input_channel});
    % grab nx6 matrix from single color input where n=number of detected
    % pmt event and each column is a pmt channel
    sc_data= sc_input{:,2:6};
    cd(currentFolder)
    %Get sample name
    filename_split = strsplit(input_info.single_color_ctrl_filename{input_channel},'.');
    filename_split = filename_split{1:end-1};
    filename_split = strsplit(filename_split,'_');
    sample_name = filename_split{end};
    figure('OuterPosition',[0*scrsize(3) 0.5*scrsize(4) 0.99*scrsize(3) 0.37*scrsize(4)]);
    compen_array = input_info.single_color_compen_array{input_channel};
        for i = 1:n_pmt_channel
                subplot(1,n_pmt_channel,i);
                    % Excluding outliers for uncompensated data
                    x_raw = sc_data(:,input_channel); 
                    un_comp_y_raw = sc_data(:,i); % +2 to account for time and smr data column
                    x_filter_ind = find(x_raw<voltage_plot_lim_higher & x_raw>voltage_plot_lim_lower);
                    un_comp_y_filter_ind = find(un_comp_y_raw<voltage_plot_lim_higher & un_comp_y_raw>voltage_plot_lim_lower);
                    [filter_ind,~] = intersect(x_filter_ind,un_comp_y_filter_ind);
                    x_filtered = x_raw(filter_ind);
                    un_comp_y_filtered = un_comp_y_raw(filter_ind);
                    % Apply compensation factor
                    comp_y_filtered = un_comp_y_filtered+x_filtered*compen_array(i);
                    % Plot
                    scatter(x_filtered,un_comp_y_filtered, 5,'b','filled',...
                        'MarkerFaceAlpha',dot_trans ,'MarkerEdgeAlpha',dot_trans ); hold on;
                    scatter(x_filtered,comp_y_filtered, 5,'r','filled'...
                        ,'MarkerFaceAlpha',dot_trans ,'MarkerEdgeAlpha',dot_trans ); hold off;
                    symlog('xy')
                    ylabel(append(channels(i), ' (mV)'))
                    xlabel(append(channels(input_channel), ' (mV)'))
                    title(sample_name)
                    legend(["Raw","Compensated"],'location',"northwest")
        end
end
end