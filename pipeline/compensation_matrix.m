%% Desciption
% compensation_matrix.m creates a 5x5 compensation matrix to correct for
% fluorescence spillover from a primary PMT channel to all other channels.
% The function takes in user-defined readout_pmt csv files from single
% color control samples and ask for manual input to set the optimal
% compensation factors. For channels that did not have a single color
% control, compensation factors from such channels to all other channels
% are set to 0, which is equivalent to setting no compensation. The
% function will save a report, as well as a 5x5 matrix as a csv file in the mother directory of
% the last user-input single color control pmt_readout file (usually the expeirment master folder)
addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');


%% Initialization
clear all;
clc;
close all;
currentFolder = pwd;
scrsize = get(0, 'Screensize');
n_pmt_channel = 5;
% Setting defaul 5x5 conpensation matrix with all zeros terms besides ones on the
% diagonal. Default matrix is equivalent to no compensation 
compen_mat = diag(ones(1,5));

% Setting up the order of channels from low wavelength to high wavelength
channels = ["Pacific Blue","FITC","PE","APC","Cy7"];

% Setting up user input info structure for report purposes
input_info.experiment_dir =[];
for i = 1:n_pmt_channel
    input_info.channel{i} = channels(i);
    input_info.single_color_ctrl_filename{i} = [];
    input_info.single_color_ctrl_filedir{i} = [];
    input_info.single_color_compen_array{i} = [];
    input_info.single_color_compen_array_indicator = [];
end

% Setting up lower and upper limit for plotting pmt values to exclude
% extreme outliers
voltage_plot_lim_lower = 1e-2;
voltage_plot_lim_higher = 1e+4;
dot_trans = 0.5;
%% Set compensation factors
output_flag = 0;

while output_flag ==0
    fprintf('\nGetting single color control pmt data...\n')
    fprintf('\nWhich pmt channel is the single color control?')
    input_channel = input('\n[1/2/3/4/5] 1-Pacific blue, 2-FITC, 3-PE, 4-APC, 5-Cy7\n');
    %Read input file
    [input_info.single_color_ctrl_filename{input_channel}, ...
        input_info.single_color_ctrl_filedir{input_channel}, ~] = uigetfile('../*.*','Select single color control File',' ');
    cd(input_info.single_color_ctrl_filedir{input_channel})
    sc_input = readtable(input_info.single_color_ctrl_filename{input_channel});
    % grab nx6 matrix from single color input where n=number of detected
    % pmt event and each column(2-6) is a pmt channel
    sc_data= sc_input{:,2:6};
    cd(currentFolder)
    %Get sample name
    filename_split = strsplit(input_info.single_color_ctrl_filename{input_channel},'.');
    filename_split = filename_split{1:end-1};
    filename_split = strsplit(filename_split,'_');
    sample_name = filename_split{end};
    
    %Plot 1x5 input pmt vs other pmts comparisions 
    close all;
    figure('OuterPosition',[0.05*scrsize(3) 0.5*scrsize(4) 0.9*scrsize(3) 0.4*scrsize(4)]);
    compen_array = zeros(5,1);
    mod_compen_flag = 0;
    while mod_compen_flag == 0
        % Plot before and after compensation
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
                    'MarkerFaceAlpha',dot_trans,'MarkerEdgeAlpha',dot_trans); hold on;
                scatter(x_filtered,comp_y_filtered, 5,'r','filled'...
                    ,'MarkerFaceAlpha',dot_trans,'MarkerEdgeAlpha',dot_trans); hold off;
                symlog('xy')
                ylabel(append(channels(i), ' (mV)'))
                xlabel(append(channels(input_channel), ' (mV)'))
                title(sample_name)
                legend(["Raw","Compensated"],'location',"northwest")
        end
        
        stay = input('Modify current compensation array? [1/0], input 1-Yes, 0-No');
        if stay == 1
            for i = 1:n_pmt_channel
                disp_pair(i) = append(channels(input_channel),' to ',channels(i)); 
            end
            fprintf('\nInput 1D compensation array: [%s,%s,%s,%s,%s]\n',disp_pair);
            compen_array = input('[_,_,_,_,_]\n');
            input_info.single_color_compen_array{input_channel} = compen_array;
            input_info.single_color_compen_array_indicator{input_channel} = disp_pair;
        elseif stay == 0
            mod_compen_flag = 1;
            compen_mat(:,input_channel) = compen_array';
            compen_mat(input_channel,input_channel)=1;
        end
    end
    
    % Exit loop or keep optimizing the matrix
    output_flag = input('\nReady to output final compensation matrix? [1/0] input 1 if yes, 0 if no\n');
end


%% Output compensation matrix and report
%  5x5 matrix format follows: 
%     column :  [PacificBlue,FITC, PE, APC, Cy7]
%     row :     [PacificBlue,FITC, PE, APC, Cy7]
col_names = {'PacificBlue','FITC','PE','APC','Cy7'};
row_names = col_names;
output_mat = array2table(compen_mat,'VariableNames',col_names,'RowName',row_names);
parts = strsplit(input_info.single_color_ctrl_filedir{input_channel}, '\');
parent_dir = string(join(parts(1:end-2),"\"));
input_info.experiment_dir = parent_dir;

cd(parent_dir)
out_file_name = ['compensation_matrix','.csv'];
writetable(output_mat,out_file_name,'WriteRowNames',true);
cd(currentFolder)

report_folder_name = '\compensation_report';
publish('compensation_matrix_report.m','outputDir',strcat(parent_dir,report_folder_name),'codeToEvaluate',...
    'compensation_matrix_report(input_info,output_mat);','showCode',false);

