%% Desciption
% apply_compensation.m takes in a pre-determined 5x5 compensation matrix
% and apply it to a choosen full-stain uncompensated pmt readout file. User
% input is required for selecting the compentation_matrix.csv file, as well
% as the full stain pmt file ('readout_pmt_uncompensated_[sample
% name].csv). An output compensated pmt file will be generated in the same directory as
% the input pmt file. A report log will also be generated and stored in the
% report folder in the fullstain sample directory. 
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

%% Get input files
fprintf('\nGetting compenstation matrix...\n')
input_info.compen_mat_filename = [];
input_info.compen_mat_filedir = [];
input_info.fullstain_filename = [];
input_info.fullstain_filedir= [];
[input_info.compen_mat_filename,input_info.compen_mat_filedir, ~] = ...
    uigetfile('../*.*','Select compenstation matrix file',' ');
cd(input_info.compen_mat_filedir)
compen_mat_input = readtable(input_info.compen_mat_filename,'ReadRowNames',true);
compen_mat = compen_mat_input{:,:};

fprintf('\nGetting full stain sample...\n')
[input_info.fullstain_filename,input_info.fullstain_filedir, ~] = ...
    uigetfile('../*.*','Select full stain pmt file',' ');
cd(input_info.fullstain_filedir)
fullstain_input = readtable(input_info.fullstain_filename);
% grab nx5 matrix from a full stain pmt input where n=number of detected
% pmt event and column 2-6 are the five pmt channels
fs_data_uncomp = fullstain_input{:,2:6};
cd(currentFolder)
%Get sample name
filename_split = strsplit(input_info.fullstain_filename,'.');
filename_split = filename_split{1:end-1};
filename_split = strsplit(filename_split,'_');
sample_name = filename_split{end};

%% Apply compensation to full stain sample

fs_data_comp = compen_mat*fs_data_uncomp';
fs_data_comp = fs_data_comp';



%% Output compensated full stain pmt data and report
% output pmt file will be the same format as the input file
output_pmt = fullstain_input{:,:};
output_pmt(:,2:6) = fs_data_comp;

cd(input_info.fullstain_filedir)
out_file_name = ['readout_pmt_compensated_' sample_name '.csv'];
dlmwrite(out_file_name, output_pmt, 'delimiter', ',', 'precision', 25);
cd(currentFolder)

% generate analysis report
report_folder_name = ['\' sample_name '_report'];
publish('apply_compensation_report.m','outputDir',strcat(input_info.fullstain_filedir,report_folder_name),'codeToEvaluate',...
    'apply_compensation_report(sample_name,input_info,compen_mat_input,fs_data_uncomp,fs_data_comp);','showCode',false);

