clear 
close all
format long;
currentFolder = pwd;
addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');

%% Define sample and output format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  User Input Required %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input UI to grab path to a single readout paired sample txt file
fprintf('\nGetting paired sample...\n')
[input_info.sample_filename, input_info.sample_dir, ~] = uigetfile('../*.*','Select Readout_paired_[sample name].txt',' ');
sample_path = [input_info.sample_dir,'\',input_info.sample_filename];
opts = detectImportOptions(sample_path,'ReadVariableNames',true,'VariableNamingRule','preserve','Delimiter','\t');
sample = readtable(sample_path,opts);
name_split = strsplit(input_info.sample_dir,'\');   
sample_name = name_split{end-1};   
sample_name= strrep(sample_name,'_',' ');

% Input UI to grab path to save FBM and metadata
fprintf('\nGetting saving location...\n')
input_info.save_dir = uigetdir(input_info.sample_dir,'Select folder to save in');

input_info.measurement_date = input('\nMeasurement_date for cell ID:','s');


% Specify base FBM and metadata table format
FBM_var_to_search = {'real_time_sec','elapsed_time_min','buoyant_mass_pg','node_deviation_hz','pmt1_mV','pmt2_mV','pmt3_mV','pmt4_mV','pmt5_mV','vol_au','pmt2smr_transit_time_ms'};
FBM_var_to_include = intersect(FBM_var_to_search,sample.Properties.VariableNames);
Meta_var_to_include = {'path','pmt1_label','pmt2_label','pmt3_label','pmt4_label','pmt5_label','fluid_density','cell_type'}; 

for i = 1:length(Meta_var_to_include)
    var = Meta_var_to_include{i};
    if strcmp(var,'path')
        input_info.(var) = sample_path;
    else
        input_info.(var) = input(strcat('\n',var," :"),'s');
    end
end

% Input UI to specify names for FBM and metadata output
assembly_name = input('Input name for this assembly (suffix annotation to be added to FBM and metadata names):\n','s');
FBM_name = ['FBM_base_',assembly_name,'.txt'];
metadata_name =['metadata_base_',assembly_name,'.txt'];



%% FBM and metadata assembly
  
% create unique cell id based on date of measurement and sample name 
cell_id = cell(height(sample),1);
for k = 1:height(sample)
    cell_id{k} = strcat([input_info.measurement_date,'_',sample_name,'_',num2str(k,'%02d')]);
end
sample.Properties.RowNames = cell_id;

% assembling base metadata sheet
pre_meta = cell2table(cell(height(sample),length(Meta_var_to_include)));
pre_meta.Properties.VariableNames = Meta_var_to_include;
pre_meta.Properties.RowNames = cell_id;

for j=1:length(Meta_var_to_include)
    var = Meta_var_to_include{j};
    pre_meta.(var) = repmat(input_info.(var),size(cell_id));
end

metadata_base = pre_meta;
FBM_base = sample(:,FBM_var_to_include);

%% generate output txt files 
% Save files
cd(input_info.save_dir)
    writetable(FBM_base,FBM_name,'Delimiter','\t','WriteRowNames',true)
    writetable(metadata_base,metadata_name,'Delimiter','\t','WriteRowNames',true)
    disp('Base FBM top rows:')
    head(FBM_base)
    disp('Base metadata top rows:')
    head(metadata_base)
cd(currentFolder)












