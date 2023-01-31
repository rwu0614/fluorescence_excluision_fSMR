clear 
close all
format long;
currentFolder = pwd;
addpath('report_functions\');
addpath('helper_functions\');
addpath('plotting_functions\');

%% Define instruction and output format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  User Input Required %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input UI to grab path to FBM_metadata_assembly_instruction sheet
fprintf('\nGetting FBM metadata assembly instruction...\n')
[input_info.instruction_filename, input_info.instruction_dir, exist_pmt] = uigetfile('../*.*','Select FBM metadata assembly instruction File',' ');

instruction_path = [input_info.instruction_dir,'\',input_info.instruction_filename];

% Specify base FBM and metadata table format
FBM_var_to_search = {'real_time_sec','elapsed_time_min','buoyant_mass_pg','node_deviation_hz','pmt1_mV','pmt2_mV','pmt3_mV','pmt4_mV','pmt5_mV','vol_au','fxm_baseline_mV','pmt2smr_transit_time_ms'};
instruct_var_to_exclude = {'path','filename'}; % metadata sheet will include all variables from the instruction sheet except variables specified here

% Input UI to specify names for FBM and metadata output
assembly_name = input('Input name for this assembly (suffix annotation to be added to FBM and metadata names):\n','s');
FBM_name = ['FBM_base_',assembly_name,'.txt'];
metadata_name =['metadata_base_',assembly_name,'.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Accessing paired data through FBM_metadata_assembly_instruction sheet
opts = detectImportOptions(instruction_path,'ReadVariableNames',true,'VariableNamingRule','preserve','Delimiter',',');
opts = setvartype(opts,'string');
instruction = readtable(instruction_path,opts);
instruct_var = instruction.Properties.VariableNames;
meta_var_to_include = instruct_var(~ismember(instruct_var, instruct_var_to_exclude));

%% reformat instruction to prep for metadata sheet
% fill all missing value with 'N/A'
instruction = fillmissing(instruction,'constant',"N/A");
for i=1:length(meta_var_to_include)
    var = meta_var_to_include{i};
    segment = instruction.(var);
    for j = 1:length(segment)
        % change all space to '_'
        segment(j) = strrep(segment(j),' ','_');
    end
    instruction.(var)= segment;
end

%% FBM and metadata assembly
merge_paired_data=[];
metadata_base =[];

for i = 1:height(instruction)
    sample_instruction = instruction(i,:);
    sample_path = [sample_instruction.path{:},'\',sample_instruction.filename{:}];
    opts = detectImportOptions(sample_path,'ReadVariableNames',true,'VariableNamingRule','preserve','Delimiter','\t');
    paired_data = readtable(sample_path,opts);
    pre_merge_data = paired_data;
    
    % create unique cell id based on date of measurement and sample name 
    cell_id = {};
    for k = 1:height(pre_merge_data)
        cell_id = [cell_id;strcat([sample_instruction.sample_collectiontime{:},'_',...
            sample_instruction.sample_name{:},'_',num2str(k,'%02d')])];
    end
    pre_merge_data.Properties.RowNames = cell_id;    
    FBM_var_to_include = intersect(FBM_var_to_search,pre_merge_data.Properties.VariableNames);
    merge_paired_data = vertcat(merge_paired_data, pre_merge_data(:,FBM_var_to_include));
    
    % assembling base metadata sheet
    pre_meta = cell2table(cell(height(pre_merge_data),length(meta_var_to_include)));
    pre_meta.Properties.VariableNames = meta_var_to_include;
    pre_meta.Properties.RowNames = cell_id;
    
    for j=1:length(meta_var_to_include)
        var = meta_var_to_include{j};
        pre_meta.(var) = repmat(string(sample_instruction.(var){:}),size(cell_id));
        %pre_meta.(var) = string(pre_meta.(var));
    end
    
    metadata_base = vertcat(metadata_base,pre_meta);
end
FBM_base = merge_paired_data;

%% generate output txt files 
% Save files in the same folder as the assembly instruction sheet
[instruction_rootdir,~,~] = fileparts(instruction_path);
cd(instruction_rootdir)
    writetable(FBM_base,FBM_name,'Delimiter','\t','WriteRowNames',true)
    writetable(metadata_base,metadata_name,'Delimiter','\t','WriteRowNames',true)
    disp('Base FBM top rows:')
    head(FBM_base)
    disp('Base metadata top rows:')
    head(metadata_base)
cd(currentFolder)












