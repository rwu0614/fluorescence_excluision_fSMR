clear all
close all 
clc
currentFolder = pwd;
addpath('plotting_functions\');
addpath('analysis_functions\');

%% Input base FBM and metadata txt files
% Ask user for base FBM txt file
fprintf('\nGetting base FBM matrix...')
[input_info.fbm_filename, input_info.fbm_dir,~] = uigetfile('../*.*','Select base FBM .txt File', ' ');
fprintf('\n%s selected for analysis\n', input_info.fbm_filename)
fbm_file_fullname = strcat(input_info.fbm_dir, input_info.fbm_filename); 
fbm_base = readtable(fbm_file_fullname,'ReadRowNames',true,'ReadVariableNames',true,'Delimiter',' ');
head(fbm_base)

% Ask user for base metadata txt file
fprintf('\nGetting base metadata sheet...')
[input_info.metadata_filename, input_info.metadata_dir,~] = uigetfile('../*.*','Select base FBM .txt File', ' ');
fprintf('\n%s selected for analysis\n', input_info.metadata_filename)
metadata_file_fullname = strcat(input_info.metadata_dir, input_info.metadata_filename); 
metadata_base = readtable(metadata_file_fullname,'ReadRowNames',true,'ReadVariableNames',true,'Delimiter',' ');
head(metadata_base)

%Convert all metadata element to string
var = metadata_base.Properties.VariableNames;
for i = var
   metadata_base.(i{:}) = string(metadata_base.(i{:}));
end
%% Initialize report table where all gating/annotation decisions are stored for each sample
% grab unique rows from the base metadata sheet
report_table = unique(metadata_base);
report_table.Row = [];


%% Gating and annotation - single cell -> live -> G1/SG2
%initializing metadata sheet with all the gating annotation
metadata_gate = metadata_base;
%initializing fbm for gating
fbm_gate = fbm_base;
%fbm_gate.pmt3overbm = fbm_gate.pmt3 ./ fbm_gate.buoyant_mass; %for viability gating

%%
% first gate - single cell gating using buoyant mass
% pg cutoffs will be consistant for any samples from a given cell type
unique_cell_type = unique(report_table.cell_type);
for i = 1:length(unique_cell_type)
    ind_to_gate = metadata_base.cell_type==unique_cell_type{i};
    cell_id_celltype_specific = metadata_base.Row(ind_to_gate);
    var_to_gate = "volume_fL"; % must be string (" "), not char (' ')
    annotation_label = "single_cell"; % must be string (" "), not char (' ')
    annotation = "single_cell"; % must be string (" "), not char (' ')
    rapid_gate = 0;
    cutoff_input=0;
    display_name = ['All ',unique_cell_type{i}];
    [metadata_gate,cutoff_out] = FBM_gating(display_name,fbm_gate,metadata_gate,cell_id_celltype_specific,...
         var_to_gate,annotation_label,annotation, cutoff_input,rapid_gate);
    rpr_ind_uniquecell = report_table.cell_type==unique_cell_type{i};
    report_table.single_cell_bm_gate_low(rpr_ind_uniquecell) = cutoff_out(1);
    report_table.single_cell_bm_gate_high(rpr_ind_uniquecell) = cutoff_out(2);
end

%%

var=["sample_name"];
annotation = ["dt40_rep3_cc"];
[~,ind_fbm_live,~] = meta_grab_cell_id(fbm_gate,metadata_gate,var,annotation);

figure()
obj2plot = fbm_gate.volume_fL(ind_fbm_live);
h1=histogram(obj2plot);
% %     h1.Normalization = 'probability';
%     h1.BinWidth = 5;
%     h1.FaceColor = [0.3010 0.7450 0.9330];

%%
% histograms that display gating of single cells for each sample and saves them

plot_path = uigetdir('../*.*','Where to save the figures?');

for i = 1:length(report_table.sample_name)
    var=["sample_name","single_cell","cell_type"];
    annotation = [report_table.sample_name(i),"single_cell",report_table.cell_type(i)];
    [~,ind_fbm_live,~] = meta_grab_cell_id(fbm_gate,metadata_gate,var,annotation);
    
    var=["sample_name","cell_type"];
    annotation=[report_table.sample_name(i),report_table.cell_type(i)];
    [~,ind_fbm_out,ind_meta_out] = meta_grab_cell_id(fbm_gate,metadata_gate,var,annotation);
    single_cell_missing = ismissing(metadata_gate.single_cell(ind_meta_out));
    ind_fbm_other = nonzeros(single_cell_missing.*ind_fbm_out);
    
    figure()
    obj2plot = fbm_gate.volume_fL(ind_fbm_live);
    h1=histogram(obj2plot);
%     h1.Normalization = 'probability';
%     h1.BinWidth = 100;
    h1.FaceColor = [0.3010 0.7450 0.9330];
    hold on
    obj2plot =fbm_gate.volume_fL(ind_fbm_other);
    h2=histogram(obj2plot);
%     h2.Normalization = 'probability';
%     h2.BinWidth = 100;
    h2.FaceColor = [0.8500 0.3250 0.0980];
    hold off
%     x_max = obj2plot(end);
%     xlim([0 10000]);
    legend('Single cells', 'Non-single cells')
    title(strrep(report_table.sample_name(i),'_',' '))
    xlabel('fSMR volume AU')
    plot_save_name = append(report_table.sample_name(i),"_gating_histo");
%     save_fig(plot_path, plot_save_name, currentFolder);
end

%%
% coulter counter calibration
cc_calibration = 1/1.1353;
fbm_gate_calibrated = fbm_gate;
fbm_gate_calibrated.volume_fL = cc_calibration.*fbm_gate_calibrated.volume_fL;

%%
% mean, median, std dev of each sample
for i = 1:length(report_table.sample_name)
    var=["sample_name","single_cell"];
    annotation = [report_table.sample_name(i),"single_cell"];
    [~,ind_fbm_live,~] = meta_grab_cell_id(fbm_gate_calibrated,metadata_gate,var,annotation);
    
    mean_sc = mean(fbm_gate_calibrated.volume_fL(ind_fbm_live));
    median_sc = median(fbm_gate_calibrated.volume_fL(ind_fbm_live));
    stddev_sc = std(fbm_gate_calibrated.volume_fL(ind_fbm_live));
    
    report_table.mean_vol_fL_single_cell(i) = mean_sc;
    report_table.median_vol_fL_single_cell(i) = median_sc;
    report_table.stddev_vol_fL_single_cell(i) = stddev_sc;
end

%% Output gating result
metadata_name = strrep(input_info.metadata_filename,'base','gated');
fbm_name = strrep(input_info.fbm_filename,'base','calibrated');
time_stamp = datestr(datetime("now"),"yyyymmddTHHMMSS");
metadata_name = strrep(metadata_name,'.txt',['_',time_stamp,'.txt']);
fbm_name = strrep(fbm_name,'.txt',['_',time_stamp,'.txt']);

report_name = strrep(metadata_name,'metadata','report');
report_name = strrep(report_name,'gated','gating');

output_dir = [input_info.metadata_dir,'\',time_stamp,'\'];
mkdir(output_dir)

cd(output_dir)
    writetable(fbm_gate_calibrated, fbm_name,'Delimiter',' ','WriteRowNames',true)
    disp('Calibrated fbm top rows:')
    head(fbm_gate_calibrated)
    writetable(metadata_gate,metadata_name,'Delimiter',' ','WriteRowNames',true)
    disp('Gated metadata top rows:')
    head(metadata_gate)
    writetable(report_table,report_name,'Delimiter',' ','WriteRowNames',true)
    disp('Report table top rows:')
    head(report_table)
cd(currentFolder)
%%
ind_to_show_sample = metadata_gate.sample_name==report_table.sample_name(rpr_ind_bio_sample(1))&...
    metadata_gate.measurement_date==report_table.measurement_date(rpr_ind_bio_sample(1));
ind_to_show_cellstate= metadata_gate.single_cell=="single_cell" & metadata_gate.Viability=="Live";
ind_to_show = ind_to_show_cellstate & ind_to_show_sample;
figure(1)
gscatter(fbm_gate.buoyant_mass(ind_to_show),fbm_gate.pmt2(ind_to_show),metadata_gate.cell_cycle(ind_to_show))
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

%%
% figure(1)
% var=["Viability","single_cell","cell_type",'pmt_3_label'];
% annotation = ["Live","single_cell","OCI-AML3","amine"];
% 
% [cell_id_out,ind_fbm_out,~] = meta_grab_cell_id(fbm_gate,metadata_gate,var,annotation);
% var_x = "buoyant_mass";
% var_y = "pmt3";
% x=log10(fbm_gate.(var_x)(ind_fbm_out));
% y=log10(fbm_gate.(var_y)(ind_fbm_out));
% [p,S] = polyfit(x,y,1); 
% [y_fit,delta] = polyval(p,x,S);
% dscatter(x,y,'logy',true,'logx',true,'SMOOTHING',10,'BINS',[3000,2000],'PLOTTYPE','scatter') 
% hold on
% plot(x,y_fit,'r-')
% plot(x,y_fit+2*delta,'m--',x,y_fit-2*delta,'m--')
% symlog('y')
% xlabel(strrep(var_x,'_',' '))
% ylabel(strrep(var_y,'_',' '))
% legend(["sample";['slope =',num2str(p(1))];"95% precision interval"],'Location','southeast')
%%
% 
% function [cell_id_out,ind_fbm_out] = meta_grab_cell_id(fbm_in,metadata_in,var,annotation)
% % Desciption: meta_grab_cell_id function serves to help user quickly grab unique cell idenfiers
% %     as well as corresponding FBM table row indecies from user-defined qualifiers on the metadata
% %     sheet by specifing the metadata variable name(s) and annotation(s) to grab.
% % Required inputs:
% %     fbm_in          : input FBM table which the user wish to grab matching unique cell ids and row indecies
% %     metadata_in     : input metadata table which the user wish to set qualifiers on
% %     var             : a single or array of strings indicating metadata variable names to qualify
% %     annontation     : *same dimension(s) as var*  a single or array of strings indicating specific annotation to grab from each metadata variable
% 
% ind_meta_interest = true;
% 
% for i = 1:length(var)
%     ind_meta_interest = ind_meta_interest & metadata_in.(var(i))==annotation(i);
% end
% 
% cell_id_out = metadata_in.Row(ind_meta_interest);
% 
% [~,ind_fbm_out]=intersect(fbm_in.Row,cell_id_out,'stable');
% 
% end


%%
function [metadata,cutoff] = FBM_gating(display_name,fbm_in,metadata_in,cell_id_to_gate,var_to_gate,annotation_label,annotation, cutoff_input,rapid_gate)
% Description: FBM_gating fucntion serves to draw 1D or square 2D gating on
% any variable(s) from the input FBM matrix. Annonation from the gating is
% defined by user and will be recorded in the output metadata sheet.
% Required inputs:
%     display_name    = display name indicating identity of the data being gated on
%     fbm             = input FBM table with unique cell ids as row names
%     metadata        = input metadata table with unique cell ids as row names
%     cell_id_to_gate = 1D array of table row indecies that user wish to include in the gating action
%     var_to_gate     = 1x1 or 2x1 array including the name(s) of variables from FBM to be gated on
%     annotation_label = string variable for the annotation label that will be added to the output metadata as a new variable 
%     annotation      = string variable for annotation of cells that pass the gating strategy
%     cutoff_input    = 0 if there's no predetermined cutoffs; else, 1x2 or 1x4 numerical array including cutoffs of the gate, 
%         format follows: for 1d gate is [lower_cutoff, higher_cutoff]; for 2d gate is [var1_lower_cutoff, var1_higher_cutoff,var2_lower_cutoff, var2_higher_cutoff]
%     rapid_gate      = 0 or 1 ; 1 - to rapidly gate out cells without graphical output or seeking user approval of the gated results
scrsize = get(0, 'Screensize');
fbm = fbm_in;
metadata = metadata_in;


[~,ind_to_gate]=intersect(fbm.Row,cell_id_to_gate,'stable');
gate_fbm = fbm(ind_to_gate,:);


% set parameteres for 1D or 2D gating
if length(var_to_gate)==1
    cutoff_initial = [0,0];
    fig_outpos = [0.2*scrsize(3) 0.2*scrsize(4) 0.6*scrsize(3) 0.5*scrsize(4)];
    cutoff_input_message =['Cutoffs? [',char(var_to_gate),' low, ',char(var_to_gate),' high]'];
elseif length(var_to_gate)==2
    cutoff_initial = [0,0,0,0];
    fig_outpos = [0.2*scrsize(3) 0.2*scrsize(4) 0.6*scrsize(3) 0.5*scrsize(4)];
    cutoff_input_message =['Cutoffs? [',char(var_to_gate(1)),' low, ',char(var_to_gate(1)),' high, ',...
        char(var_to_gate(2)),' low, ',char(var_to_gate(2)),' high]'];
end


if cutoff_input == 0
    cutoff = cutoff_initial;
else
    cutoff = cutoff_input;
end


if rapid_gate == 1
    % Apply gating
        if length(var_to_gate)==1
            % 1D gate
            row_ind_gate_fbm = gate_fbm.(var_to_gate)>cutoff(1)&gate_fbm.(var_to_gate)<cutoff(2);
            cell_id_gate = gate_fbm.Row(row_ind_gate_fbm);
        elseif length(var_to_gate)==2
            % 2D gate
            TF1 = gate_fbm.(var_to_gate(1))>cutoff(1);
            TF2 = gate_fbm.(var_to_gate(1))<cutoff(2);
            TF3 = gate_fbm.(var_to_gate(2))>cutoff(3);
            TF4 = gate_fbm.(var_to_gate(2))<cutoff(4);
            % combine them
            row_ind_gate_fbm = TF1 & TF2 & TF3& TF4;
        end  
    cell_id_gate = gate_fbm.Row(row_ind_gate_fbm);
else
    pass_flag = 0;
    while pass_flag ~= 1
        figure('OuterPosition',fig_outpos);
        fig = gcf;
        % plot pre-gated data
        if length(var_to_gate)==1
            subplot(1,2,1)
                h_gate = histogram(gate_fbm.(var_to_gate));
                high_exclude_pct = 99.9;
                plot_lim_lower = prctile(gate_fbm.(var_to_gate),100-high_exclude_pct);
                plot_lim_higher = prctile(gate_fbm.(var_to_gate),high_exclude_pct);
                xlim([plot_lim_lower ,plot_lim_higher])
                xlabel(strrep(var_to_gate,'_',' '))
                title([strrep(display_name,'_',' '),' Pre-gating'])
        elseif length(var_to_gate)==2
            subplot(1,2,1)
                dscatter(gate_fbm.(var_to_gate(1)),gate_fbm.(var_to_gate(2)),'logy',true,'logx',true,'SMOOTHING',10,'BINS',[3000,2000],'PLOTTYPE','scatter')
                %scatter(gate_fbm.(var_to_gate(1)),gate_fbm.(var_to_gate(2)),5,'filled')
                high_exclude_pct = 99.9;
                xlim([prctile(gate_fbm.(var_to_gate(1)),100-high_exclude_pct) ,prctile(gate_fbm.(var_to_gate(1)),high_exclude_pct)])
                ylim([prctile(gate_fbm.(var_to_gate(2)),100-high_exclude_pct) ,prctile(gate_fbm.(var_to_gate(2)),high_exclude_pct)])
                xlabel(strrep(var_to_gate(1),'_',' '))
                ylabel(strrep(var_to_gate(2),'_',' '))
                title([strrep(display_name,'_',' '),' Pre-gating'])
                symlog()
        end
        
        % seek cutoffs from user
        cutoff = input(cutoff_input_message);

        % Apply gating
        if length(var_to_gate)==1
            % 1D gate
            row_ind_gate_fbm = gate_fbm.(var_to_gate)>cutoff(1)&gate_fbm.(var_to_gate)<cutoff(2);
            cell_id_gate = gate_fbm.Row(row_ind_gate_fbm);
        elseif length(var_to_gate)==2
            % 2D gate
            TF1 = gate_fbm.(var_to_gate(1))>cutoff(1);
            TF2 = gate_fbm.(var_to_gate(1))<cutoff(2);
            TF3 = gate_fbm.(var_to_gate(2))>cutoff(3);
            TF4 = gate_fbm.(var_to_gate(2))<cutoff(4);
            % combine them
            row_ind_gate_fbm = TF1 & TF2 & TF3& TF4;
        end  
        % plot gating results
        if length(var_to_gate)==1
            subplot(1,2,2)
                h_gated_in = histogram(gate_fbm.(var_to_gate)(row_ind_gate_fbm));
                hold on
                h_gated_out = histogram(gate_fbm.(var_to_gate)(~row_ind_gate_fbm));
                h_gated_in.BinEdges = h_gate.BinEdges;
                h_gated_out.BinEdges = h_gate.BinEdges;
                xlim([plot_lim_lower ,plot_lim_higher])
                xlabel(strrep(var_to_gate,'_',' '))
                title([strrep(display_name,'_',' '),' ',strrep(char(annotation_label),'_',' '),'-gated'])
                legend([strrep(annotation,'_',' '),'other'])
        elseif length(var_to_gate)==2
            subplot(1,2,2)
                scatter(gate_fbm.(var_to_gate(1))(row_ind_gate_fbm),gate_fbm.(var_to_gate(2))(row_ind_gate_fbm),5,'filled')
                scatter(gate_fbm.(var_to_gate(1))(~row_ind_gate_fbm),gate_fbm.(var_to_gate(2))(~row_ind_gate_fbm),5,'filled')
                hold on
                xlim([prctile(gate_fbm.(var_to_gate(1)),100-high_exclude_pct) ,prctile(gate_fbm.(var_to_gate(1)),high_exclude_pct)])
                ylim([prctile(gate_fbm.(var_to_gate(2)),100-high_exclude_pct) ,prctile(gate_fbm.(var_to_gate(2)),high_exclude_pct)])
                xlabel(strrep(var_to_gate(1),'_',' '))
                ylabel(strrep(var_to_gate(2),'_',' '))
                title([strrep(display_name,'_',' '),' ',strrep(char(annotation_label),'_',' '),'-gated'])
                legend([strrep(annotation,'_',' '),'other'])
                symlog()
            
        end    
        pass_flag = input('Pass? 1 - pass, 0 - reset cutoffs');
        close(fig)
    end
end
% annotate the metadata
cell_id_gate = gate_fbm.Row(row_ind_gate_fbm);
[~,row_ind_gate_meta]=intersect(metadata.Row,cell_id_gate,'stable');
metadata.(annotation_label)(row_ind_gate_meta)=annotation;

end

%%
function save_fig(save_path,title_name,current_folder)
    cd(save_path)
    print(gcf,[char(title_name),'.png'],'-dpng','-r800');
    cd(current_folder)
    close all
end





















