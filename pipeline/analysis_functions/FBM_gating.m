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
                hold on
                scatter(gate_fbm.(var_to_gate(1))(~row_ind_gate_fbm),gate_fbm.(var_to_gate(2))(~row_ind_gate_fbm),5,'filled')
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