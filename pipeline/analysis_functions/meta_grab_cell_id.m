function [cell_id_fbm,ind_fbm_out,ind_meta_out] = meta_grab_cell_id(fbm_in,metadata_in,var,annotation)
% Desciption: meta_grab_cell_id function serves to help user quickly grab unique cell idenfiers
%     as well as corresponding FBM table row indecies from user-defined qualifiers on the metadata
%     sheet by specifing the metadata variable name(s) and annotation(s) to grab.
% Required inputs:
%     fbm_in          : input FBM table which the user wish to grab matching unique cell ids and row indecies
%     metadata_in     : input metadata table which the user wish to set qualifiers on
%     var             : a single or array of strings indicating metadata variable names to qualify
%     annontation     : *same dimension(s) as var*  a single or array of strings indicating specific annotation to grab from each metadata variable

ind_meta_interest = true;

for i = 1:length(var)
    ind_meta_interest = ind_meta_interest & metadata_in.(var(i))==annotation(i);
end

cell_id_out = metadata_in.Row(ind_meta_interest);

[cell_id_fbm,ind_fbm_out]=intersect(fbm_in.Row,cell_id_out,'stable');

[~,ind_meta_out] = intersect(metadata_in.Row,cell_id_fbm,'stable');

if length(cell_id_fbm)~=length(cell_id_out)
    warning('Warning: input FBM matrix does not contain all the cells that satisfy input metadata requirements')
    warning_count = [num2str(length(cell_id_out)-length(cell_id_fbm)),' out of ',num2str(length(cell_id_out)),' cells are dropped'];
    warning(warning_count)
end
end