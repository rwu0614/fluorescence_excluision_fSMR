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
fbm_base = readtable(fbm_file_fullname,'ReadRowNames',true,'ReadVariableNames',true,'Delimiter','\t');
head(fbm_base)

% Ask user for base metadata txt file
fprintf('\nGetting base metadata sheet...')
[input_info.metadata_filename, input_info.metadata_dir,~] = uigetfile('../*.*','Select base metadata .txt File', ' ');
fprintf('\n%s selected for analysis\n', input_info.metadata_filename)
metadata_file_fullname = strcat(input_info.metadata_dir, input_info.metadata_filename); 
metadata_base = readtable(metadata_file_fullname,'ReadRowNames',true,'ReadVariableNames',true,'Delimiter','\t');
head(metadata_base)


%Convert all metadata element to string
var = metadata_base.Properties.VariableNames;
for i = var
   metadata_base.(i{:}) = string(metadata_base.(i{:}));
end

%initializing metadata sheet with all the gating annotation
metadata_gate = metadata_base;
%initializing fbm for gating
fbm_gate = fbm_base;
%%
sample_median_vol = 1100; 
        
vol_calibration_factor = sample_median_vol/median(fbm_gate.vol_au);
fbm_gate.volume_fL = fbm_gate.vol_au*vol_calibration_factor;
fluid_den = unique(metadata_base.fluid_density);
for i = 1:length(fluid_den)
    var_grab = ["fluid_density"];
    annotation_grab = [fluid_den(i)];
    [~,ind_fbm_sample_fluid_temp,~] = meta_grab_cell_id(fbm_gate,metadata_gate,var_grab,annotation_grab);
    fbm_gate.density_gcm3(ind_fbm_sample_fluid_temp) = fbm_gate.buoyant_mass_pg(ind_fbm_sample_fluid_temp)./fbm_gate.volume_fL(ind_fbm_sample_fluid_temp)+str2num(fluid_den(i));
    fbm_gate.buoyant_mass_pg_inwater(ind_fbm_sample_fluid_temp) = fbm_gate.buoyant_mass_pg(ind_fbm_sample_fluid_temp) +(str2num(fluid_den(i))-1)*fbm_gate.volume_fL(ind_fbm_sample_fluid_temp);
end
%%
figure(1)
yyaxis left
scatter(fbm_gate.buoyant_mass_pg,fbm_gate.volume_fL,5,'filled','MarkerFaceAlpha',0.3)
ylim([400,2300])
ylabel('Volume (fL)')
hold on
yyaxis right
scatter(fbm_gate.buoyant_mass_pg,fbm_gate.density_gcm3,5,'filled','MarkerFaceAlpha',0.3)
ylim([0.95,1.15])
ylabel('Density (g/mL)')
hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',-90,'VerticalAlignment','bottom')
xlim([30,100])
xlabel('Buoyant mass (pg)')
legend(['n = ',num2str(length(fbm_gate.buoyant_mass_pg))])
%%
% SNACS = NV - m(V_ref-V)
fbm_gate.nv = fbm_gate.node_deviation_hz./fbm_gate.volume_fL;

%%
figure(1)
scatter(fbm_gate.volume_fL,fbm_gate.nv)
xlabel('Volume fL')
ylabel('Node deviation / Volume')
hold on

vol_cutoff_forfit = [776,3148];
nv_cutoff_forfit = [-10^-3,2*10^-3];

for_fit_ind = find(fbm_gate.volume_fL>vol_cutoff_forfit(1)&...
    fbm_gate.volume_fL<vol_cutoff_forfit(2)&...
    fbm_gate.nv>nv_cutoff_forfit(1)&...
    fbm_gate.nv<nv_cutoff_forfit(2));

scatter(fbm_gate.volume_fL(for_fit_ind),fbm_gate.nv(for_fit_ind))

temp_stats = regstats(fbm_gate.nv(for_fit_ind),fbm_gate.volume_fL(for_fit_ind),'linear');

m = temp_stats.beta(2);

fbm_gate.snacs = fbm_gate.nv + m*(median(fbm_gate.volume_fL)-fbm_gate.volume_fL);

%%
figure(3)
scatter(fbm_gate.volume_fL(for_fit_ind),fbm_gate.nv(for_fit_ind),5,fbm_gate.snacs(for_fit_ind))
figure(4)
scatter(fbm_gate.volume_fL,fbm_gate.snacs,5)

%%
figure(5)
scatter(fbm_gate.density_gcm3(for_fit_ind),fbm_gate.snacs(for_fit_ind),5,fbm_gate.volume_fL(for_fit_ind))


