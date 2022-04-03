%% Read excel files (Readout_paird)
clear all
close all
clc
currentFolder=pwd;

%read file
[ReadoutPaired.file_1, ReadoutPaired.path_1] = uigetfile({'*.*','Select the data file.'});
%fid=fopen(strcat(Coulter.path_1,Coulter.file_1),'r');
ReadoutPaired_table = readtable(strcat(ReadoutPaired.path_1, ReadoutPaired.file_1));
ReadoutPaired_matrix = ReadoutPaired_table{:,:};
sample_name= extractBefore(ReadoutPaired.file_1,'.csv');  

BM=ReadoutPaired_matrix(:,2); % buoyant mass
FITC=ReadoutPaired_matrix(:, 4); % FITC
PE=ReadoutPaired_matrix(:,5);  %PE
APC=ReadoutPaired_matrix(:,6); %APC
ND=ReadoutPaired_matrix(:,9); %ND

%% Read columns
PMT_to_fL_conversion_factor=10.6148;
ReadoutPaired_matrix(:,10)= ReadoutPaired_matrix(:, 4)*PMT_to_fL_conversion_factor;% volume
ReadoutPaired_matrix(:,11)= 1.005 + ReadoutPaired_matrix(:,2)./ReadoutPaired_matrix(:,10); % density
Volume=ReadoutPaired_matrix(:,10); % cell volume
Density=ReadoutPaired_matrix(:,11); % cell density

%% Plotting density vs volume
figure(1)
scatter(Volume, Density,12, BM,'fill')
xlim([0, 1000])
ylim([1.01, 1.07])
set(gca,'FontSize',12)
xlabel('Cell volume (fL)','FontSize',18) 
ylabel('Cell density (pg/mL)','FontSize',18)
colorbar

%% Save
PairedData = array2table(ReadoutPaired_matrix, 'VariableNames', ["Time", "Mass","PacificBlue","FITC","PE","APC","Cy7","Transit_time","ND","Volume","Density"]);
cd(ReadoutPaired.path_1)
out_file_name = [sample_name '_processed.xlsx'];
writetable(PairedData, out_file_name)
cd(currentFolder)

