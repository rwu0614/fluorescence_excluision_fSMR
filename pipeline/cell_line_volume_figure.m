%%
% import data (both report tables of fSMR and coulter counter data

%fSMR
fprintf('\nGetting gating report file...')
[input_info.report_filename_fSMR, input_info.report_dir_fSMR,~] = uigetfile('../*.*','Select gating report .txt File', ' ');
fprintf('\n%s selected for analysis\n', input_info.report_filename_fSMR)
report_file_fSMR_fullname = strcat(input_info.report_dir_fSMR, input_info.report_filename_fSMR); 
report_table_fSMR = readtable(report_file_fSMR_fullname,'ReadRowNames',false,'ReadVariableNames',true,'Delimiter',' ');
head(report_table_fSMR)

%Coulter counter
fprintf('\nGetting gating report file...')
[input_info.report_filename_CC, input_info.report_dir_CC,~] = uigetfile('../*.*','Select gating report .txt File', ' ');
fprintf('\n%s selected for analysis\n', input_info.report_filename_CC)
report_file_CC_fullname = strcat(input_info.report_dir_CC, input_info.report_filename_CC); 
report_table_CC = readtable(report_file_CC_fullname,'ReadRowNames',false,'ReadVariableNames',true,'Delimiter',' ');
head(report_table_CC)

plot_path = uigetdir('../*.*','Where to save the fitting data?');

%%

idx = ismember(report_table_CC.cell_type, report_table_fSMR.cell_type);

idx_beads = find(idx==0);

report_table_CC(idx_beads,:) = [];

%%
% scatter(report_table_CC.median_vol_fL_single_cell, report_table_fSMR.median_vol_au_single_cell, 10)
% hold on

x = report_table_CC.median_vol_fL_single_cell;
y = report_table_fSMR.median_vol_au_single_cell;
gscatter(report_table_CC.median_vol_fL_single_cell, report_table_fSMR.median_vol_au_single_cell, report_table_fSMR.cell_type,[0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330],'.',20)
Fit = polyfit(report_table_CC.median_vol_fL_single_cell,report_table_fSMR.median_vol_au_single_cell,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
hold on
fit_x = linspace(750,4300);
z = (Fit(1).*fit_x) + Fit(2);
plot(fit_x, z,'k--')
% errorbar(x(),y,yneg,ypos,xneg,xpos,'LineStyle','None')
legend('patu-8902','BAF3','DT40','L1210','sHela','THP1',['y = 0.077x + 12.663' newline 'R^2 = 0.99']);
xlabel('Coulter counter volume (fL)'); ylabel('fSMR volume (AU)');
title('Coulter counter vs fSMR volumes');
currentFolder = pwd;
% save_fig(plot_path,'fSMR_vs_cc_vol_all_reps_trendline',currentFolder)

%%
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330]};

yneg = report_table_fSMR.stddev_vol_au_single_cell;
ypos = yneg;
xneg = report_table_CC.stddev_vol_fL_single_cell;
xpos = xneg;
categories = unique(report_table_fSMR.cell_type);


% gscatter(x,y,report_table_fSMR.cell_type,[0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330],'.')
% hold on
for k = 1:length(categories)
    idx = ismember(report_table_CC.cell_type, categories(k));
    idx_cell_type = find(idx==1);
    for i = 1:length(idx_cell_type)
        errorbar(x(idx_cell_type(i)),y(idx_cell_type(i)),yneg(idx_cell_type(i)),ypos(idx_cell_type(i)),xneg(idx_cell_type(i)),xpos(idx_cell_type(i)),'.','MarkerSize',12,'LineWidth',1,'Color',colors{k})
        hold on
    end
end

xlabel('Coulter counter volume (fL)'); ylabel('fSMR volume (AU)');
title('Coulter counter vs fSMR volumes');

save_fig(plot_path,'fSMR_vs_cc_vol_all_stddev',currentFolder)


%%
% averages
cell_type = unique(report_table_fSMR.cell_type,'stable');
averages = table(cell_type);
k = 1;
for i = 1:3:length(report_table_fSMR.sample_name)
    avg_fSMR = mean([report_table_fSMR.median_vol_au_single_cell(i) report_table_fSMR.median_vol_au_single_cell(i+1) report_table_fSMR.median_vol_au_single_cell(i+2)]);
    avg_CC = mean([report_table_CC.median_vol_fL_single_cell(i) report_table_CC.median_vol_fL_single_cell(i+1) report_table_CC.median_vol_fL_single_cell(i+2)]);
    stderror_fSMR = std([report_table_fSMR.median_vol_au_single_cell(i) report_table_fSMR.median_vol_au_single_cell(i+1) report_table_fSMR.median_vol_au_single_cell(i+2)])./sqrt(3);
    stderror_CC = std([report_table_CC.median_vol_fL_single_cell(i) report_table_CC.median_vol_fL_single_cell(i+1) report_table_CC.median_vol_fL_single_cell(i+2)])./sqrt(3);
    averages.avg_vol_fSMR_AU(k) = avg_fSMR;
    averages.avg_vol_CC_fL(k) = avg_CC;
    averages.stderror_fSMR(k) = stderror_fSMR;
    averages.stderror_CC(k) = stderror_CC;
    k=k+1;
end
%%
x = averages.avg_vol_CC_fL;
y = averages.avg_vol_fSMR_AU;
yneg = averages.stderror_fSMR;
ypos = yneg;
xneg = averages.stderror_CC;
xpos = xneg;
for i = 1:6
    errorbar(x(i), y(i),yneg(i),ypos(i),xneg(i),xpos(i),'.','MarkerSize',12,'LineWidth',1)
    hold on
end
hold off
% Fit = polyfit(report_table_CC.median_vol_fL_single_cell,report_table_fSMR.median_vol_au_single_cell,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
% hold on
% z = Fit(1).*report_table_CC.median_vol_fL_single_cell + Fit(2);
% plot(report_table_CC.median_vol_fL_single_cell, z)
% errorbar(x(),y,yneg,ypos,xneg,xpos,'LineStyle','None')
legend('patu-8902','BAF3','DT40','L1210','sHela','THP1','Trendline','Location','southeast');
xlabel('Coulter counter volume (fL)'); ylabel('fSMR volume (AU)');
title('Average coulter counter vs fSMR volumes');

save_fig(plot_path,'avg_fSMR_vs_cc_vol_2',currentFolder)

%%
function save_fig(save_path,title_name,current_folder)
    cd(save_path)
    print(gcf,[char(title_name),'.png'],'-dpng','-r800');
    cd(current_folder)
    close all
end