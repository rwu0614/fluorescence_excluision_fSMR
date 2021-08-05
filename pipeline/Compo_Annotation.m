clear all;
clc;
close all;
currentFolder = pwd;
addpath('report_functions\');
addpath('helAPCr_functions\');
addpath('plotting_functions\');

%%
%%==================Initializing=====================%%

%initialize dataframe parameters
cell_id = [];
bm = [];
pmt1 = [];
pmt2 = [];
sample_id = [];
treatment = [];
dose = [];
treatment_time = [];
annotation = [];
group_apop = [];
%create FBM matrix (fluorescence-linked buoyant mass matrix)
FBM_raw = table(cell_id, bm, pmt1, pmt2, sample_id, treatment,...
    dose, treatment_time, annotation, group_apop);


%%
ask_input = input('Want to input raw (1) or saved .mat file (2)?');

if ask_input == 2
    fprintf('\nGetting paired data...\n')
        [paired_filename, paired_dir, exist_paired] = uigetfile('../*.*','Select paired File',' ');
        cd(paired_dir)
        FBM_raw = load(paired_filename);
         FBM_raw = FBM_raw.FBM_raw;
%         FBM_anno = FBM_raw.FBM_anno ;
        cd(currentFolder)
        filename_split = strsplit(paired_filename,'.');
        filename_split = filename_split{1:end-1};
        filename_split = strsplit(filename_split,'_');
        sample_name = filename_split{end};
elseif ask_input == 1

    data_input_flag =0;
    while data_input_flag ==0

        fprintf('\nGetting paired data...\n')
        [paired_filename, paired_dir, exist_paired] = uigetfile('../*.*','Select paired File',' ');
        if  (exist_paired ~= 0)
            filename_split = strsplit(paired_filename,'.');
            filename_split = filename_split{1:end-1};
            filename_split = strsplit(filename_split,'_');
            sample_name = filename_split{end};
            time_st = strsplit(paired_dir,'\');
            sample_name = apAPCnd(sample_name,' ',time_st{end-1});
            sample_split = strsplit(sample_name,' ');

            % grab and initialize apAPCnding parameters
            cd(paired_dir)
            paired_input = readtable(paired_filename);
            cd(currentFolder)
            input_size = size(paired_input);
            cell_id_apAPCnd = strings(input_size(1),1);
            bm_apAPCnd = paired_input{:,2};
            pmt1_apAPCnd = paired_input{:,3};
            pmt2_apAPCnd = paired_input{:,5};
            sample_id_apAPCnd = strings(input_size(1),1);
            treatment_apAPCnd = strings(input_size(1),1);
            dose_apAPCnd = strings(input_size(1),1);
            treatment_time_apAPCnd = strings(input_size(1),1);
            annotation_apAPCnd = strings(input_size(1),1);
            group_apop_apAPCnd = strings(input_size(1),1);

            sample_id_apAPCnd(:) = string(sample_name);
            cell_id_apAPCnd = strcat(sample_id_apAPCnd,'_',string(1:1:input_size(1))');
            if strcmp(sample_split{1},'gem')
                treatment_apAPCnd(:) = 'Gemcitabine';
            elseif strcmp(sample_split{1},'tra')
                treatment_apAPCnd(:) = 'Trametinib';
            elseif strcmp(sample_split{1},'dmso')
                treatment_apAPCnd(:) = 'DMSO';
            end

            if strcmp(sample_split{1},'dmso')
                dose_apAPCnd(:) = 'N/A';
                treatment_time_apAPCnd(:) = sample_split{end};
            else
                dose_apAPCnd(:) = sample_split(2);
                treatment_time_apAPCnd(:) = sample_split{end};
            end

            FBM_apAPCnd = table(cell_id_apAPCnd, bm_apAPCnd, pmt1_apAPCnd, pmt2_apAPCnd,...
            sample_id_apAPCnd, treatment_apAPCnd, dose_apAPCnd, treatment_time_apAPCnd,...
            annotation_apAPCnd, group_apop_apAPCnd,'VariableNames',{'cell_id','bm',...
            'pmt1','pmt2','sample_id','treatment','dose','treatment_time','annotation',...
            'group_apop'});

            FBM_raw = [FBM_raw;FBM_apAPCnd];
            fprintf(['\nInputing ',sample_name,'\n'])
        else
            fprintf('Finishing paired data input...')
            data_input_flag = 1;
        end

    end
end
%%
% readout_paired = [realtime, smr, pmt1, pmt2, pmt3]
sudo_FBM = FBM_raw;
% SAPCcify you conditions
TF1 = sudo_FBM.pmt1< 10^-4 ;
TF2 = sudo_FBM.pmt2< 10^-4 ;
% combine them
TFall = TF1 | TF2 ;
% remove
sudo_FBM(TFall,:) = [];

paired_pmt = [sudo_FBM.pmt1,sudo_FBM.pmt2]'; % row 1 pmt 1, row2 pmt2

%composition matrix
compo_mat = [1,-0.02;-0.0028,1];
%compo_mat = [1,-0.0006;-0.08,1]; % from 8902 single color controls
compo_pmt = compo_mat*paired_pmt;

scrsize = get(0, 'Screensize');

figure (1)
subplot(1,2,1);
     scatter(paired_pmt(2,:),paired_pmt(1,:),5,'r','filled')
    axis manual
%     plot(paired_pmt(2,:),paired_pmt(1,:),'-mo','LineStyle','none',...
%     'MarkerEdgeColor','r',...
%     'MarkerFaceColor','r',...
%     'MarkerSize',2)
    title(sample_name)
    xlabel("APC (Annexin V) mV")
    ylabel("FITC (Calcien AM) mV")
    %set(gca, 'YScale', 'log')
    %set(gca, 'XScale', 'log')
    legend(['n' '=' int2str(length(paired_pmt))],'location',"southeast")
    symlog('xy')
%     colorbar;
%     caxis([0 300]);
    
subplot(1,2,2);
    scatter(compo_pmt(2,:),compo_pmt(1,:),3,sudo_FBM.bm,'filled')
    title('All cells from all timepoints')
    xlabel("APC (Annexin V) mV")
    ylabel("FITC (Calcien AM) mV")
    legend(['n' '=' int2str(length(paired_pmt))],'location',"southeast")
    symlog('xy',0.8)
    cmap = crameri('roma');
    set(gca,'ColorScale','log')
    colormap(cmap);
    c=colorbar;
    caxis([5 500]);
    set(get(c,'title'),'string','Buoyant mass(pg)','Rotation',0);
%%
 change_time =["6hr","12hr","24hr","48hr"]; 
% change_time =["6hr","12hr"];   
for i = 1:length(change_time)
    change_hztopg_FBM =  sudo_FBM(sudo_FBM.treatment_time==change_time(i),:);
     change_hztopg_FBM.bm = change_hztopg_FBM.bm.*(0.59553/0.8);
%      change_hztopg_FBM.bm = change_hztopg_FBM.bm.*(0.7/0.8);
    sudo_FBM(sudo_FBM.treatment_time==change_time(i),:)=change_hztopg_FBM;
end
%%
sudo_FBM.pmt1= compo_pmt(1,:)';
sudo_FBM.pmt2= compo_pmt(2,:)';

figure (2)
sub1 = subplot(1,2,1);
%cmap = crameri('vanimo');
     cmap = viridis(10);
    scatter(sudo_FBM.bm,sudo_FBM.pmt2,3,sudo_FBM.pmt1,'filled',...
    'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
    title('Calcien AM projection on all cells')
    xlabel("Buoyant mass (pg)")
    ylabel("APC (Annexin V) mV")
    symlog('xy')
    colormap(sub1, cmap);
    c=colorbar;
    set(gca,'ColorScale','log')
    caxis([1 100]);
    %xlim([0,400]);
    legend(['n' '=' int2str(length(paired_pmt))],'location',"southeast")
    set(get(c,'title'),'string','FITC (Calcien AM) mV','Rotation',0);
sub2 = subplot(1,2,2);
%cmap = crameri('vanimo');
    %cmap = crameri('buda');
    cmap = tab20(5);
    scatter(sudo_FBM.bm,sudo_FBM.pmt1,3,sudo_FBM.pmt2,'filled',...
    'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
    title('Annexin V projection on all cells')
    xlabel("Buoyant mass (pg)")
    ylabel("FITC (Calcien AM) mV")
    symlog('xy')
    colormap(sub2,cmap);
    c=colorbar;
    set(gca,'ColorScale','log')
    caxis([2 100]);
    %xlim([0,400]);
    legend(['n' '=' int2str(length(paired_pmt))],'location',"southeast")
    set(get(c,'title'),'string','APC (Annexin V) mV','Rotation',0);
%%
figure (3)
dscatter(sudo_FBM.bm,sudo_FBM.pmt2,'logy',true,'logx',true,'SMOOTHING',10,'BINS',[3000,2000],'PLOTTYPE','scatter')
%dscatter(sudo_FBM.bm,sudo_FBM.pmt2,'SMOOTHING',10,'BINS',[3000,2000],'PLOTTYAPC','scatter')
title('Density projection on all cells')
xlabel("Buoyant mass (pg)")
ylabel("APC (Annexin V) mV")
symlog('xy')
caxis([0.01 0.18]);
c=colorbar;
legend(['n' '=' int2str(length(paired_pmt))],'location',"northeast")
set(get(c,'title'),'string','Density','Rotation',0);
%%
figure (4)
dscatter(sudo_FBM.pmt2,sudo_FBM.pmt1,'logy',true,'logx',true,'SMOOTHING',10,'BINS',[3000,2000],'PLOTTYPE','scatter')
title('All cells from all timepoints')
xlabel("APC (Annexin V) mV")
ylabel("FITC (Calcien AM) mV")
symlog('xy')
caxis([0.01 0.21]);
c=colorbar;
%legend(['n' '=' int2str(length(paired_pmt))],'location',"northeast")
set(get(c,'title'),'string','Density','Rotation',0);

%%
figure (5)
dscatter(sudo_FBM.bm,sudo_FBM.pmt1,'logy',true,'logx',true,'SMOOTHING',10,'BINS',[3000,2000],'PLOTTYPE','scatter')
title(sample_name)
xlabel("Buoyant mass (pg)")
ylabel("FITC (Calcien AM) mV")
symlog('xy')
caxis([0.01 0.15]);
colorbar;
legend(['n' '=' int2str(length(paired_pmt))],'location',"northeast")

%%
% Time point batch effect evaluation
figure(6)
timepoints = ["6hr","12hr","24hr","48hr"];
for i = 1:length(timepoints) 
    temp_FBM = sudo_FBM(sudo_FBM.treatment_time == timepoints(i),:);
    subplot(1,4,i);
        h = gscatter(temp_FBM.pmt2,temp_FBM.pmt1,string(temp_FBM.treatment),hsv(3),'.',2,'of');
        title(timepoints(i))
        xlabel("APC (Annexin V) mV")
        ylabel("FITC (Calcien AM) mV")
        xlim([2*10^-2,2*10^3])
        ylim([2*10^-1,0.5*10^4])
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')

end
%legend('location','northeastoutside','FontSize',14)
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')

%%
% All samples dying population drifts
scrsize = get(0, 'Screensize');
figure('OuterPosition',[0*scrsize(3) 0.01*scrsize(4) 0.32*scrsize(3) 1*scrsize(4)]);

timepoints = ["6hr","24hr","48hr"];
dosage = ["N/A","50uM","10uM","100nM","32nM"];
%dosage = ["N/A"];

ha = tight_subplot(5,3,[.06 .08],[.07 .02],[.08 .05]);

p=1;
for d = 1:length(dosage)
    sample_FBM = sudo_FBM(sudo_FBM.dose == dosage(d),:);
    for i = 1:length(timepoints) 
        temp_FBM = sample_FBM(sample_FBM.treatment_time == timepoints(i),:);
        axes(ha(p));
            %scatter(temp_FBM.bm,temp_FBM.pmt2,1,'r')
            dscatter(temp_FBM.bm,temp_FBM.pmt2,'logy',true,'logx',true,'SMOOTHING',10,'BINS',[3000,2000],'PLOTTYPE','scatter')
            title(temp_FBM.sample_id(1))
            xlabel("Buoyant mass (pg)")
            ylabel("APC (Annexin V) mV")
            xlim([2,1000])
            ylim([0,2000])
            symlog('xy')
            
        p=p+1;
    end
end
%legend('location','northeastoutside','FontSize',14)
%%
% All samples dying population drifts
figure('OuterPosition',[0*scrsize(3) 0.01*scrsize(4) 0.32*scrsize(3) 1*scrsize(4)]);
ha = tight_subplot(5,3,[.06 .08],[.07 .02],[.08 .05]);
timepoints = ["6hr","24hr","48hr"];
dosage = ["N/A","50uM","10uM","100nM","32nM"];
%dosage = ["N/A"];
p=1;
for d = 1:length(dosage)
    sample_FBM = sudo_FBM(sudo_FBM.dose == dosage(d),:);
    for i = 1:length(timepoints) 
        temp_FBM = sample_FBM(sample_FBM.treatment_time == timepoints(i),:);
        axes(ha(p));
            %scatter(temp_FBM.bm,temp_FBM.pmt2,1,'r')
            dscatter(temp_FBM.bm,temp_FBM.pmt1,'logy',true,'logx',true,'SMOOTHING',10,'BINS',[3000,2000],'PLOTTYPE','scatter')
            title(temp_FBM.sample_id(1))
            xlabel("Buoyant mass (pg)")
            ylabel("FITC (Calcien AM) mV")
            xlim([5,1000])
            ylim([0,5000])
            symlog('xy')
            
        p=p+1;
    end
end
%legend('location','northeastoutside','FontSize',14)
%%
%%==================annotation====================%%
FBM_anno = sudo_FBM;
% 
% figure(10)
% sample_id_list = unique(sudo_FBM.sample_id);
% i=1;
% while i <= length(sample_id_list)
%     sample_FBM = sudo_FBM(sudo_FBM.sample_id == sample_id_list(i),:);
%     sub1 = subplot(1,2,1);
%         cmap = viridis(5);
%         scatter(sample_FBM.bm,sample_FBM.pmt2,10,sample_FBM.pmt1,'filled')
%         title(sample_id_list(i))
%         xlabel("Buoyant mass (pg)")
%         ylabel("APC (Annexin V) mV")
%         %symlog('xy')
%         set(gca, 'YScale', 'log')
%         set(gca, 'XScale', 'log')
%         colormap( cmap);
%         colorbar;
%         set(gca,'ColorScale','log')
%         caxis([1 10]);
%         legend(['n' '=' int2str(length(paired_pmt))],'location',"southeast")
%     thresh = input('[lower_annex_thresh,upper_annex_thresh,upper_cal_thresh,lower_cal_thresh]:\n');
%     lower_annex_thresh = thresh(1);
%     upper_annex_thresh = thresh(2);
%     upper_cal_thresh = thresh(3);
%     lower_cal_thresh = thresh(4);
% %     dot1 = input('lower x-y cut [x1,y1]:\n');x1=dot1(1);y1=dot1(2);
% %     dot2 = input('upper x-y cut [x2,y2]:\n');x2=dot2(1);y2=dot2(2);
% %     p1=(y2-y1)/(x2-x1);
% %     p2= y1 -(y2-y1)/(x2-x1)*x1;
%     gate_FBM = sample_FBM;
%     gate_FBM.annotation(:) = 'not_live';
%         % SAPCcify you conditions
%         TF1 = gate_FBM.pmt2<upper_cal_thresh;
%         TF2 = gate_FBM.bm<upper_annex_thresh;
%         TF3 = gate_FBM.bm>lower_annex_thresh;
%         TF4 = gate_FBM.pmt1>lower_cal_thresh;
%         %TF4 = gate_FBM.pmt2<gate_FBM.pmt1*p1+p2;
%         % combine them
%         TFall = TF1 & TF2 & TF3& TF4;
%         % remove
%         gate_FBM.annotation(TFall) = 'live';
%     sub2 = subplot(1,2,2);
%         gscatter(gate_FBM.bm,gate_FBM.pmt2,gate_FBM.annotation,'rg','.',10,'off')
%         title(sample_id_list(i))
%         ylabel("APC (Annexin V) mV")
%         xlabel("Buoyant mass (pg)")
%         set(gca, 'YScale', 'log')
%         set(gca, 'XScale', 'log')
%     
%         Anno = FBM_anno.sample_id == sample_id_list(i);
%         FBM_anno.annotation(Anno) = gate_FBM.annotation;
%         pass = input('pass?(yes-1,no-2)');
%         if pass==1
%             i=i+1;
%         end
% end
% %legend('location','northeastoutside','FontSize',14)
sample_FBM = FBM_anno;

    lower_annex_thresh = -100;
    upper_annex_thresh = 60;
    upper_cal_thresh = 5000;
    lower_cal_thresh = 70;
    
    gate_FBM = sample_FBM;
    gate_FBM.annotation(:) = 'not_live';
        % SAPCcify you conditions
        TF1 = gate_FBM.pmt1<upper_cal_thresh;
        TF2 = gate_FBM.pmt2<upper_annex_thresh;
        TF3 = gate_FBM.pmt2>lower_annex_thresh;
        TF4 = gate_FBM.pmt1>lower_cal_thresh;
        %TF4 = gate_FBM.pmt2<gate_FBM.pmt1*p1+p2;
        % combine them
        TFall = TF1 & TF2 & TF3& TF4;
        % remove
        gate_FBM.annotation(TFall) = 'live';
        FBM_anno = gate_FBM;

%%
gscatter(FBM_anno.bm,FBM_anno.pmt2,FBM_anno.annotation,'rg','.',2,'off')
        title(sample_id_list(i))
        ylabel("APC (Annexin V) mV")
        xlabel("Buoyant mass (pg)")
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
%%
% annotation check
FBM_anno= sortrows(FBM_anno,{'annotation'},{'ascend'});

figure(9)
timepoints = ["6hr","12hr","24hr","48hr"];
%dosage = ["N/A","50uM","10uM","100nM","32nM"];
dosage = ["N/A"];
p=1;
for d = 1:length(dosage)
    sample_FBM = FBM_anno(FBM_anno.dose == dosage(d),:);
    for i = 1:length(timepoints) 
        temp_FBM = sample_FBM(sample_FBM.treatment_time == timepoints(i),:);
        subplot(1,6,p);
        gscatter(temp_FBM.bm,temp_FBM.pmt2,temp_FBM.annotation,tab10(2),'.',1,'off','siz',15)
        title(temp_FBM.sample_id(1))
        ylabel("APC (Annexin V) mV")
        xlabel("Buoyant mass (pg)")
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
        xlim([5,1000])
        ylim([0,2000])
        p=p+1;
    end
end

%%
timepoints = ["6hr","12hr","24hr","36hr","48hr","72hr"];
dosage = ["N/A","50uM","10uM","100nM","32nM"];
mean_bm_time_dmso = zeros(length(timepoints),1);
mean_bm_time_gem_50uM = zeros(length(timepoints),1);
mean_bm_time_gem_10uM = zeros(length(timepoints),1);
mean_bm_time_tra_100nM = zeros(length(timepoints),1);
mean_bm_time_tra_32nM = zeros(length(timepoints),1);

    for i = 1:length(timepoints) 
        temp_FBM = FBM_live(FBM_live.treatment_time == timepoints(i),:);
        
        mean_bm_time_dmso(i) = mean(temp_FBM(temp_FBM.dose == dosage(1),:).bm);
        mean_bm_time_gem_50uM(i) = mean(temp_FBM(temp_FBM.dose == dosage(2),:).bm);
        mean_bm_time_gem_10uM(i) = mean(temp_FBM(temp_FBM.dose == dosage(3),:).bm);
        mean_bm_time_tra_100nM(i) = mean(temp_FBM(temp_FBM.dose == dosage(4),:).bm);
        mean_bm_time_tra_32nM(i) = mean(temp_FBM(temp_FBM.dose == dosage(5),:).bm);
    end

figure(10)
time = [6,24,48];
plot(time,mean_bm_time_dmso,'--o','color',[0.466,0.674,0.188],'MarkerEdgeColor',[0.466,0.674,0.188],...
    'MarkerFaceColor',[0.466,0.674,0.188])
legend(["DMSO","Gem 50uM","Gem 10uM","Tra 100nM","Tra 32nM"])
hold on
plot(time,mean_bm_time_gem_50uM,'-o','color',[0.929,0.694,0.125],'MarkerEdgeColor',[0.929,0.694,0.125],...
    'MarkerFaceColor',[0.929,0.694,0.125])
legend(["DMSO","Gem 50uM","Gem 10uM","Tra 100nM","Tra 32nM"])
hold on
plot(time,mean_bm_time_gem_10uM,'-o','color',[0.829,0.694,0.125],'MarkerEdgeColor',[0.829,0.694,0.125],...
    'MarkerFaceColor',[0.829,0.694,0.125])
legend(["DMSO","Gem 50uM","Gem 10uM","Tra 100nM","Tra 32nM"])
hold on
plot(time,mean_bm_time_tra_100nM,'-o','color',[0.494,0.184,0.556],'MarkerEdgeColor',[0.494,0.184,0.556],...
    'MarkerFaceColor',[0.494,0.184,0.556])
legend(["DMSO","Gem 50uM","Gem 10uM","Tra 100nM","Tra 32nM"])
hold on
plot(time,mean_bm_time_tra_32nM,'-o','color',[0.594,0.184,0.556],'MarkerEdgeColor',[0.594,0.184,0.556],...
    'MarkerFaceColor',[0.594,0.184,0.556])
hold off
legend(["DMSO","Gem 50uM","Gem 10uM","Tra 100nM","Tra 32nM"],'location','southeast')
xlabel('Hours post-treatment')
ylabel('Buoyant mass (pg)')
title('Live cell populations')


%%

timepoints = ["6hr","24hr","48hr"];

hellinger_bm_time_gem_50uM = zeros(length(timepoints),1);
hellinger_bm_time_gem_10uM = zeros(length(timepoints),1);
hellinger_bm_time_tra_100nM = zeros(length(timepoints),1);
hellinger_bm_time_tra_32nM = zeros(length(timepoints),1);

    for i = 1:length(timepoints) 
        hellinger_bm_time_gem_50uM(i) = helingerall.Gem50(i);
        hellinger_bm_time_gem_10uM(i) = helingerall.Gem10(i);
        hellinger_bm_time_tra_100nM(i) = helingerall.Tra100(i);
        hellinger_bm_time_tra_32nM(i) = helingerall.Tra32(i);
    end

figure(12)
time = [6,24,48];
plot(time,hellinger_bm_time_gem_50uM,'-o','color',[0.929,0.694,0.125],'MarkerEdgeColor',[0.929,0.694,0.125],...
    'MarkerFaceColor',[0.929,0.694,0.125])
legend(["DMSO","Gem 50uM","Gem 10uM","Tra 100nM","Tra 32nM"])
hold on
plot(time,hellinger_bm_time_gem_10uM,'-o','color',[0.829,0.694,0.125],'MarkerEdgeColor',[0.829,0.694,0.125],...
    'MarkerFaceColor',[0.829,0.694,0.125])
legend(["DMSO","Gem 50uM","Gem 10uM","Tra 100nM","Tra 32nM"])
hold on
plot(time,hellinger_bm_time_tra_100nM,'-o','color',[0.494,0.184,0.556],'MarkerEdgeColor',[0.494,0.184,0.556],...
    'MarkerFaceColor',[0.494,0.184,0.556])
legend(["DMSO","Gem 50uM","Gem 10uM","Tra 100nM","Tra 32nM"])
hold on
plot(time,hellinger_bm_time_tra_32nM,'-o','color',[0.594,0.184,0.556],'MarkerEdgeColor',[0.594,0.184,0.556],...
    'MarkerFaceColor',[0.594,0.184,0.556])
hold off
legend(["Gem 50uM","Gem 10uM","Tra 100nM","Tra 32nM"],'location','northeast')
xlabel('Hours post-treatment')
ylabel('Hellinger distance to DMSO')
title('All cells')

%%

timepoints = ["6hr","12hr","24hr","36hr","48hr","72hr"];

hellinger_bm_time_gem_50uM = zeros(length(timepoints),1);
hellinger_bm_time_gem_10uM = zeros(length(timepoints),1);
hellinger_bm_time_tra_100nM = zeros(length(timepoints),1);
hellinger_bm_time_tra_32nM = zeros(length(timepoints),1);

    for i = 1:length(timepoints) 
        hellinger_bm_time_gem_50uM(i) = helingerall3.V1(i);
        hellinger_bm_time_gem_10uM(i) = helingerall3.V2(i);
        hellinger_bm_time_tra_100nM(i) = helingerall3.V3(i);
        hellinger_bm_time_tra_32nM(i) = helingerall3.V4(i);
    end

figure(13)
time = [6,12,24,36,48,72];
plot(time,hellinger_bm_time_gem_50uM,'-o','color',[0.929,0.694,0.125],'MarkerEdgeColor',[0.929,0.694,0.125],...
    'MarkerFaceColor',[0.929,0.694,0.125])
legend(["DMSO","Gem 50uM","Gem 10uM","Tra 100nM","Tra 32nM"])
hold on
plot(time,hellinger_bm_time_gem_10uM,'-o','color',[0.829,0.694,0.125],'MarkerEdgeColor',[0.829,0.694,0.125],...
    'MarkerFaceColor',[0.829,0.694,0.125])
legend(["DMSO","Gem 50uM","Gem 10uM","Tra 100nM","Tra 32nM"])
hold on
plot(time,hellinger_bm_time_tra_100nM,'-o','color',[0.494,0.184,0.556],'MarkerEdgeColor',[0.494,0.184,0.556],...
    'MarkerFaceColor',[0.494,0.184,0.556])
legend(["DMSO","Gem 50uM","Gem 10uM","Tra 100nM","Tra 32nM"])
hold on
plot(time,hellinger_bm_time_tra_32nM,'-o','color',[0.594,0.184,0.556],'MarkerEdgeColor',[0.594,0.184,0.556],...
    'MarkerFaceColor',[0.594,0.184,0.556])
hold off
legend(["Gem 50uM","Gem 10uM","Tra 100nM","Tra 32nM"],'location','northeast')
xlabel('Hours post-treatment')
ylabel('Hellinger distance to DMSO')
title('All cells')


%%
AllPoints = [log10(sudo_FBM.bm),log10(sudo_FBM.pmt2),log10(sudo_FBM.pmt1)];
density = [];
R = 0.1; % deAPCnds on your data
for idx=1:length(AllPoints)
Distances = sqrt( sum( (AllPoints-AllPoints(idx,:)).^2 ,2) );
Ninside   = length( find(Distances<=R) );
density(idx) = Ninside/(4*pi*R.^3/3);
end
%%
figure('units','pixels','position',[0 0 1300 2000]) 
ind_keep = find(density<1*10^5);
scatter3(sudo_FBM.bm(ind_keep),sudo_FBM.pmt2(ind_keep),sudo_FBM.pmt1(ind_keep),7*ones(length(sudo_FBM.pmt1(ind_keep)),1),density(ind_keep),...
    'filled','MarkerFaceAlpha',.4)
% scatter3(sudo_FBM.bm(ind_keep),sudo_FBM.pmt2(ind_keep),sudo_FBM.pmt1(ind_keep),4)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'ZScale', 'log')
    cmap = crameri('bilbao');
    set(gca,'ColorScale','log')
    colormap(cmap);
    c=colorbar;
    caxis([200 10000]);
%     symlog('xyz',0.00001)
    ylabel("Inner-membrane flipping (Annexin V) mV")
    xlabel("Buoyant mass (pg)")
    zlabel("Cellular viability (Calcien AM) mV")
    title('Model NIBRX-1362 (52339 single cells)')
    xlim([2,1000])
    ylim([0.1,4000])
    zlim([0.1,10000])
    set(get(c,'title'),'string','Density','Rotation',0);

%%
figure('units','pixels','position',[0 0 1300 2000]) 
ind_keep = find(density<1*10^5);
scatter3(sudo_FBM.bm(ind_keep),sudo_FBM.pmt2(ind_keep),sudo_FBM.pmt1(ind_keep),7*ones(length(sudo_FBM.pmt1(ind_keep)),1),density(ind_keep),...
    'filled','MarkerFaceAlpha',.4)
% scatter3(sudo_FBM.bm(ind_keep),sudo_FBM.pmt2(ind_keep),sudo_FBM.pmt1(ind_keep),4)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'ZScale', 'log')
    cmap = crameri('bilbao');
    set(gca,'ColorScale','log')
    colormap(cmap);
    c=colorbar;
    caxis([200 10000]);
%     symlog('xyz',0.00001)
    ylabel("Inner-membrane flipping (Annexin V) mV")
    xlabel("Buoyant mass (pg)")
    zlabel("Cellular viability (Calcien AM) mV")
    title('Model NIBRX-1362 (52339 single cells)')
    xlim([2,1000])
    ylim([0.1,4000])
    zlim([0.1,10000])
    set(get(c,'title'),'string','Density','Rotation',0);
    
%%
figure('units','pixels','position',[0 0 1300 2000]) 
ind_keep = find(density<1*10^5);
plot_FBM = FBM_anno(ind_keep,:);
color = tab10(2);
temp_FBM = plot_FBM(plot_FBM.annotation ==  "not_live",:);
scatter3(temp_FBM.bm,temp_FBM.pmt2,temp_FBM.pmt1,7*ones(length(temp_FBM.pmt1),1),color(1,:),...
    'filled','MarkerFaceAlpha',.4);hold on
temp_FBM = plot_FBM(plot_FBM.annotation ==  "live",:);
scatter3(temp_FBM.bm,temp_FBM.pmt2,temp_FBM.pmt1,7*ones(length(temp_FBM.pmt1),1),color(2,:),...
    'filled','MarkerFaceAlpha',.4);hold off
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'ZScale', 'log')
%     cmap = crameri('bilbao');
%     set(gca,'ColorScale','log')
%     colormap(cmap);
%     c=colorbar;
%     caxis([200 10000]);
%     symlog('xyz',0.00001)
    ylabel("Inner-membrane flipping (Annexin V) mV")
    xlabel("Buoyant mass (pg)")
    zlabel("Cellular viability (Calcien AM) mV")
    title('Model NIBRX-1362 (52339 single cells)')
    xlim([2,1000])
    ylim([0.1,4000])
    zlim([0.1,10000])
    legend(["Non-live cells","Live cells"],'location','eastOutside')  
%% make movie for the staircase
 OptionZ.FrameRate=100;OptionZ.Duration=45;OptionZ.APCriodic=true;CaptureFigVid([0,90;-40,10;-110,5;-190,30;-270,0;-360,0], 'WellMadeVid',OptionZ)