function Peak_shape = P2_peak_shape(t,Data,Peak,channelID_max_amp,plotrange,maxrange,disp_progress,plotspace,channel_name,input_info)
% P2_peak_shape.m serves as a subfunction for P1_peakanalysis_pmt, to
% specificly analyze the shape factors of a pmt peak. The goal is to pull
% out features from the individual pmt peak that may link to the morphology
% of the paricle being measured.

currentFolder = pwd;
%% Initialization

% define region of interest, choose between plotrange and maxrange
region = plotrange;

%============================== WARNING! ===============================%
% Set median and moving average filter length.
med_filt_length = 5;
moving_average_window_size = 5;
%==== WARNING: filters are not the inherited from main P1 function======%
    
% Apply Data filtering
Data.filtered_med = cellfun(@(x) medfilt1(x,med_filt_length),Data.normalized,...
        'UniformOutput',false);

%setting up parameter for moving average filter
b = (1/moving_average_window_size)*ones(1,moving_average_window_size); a=1; 

Data.filtered_med_ave = cellfun(@(x) filter(b,a,x), Data.filtered_med,...
    'UniformOutput',false);

% Get 1st,2nd and 3rd derivatives for non-filtered, median filtered, and
% median+moving-average filtered data 
Peak_shape.no_fil = Data.normalized{channelID_max_amp}(region);
Peak_shape.no_fil_diff1 = diff(Peak_shape.no_fil);
Peak_shape.no_fil_diff2 = diff(Peak_shape.no_fil_diff1);
Peak_shape.no_fil_diff3 = diff(Peak_shape.no_fil_diff2);

Peak_shape.med_fil = Data.filtered_med{channelID_max_amp}(region);
Peak_shape.med_fil_diff1 = diff(Peak_shape.med_fil);
Peak_shape.med_fil_diff2 = diff(Peak_shape.med_fil_diff1);
Peak_shape.med_fil_diff3 = diff(Peak_shape.med_fil_diff2);

Peak_shape.med_ave_fil =Data.filtered_med_ave{channelID_max_amp}(region);
Peak_shape.med_ave_fil_diff1 = diff(Peak_shape.med_ave_fil);
Peak_shape.med_ave_fil_diff2 = diff(Peak_shape.med_ave_fil_diff1);
Peak_shape.med_ave_fil_diff3 = diff(Peak_shape.med_ave_fil_diff2);

%% plot first, second and third derivative of the signal
if disp_progress == 1
    subplot(plotspace.n_col,plotspace.n_row,plotspace.peakshape(1));
        plot(region,Peak_shape.no_fil,'-b');hold on;
        plot(region,Peak_shape.med_fil,'-g');hold on;
        plot(region,Peak_shape.med_ave_fil,'-r');hold off;
        ylabel('PMT voltage (V)')
        title(append(channel_name(channelID_max_amp),' peak view'))
    subplot(plotspace.n_col,plotspace.n_row,plotspace.peakshape(2));
        plot(region(2:end),Peak_shape.no_fil_diff1,'-b');hold on;
        plot(region(2:end),Peak_shape.med_fil_diff1,'-g');hold on;
        plot(region(2:end),Peak_shape.med_ave_fil_diff1,'-r');hold off;
        ylabel('First derivative')
        title(append(channel_name(channelID_max_amp),' peak view'))
    subplot(plotspace.n_col,plotspace.n_row,plotspace.peakshape(3));
        plot(region(3:end),Peak_shape.no_fil_diff2,'-b');hold on;
        plot(region(3:end),Peak_shape.med_fil_diff2,'-g');hold on;
        plot(region(3:end),Peak_shape.med_ave_fil_diff2,'-r');hold off;
        ylabel('Second derivative')
        title(append(channel_name(channelID_max_amp),' peak view'))
    subplot(plotspace.n_col,plotspace.n_row,plotspace.peakshape(4));
        plot(region(4:end),Peak_shape.no_fil_diff3,'-b');hold on;
        plot(region(4:end),Peak_shape.med_fil_diff3,'-g');hold on;
        plot(region(4:end),Peak_shape.med_ave_fil_diff3,'-r');hold off;
        ylabel('Third derivative')
        title(append(channel_name(channelID_max_amp),' peak view'))
        axP = get(gca,'Position');
        legend(["Raw nomalized","Median filter",...
            "Median+moving average filter"]...
        ,'Location','SouthOutside');
        set(gca, 'Position', axP);
end

%% Output generation
%output matrix will be appeneded by rows to an existing matrix, with a
%spacer row in the following format:
%         [computer-realtime;
%         raw-normalized pmt(v);
%         med-ave filered(v);
%         med-ave filter 1st derivative;
%         med-ave filter 2nd derivative;
%         med-ave filter 3rd derivative]
output_mat =[t(region(4:end))';Peak_shape.no_fil(4:end)';Peak_shape.med_ave_fil(4:end)';...
    Peak_shape.med_ave_fil_diff1(3:end)';Peak_shape.med_ave_fil_diff2(2:end)';...
    Peak_shape.med_ave_fil_diff3'];

out_type_ID = input('Want to save as singlet[1], or non-singlet[2]? input [1/2]\n');
out_type = ["example_singlet_","example_non-singlet_"];

if ~isempty(out_type_ID)
    if out_type_ID == 1 ||out_type_ID == 2
        name_split = strsplit(input_info.time_dir,'\');
        sample_name = name_split{end-1};
        sample_name= strrep(sample_name,'_',' ');
        cd(input_info.pmt_dir(1))
        out_file_name = [out_type{out_type_ID} sample_name '.csv'];
        dlmwrite(out_file_name, output_mat,'-append', 'delimiter', ',', 'precision', 25,'roffset',1);
        cd(currentFolder)
    end
end

end











