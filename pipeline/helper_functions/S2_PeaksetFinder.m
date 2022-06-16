function [peaks] = S2_PeaksetFinder(xdata, ydata, offset_input, baseparams, analysismode)

%%% Original by Sungmin Son
%%% Extensively modified by Nikita Khlystov
%%% Latest edit 12/16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    S2 PART I: Find primary and secondary peaks   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% minpkht_thres = min(ydata) + 0.5*(max(ydata) - min(ydata));               

idx = find(abs(diff(ydata))<baseparams(5));                                                  

% ignore the flat points found over the anti-node
mf_ydata_thres=medfilt1(ydata(idx), baseparams(6));
idx_f=find(abs(ydata(idx)-mf_ydata_thres)< baseparams(7));
idx=idx(idx_f);
%ydata_thres(idx) = smooth(ydata(idx), 2);
% filter extreme outliers.
% in density measurement flat part forms in outliers

idx=[idx; length(ydata)];

% if no flat part is found
if length(idx)<2
    peaks = [];
    %disp('no flat part found');
    return
end

ydata_thres=interp1(xdata(idx), ydata(idx), xdata);

% locally smooth ydata again
ydata=smooth(ydata, 3);
offset_input = (max(ydata) - min(ydata))*0.5;
minpkht_thres=ydata_thres-offset_input;

if analysismode == 0
    hold off;
    figure(1);
    subplot(2,2,3); plot(xdata, ydata, '-');
    hold on;
    subplot(2,2,3); plot(xdata(idx), ydata(idx), '.g'); 
    subplot(2,2,3); plot(xdata, minpkht_thres, '-k');
    input('go');
end

repeatflag = 1;
while(repeatflag == 1)
    pkidx = find(ydata < minpkht_thres);                                                  % find all freq values that pass this threshold - these are the main peaks
    if ~isempty(pkidx)
        idx_ends = find(abs(diff(pkidx)) > 1);                                                  % find all breaks in y_diff indices; each segment is a peak of interest, with end on idx_end and start on (idx_end + 1)
        idx_ends = [0 idx_ends' length(pkidx)];                                                   % cap the beginning and end of idx with markers
        peaks = zeros(1,length(idx_ends) - 1);                                            % set empty matrix for loop efficiency
        for i = 1:length(idx_ends) - 1
            data_segment = ydata(pkidx(idx_ends(i)+1:idx_ends(i + 1)));                        % focus on the piece of y_diff between two consecutive idx_end markers
            segment_idx_max = min(find(data_segment == min(data_segment)));                      % find index of the maximum deviation, i.e. the apex within the piece
            global_idx_max = pkidx(idx_ends(i) + 1) + segment_idx_max - 1;
            peaks(i) = global_idx_max;                                                       % convert the apex index to the global index within y_diff      
        end
        
        if analysismode == 0
        figure(1);    
        hold off                                                                                    %%% the entire frequency data
        subplot(2,2,3); plot(xdata, ydata, '-');
        hold on
        subplot(2,2,3); plot(xdata(idx), ydata(idx), '.g');
        subplot(2,2,3); plot(xdata, minpkht_thres, '-k');
        subplot(2,2,3); plot(peaks, ydata(peaks), '.r');                                      % peak location is only estimated - due to smoothing, actual apex is a bit different; refinement is below
        set(gca,'XLim',[0 length(xdata)]);
        set(gca,'YLim',[min(ydata) - 10 max(ydata) + 10]);
        
        
            %input('go?');
        end
        
        goflag = 0;
        if mod(length(peaks),3) ~= 0 
            %disp('Warning: the number of peaks found is not a multiple of three. Check plot.')
%             peaks=[];
            checkthres = 0;%input('Enter 1 if you want to adjust threshold, enter 0 to skip:   ');
            if length(peaks)==1
                %disp('Skipping this peak set. There is only one peak..');
                peaks=[];
                repeatflag=0;
            end
            
            if checkthres == 1
               % fprintf('Previous threshold multiplier:   %3.0f \n', minpkht_thres);
                minpkht_thres = input('Input new threshold:   ');
            else  %%%%% add elseif statement to try anyways
%                 peaks = [];
                goflag = 1;
            end
            
        else
            goflag = 1;
        end

        if goflag == 1
            repeatflag = 0;
        end
        
    elseif isempty(pkidx)
%         figure(1);
%         hold off                                                                                   
%         subplot(2,2,3); plot(xdata, ydata, '-');
%         hold on
%         set(gca,'XLim',[0 length(xdata)]);
%         set(gca,'YLim',[min(ydata) - 10 max(ydata) + 10]);
        %disp('No peaks were found. Please adjust threshold baseline.')
%         fprintf('Previous threshold multiplier:    %3.0f \n', minpkht_thres)
%         minpkht_thres = input('Input new threshold multiplier:    ');
        peaks = [];
        repeatflag=0;
    end
    
end
return
end
