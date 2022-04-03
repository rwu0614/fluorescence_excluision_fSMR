function [left_base, right_base, edgeidx] = S2_BaselineFinder(xdata, ydata, peaks, baseparams,peakdist, analysismode)

%%% Written by Nikita Khlystov
%%% Latest edit: 12/16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    S2 PART II: Find baselines    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%New way of selecting baseline motivated by small size
edgethres = baseparams(3)*abs(max(ydata)-min(ydata))/100; %Accounts for edge threshold to be proportional to cell size
winsize = baseparams(4);


%New way of accounting for transit time
edgethres = edgethres/(peakdist/200); %assuming 400 data points per transit gives peak_dist of 200
winsize = round(winsize*peakdist/200);

% locally smooth ydata again just for baseline finding...
sgolay_length_idx_base = 7; %this will make default length of 5 for 100 idx index
sgolay_length_base = 2*round(sgolay_length_idx_base*(peakdist/200))+1;
if sgolay_length_base >5
ydata= sgolayfilt(ydata, 3, sgolay_length_base);
else
ydata= sgolayfilt(ydata, 3, 5);
end

exitflag = 0;
while exitflag == 0
    leftedgesearch = max(1,(peaks(1) - winsize)):1:max(1,peaks(1)-round(peakdist/6));
    rightedgesearch =min(length(ydata), peaks(end)+round(peakdist/6)):1:min((peaks(end) + winsize) , length(ydata));
    
    if (isempty(find((diff(ydata(leftedgesearch))) < -edgethres)) || isempty(find((diff(ydata(rightedgesearch))) > edgethres)))
        disp('Skipping this peakset... no proper edge found')
        left_base = [];
        right_base = [];
        edgeidx = []; exitflag=1;
        return
    end
    
    idx_idx1 = find(diff(ydata(leftedgesearch)) > -edgethres, 1, 'last');
    idx_idx2 = find(diff(ydata(rightedgesearch)) < edgethres, 1);
    
    if (isempty(idx_idx1) || isempty(idx_idx2))
        left_base=[];
        right_base=[];
    else
    edgeidx(1) = leftedgesearch(idx_idx1);
    edgeidx(2) = rightedgesearch(idx_idx2);
    

    sidelength=round(peakdist); % change from 0.25*peakdist
    offset_length=20;   % change from 0
    
    left_base = max(edgeidx(1)-sidelength-offset_length,1):1:edgeidx(1)-offset_length;
    right_base = edgeidx(2)+offset_length:1:min(edgeidx(2)+sidelength+offset_length, length(ydata));
    
    if analysismode==0
        figure(1);
        freqdata = ydata - median(ydata(left_base));
        drawnow;
        hold off
        subplot(2,2,3); plot(xdata, freqdata, '-');
        hold on
        subplot(2,2,3); plot(peaks, freqdata(peaks), '.r');
        subplot(2,2,3); plot(edgeidx(1),ydata(edgeidx(1))- median(ydata(left_base)),'ok',edgeidx(2),ydata(edgeidx(2)) -median(ydata(left_base)),'ok')
        subplot(2,2,3); plot(left_base, freqdata(left_base), '.g'); plot(right_base, freqdata(right_base), '.g');
        subplot(2,2,3); plot(leftedgesearch, freqdata(leftedgesearch), '.c'); plot(rightedgesearch, freqdata(rightedgesearch), '.c');
        
        set(gca,'XLim',[min(left_base) max(right_base)]);
        set(gca,'YLim',[1.1*min(freqdata)  1.1*max(freqdata)]);
        
        
        input('go')
        
    end
    end
    
    if(isempty(left_base) || isempty(right_base))
        disp('Skipping this peakset...')
        left_base = [];
        right_base = [];
        edgeidx = [];
        return
        
        %             fprintf('No valid baselines could be found because selection criteria is too strict. Please adjust parameters. \n')
        %             fprintf('Previous stdevmultiplier:      %1.3f \n', stdevmultiplier)
        %             fprintf('Previous diffmultiplier:       %1.3f \n', diffmultiplier)
        %             fprintf('Previous edgethres:         %1.1f \n', edgethres)
        %             stdevmultiplier = input('Input new stdevmultiplier:     ');
        %             diffmultiplier = input('Input new diffmultiplier:      ');
        %             edgethres = input('Input new edgethres:        ');
        %             if stdevmultiplier == 0  || diffmultiplier == 0
        %                 disp('Skipping this peakset...')
        %                 left_base = [];
        %                 right_base = [];
        %                 edgeidx = [];
        %                 return
        %             end
    else
        exitflag = 1;
    end
    
end
return
end
