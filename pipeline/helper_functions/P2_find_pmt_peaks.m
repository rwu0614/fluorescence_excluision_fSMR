function [Peak, tempstart, tempend] = P2_find_pmt_peaks(peak,uni_peak_window_ext)

% finding the first index for peak analysis across all pmt channels
non_empty_ch = find(~cellfun(@isempty, peak.start)); %find channels that have non empty peak.start
minstart = min([peak.start{non_empty_ch}]); %find smallest index in peak.start in all channels

min_ch_search=false(length(non_empty_ch),1);
for i=1:length(non_empty_ch)
   min_ch_search(i) = ismember(minstart, peak.start{non_empty_ch(i)});
end

min_ch = non_empty_ch(min_ch_search(:)); %PMT channel with the smallest index for its first peak
other_ch = setxor(min_ch(1), non_empty_ch); % The other pmt channels that have unique peaks

if isempty(other_ch)==1 %if there is only one channel remaining to be analyzed
   tempstart = peak.start{min_ch(1)}(1)-uni_peak_window_ext;
   tempend = peak.end{min_ch(1)}(1)+uni_peak_window_ext;

else %there are other channel(s) remaining to be analyzed
    unionrange1to = cell(length(other_ch),1);
    unionrange = [];
    for i = 1:length(other_ch)
        unionrange1to{i} = union([peak.start{min_ch(1)}(1):peak.end{min_ch(1)}(1)], [peak.start{other_ch(i)}(1): peak.end{other_ch(i)}(1)]);
         if ~any(diff(unionrange1to{i})>1) %no breakage between min_ch1 and other_ch_i
              peak.start{other_ch(i)}(1)=[];
              peak.end{other_ch(i)}(1)=[];
         end
        unionrange = union(unionrange,unionrange1to{i});
    end
     
     if ~any(diff(unionrange)>1)
         tempstart = min(unionrange)-uni_peak_window_ext;
         tempend = max(unionrange)+uni_peak_window_ext;
     else
         tempstart = min(unionrange)-uni_peak_window_ext;
         tempend = unionrange(min(find(diff(unionrange)>1,1)));
     end
    
end
peak.start{min_ch(1)}(1)=[];
peak.end{min_ch(1)}(1)=[]; %delete peak region for next iteration
Peak = peak;
end
