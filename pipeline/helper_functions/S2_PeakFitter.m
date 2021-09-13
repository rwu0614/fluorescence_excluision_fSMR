function [pkidx_poly, pkht_poly, apkidx_poly, apkht_poly, baselineslope, htdiff_poly, ahtdiff_poly] = S2_PeakFitter(xdata, ydata, baseline, peaks, peakwidth, dispprogress)

%%% Original by Sungmin Son
%%% Modified by Nikita Khlystov
%%% Latest edit on 12/16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    S2 PART III: Imperically fit peaks   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

xbasedata = baseline;
ybasedata = ydata(baseline)';

antipeaks = zeros(1,length(peaks)-1);
for i = 1:(length(peaks)-1)
    antiy_max_idx=find(ydata == max(ydata(peaks(i):peaks(i+1))));
    antipeaks(i) = antiy_max_idx(1);
end

linbaseline = polyfit(xbasedata, ybasedata, 1);
baselineslope = linbaseline(1);
baselinefreq = (polyval(linbaseline, xdata))';
median(ydata);

freqdata = ydata' - baselinefreq; 

freqbasedata = ybasedata - baselinefreq(xbasedata);
% N=length(freqdata);
% % if N<=101
freqsmth = sgolayfilt(freqdata, 3, 11);
% % else
% %     freqsmth = sgolayfilt(freqdata, 3, 101);
% % end
if dispprogress ==1
   
% subplot(2,2,3); 
%  hold off; cla
% plot(xdata, ydata, '-'); hold on; plot(xdata, baselinefreq); 
% input('go');

figure(1);
subplot(2,2,3); hold off; plot(xdata, freqdata, '-')
hold on
subplot(2,2,3); plot(xbasedata, freqbasedata, '.g');
% subplot(2,2,3); plot(xdata, median(ydata)-baselinefreq, '--k');
subplot(2,2,3); plot(xdata, baselinefreq - baselinefreq, '-c');
set(gca,'XLim',[0 length(xdata)]);
set(gca,'YLim',[1.1*min(freqdata) max(freqdata)-0.1*min(freqdata)]);
end

hahtwd_ll = find(freqsmth <= freqsmth(peaks(1))/2, 1);
hahtwd_rr = find(freqsmth <= freqsmth(peaks(end))/2, 1, 'last');

if dispprogress ==1
subplot(2,2,3); plot(xdata(peaks), freqdata(peaks), '.r');                                        % peak minima are in red
subplot(2,2,3); plot(xdata(antipeaks), freqdata(antipeaks), '.r');                                % antipeak maxima are in red
subplot(2,2,3); plot(xdata(hahtwd_ll), freqdata(hahtwd_ll), '*g');                                % left half-peak height as a green star
subplot(2,2,3); plot(xdata(hahtwd_rr), freqdata(hahtwd_rr), '*g');                                % right half-peak height as a green star
end

pk_fit_wd = round(peakwidth/30);                                                                  % width of segment of peak to be fitted
%% output half-height indices as peakwidth
pkidx_poly = zeros(1,length(peaks));
apkidx_poly = zeros(1,length(antipeaks));
pkht_poly = zeros(1,length(peaks));
apkht_poly = zeros(1,length(antipeaks));
for i = 1:length(peaks)
    pk_fit_segx = (max((xdata(peaks(i)) - pk_fit_wd),1):min((xdata(peaks(i)) + pk_fit_wd),length(xdata)))';                 % actual segment of peak to be fitted
    pk_fit_segy = freqdata(pk_fit_segx);
    polypkeq = polyfit([1:length(pk_fit_segx)], pk_fit_segy, 4);                                  % perform quarternary polynomial fit
    peakfit = polyval(polypkeq, [1:length(pk_fit_segx)]);                                            % evaluate peak fit 
    [pky, pkx] = min(peakfit);                                                                    % find peak min value and location
    pkidx_poly(i) = pkx + xdata(peaks(i)) - pk_fit_wd - 1;                                        % adjust the position of peak apex
    pkht_poly(i) = -1*pky;                                                                        % define peak height as vertical distance from fitted peak minimum and 
    if dispprogress ==1

    subplot(2,2,3); plot(pk_fit_segx, pk_fit_segy, '*g') 
    subplot(2,2,3); plot(pk_fit_segx, peakfit, 'r')                 
    subplot(2,2,3); plot(pkidx_poly(i), pky, 'or')    
    end
end

for i = 1:length(antipeaks)
    antipk_fit_segx = (max((xdata(antipeaks(i)) - pk_fit_wd),1):min(length(xdata),(xdata(antipeaks(i)) + pk_fit_wd)))';     % actual segment of peak to be fitted
    antipk_fit_segy = freqdata(antipk_fit_segx);
    polyantipkeq = polyfit((-pk_fit_wd:pk_fit_wd), antipk_fit_segy, 4);                          % perform quarternary polynomial fit
    antipeakfit = polyval(polyantipkeq, -pk_fit_wd:pk_fit_wd);                                    % evaluate peak fit 
    [apky, apkx] = max(antipeakfit);                                                              % find peak min value and location
    apkidx_poly(i) = xdata(antipeaks(i)) - pk_fit_wd + apkx - 1;                                  % adjust the position of peak apex
    apkht_poly(i) = apky;      % define peak height as vertical distance from fitted peak minimum and 
%     
%     antipeakfit = freqsmth(antipk_fit_segx);
%     [apky, apkx] = max(antipeakfit);
    apkidx_poly(i) = xdata(antipeaks(i)) - pk_fit_wd + apkx - 1;                                  % adjust the position of peak apex
    apkht_poly(i) = apky; 

if dispprogress ==1
    subplot(2,2,3); plot(antipk_fit_segx, antipeakfit, 'r')          
    subplot(2,2,3); plot(apkidx_poly(i), apky, 'or')       
    end
end

% htdiff_poly = 100*abs(diff(pkht_poly([1 end])))/mean(pkht_poly([1 end]));
htdiff_poly = mean(apkht_poly([1 end]));
% ahtdiff_poly = 100*abs(diff(apkht_poly([1 end])))/mean(apkht_poly([1 end]));
ahtdiff_poly = hahtwd_rr - hahtwd_ll;   % FWHM

hold off

return

end