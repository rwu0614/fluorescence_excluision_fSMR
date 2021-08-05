global datafull
global datasmr_good
S3_Merge;
noise = 15;
coeff_peak_bal = 0.1; %the peak imbalance in percentage relative to mean mass
coeff_node_bal = noise*2;
coeff_node_peak = 0.2; %node that is bigger than 20% is ditched.

idx_discard = find(abs((datasmr(:,6)-datasmr(:,8))./datasmr(:,3)) > coeff_peak_bal | abs(datasmr(:,9)-datasmr(:,10)) > coeff_node_bal | abs(datasmr(:,11)./datasmr(:,3))>coeff_node_peak);

datasmr_good = datasmr;
datasmr_good(idx_discard,:)=[];