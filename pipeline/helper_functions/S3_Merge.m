close all;
global datasmr
t=0;
m=0;
tm=[];
mm=[];
bm=[];
t1=[];
t2=[];
m1=[];
m2=[];
m3=[];
b1=[];
b2=[];
w=[];
ndm=[];
nd1=[];
nd2=[];
bs=[];
sectnum=[];
bd=[];
pkorder=[];
vs=[];

% filter out zero columns
idx_nonzero=find(datafull(12,:)>0);
datafull_processed=datafull(:,idx_nonzero);

flag=0;idx=0;
idx=find(diff(datafull_processed(12,:))~=0);
idx=[0 idx];

for i=1:length(idx)
    if i==length(idx)
        temp_idx=idx(i)+1:length(datafull_processed);
    else
    temp_idx=idx(i)+1:idx(i+1);
    end
    
    t1(i)=datafull_processed(1, temp_idx(1));
    t2(i)=datafull_processed(1, temp_idx(end));
    m1(i)=datafull_processed(2, temp_idx(1)); %left peak height
    m2(i)=datafull_processed(2, temp_idx(2));  %middle peak height (path dependent)
    m3(i)=datafull_processed(2, temp_idx(end)); %right peak height
%     m3(i)=mean([datafull(8, temp_idx(1)) datafull(8, temp_idx(3))]);
    b1(i)=mean([datafull_processed(4, temp_idx(1)) datafull_processed(5, temp_idx(1))]);
    b2(i)=mean([datafull_processed(4, temp_idx(end)) datafull_processed(5, temp_idx(end))]);
    tm(i)=mean([t1(i), t2(i)]);
    mm(i)=mean([m1(i), m3(i)]);
    bm(i)=mean([b1(i), b2(i)]);
    w(i)=datafull_processed(9, temp_idx(1)); %FWHM
    nd1(i)=datafull_processed(8, temp_idx(1));
    nd2(i)=datafull_processed(8, temp_idx(2));
    ndm(i)=mean([nd1(i), nd2(i)]);
    bs(i)=datafull_processed(7, temp_idx(1));
    vs(i)=datafull(13, temp_idx(2));
    sectnum(i)=datafull_processed(10, temp_idx(1)); %sectnum 
    bd(i)=datafull_processed(6, temp_idx(1)); % added by JK (transit time in ms)
    pkorder(i)=datafull_processed(12,temp_idx(1));
    
%     if(abs(datafull(7,temp_idx(1)))>0.001)
%     if(abs(datafull(6,temp_idx(1)))>800)
%     if(abs(datafull(7, temp_idx(1)))<0.0015)
%         flag(i)=0;
%     else
%         flag(i)=1;
%     end
%     plot(t1(i), m1(i), '.c');
%     plot(t2(i), m2(i), '.r');
%     hold on;
%      input('?');
end


datasmr=[tm' tm'/60 mm' bm' bs' m1' m2' m3' nd1' nd2' ndm' w' bd' vs' sectnum' tm'/3600 mm'/2 pkorder' ndm'./mm'];
% idx=find(abs(data(:,10))<10);
% data=data(idx,:);
datafull_processed = datafull_processed';

