clear; close all
R = makeHeader_SubCort_Cort_Networks_SFLAP();
p = 0;
cmap(:,:,1) = [0 0 1; 0 0.1 0.8];
cmap(:,:,3) = [1 0 0; 0.8 0 0.1];
for cond = 1%:2
    
    for sub = 1;%:numel(R.subname)
        for side = 1:2
            for band = [1 3]; %1:numel(R.bandname)
                p = p+1;
                load([R.datapathr R.subname{sub} '\ftdata\networks\threshnetwork_' R.ipsicon '_' R.siden{side} '_' R.bandname{band} '_' R.condname{cond}],'threshnetwork')
                %                 a(1) =subplot(3,1,1)
                %                 plot(threshnetwork.Xdata.time{1},log10(threshnetwork.amp(:,1))); hold on
                %                 a(2)= subplot(3,1,2)
                %                 plot(threshnetwork.Xdata.time{1},log10(threshnetwork.amp(:,2))); hold on
                
                a(1) =subplot(3,1,2)
                b1(p) = plot(threshnetwork.PLV_tvec,threshnetwork.PLV,'color',cmap(side,:,band)); hold on
                y1(:,band,side) = threshnetwork.PLV;
                a(2) =subplot(3,1,1)
                tvec = linspace(0,size(threshnetwork.amp,1)/1024,size(threshnetwork.amp,1));
                y =normaliseV(threshnetwork.amp(:,1)');
                y = y + abs(min(y));
                y = movmean(y,2*1024);
                b2(p) = plot(tvec,y,'color',cmap(side,:,band)); hold on
                y2(:,band,side) = y;
                %                 linkaxes(a,'x'); %xlim([80 82])
                
                %                 threshnetwork.PLV = PLV;
                %                 threshnetwork.PLV_tvec = PLV_tvec;
                %                 threshnetwork.stn_lb_frq = stn_lb_frq;
                %                 threshnetwork.maxfrq = maxfrq;
                %                 threshnetwork.betaS = betaS;
                %                 threshnetwork.amp = amp;
                %                 threshnetwork.Xdata = Xdata;
                %                 threshnetwork.consecSegs = consecSegs;
                %                 mkdir([R.datapathr R.subname{sub} '\ftdata\networks\'])
                %                 save([R.datapathr R.subname{sub} '\ftdata\networks\threshnetwork_' R.ipsicon '_' R.siden{side} '_' R.bandname{band} '_' R.condname{cond}],'threshnetwork')
            end
        end
    end
end



% Amp Env Corrcoef
winLength =2;
x1 = y2{1}; x2 = y2{3};
sfrq = 1/median(diff(tvec));
[slidex1,~] = slideWindow(x1, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
[slidex2,sind] = slideWindow(x2, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
for i = 1:size(slidex1,2)-1
    sR = corrcoef(slidex1(:,i),slidex2(:,i));
    R_alpha_beta(i) = sR(2);
    tvecAB(i) = tvec(median(sind(:,i)));
end
a(3) =subplot(3,1,3)
b3(1) = plot(tvecAB,R_alpha_beta,'b'); hold on

clear R_alpha_beta tvecAB
for rep = 1:100
    shuffX1 =   x1(randperm(length(x1)));
    shuffX2 =   x2(randperm(length(x2)));
    [slidex1,~] = slideWindow(x1, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
    [slidex2,sind] = slideWindow(x2, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
    
    for i = 1:size(slidex1,2)-1
        sR = corrcoef(slidex1(:,i),slidex2(:,i));
        R_alpha_beta(i,rep) = sR(2);
        tvecAB(i,rep) = tvec(min(sind(:,i)));
    end
    rep
end
upperBnd = prctile(R_alpha_beta(:),95); upperBnd2 = prctile(R_alpha_beta(:),99);

lowerBnd = prctile(R_alpha_beta(:),5); lowerBnd2 = prctile(R_alpha_beta(:),1);
b3(2) = plot(tvecAB(:,1),repmat(upperBnd,1,size(tvecAB,1)),'b--'); hold on
b3(3) = plot(tvecAB(:,1),repmat(lowerBnd,1,size(tvecAB,1)),'b--');

clear R_alpha_beta tvecAB
%% PLV Corrcoef
winLength =2;
x1 = y1{1}; x2 = y1{3};
sfrq = 1/median(diff(threshnetwork.PLV_tvec));
[slidex1,~] = slideWindow(x1, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
[slidex2,sind] = slideWindow(x2, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
for i = 1:size(slidex1,2)-1
    sR = corrcoef(slidex1(:,i),slidex2(:,i));
    R_alpha_beta(i) = sR(2);
    tvecAB(i) = threshnetwork.PLV_tvec(median(sind(:,i)));
end
a(3) =subplot(3,1,3)
b3(4) = plot(tvecAB,R_alpha_beta,'k'); hold on

clear R_alpha_beta tvecAB
for rep = 1:100
    shuffX1 =   x1(randperm(length(x1)));
    shuffX2 =   x2(randperm(length(x2)));
    [slidex1,~] = slideWindow(x1, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
    [slidex2,sind] = slideWindow(x2, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
    
    for i = 1:size(slidex1,2)-1
        sR = corrcoef(slidex1(:,i),slidex2(:,i));
        R_alpha_beta(i,rep) = sR(2);
        tvecAB(i,rep) = threshnetwork.PLV_tvec(min(sind(:,i)));
    end
    rep
end
upperBnd = prctile(R_alpha_beta(:),95); upperBnd2 = prctile(R_alpha_beta(:),99);

lowerBnd = prctile(R_alpha_beta(:),5); lowerBnd2 = prctile(R_alpha_beta(:),1);
b3(5) = plot(tvecAB(:,1),repmat(upperBnd,1,size(tvecAB,1)),'k--'); hold on
b3(6) = plot(tvecAB(:,1),repmat(lowerBnd,1,size(tvecAB,1)),'k--');
% b2(4) = plot(tvecAB(:,1),repmat(upperBnd2,1,size(tvecAB,1)),'k--'); hold on
% b2(5) = plot(tvecAB(:,1),repmat(lowerBnd2,1,size(tvecAB,1)),'k--');


linkaxes(a,'x'); %xlim([80 82])

subplot(3,1,2)
xlabel('Frame Time (s)'); ylabel('Within-Frame PLV');
legend(b1,{'Alpha STG Source','Beta Premotor Source'})
subplot(3,1,1)
xlabel('Frame Time (s)'); ylabel('Amplitude Envelope');
legend(b2,{'Alpha STN/STG Network','PLV STN/Premotor Network'})
subplot(3,1,3)
xlabel('Frame Time (s)'); ylabel('Correlation Coeff.');
legend(b3,{'Amp Env Correlation R','95th Percentile','5th Percentile','PLV Correlation R','95th Percentile','5th Percentile'}); %,'99th Percentile','1st Percentile'})

set(gcf,'Position',[4141 -133 1302 821])
