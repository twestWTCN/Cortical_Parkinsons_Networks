clear; close all
load('Xdata.mat')
load('surdata.mat')
Odata= Xdata;
% [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Odata,37,16,0.5);
% qstable = find(abs(dphi_12_dt')<0.01); % The points that are below threshold
% consecSegs = SplitVec(qstable,'consecutive');
% % lengths
% segL_ddt = cellfun('length',consecSegs);
% histogram(log(segL_ddt/256),'FaceColor','b'); hold on



surfft1 = fft(surdata(1,:));
Amp1 = abs(surfft1);
Phi1 = angle(surfft1);

surfft2 = fft(surdata(2,:));
Amp2 = abs(surfft2);
Phi2 = angle(surfft2);


N = 6;
eps = logspace(-3,-2,N);
eps = 0.27; segL_eos = 1; WinN = 1.5;
WinN = [0.5 1 1.5 2]; %[0.05 0.1 0.25 0.5 2.5 3 5];
for i = 1:numel(WinN)
    for ha = 1:5
        figure
        [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Odata,37,16,1);
        WinSize =WinN(i)*1024;
        [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,0.95);
        SW_sampr = max(diff(PLV_tvec))/1024;
        figure(1)
        histogram(PLV,'FaceColor','b','Normalization','pdf'); hold on
        qstable = find(PLV>eps); % The points that are below threshold
        consecSegs = SplitVec(qstable,'consecutive'); PLVO = PLV;
        % lengths
        segL_ddt = cellfun('length',consecSegs);
        figure(2)
        histogram(segL_ddt/SW_sampr,'FaceColor','b','Normalization','pdf'); hold on
        figure(3)
        plot(PLV_tvec/1024,PLV,'b'); hold on; plot(PLV_tvec/1024,repmat(eps,1,numel(PLV_tvec)))
        
        
        Phishuff = Phi1(randperm(length(Phi1)));
        phasescm1 = randn(1,length(Phishuff)); %real(ifft(Amp1.*exp(sqrt(-1)*Phishuff)));
        Phishuff = Phi2(randperm(length(Phi2)));
        phasescm2 = randn(1,length(Phishuff)); %real(ifft(Amp2.*exp(sqrt(-1)*Phishuff)));
        
        Xdata.trial{1} = [phasescm1; phasescm2];
        
        [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Xdata,37,16,1);
        WinSize = WinN(i)*1024;
        [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,0.95);
        figure(1)
        histogram(PLV,'FaceColor','r','Normalization','pdf'); hold on
        
        
        qstable = find(PLV>eps); % The points that are below threshold
        consecSegs = SplitVec(qstable,'consecutive');
        segL_ddt = cellfun('length',consecSegs);
        figure(2)
        histogram(segL_ddt/1024,'FaceColor','r','Normalization','pdf'); hold on
        figure(3)
        plot(PLV_tvec/1024,PLV,'r'); hold on; plot(PLV_tvec/1024,repmat(eps,1,numel(PLV_tvec)))
        
        epsi(i,ha) = prctile(PLV,95);
        thresh(i,ha) = prctile(segL_ddt,95);
        [h pval(i,ha)] = ttest(PLV,PLVO);
        close all
    end
end
subplot(3,1,1)
plot(WinN,mean(pval,2)); xlabel('Window Size'); ylabel('PVal testing Xdata vs Surrogate');
subplot(3,1,2)
plot(WinN,mean(thresh,2)); xlabel('Window Size'); ylabel('Segment Length 95%');
subplot(3,1,3)
plot(WinN,mean(epsi,2)); xlabel('Window Size'); ylabel('Segment Length 95%');
savefigure_v2([R.datapathr '\results\seganalysis\PLI'],['PLI_surrogate'],[],[],[]); 
%close all
