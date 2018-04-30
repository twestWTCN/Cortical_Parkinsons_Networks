function [SRPeps Ampeps SNReps] = phase_amp_surrComp(Xdata,R,frq,NR)
fsamp = R.pp.cont.full.fs;
surdata = Xdata.trial{1};
surfft(:,1) = fft(surdata(1,:));
surfft(:,2) = fft(surdata(2,:));
Amp = abs(surfft)';
Phi = angle(surfft)';
% progressStepSize = 5;
% gcp;
% ppm = ParforProgMon('Power Estimation: ', NR, progressStepSize, 800, 300);


parfor N = 1:NR
    phasescm =[];
    Phishuff = Phi(randperm(length(Phi(1,:))));
    phasescm(1,:) = abs(ifft(Amp(1,:).*exp(sqrt(-1)*Phishuff)));
    Phishuff = Phi(randperm(length(Phi(2,:))));
    phasescm(2,:) = abs(ifft(Amp(2,:).*exp(sqrt(-1)*Phishuff)));
    
    Rdata = [];
    Rdata.trial{1} = phasescm;
    Rdata.time = Xdata.time;
    [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Rdata,frq,R.PA.stn_lb_frq,R.PA.bwid(2),fsamp);
    ampbank(:,N) = median(amp);
    SRPbank(:,N) = median(abs(dphi_12_dt));
%     ppm.increment();
end
% ppm.delete();
SRPeps = prctile(SRPbank,5);
Ampeps = prctile(ampbank,50,2);
SNReps = prctile(ampbank,5,2);