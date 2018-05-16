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
    Phishuff = Phi(1,randperm(length(Phi(1,:))));
    phasescm(1,:) = abs(ifft(Amp(1,:).*exp(sqrt(-1)*Phishuff)));
    Phishuff = Phi(2,randperm(length(Phi(2,:))));
    phasescm(2,:) = abs(ifft(Amp(2,:).*exp(sqrt(-1)*Phishuff)));
    Rdata = [];
    Rdata.trial{1} = phasescm;
    Rdata.time = Xdata.time;
    [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Rdata,frq,R.PA.stn_lb_frq,R.PA.bwid(2),fsamp);
    SRPbank(:,N) = abs(dphi_12_dt);
    
    ampscm =[];
    Ampshuff = Amp(1,randperm(length(Phi(1,:))));
    ampscm(1,:) = abs(ifft(Ampshuff.*exp(sqrt(-1)*Phi(1,:))));
    Ampshuff = Amp(2,randperm(length(Phi(2,:))));
    ampscm(2,:) = abs(ifft(Ampshuff.*exp(sqrt(-1)*Phi(2,:))));
    Rdata = [];
    Rdata.trial{1} = ampscm;
    Rdata.time = Xdata.time;
    [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Rdata,frq,R.PA.stn_lb_frq,R.PA.bwid(2),fsamp);
    ampbank(:,:,N) = amp;
%     ppm.increment();
end
% ppm.delete();
SRPeps = prctile(SRPbank(:),R.PA.SRPeps_prctile)*2;
Ampeps = prctile(reshape(permute(ampbank,[1 3 2]),[],3),50,1);
SNReps = prctile(reshape(permute(ampbank,[1 3 2]),[],3),R.PA.SNReps_prctile,1);