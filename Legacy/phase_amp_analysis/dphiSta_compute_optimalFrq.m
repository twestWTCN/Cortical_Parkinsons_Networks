function [maxfrq maxPLV] = dphiSta_compute_optimalFrq(Xdata,R,band)
frqlist = R.PA.frqrange{band};
fsamp = R.pp.cont.full.fs;
surdata = Xdata.trial{1};
surfft = fft(surdata);
Amp = abs(surfft);
Phi = angle(surfft);
for frqn = 1:numel(frqlist)
    [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Xdata,frqlist(frqn),R.PA.stn_lb_frq,R.PA.bwid(band),fsamp);
%     qstable = find(abs(dphi_12_dt')<0.005)
    dataPLV = mean(dphi_12_dt);
    
    phasescm =[];
    Phishuff = Phi(randperm(length(Phi(1,:))));
    phasescm(1,:) = abs(ifft(Amp(1,:).*exp(sqrt(-1)*Phishuff)));
    Phishuff = Phi(randperm(length(Phi(2,:))));
    phasescm(2,:) = abs(ifft(Amp(2,:).*exp(sqrt(-1)*Phishuff)));
    
    Rdata = [];
    Rdata.trial{1} = phasescm;
    [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Rdata,frqlist(frqn),R.PA.stn_lb_frq,R.PA.bwid(2),fsamp);
%     qstable = find(abs(dphi_12_dt')<0.005)
    surrPLV = mean(dphi_12_dt);
    diffPLV(frqn) = dataPLV-surrPLV;
end
maxfrq = frqlist(diffPLV==max(diffPLV));
maxPLV = max(diffPLV);
