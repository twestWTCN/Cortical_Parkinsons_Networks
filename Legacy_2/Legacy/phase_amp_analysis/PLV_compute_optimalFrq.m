function [maxfrq maxPLV] = PLV_compute_optimalFrq(Xdata,R,band)
frqlist = R.PA.frqrange{band};
fsamp = R.pp.cont.full.fs;
surdata = Xdata.trial{1};
surfft = fft(surdata);
Amp = abs(surfft);
Phi = angle(surfft);
frqlist = 5:.5:49
for frqn = 1:numel(frqlist)
    [amp phi dphi_12 , ~, betaS] = comp_instant_angle_phase(Xdata,frqlist(frqn),R.PA.stn_lb_frq,R.PA.bwid(band),fsamp);
    %Epoched
    WinSize = R.PA.slidingwindow*fsamp;
    [dphi_12,sind] = slideWindow(dphi_12, floor(WinSize*5), 0);
    dphi_12 = wrapToPi(dphi_12);
    PLV(frqn) = mean(abs(mean(exp(1i*dphi_12),1)));
    PLI(frqn) = mean(abs(mean(sign(dphi_12),1)));
    % SLiding window
%     WinSize = R.PA.slidingwindow*fsamp;
%     [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,R.PA.WinOver);
%     dataPLV = mean(PLV);
%     
%     phasescm =[];
%     Phishuff = Phi(randperm(length(Phi(1,:))));
%     phasescm(1,:) = abs(ifft(Amp(1,:).*exp(sqrt(-1)*Phishuff)));
%     Phishuff = Phi(randperm(length(Phi(2,:))));
%     phasescm(2,:) = abs(ifft(Amp(2,:).*exp(sqrt(-1)*Phishuff)));
%     
%     Rdata = [];
%     Rdata.trial{1} = phasescm;
%     Rdata.time = Xdata.time;
%     [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Rdata,frqlist(frqn),R.PA.stn_lb_frq,R.PA.bwid(band),fsamp);
%     WinSize =R.PA.slidingwindow*fsamp;
%     [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,R.PA.WinOver);
%     surrPLV = mean(PLV);
%     diffPLV(frqn) = dataPLV-surrPLV;
end
maxfrq = frqlist(diffPLV==max(diffPLV));
maxPLV = max(diffPLV);
