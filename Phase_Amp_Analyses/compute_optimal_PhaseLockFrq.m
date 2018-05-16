function [maxfrq maxPLV] = compute_optimal_PhaseLockFrq(R,Xdata,band,stn_lb_frq)
frqlist = R.PA.frqrange{band};
for frqn = 1:numel(frqlist)
    [amp phi dphi_12 , ~, ~] = comp_instant_angle_phase(Xdata,frqlist(frqn),stn_lb_frq,R.PA.bwid,band);
    %Epoched
    WinSize = R.PA.slidingwindow*Xdata.fsample;
    [dphi_12,sind] = slideWindow(dphi_12, floor(WinSize*5), 0);
    dphi_12 = wrapToPi(dphi_12);
    PLV(frqn) = mean(abs(mean(exp(1i*dphi_12),1)));
    PLI(frqn) = mean(abs(mean(sign(dphi_12),1)));
end
crit = eval(R.PA.optimalPLFrqMeth);
[~,loc] = findpeaks(crit);
[~,i] = min(abs(loc - median(1:numel(frqlist)))); % Find the peak closest to the centre of the band
maxfrq = frqlist(i);
maxPLV = crit(i);
