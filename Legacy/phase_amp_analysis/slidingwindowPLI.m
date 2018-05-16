function [PLI PLI_tvec] = slidingwindowPLI(period,phi,time)
[slide_dphi_12,sind] = slideWindow(phi(:,1)-phi(:,2), floor(period), floor(period*0.75));
% PLV = abs(mean(exp(1i*slide_dphi_12),1));
% PLV_tvec = (round(median(sind,1)));

PLI = abs(mean(sign(slide_dphi_12),1));  % PLI
PLI_tvec = time(round(median(sind,1)));