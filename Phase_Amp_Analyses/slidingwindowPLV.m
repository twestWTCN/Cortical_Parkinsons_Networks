function [PLV PLV_tvec] = slidingwindowPLV(period,phi,overlap)
[slide_dphi_12,sind] = slideWindow(phi(:,1)-phi(:,2), floor(period), floor(period*overlap));
PLV = abs(mean(exp(1i*slide_dphi_12),1));

% PLV_tvec = (round(median(sind,1)));
% PLV = movmean(PLV,8);

% PLVi = abs(sum(exp(1i*slide_dphi_12),1));
% PLV_tvec = (round(median(sind,1)));

% PLI = abs(mean(sign(slide_dphi_12),1));  % PLI
% PLV= abs(mean(abs(slide_dphi_12).*sign(slide_dphi_12),1))./ mean(abs(slide_dphi_12),1);  %wPLI
PLV_tvec = (round(median(sind,1)));