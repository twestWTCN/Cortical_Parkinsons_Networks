% function sigcohlist = findsubSTNCTXsigcoh
clear
I = i;
% clear;
close all;
datapathr = 'C:\Users\Tim\Documents\Work\Cortical_Networks\Data\';
% subname = {'JN','MC','SW','DF','JB','MW','DP','DS','JA'};
subname = {'JN','MC','SW','DF','JB','MW','DP','DS','JA'}
ipsicon = 'ipsi';

for sub = 1; %1:numel(subname)
    load([datapathr subname{sub} '\ftdata\ROI_analy\ROIvoxel_power_' ipsicon],'powsave','frqsave')
    load([datapathr subname{sub} '\ftdata\ROI_analy\ROIvoxel_coh_' ipsicon],'cohsave','frqsave')
    
    
    
end
