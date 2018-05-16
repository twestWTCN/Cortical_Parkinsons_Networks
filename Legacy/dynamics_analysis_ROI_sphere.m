close all; clear 
i = 1; subname = {'LM'};
ref_chan = 'STN_R23'; nr = 1;
datapathr = 'C:\Users\Tim\Documents\Work\Cortical_Networks\Data\';

load([datapathr subname{i} '\ftdata\ROI_source_granger_' ref_chan '_' num2str(nr)])
load([datapathr subname{i} '\ftdata\ROI_source_coh_' ref_chan '_' num2str(nr)])
load([datapathr subname{i} '\ftdata\virtual_sources_' num2str(nr) '_ROI'],'source_sens','vtime','ROI_list')

tvec = 0:round(max(vtime{1}),2):(round(max(vtime{1}),2)*(numel(vtime)-1));
cX = mean(cASym,1);
gX = mean(gASym,1);
figure
A(1) = plot(tvec,cX,'b-','LineWidth',3); hold on; 
A(2) = plot(tvec,gX,'r-','LineWidth',3);
xlabel('Time'); ylabel('Granger/coherence');
hold on
sigC = 1;
A(3:4) = plot(repmat(tvec,2,1)',repmat([-sigC*std(cX); sigC*std(cX)],1,size(tvec,2))','b--','LineWidth',1);
A(5:6) = plot(repmat(tvec,2,1)',repmat([-sigC*std(gX); sigC*std(gX)],1,size(tvec,2))','r--','LineWidth',1);
legend(A([1 3 2 5]),'signed iCoherence','iCoherence std','Granger','Granger std');
set(gcf,'Position',[37 181.5 918 550])

% Dwell Time
x = abs(cX)>sigC*std(cX); 
x = sprintf('%d',x);
t1=textscan(x,'%s','delimiter','0','multipleDelimsAsOne',1);
d = t1{:};
for k = 1:length(d)
      Z1(k) = length(d{k});
end

x = abs(gX)>sigC*std(gX); 
x = sprintf('%d',x);
t1=textscan(x,'%s','delimiter','0','multipleDelimsAsOne',1);
d = t1{:};
for k = 1:length(d)
      Z2(k) = length(d{k});
end
figure
h(1) = histogram(Z1,1:max([Z1 Z2]+1),'FaceColor','b','FaceAlpha',0.6);
hold on
h(2) = histogram(Z2,1:max([Z1 Z2]+1),'FaceColor','r','FaceAlpha',0.6);
xlabel(sprintf('number of bins (of %.1f seconds)',round(max(vtime{1}),2)));
ylabel('frequency');
legend(h,'icoherence','granger')
set(gcf,'Position',[982 182.5 605 550])
