function [] = spectralplots_groups(R,powsubgrand,cohsubgrand,frq,titular,statflag)
if nargin<6
    statflag = 0;
end
grouperN = {'Left','Right','Both'};
grouperI = {1,2,1:2};
        cmap = linspecer(2);

for i = 3 %1:3
figure('Name',['Group Results ' grouperN{i}])
set(gcf,'Position',[642   263   915   239])

subplot(1,3,1)
OFF = [powsubgrand{1,1,grouperI{i},:}];
ON = [powsubgrand{2,1,grouperI{i},:}];
ax(1) = boundedline(frq,mean(OFF,2),std(OFF,0,2)./sqrt(size(OFF,2)),'cmap',cmap(1,:),'alpha','transparency',0.8); %std(ON,0,2)/sqrt(size(ON,2))
hold on
ax(2) = boundedline(frq,mean(ON,2),std(ON,0,2)./sqrt(size(ON,2)),'cmap',cmap(2,:),'alpha','transparency',0.8);
% set(gca,'xscale','log');
% set(gca,'yscale','log')
xlim([6 45]);%ylim([0.005 0.05]);
xlabel('Frequency (Hz)','FontSize',10);ylabel('Normalized Power','FontSize',10); title(titular{1},'FontSize',10)
ax(1).LineWidth = 1.5; ax(2).LineWidth = 1.5;
legend(ax,R.condname); grid on

if statflag==1
    [specstat] = spectralclusterstats230518(OFF',ON',frq,'ll',100);
     clustat = plotFreqCluster(specstat)
end

subplot(1,3,2)
OFF = [powsubgrand{1,2,grouperI{i},:}];
ON = [powsubgrand{2,2,grouperI{i},:}];
ax(1) = boundedline(frq,mean(OFF,2),std(OFF,0,2)./sqrt(size(OFF,2)),'cmap',cmap(1,:),'alpha','transparency',0.8);
hold on
ax(2) = boundedline(frq,mean(ON,2),std(ON,0,2)./sqrt(size(ON,2)),'cmap',cmap(2,:),'alpha','transparency',0.8);
% set(gca,'xscale','log'); 
% set(gca,'yscale','log')
xlim([6 45]);%ylim([0.0085 0.05]);
xlabel('Frequency (Hz)','FontSize',10);ylabel('Normalized Power','FontSize',10); title(titular{2},'FontSize',10)
ax(1).LineWidth = 1.5; ax(2).LineWidth = 1.5;
legend(ax,R.condname); grid on

if statflag==1
    [specstat] = spectralclusterstats230518(OFF',ON',frq,'ll',500)
     clustat = plotFreqCluster(specstat)
end

subplot(1,3,3)
OFF = [cohsubgrand{1,grouperI{i},:}];
ON = [cohsubgrand{2,grouperI{i},:}];
ax(1) = boundedline(frq,mean(OFF,2),std(OFF,0,2)./sqrt(size(OFF,2)),'cmap',cmap(1,:),'alpha','transparency',0.8);
hold on
ax(2) = boundedline(frq',mean(ON,2),std(ON,0,2)./sqrt(size(ON,2)),'cmap',cmap(2,:),'alpha','transparency',0.8);
xlim([6 45]); %ylim([0 0.35]); 
xlabel('Frequency (Hz)','FontSize',10);ylabel('Coherence','FontSize',10); title([titular{1} '/' titular{2} ' coherence'],'FontSize',10)
ax(1).LineWidth = 1.5; ax(2).LineWidth = 1.5;
legend(ax,R.condname);grid on

if statflag==1
    [specstat] = spectralclusterstats230518(OFF',ON',frq,'ll',100);
     clustat = plotFreqCluster(specstat)
end

% set(gca,'xscale','log')
% close all
end