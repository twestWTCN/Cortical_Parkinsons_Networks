function [] = spectralplots_groups(datapathr,powsubgrand,cohsubgrand,frq,titular,bandname)
grouperN = {'Left','Right','Both'};
grouperI = {1,2,1:2}
        cmap = linspecer(2);

for i = 1:3
figure('Name',['Group Results ' grouperN{i}])
subplot(1,3,1)
ON = [powsubgrand{1,1,grouperI{i},:}];
OFF = [powsubgrand{2,1,grouperI{i},:}];
ax(1) = boundedline(frq,mean(ON,2),std(ON,0,2),'cmap',cmap(1,:),'alpha','transparency',0.45); %std(ON,0,2)/sqrt(size(ON,2))
hold on
ax(2) = boundedline(frq,mean(OFF,2),std(OFF,0,2),'cmap',cmap(2,:),'alpha','transparency',0.45);
% set(gca,'xscale','log');
%set(gca,'yscale','log')
xlim([6 45]);ylim([0.005 0.05]);xlabel('Frequency (Hz)');ylabel('Normalized Power'); title(titular{1})
ax(1).LineWidth = 1.5; ax(2).LineWidth = 1.5;
legend(ax,{'ON','OFF'}); grid on

subplot(1,3,2)
ON = [powsubgrand{1,2,grouperI{i},:}];
OFF = [powsubgrand{2,2,grouperI{i},:}];
ax(1) = boundedline(frq,mean(ON,2),std(ON,0,2),'cmap',cmap(1,:),'alpha','transparency',0.45);
hold on
ax(2) = boundedline(frq,mean(OFF,2),std(OFF,0,2),'cmap',cmap(2,:),'alpha','transparency',0.45);
% set(gca,'xscale','log'); 
%set(gca,'yscale','log')
xlim([6 45]);ylim([0.0085 0.025]);xlabel('Frequency (Hz)');ylabel('Normalized Power'); title(titular{2})
ax(1).LineWidth = 1.5; ax(2).LineWidth = 1.5;
legend(ax,{'ON','OFF'}); grid on

subplot(1,3,3)
ON = [cohsubgrand{1,grouperI{i},:}];
OFF = [cohsubgrand{2,grouperI{i},:}];
ax(1) = boundedline(frq,mean(ON,2),std(ON,0,2),'cmap',cmap(1,:),'alpha','transparency',0.45);
hold on
ax(2) = boundedline(frq',mean(OFF,2),std(OFF,0,2),'cmap',cmap(2,:),'alpha','transparency',0.45);
xlim([6 45]);ylim([0 0.35]); xlabel('Frequency (Hz)');ylabel('WPLI'); title('STN/CTX WPLI')
ax(1).LineWidth = 1.5; ax(2).LineWidth = 1.5;
legend(ax,{'ON','OFF'})
% set(gca,'xscale','log')
grid on
set(gcf,'Position',[314         551        1208         307])
savefigure_v2([datapathr 'results\spectral\'],[bandname '_STN_Source_Power_analysis_GroupAverage_' grouperN{i}],[],[],[]);
close all
end