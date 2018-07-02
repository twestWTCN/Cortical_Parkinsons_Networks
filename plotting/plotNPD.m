function [] = plotNPD(Hz,npdspctrm,R,ylimz,statflag)
if nargin<4
    ylimz = [0 0.1];
end
if nargin<5
    statflag= 0;
end
subplot(1,4,1);
ax(1) = boundedline(Hz,nanmean(npdspctrm{1,1,1},2),nanstd(npdspctrm{1,1,1},0,2)/sqrt(size(npdspctrm{1,1,1},2)),'cmap',R.condcmap(1,:),'alpha','transparency',0.4);
ax(2) = boundedline(Hz,nanmean(npdspctrm{2,1,1},2),nanstd(npdspctrm{2,1,1},0,2)/sqrt(size(npdspctrm{2,1,1},2)),'cmap',R.condcmap(2,:),'alpha','transparency',0.4);
title('Zero','FontSize',10); xlim([6 45]); ylim(ylimz);
xlabel('Frequency (Hz)','FontSize',10); ylabel('NPD','FontSize',10);

if statflag==1
    OFF = npdspctrm{1,1,1}'; ON = npdspctrm{2,1,1}';
    OFF = OFF(~any(isnan(OFF),2),:); ON = ON(~any(isnan(ON),2),:);
    [specstat] = spectralclusterstats230518(OFF,ON,Hz,'ll',100);
    clustat = plotFreqCluster(specstat)
end
legend(ax,R.condname); grid on

subplot(1,4,2);
ax(1) = boundedline(Hz,nanmean(npdspctrm{1,1,2},2),nanstd(npdspctrm{1,1,2},0,2)/sqrt(size(npdspctrm{1,1,2},2)),'cmap',R.condcmap(1,:),'alpha','transparency',0.8);
ax(2) = boundedline(Hz,nanmean(npdspctrm{2,1,2},2),nanstd(npdspctrm{2,1,2},0,2)/sqrt(size(npdspctrm{2,1,2},2)),'cmap',R.condcmap(2,:),'alpha','transparency',0.8);
title('Forward','FontSize',10); xlim([6 45]); ylim(ylimz); legend(ax,R.condname); grid on
xlabel('Frequency (Hz)','FontSize',10);ylabel('NPD','FontSize',10);

if statflag==1
    OFF = npdspctrm{1,1,2}'; ON = npdspctrm{2,1,2}';
    OFF = OFF(~any(isnan(OFF),2),:); ON = ON(~any(isnan(ON),2),:);
    [specstat] = spectralclusterstats230518(OFF,ON,Hz,'ll',100);
    clustat = plotFreqCluster(specstat)
end

subplot(1,4,3);
ax(1) = boundedline(Hz,nanmean(npdspctrm{1,1,3},2),nanstd(npdspctrm{1,1,3},0,2)/sqrt(size(npdspctrm{1,1,3},2)),'cmap',R.condcmap(1,:),'alpha','transparency',0.8);
ax(2) = boundedline(Hz,nanmean(npdspctrm{2,1,3},2),nanstd(npdspctrm{2,1,3},0,2)/sqrt(size(npdspctrm{2,1,3},2)),'cmap',R.condcmap(2,:),'alpha','transparency',0.8);
title('Reverse','FontSize',10); xlim([6 45]); ylim(ylimz); 
xlabel('Frequency (Hz)','FontSize',10);ylabel('NPD','FontSize',10);

if statflag==1
    OFF = npdspctrm{1,1,3}'; ON = npdspctrm{2,1,3}';
    OFF = OFF(~any(isnan(OFF),2),:); ON = ON(~any(isnan(ON),2),:);
    [specstat] = spectralclusterstats230518(OFF,ON,Hz,'ll',100);
    clustat = plotFreqCluster(specstat)
end
legend(ax,R.condname); grid on

subplot(1,4,4);
ax(1) = boundedline(Hz,nanmean(npdspctrm{1,1,4},2),nanstd(npdspctrm{1,1,4},0,2)/sqrt(size(npdspctrm{1,1,4},2)),'cmap',R.condcmap(1,:),'alpha','transparency',0.8);
ax(2) = boundedline(Hz,nanmean(npdspctrm{2,1,4},2),nanstd(npdspctrm{2,1,4},0,2)/sqrt(size(npdspctrm{2,1,4},2)),'cmap',R.condcmap(2,:),'alpha','transparency',0.8);
title('Sum Coherence','FontSize',10); xlim([6 45]); ylim([ylimz(1) ylimz(2)+0.1]);
xlabel('Frequency (Hz)','FontSize',10);ylabel('NPD','FontSize',10);
if statflag==1
    OFF = npdspctrm{1,1,4}'; ON = npdspctrm{2,1,4}';
    OFF = OFF(~any(isnan(OFF),2),:); ON = ON(~any(isnan(ON),2),:);
    [specstat] = spectralclusterstats230518(OFF,ON,Hz,'ll',500);
    clustat = plotFreqCluster(specstat)
end
legend(ax,R.condname); grid on
% % subplot(1,4,1);
% %
% % ax(1) = plot(Hz,nanmean(npdspctrm{1,1,1},2),'color',R.condcmap(1,:)); hold on
% % ax(2) = plot(Hz,nanmean(npdspctrm{2,1,1},2),'color',R.condcmap(2,:));
% % title('Zero'); xlim([6 45]); ylim([0 0.8]); legend(ax,R.condname); grid on; ylabel('NPD'); xlabel('Hz')
% %
% % subplot(1,4,2);
% % ax(1) = plot(Hz,nanmean(npdspctrm{1,1,2},2),'color',R.condcmap(1,:)); hold on
% % ax(2) = plot(Hz,nanmean(npdspctrm{2,1,2},2),'color',R.condcmap(2,:));
% % title('Forward'); xlim([6 45]); ylim([0 0.8]); legend(ax,R.condname); grid on; ylabel('NPD'); xlabel('Hz')
% %
% % subplot(1,4,3);
% % ax(1) = plot(Hz,nanmean(npdspctrm{1,1,3},2),'color',R.condcmap(1,:)); hold on
% % ax(2) = plot(Hz,nanmean(npdspctrm{2,1,3},2),'color',R.condcmap(2,:));
% % title('Reverse'); xlim([6 45]); ylim([0 0.8]); legend(ax,R.condname); grid on; ylabel('NPD'); xlabel('Hz')
% %
% % subplot(1,4,4);
% % ax(1) = plot(Hz,nanmean(npdspctrm{1,1,4},2),'color',R.condcmap(1,:)); hold on
% % ax(2) = plot(Hz,nanmean(npdspctrm{2,1,4},2),'color',R.condcmap(2,:));
% title('Sum Coherence'); xlim([6 45]); ylim(ylimz); legend(ax,R.condname); grid on; ylabel('NPD'); xlabel('Hz')
set(gcf,'Position',[592         594        1266         236])