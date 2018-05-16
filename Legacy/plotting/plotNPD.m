function [] = plotNPD(Hz,npdspctrm,R)
cmap = linspecer(2);

% subplot(1,4,1);
% ax(1) = boundedline(Hz,mean(npdspctrm{1,1,1},2),std(npdspctrm{1,1,1},0,2)/sqrt(size(npdspctrm{1,1,1},2)),'cmap',R.condcmap(1,:),'alpha','transparency',0.45);
% ax(2) = boundedline(Hz,mean(npdspctrm{2,1,1},2),std(npdspctrm{2,1,1},0,2)/sqrt(size(npdspctrm{2,1,1},2)),'cmap',R.condcmap(2,:),'alpha','transparency',0.45);
% title('Zero'); xlim([2 48]); ylim([0 0.1]); legend(ax,R.condname); grid on
% 
% subplot(1,4,2);
% ax(1) = boundedline(Hz,mean(npdspctrm{1,1,2},2),std(npdspctrm{1,1,2},0,2)/sqrt(size(npdspctrm{1,1,2},2)),'cmap',R.condcmap(1,:),'alpha','transparency',0.45);
% ax(2) = boundedline(Hz,mean(npdspctrm{2,1,2},2),std(npdspctrm{2,1,2},0,2)/sqrt(size(npdspctrm{2,1,2},2)),'cmap',R.condcmap(2,:),'alpha','transparency',0.45);
% title('Forward'); xlim([2 48]); ylim([0 0.1]); legend(ax,R.condname); grid on
% 
% subplot(1,4,3);
% ax(1) = boundedline(Hz,mean(npdspctrm{1,1,3},2),std(npdspctrm{1,1,3},0,2)/sqrt(size(npdspctrm{1,1,3},2)),'cmap',R.condcmap(1,:),'alpha','transparency',0.45);
% ax(2) = boundedline(Hz,mean(npdspctrm{2,1,3},2),std(npdspctrm{2,1,3},0,2)/sqrt(size(npdspctrm{2,1,3},2)),'cmap',R.condcmap(2,:),'alpha','transparency',0.45);
% title('Reverse'); xlim([2 48]); ylim([0 0.1]); legend(ax,R.condname); grid on
% 
% subplot(1,4,4);
% ax(1) = boundedline(Hz,mean(npdspctrm{1,1,4},2),std(npdspctrm{1,1,4},0,2)/sqrt(size(npdspctrm{1,1,4},2)),'cmap',R.condcmap(1,:),'alpha','transparency',0.45);
% ax(2) = boundedline(Hz,mean(npdspctrm{2,1,4},2),std(npdspctrm{2,1,4},0,2)/sqrt(size(npdspctrm{2,1,4},2)),'cmap',R.condcmap(2,:),'alpha','transparency',0.45);
% title('Sum Coherence'); xlim([2 48]); ylim([0 0.1]); legend(ax,R.condname); grid on

subplot(1,4,1);

ax(1) = plot(Hz,mean(npdspctrm{1,1,1},2),'color',R.condcmap(1,:)); hold on
ax(2) = plot(Hz,mean(npdspctrm{2,1,1},2),'color',R.condcmap(2,:));
title('Zero'); xlim([6 45]); ylim([0 0.1]); legend(ax,R.condname); grid on; ylabel('NPD'); xlabel('Hz')

subplot(1,4,2);
ax(1) = plot(Hz,mean(npdspctrm{1,1,2},2),'color',R.condcmap(1,:)); hold on
ax(2) = plot(Hz,mean(npdspctrm{2,1,2},2),'color',R.condcmap(2,:));
title('Forward'); xlim([6 45]); ylim([0 0.1]); legend(ax,R.condname); grid on; ylabel('NPD'); xlabel('Hz')

subplot(1,4,3);
ax(1) = plot(Hz,mean(npdspctrm{1,1,3},2),'color',R.condcmap(1,:)); hold on
ax(2) = plot(Hz,mean(npdspctrm{2,1,3},2),'color',R.condcmap(2,:));
title('Reverse'); xlim([6 45]); ylim([0 0.1]); legend(ax,R.condname); grid on; ylabel('NPD'); xlabel('Hz')

subplot(1,4,4);
ax(1) = plot(Hz,mean(npdspctrm{1,1,4},2),'color',R.condcmap(1,:)); hold on
ax(2) = plot(Hz,mean(npdspctrm{2,1,4},2),'color',R.condcmap(2,:));
title('Sum Coherence'); xlim([6 45]); ylim([0 0.1]); legend(ax,R.condname); grid on; ylabel('NPD'); xlabel('Hz')
set(gcf,'Position',[400         698        1311         284])