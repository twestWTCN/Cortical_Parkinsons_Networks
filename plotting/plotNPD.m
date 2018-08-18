function [] = plotNPD(Hz,npdspctrm,R,ylimz,statflag)
if nargin<4
    ylimz = [0 0.1];
end
if nargin<5
    statflag= 0;
end
titlist = {'Zero','SMA \rightarrow STN','STN \rightarrow SMA','Sum Coherence'};
figure('Name',['Group NPD Results'])
set(gcf,'Position',[642   263   915   239])
for i = 1:3
    subplot(1,3,i);
    for cond = 1:size(npdspctrm,1)
        X = [npdspctrm{cond,:,:,i}];
        ax(cond) = boundedline(Hz,nanmean(X,2),nanstd(X,0,2)/sqrt(size(X,2)),'cmap',R.condcmap(cond,:),'alpha','transparency',0.8);
    end
    ax(1).LineWidth = 1.5; ax(2).LineWidth = 1.5;
    
    title(titlist{i},'FontSize',10); xlim([6 45]); ylim(ylimz);
    xlabel('Frequency (Hz)','FontSize',10); ylabel('NPD','FontSize',10);
    
    if statflag==1
        OFF = [npdspctrm{1,:,:,i}]'; ON = [npdspctrm{2,:,:,1}]';
        OFF = OFF(~any(isnan(OFF),2),:); ON = ON(~any(isnan(ON),2),:);
        [specstat] = spectralclusterstats230518(OFF,ON,Hz,'ll',100);
        clustat = plotFreqCluster(specstat)
    end
    legend(ax,R.condname); grid on
end


%%% Script Grave
% % % subplot(1,4,1);
% % %
% % % ax(1) = plot(Hz,nanmean(npdspctrm{1,1,1},2),'color',R.condcmap(1,:)); hold on
% % % ax(2) = plot(Hz,nanmean(npdspctrm{2,1,1},2),'color',R.condcmap(2,:));
% % % title('Zero'); xlim([6 45]); ylim([0 0.8]); legend(ax,R.condname); grid on; ylabel('NPD'); xlabel('Hz')
% % %
% % % subplot(1,4,2);
% % % ax(1) = plot(Hz,nanmean(npdspctrm{1,1,2},2),'color',R.condcmap(1,:)); hold on
% % % ax(2) = plot(Hz,nanmean(npdspctrm{2,1,2},2),'color',R.condcmap(2,:));
% % % title('Forward'); xlim([6 45]); ylim([0 0.8]); legend(ax,R.condname); grid on; ylabel('NPD'); xlabel('Hz')
% % %
% % % subplot(1,4,3);
% % % ax(1) = plot(Hz,nanmean(npdspctrm{1,1,3},2),'color',R.condcmap(1,:)); hold on
% % % ax(2) = plot(Hz,nanmean(npdspctrm{2,1,3},2),'color',R.condcmap(2,:));
% % % title('Reverse'); xlim([6 45]); ylim([0 0.8]); legend(ax,R.condname); grid on; ylabel('NPD'); xlabel('Hz')
% % %
% % % subplot(1,4,4);
% % % ax(1) = plot(Hz,nanmean(npdspctrm{1,1,4},2),'color',R.condcmap(1,:)); hold on
% % % ax(2) = plot(Hz,nanmean(npdspctrm{2,1,4},2),'color',R.condcmap(2,:));
% % title('Sum Coherence'); xlim([6 45]); ylim(ylimz); legend(ax,R.condname); grid on; ylabel('NPD'); xlabel('Hz')
%