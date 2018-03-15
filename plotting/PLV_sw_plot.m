function [] = PLV_sw_plot(Xdata,betaS,amp,phi,amp_sw,seg_ddt,PLV,PLV_tvec,consecSegs,R)
        ax(1) = subplot(4,1,1);
        cmap = linspecer(2);
        plot(Xdata.time{1},betaS.trial{1}(1,:),'color',cmap(1,:));hold on
        plot(Xdata.time{1},betaS.trial{1}(2,:),'color',cmap(2,:))
        plot(Xdata.time{1},amp(:,1),'color',cmap(1,:)); plot(Xdata.time{1},amp(:,2),'color',cmap(2,:));
        xmed1 = median(amp(:,1)); xmed2 = median(amp(:,2));
        plot(Xdata.time{1},repmat(xmed1,1,size(Xdata.time{1},2)),'--','color',cmap(1,:));
        plot(Xdata.time{1},repmat(xmed2,1,size(Xdata.time{1},2)),'--','color',cmap(2,:));
        ylabel('\beta activity','FontSize',14,'FontWeight','bold'); title('Example Analysis')

        ax(2) = subplot(4,1,2);
        plot(Xdata.time{1},phi(:,1),Xdata.time{1},phi(:,2),'linestyle','none','marker','.');%xlim([60 70])
        ylabel('\phi_{1/2}','FontSize',14,'FontWeight','bold');
        
        ax(3) = subplot(4,1,3);
        plot(PLV_tvec,amp_sw); %xlim([60 70]);
        tvec = nan(size(Xdata.time{1}(2:end))); tvec(seg_ddt) = Xdata.time{1}(seg_ddt);
        ylabel('Amplitudes (dB)','Interpreter','latex','FontSize',14,'FontWeight','bold')
        
        ax(4) = subplot(4,1,4);
%         PLV_tvec = PLV_tvec.*R.pp.cont.full.fs;
        plot(PLV_tvec,PLV); ylim([0 1]); % xlim([60 70]);
        tvec = nan(size(PLV_tvec)); tvec([consecSegs{:}]) = PLV_tvec([consecSegs{:}]);
        yvec = nan(size(PLV)); yvec([consecSegs{:}]) = PLV([consecSegs{:}]);
        hold on; plot(tvec,yvec,'LineWidth',2)
        ylabel('PLV/PLI','FontSize',14,'FontWeight','bold'); xlabel('Time (s)','FontSize',14,'FontWeight','bold'); legend({'PLV','PLI'})
        hold on; %plot([0 Xdata.time{1}(end)],[wpli_ci wpli_ci],'k--')
        linkaxes(ax,'x'); %xlim([80 82])
