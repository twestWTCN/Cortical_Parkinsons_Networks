function plot_example_phaseanalysis_trace(Xdata,amp,phi,dphi_12_dt,seg_ddt,ddphi_ci) %,PLI_tvec,PLI,consecSegs)     
        ax(1) = subplot(4,1,1);
        cmap = linspecer(2);
        plot(Xdata.time{1},normaliseV(Xdata.trial{1}(1,:)),'color',cmap(1,:));hold on
        plot(Xdata.time{1},normaliseV(Xdata.trial{1}(2,:)),'color',cmap(2,:))
        plot(Xdata.time{1},normaliseV(amp(:,1)'),'color',cmap(1,:)); plot(Xdata.time{1},normaliseV(amp(:,2)'),'color',cmap(2,:));
%         xmed1 = median(amp(:,1)); xmed2 = median(amp(:,2));
%         plot(Xdata.time{1},repmat(xmed1,1,size(Xdata.time{1},2)),'--','color',cmap(1,:));
%         plot(Xdata.time{1},repmat(xmed2,1,size(Xdata.time{1},2)),'--','color',cmap(2,:));
        
        %xlim([60 70])
        ylabel('\beta activity','FontSize',14,'FontWeight','bold'); title('Example Analysis')


        ax(2) = subplot(4,1,2);
        plot(Xdata.time{1},phi(:,1),Xdata.time{1},phi(:,2),'linestyle','none','marker','.');%xlim([60 70])
        ylabel('\phi_{1/2}','FontSize',14,'FontWeight','bold');


        % ddt Time Series
        ax(3) = subplot(4,1,3);
        plot(Xdata.time{1}(2:end),dphi_12_dt); %xlim([60 70]);

        tvec = nan(size(Xdata.time{1}(2:end))); tvec(seg_ddt) = Xdata.time{1}(seg_ddt);
        yvec = nan(size(dphi_12_dt)); yvec(seg_ddt) = dphi_12_dt(seg_ddt);
        hold on; plot(tvec,yvec,'LineWidth',2)
        hold on; plot([0 Xdata.time{1}(end)],[ddphi_ci ddphi_ci],'k--');
        plot([0 Xdata.time{1}(end)],[-ddphi_ci -ddphi_ci],'k--'); ylim([-0.1 0.1])

        ylabel('$$ d \frac{(\phi_1 - \phi_2)}{dt} $$','Interpreter','latex','FontSize',14,'FontWeight','bold')

        %         % PA Time Series
        %         ax(3) = subplot(4,1,3);
        %         plot(Xdata.time{1}(:),wrapToPi(unwrap(phi1-phi2))); %xlim([60 70]);
        %         ylabel('\phi_1 - \phi_2')

        % PLI Time Series
%         ax(4) = subplot(4,1,4);
%         %             plot(PLV_tvec,PLV);hold on
%         PLV_tvec = PLV_tvec.*R.pp.cont.full.fs;
%         plot(PLV_tvec,PLV); ylim([0 1]); % xlim([60 70]);
% 
%         tvec = nan(size(PLV_tvec)); tvec([consecSegs{:}]) = PLV_tvec([consecSegs{:}]);
%         yvec = nan(size(PLV)); yvec([consecSegs{:}]) = PLV([consecSegs{:}]);
%         hold on; plot(tvec,yvec,'LineWidth',2)
%         ylabel('PLV/PLI','FontSize',14,'FontWeight','bold'); xlabel('Time (s)','FontSize',14,'FontWeight','bold'); legend({'PLV','PLI'})
%         hold on; %plot([0 Xdata.time{1}(end)],[wpli_ci wpli_ci],'k--')
%         %             hold on; plot([0 Xdata.time{1}(end)],[plv_ci plv_ci],'k--')
        linkaxes(ax,'x');
%         xlim([80 82])
        set(gcf,'Position',[287 72 1446 932]);
%         xlim([80 85])
