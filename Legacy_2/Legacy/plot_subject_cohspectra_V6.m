function plot_subject_cohspectra_V6(R)
%%%
% This function first computes coherences for all vchans and then selects a
% v-chan/STN ref pair that has maximal WPLI. The power and coherences are
% plot on a subject level and the distribution of cortical locations are
% plot last. Barplots are also computed for differences in these measures
% and UPDRS correlations can also be done here.
% TO DO:
% SElect trial repeat with max WPLI
%%%
load([R.datapathr '\UPDRS_Scores.mat'])
for band = [1 3]
    for sub = 1:numel(R.subname)
        for side = 1:2
            for cond = 1:2
                load([R.datapathr R.subname{sub} '\ftdata\virtualV6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bandname{band}],'vc_clean')
                if cond == 1
                    vchansave_ON = vc_clean;
                else
                    vchansave_OFF = vc_clean;
                end
            end
            
            if (vchansave_ON.specanaly.flag || vchansave_OFF.specanaly.flag) == 0
                %                 locbank(:,1,side,sub) = vchansave(id(1)).loc; vchansave_ON = vchansave(id(1));
                %                 locbank(:,2,side,sub) = vchansave(id(2)).loc; vchansave_OFF = vchansave(id(2));
                
                [npdspctrm, Hz] = NPD_cortical_STN(vchansave_ON,vchansave_OFF,R,1,[2 2],8.5);
                savefigure_v2([R.datapathr R.subname{sub} '\images\spectral\'],['NPD_STN_Source_analysis_' R.subname{sub} '_' R.siden{side} '_'  R.bandname{band}],[],[],[]);
                close all
                for cond =1:2
                    npdspctrm_group{cond,side,1}(:,sub) = npdspctrm{cond,1,1}(:,1);
                    npdspctrm_group{cond,side,2}(:,sub) = npdspctrm{cond,1,2}(:,1);
                    npdspctrm_group{cond,side,3}(:,sub) = npdspctrm{cond,1,3}(:,1);
                    npdspctrm_group{cond,side,4}(:,sub) = npdspctrm{cond,1,4}(:,1);
                    x = npdspctrm{cond,1,4}(:,1);
                end
                figure('Name',[R.subname{sub} ' ' R.siden{side}])
                
                cmap = linspecer(2);
                for ch = 1:2
                    ON = vchansave_ON.specanaly.normpow;
                    OFF = vchansave_OFF.specanaly.normpow;
                    ON = squeeze(ON(ch,:,:)); if size(ON,1)<size(ON,2); ON = ON'; end
                    OFF = squeeze(OFF(ch,:,:));if size(OFF,1)<size(OFF,2); OFF = OFF'; end
                    
                    subplot(1,3,ch)
                    frq = vchansave_ON.specanaly.frq; powsubgrand{1,ch,side,sub} = ON(:,1);
                    ax(1) = boundedline(frq,mean(ON,2),std(ON,0,2)/sqrt(size(ON,2)),'cmap',cmap(1,:),'alpha','transparency',0.45);
                    set(gca,'xscale','log'); set(gca,'yscale','log')
                    xlim([6 45]);%ylim([0.005 0.05])
                    hold on
                    %                             [xCalc,yCalc] = linregress(log10(frq)',log10(ON));
                    %                             plot(xCalc(:,2),yCalc,'--','color',cmap(1,:),'linewidth',2);
                    frq =  vchansave_ON.specanaly.frq; powsubgrand{2,ch,side,sub} = OFF(:,1);
                    ax(2) = boundedline(frq,mean(OFF,2),std(OFF,0,2)/sqrt(size(OFF,2)),'cmap',cmap(2,:),'alpha','transparency',0.45);
                    set(gca,'xscale','log'); set(gca,'yscale','log')
                    xlim([6 45]);%ylim([0.005 0.05]);
                    xlabel('Frequency (Hz)');ylabel('log Normalized Power'); title(R.titular{ch})
                    legend(ax,{'ON','OFF'});
                    hold on
                    %                            [xCalc,yCalc] = linregress(frq',mean(OFF));
                    %                             plot(xCalc(:,2),yCalc,'--','color',cmap(2,:),'linewidth',2);
                    
                    grid on
                end
                clear ax
                
                ON =  vchansave_ON.specanaly.icoh'; cohsubgrand{1,side,sub} = ON(:,1);
                OFF = vchansave_OFF.specanaly.icoh'; cohsubgrand{2,side,sub} = OFF(:,1);
                
                subplot(1,3,3)
                frq = vchansave_ON.specanaly.frq'; frqsave =vchansave_ON.specanaly.frq';
                ax(1) = plot(frq,mean(ON,2),'color',cmap(1,:))
                %             ax(1) = boundedline(frq,mean(ON,2),0,'cmap',cmap(1,:),'alpha','transparency',0.45);
                xlim([4 45]);
                hold on
                %         [xCalc,yCalc] = linregress(frq',mean(ON,2));
                %         plot(xCalc(:,2),yCalc,'--','color',cmap(1,:),'linewidth',2);
                frq = vchansave_ON.specanaly.frq';
                ax(2) = plot(frq,mean(OFF,2),'color',cmap(2,:));
                %             ax(2) = boundedline(frq,mean(OFF,2),std(OFF,0,2)/sqrt(size(OFF,2)),'cmap',cmap(2,:),'alpha','transparency',0.45);
                xlim([6 45]); xlabel('Frequency (Hz)');ylabel('Coherence'); title('STN/CTX Coherence')
                legend(ax,{'ON','OFF'})
                hold on
                %         [xCalc,yCalc] = linregress(frq',mean(OFF,2));
                %         plot(xCalc(:,2),yCalc,'--','color',cmap(2,:),'linewidth',2);
                %         set(gca,'xscale','log')
                grid on
                set(gcf,'Position',[314         551        1208         307]); shg
            else
                ON = NaN; OFF = NaN; frq = NaN;
            end
            
            pind = find(strcmp(R.subname{sub},patient));
            %         if side == 1; hemiI = 2; elseif side == 2; hemiI = 1; end
            updrsSave{1,side,sub} = eval([R.siden{side}(1) '_akinesia_' R.condname{2} '(' num2str(pind) ')']);
            updrsSave{2,side,sub} = eval(['total_' R.condname{2} '(' num2str(pind) ')']);
            [maxcoh1 ind] = max(mean(ON(frq>R.bandef(band,1) & frq<R.bandef(band,2),:),1));
            highbetacoh{1,side,sub} = [maxcoh1 ind maxcoh1>0.05];
            [maxcoh2 ind] = max(mean(OFF(frq>R.bandef(band,1) & frq<R.bandef(band,2),:),1));
            highbetacoh{2,side,sub} = [maxcoh2 ind maxcoh2>0.05];
            highbetacoh{3,side,sub} = maxcoh1>0.05 && maxcoh2>0.05;
            %             refidSave{side,sub} = refid;
            savefigure_v2([R.datapathr R.subname{sub} '\images\spectral\'],['STN_Source_Power_analysis_' R.subname{sub} '_' R.siden{side}  '_'  R.bandname{band}],[],[],'-r200');
            close all
        end
    end
%     figure
%     load('source_template.mat')
%     A = min(source.pos); B = max(source.pos);
%     X = A(1):.75:B(1);Y = A(2):.75:B(2);Z = A(3):.75:B(3);
%     [Xz Yz Zz] = meshgrid(Y,X,Z);
%     mask =  reshape(source.avg.pow,source.dim)>0;
%     fv = isosurface(double(mask),0); %Yz,Xz,Zz,
%     patch(fv,'FaceColor',[.1 .1 .1],'EdgeColor',[0 0 0],'FaceAlpha',0.05);
%     hold on
%     xyz = squeeze(locbank(:,1,1,:));
%     clear a
%     a(1) = scatter3(xyz(2,:),xyz(1,:),xyz(3,:),100,'b','x','LineWidth',2); hold on
%     xyz = squeeze(locbank(:,1,2,:));
%     a(2) = scatter3(xyz(2,:),xyz(1,:),xyz(3,:),100,'b','o','LineWidth',2);
%     xyz = squeeze(locbank(:,2,1,:));
%     a(3) = scatter3(xyz(2,:),xyz(1,:),xyz(3,:),100,'r','x','LineWidth',2); hold on
%     xyz = squeeze(locbank(:,2,2,:));
%     a(4) = scatter3(xyz(2,:),xyz(1,:),xyz(3,:),100,'r','o','LineWidth',2);
%     legend(a,{'Left ON','Right ON','Left OFF','Right OFF'})
    
    % figure
    barplot_coh_groups(R.datapathr,highbetacoh)
    savefigure_v2([R.datapathr 'results\spectral\'],['STN_Source_Power_analysis_GroupAverage_boxplots_'  R.bandname{band}],[],[],'-r100'); close all
    
    spectralplots_groups(R.datapathr,powsubgrand,cohsubgrand,frqsave,R.titular,R.bandname{band})
    
    save([R.datapathr 'results\spectral\UPDRS_Connect_groupdata.m'],'updrsSave','highbetacoh')
    a = vertcat(updrsSave{2,:,:});
    b = vertcat(highbetacoh{2,:,:}); b = b(:,1);
    % a(b<0.075) = []; b(b<0.075) = [];
    [r p] = corrcoef(a,b)
    figure; scatter(a,b,75,'ro','filled');
    [yCalc ba Rsq] = linregress(a,b);
    hold on; plot(a,yCalc,'k-','LineWidth',2)
    xlabel('OFF Hemi Aknesia Score'); ylabel(['Maximum ' R.bandname{band} ' WPLI'])
    annotation(gcf,'textbox',...
        [0.16975 0.703571428571429 0.268642857142857 0.165476190476191],...
        'String',{sprintf('y = %0.3f + %0.3fx',ba),sprintf('R^2 = %0.2f',Rsq),sprintf('P = %0.3f',p(2))},...
        'LineStyle','none',...
        'FontWeight','bold',...
        'FontSize',12,...
        'FitBoxToText','off'); grid on
    savefigure_v2([R.datapathr 'results\spectral\'],['STN_WPLI_H_UPDRS_Corr_'  R.bandname{band}],[],[],'-r100'); close all
    
end
