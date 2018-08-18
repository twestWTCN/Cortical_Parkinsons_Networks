function plot_subject_cohspectra_V6(R)
if nargin<1
    R = makeHeader_SubCort_Cort_Networks();
end
load([R.datapathr 'UPDRS_Scores'])
for breg = 2:length(R.bregname)
    for sub = 13:length(R.subname)
        for side = 1:2
            for cond = 1:length(R.condname)
                load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                if strncmp('OFF',R.condname{cond},2)
                    vchansave_OFF = vc_clean;
                else
                    vchansave_ON = vc_clean;
                end
            end
            
            if (vchansave_ON.specanaly.flag || vchansave_OFF.specanaly.flag) == 0
                %                 locbank(:,1,side,sub) = vchansave(id(1)).loc; vchansave_ON = vchansave(id(1));
                %                 locbank(:,2,side,sub) = vchansave(id(2)).loc; vchansave_OFF = vchansave(id(2));
                
                [npdspctrm, Hz] = NPD_cortical_STN(vchansave_OFF,vchansave_ON,R,1,[2 2],log2(R.pp.cont.thin.fs));
                set(gcf,'Position',[197         539        1633         310]); shg
                %                 mkdir([R.datapathr R.subname{sub} '\images\spectral\'])
                %                 savefigure_v2([R.datapathr R.subname{sub} '\images\spectral\'],['NPD_STN_Source_analysis_' R.subname{sub} '_' R.siden{side} '_'  R.bregname{breg}],[1],[],[]);
                %                 close all
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
                    %% ON
                    frq = vchansave_ON.specanaly.frq; powsubgrand{2,ch,side,sub} = ON(:,1);
                    ax(2) = plot((frq),(mean(ON,2)),'color',cmap(2,:));
                    %                     %                     ax(1) = boundedline(frq,mean(ON,2),std(ON,0,2)/sqrt(size(ON,2)),'cmap',cmap(1,:),'alpha','transparency',0.45);
                    %                     %                     set(gca,'xscale','log'); set(gca,'yscale','log')
                    xlim(([6 45]));%ylim([0.005 0.05])
                    hold on
                    %                     %                             [xCalc,yCalc] = linregress(log10(frq)',log10(ON));
                    %                     %                             plot(xCalc(:,2),yCalc,'--','color',cmap(1,:),'linewidth',2);
                    %% OFF
                    frq =  vchansave_OFF.specanaly.frq; powsubgrand{1,ch,side,sub} = OFF(:,1);
                    ax(1) = plot((frq),(mean(OFF,2)),'color',cmap(1,:));
                    %                     ax(2) = boundedline(frq,mean(OFF,2),std(OFF,0,2)/sqrt(size(OFF,2)),'cmap',cmap(2,:),'alpha','transparency',0.45);
                    %                     set(gca,'xscale','log'); set(gca,'yscale','log')
                    xlim(([6 45]));%ylim([0.005 0.05]);
                    if ch == 1; tit = [R.siden{side} ' ' R.bregname{breg}]; else tit = [R.siden{side} ' STN']; end
                    xlabel('Frequency (Hz)');ylabel('log Normalized Power'); title(tit)
                    legend(ax,R.condname);
                    hold on
                    %                            [xCalc,yCalc] = linregress(frq',mean(OFF));
                    %                             plot(xCalc(:,2),yCalc,'--','color',cmap(2,:),'linewidth',2);
                    grid on
                end
                clear ax
                
                OFF = vchansave_OFF.specanaly.coh'; cohsubgrand{1,side,sub} = OFF(:,1);
                ON =  vchansave_ON.specanaly.coh'; cohsubgrand{2,side,sub} = ON(:,1);
                subplot(1,3,3)
                %% Plot ON
                frq = vchansave_ON.specanaly.frq'; frqsave =vchansave_ON.specanaly.frq';
                ax(2) = plot(frq,mean(ON,2),'color',cmap(2,:));
                %             ax(1) = boundedline(frq,mean(ON,2),0,'cmap',cmap(1,:),'alpha','transparency',0.45);
                xlim([4 45]);
                hold on
                %         [xCalc,yCalc] = linregress(frq',mean(ON,2));
                %         plot(xCalc(:,2),yCalc,'--','color',cmap(1,:),'linewidth',2);
                %% Plot OFF
                frq = vchansave_ON.specanaly.frq';
                ax(1) = plot(frq,mean(OFF,2),'color',cmap(1,:));
                %             ax(2) = boundedline(frq,mean(OFF,2),std(OFF,0,2)/sqrt(size(OFF,2)),'cmap',cmap(2,:),'alpha','transparency',0.45);
                xlim([6 45]); xlabel('Frequency (Hz)');ylabel('coherence'); title(['STN/' R.bregname{breg} ' coherence'])
                legend(ax,R.condname)
                hold on
                %         [xCalc,yCalc] = linregress(frq',mean(OFF,2));
                %         plot(xCalc(:,2),yCalc,'--','color',cmap(2,:),'linewidth',2);
                %         set(gca,'xscale','log')
                grid on
                set(gcf,'Position',[314         551        1208         307]); shg
            else
                ON = NaN; OFF = NaN; frq = NaN;
                for cond =1:2
                    npdspctrm_group{cond,side,1}(:,sub) = nan(size(npdspctrm{cond,1,1}(:,1)));
                    npdspctrm_group{cond,side,2}(:,sub) = nan(size(npdspctrm{cond,1,2}(:,1)));
                    npdspctrm_group{cond,side,3}(:,sub) = nan(size(npdspctrm{cond,1,3}(:,1)));
                    npdspctrm_group{cond,side,4}(:,sub) = nan(size(npdspctrm{cond,1,4}(:,1)));
                end
            end
%             % If testing NPD
%             OFF = npdspctrm_group{1,side,3}(:,sub);  ON = npdspctrm_group{2,side,3}(:,sub);
%             frq = Hz;

            pind = find(strcmp(R.subname{sub},patient));
            if side == 1; hemiI = 2; elseif side == 2; hemiI = 1; end
            updrsSave{1,side,sub} = eval([R.siden{hemiI}(1) '_akinesia_' R.condname{1} '(' num2str(pind) ')']);
            updrsSave{2,side,sub} = eval(['total_' R.condname{1} '(' num2str(pind) ')']);
            updrsSave{3,side,sub} = eval(['total_' R.condname{2} '(' num2str(pind) ')'])-eval(['total_' R.condname{1} '(' num2str(pind) ')']);
            updrsSave{4,side,sub} = eval([R.siden{hemiI}(1) '_akinesia_' R.condname{2} '(' num2str(pind) ')'])-eval([R.siden{hemiI}(1) '_akinesia_' R.condname{1} '(' num2str(pind) ')']);
            [maxcoh1 ind] = max(mean(OFF(frq>R.bandef(R.bregband{breg},1) & frq<R.bandef(R.bregband{breg},2),:),1));
            highbetacoh{1,side,sub} = [maxcoh1 ind maxcoh1>0.05];
            [maxcoh2 ind] = max(mean(ON(frq>R.bandef(R.bregband{breg},1) & frq<R.bandef(R.bregband{breg},2),:),1));
            highbetacoh{2,side,sub} = [maxcoh2 ind maxcoh2>0.05];
            highbetacoh{3,side,sub} = maxcoh1>0.05 && maxcoh2>0.05;
            highbetacoh{3,side,sub} = maxcoh2 - maxcoh1;
            %             refidSave{side,sub} = refid;
            title(['OFF ' num2str(updrsSave{2,side,sub}) ' OFF-ON: ' num2str(updrsSave{3,side,sub})])
                        savefigure_v2([R.datapathr '\images\spectral\' R.subname{sub} '\'],['STN_Source_Power_analysis_' R.subname{sub} '_' R.siden{side}  '_'  R.bregname{breg}],[1:2],[],'-r200');
            close all
        end
    end
    
    % figure
    barplot_coh_groups(R.datapathr,highbetacoh)
    savefigure_v2([R.datapathr 'results\spectral\'],['STN_Source_Power_analysis_GroupAverage_boxplots_'  R.bregname{breg}],[],[],'-r100'); close all
    switch breg; case 1; ylimz = [0 0.05]; case 2; ylimz = [0 0.1]; end
    plotNPD(Hz,npdspctrm_group,R,ylimz,1)
    
    spectralplots_groups(R,powsubgrand,cohsubgrand,frqsave,{R.bregname{breg},'STN'},1)
    close all
    save([R.datapathr 'results\spectral\UPDRS_Connect_groupdata.m'],'updrsSave','highbetacoh')
    
    a = vertcat(updrsSave{1,:,:});
    b = vertcat(highbetacoh{1,:,:}); b = b(:,1);
    a(isnan(b)) = []; b(isnan(b)) = [];
    xlab = 'OFF Total Score'; ylab = ['Maximum OFF ' R.bregname{breg} '/STN ' R.bandname{R.bregband{breg}} ' Coh.'];
    figure; linplot(a,b,xlab,ylab,'ro');
    
    a = vertcat(updrsSave{2,:,:});
    b = vertcat(highbetacoh{1,:,:}); b = b(:,1);
    a(isnan(b)) = []; b(isnan(b)) = [];
    xlab = 'OFF Hemi Aknesia Score'; ylab = ['Maximum OFF ' R.bregname{breg} '/STN ' R.bandname{R.bregband{breg}} ' Coh.'];
    figure; linplot(a,b,xlab,ylab,'ro');
    
    a = vertcat(updrsSave{3,:,:});
    b = vertcat(highbetacoh{3,:,:}); b = b(:,1);
    a(isnan(b)) = []; b(isnan(b)) = [];
    xlab = 'ON-OFF Total Score'; ylab = ['ON-OFF ' R.bregname{breg} '/STN ' R.bandname{R.bregband{breg}} ' Coh.'];
    figure; linplot(a,b,xlab,ylab,'bo');
    
    a = vertcat(updrsSave{4,:,:});
    b = vertcat(highbetacoh{3,:,:}); b = b(:,1);
    a(isnan(b)) = []; b(isnan(b)) = [];
    xlab = 'ON-OFF Hemi Aknesia Score'; ylab = ['ON-OFF ' R.bregname{breg} '/STN ' R.bandname{R.bregband{breg}} ' Coh.'];
    figure; linplot(a,b,xlab,ylab,'bo');    
    
    savefigure_v2([R.datapathr 'results\spectral\'],['STN_WPLI_H_UPDRS_Corr_'  R.bregname{breg}],[],[],'-r100'); close all
end


%% For Location Plotting

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
