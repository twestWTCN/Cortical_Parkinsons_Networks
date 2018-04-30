function plot_subject_cohspectra(R)
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
        if band == 1
            delete([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon '.mat']);
            eval(['!del /q ' R.datapathr R.subname{sub} '\ftdata\ROI_analy\'])
        end
        if 1% exist([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon '_'  R.bandname{band} '.mat']) == 0
            [idbank frqbank stn_lb_frqbank] = find_voxel_pow_coh_v4(R,sub,band);
            load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_power_' R.ipsicon '_'  R.bandname{band} '.mat'],'powsave')
            load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_coh_' R.ipsicon '_'  R.bandname{band}],'cohsave','frqsave')
            load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon '_'  R.bandname{band}],'frqbank','idbank','stn_lb_frqbank','refsave')
        else
            load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_power_' R.ipsicon '_'  R.bandname{band}],'powsave')
            load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_coh_' R.ipsicon '_'  R.bandname{band}],'cohsave','frqsave')
            load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon '_'  R.bandname{band}],'frqbank','idbank','stn_lb_frqbank','refsave')
        end
        [~,~,nrep(1),~] = data_fileguide(R.subname{sub},0);
        [~,~,nrep(2),~] = data_fileguide(R.subname{sub},1);
        for side = 1:2
            refid(1) = refsave{nrep(1),side,1};
            refid(2) = refsave{nrep(2),side,2};

            id(1) = idbank(nrep(1),side,1);
            id(2) = idbank(nrep(2),side,2);
            
            load([R.datapathr R.subname{sub} '\ftdata\virtual_sources_' num2str(nrep(1)) '_ROI_' R.condname{1} '_' R.siden{side} '_' R.ipsicon '_'  R.bandname{band}]);
            locbank(:,1,side,sub) = vchansave(id(1)).loc; vchansave_ON = vchansave(id(1));
            load([R.datapathr R.subname{sub} '\ftdata\virtual_sources_' num2str(nrep(2)) '_ROI_' R.condname{2} '_' R.siden{side} '_' R.ipsicon '_'  R.bandname{band}]);
            locbank(:,2,side,sub) = vchansave(id(2)).loc; vchansave_OFF = vchansave(id(2));
            
            [npdspctrm, Hz] = NPD_cortical_STN(vchansave_ON,vchansave_OFF,R,1,refid);
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
                ON = powsave{nrep(1),side,1};
                OFF = powsave{nrep(2),side,2};
                %             id(1) = idbank(end,side,1);
                %             id(2) = idbank(end,side,1)
                ON = squeeze(ON(ch,:,:)); if size(ON,1)<size(ON,2); ON = ON'; end
                OFF = squeeze(OFF(ch,:,:));if size(OFF,1)<size(OFF,2); OFF = OFF'; end
                
                subplot(1,3,ch)
                frq = frqsave{nrep(1),side,1}(:,1)'; powsubgrand{1,ch,side,sub} = ON(:,id(1));
                ax(1) = boundedline(frq,mean(ON,2),std(ON,0,2)/sqrt(size(ON,2)),'cmap',cmap(1,:),'alpha','transparency',0.45);
                set(gca,'xscale','log'); set(gca,'yscale','log')
                xlim([6 45]);%ylim([0.005 0.05])
                hold on
                %             [xCalc,yCalc] = linregress(log10(frq)',log10(mean(ON,2)));
                %             plot(xCalc(:,2),yCalc,'--','color',cmap(1,:),'linewidth',2);
                frq = frqsave{nrep(2),side,2}(:,1)'; powsubgrand{2,ch,side,sub} = OFF(:,id(2));
                ax(2) = boundedline(frq,mean(OFF,2),std(OFF,0,2)/sqrt(size(OFF,2)),'cmap',cmap(2,:),'alpha','transparency',0.45);
                set(gca,'xscale','log'); set(gca,'yscale','log')
                xlim([6 45]);%ylim([0.005 0.05]);
                xlabel('Frequency (Hz)');ylabel('log Normalized Power'); title(R.titular{ch})
                legend(ax,{'ON','OFF'});
                hold on
                %           x  [xCalc,yCalc] = linregress(frq',mean(OFF,2));
                %             plot(xCalc(:,2),yCalc,'--','color',cmap(2,:),'linewidth',2);
                
                grid on
            end
            
            ON =  cohsave{nrep(1),side,1}; cohsubgrand{1,side,sub} = ON(:,id(1));
            OFF = cohsave{nrep(2),side,2}; cohsubgrand{2,side,sub} = OFF(:,id(2));
            
            subplot(1,3,3)
            frq = frqsave{nrep(1),side,1}(:,1)';
            ax(1) = boundedline(frq,mean(ON,2),std(ON,0,2)/sqrt(size(ON,2)),'cmap',cmap(1,:),'alpha','transparency',0.45);
            xlim([4 45]);
            hold on
            %         [xCalc,yCalc] = linregress(frq',mean(ON,2));
            %         plot(xCalc(:,2),yCalc,'--','color',cmap(1,:),'linewidth',2);
            frq = frqsave{nrep(2),side,2}(:,1)';
            ax(2) = boundedline(frq,mean(OFF,2),std(OFF,0,2)/sqrt(size(OFF,2)),'cmap',cmap(2,:),'alpha','transparency',0.45);
            xlim([6 45]); xlabel('Frequency (Hz)');ylabel('Coherence'); title('STN/CTX Coherence')
            legend(ax,{'ON','OFF'})
            hold on
            %         [xCalc,yCalc] = linregress(frq',mean(OFF,2));
            %         plot(xCalc(:,2),yCalc,'--','color',cmap(2,:),'linewidth',2);
            %         set(gca,'xscale','log')
            grid on
            set(gcf,'Position',[314         551        1208         307]); shg
            
            
            pind = find(strcmp(R.subname{sub},patient));
            %         if side == 1; hemiI = 2; elseif side == 2; hemiI = 1; end
            updrsSave{1,side,sub} = eval([R.siden{side}(1) '_akinesia_' R.condname{2} '(' num2str(pind) ')']);
            updrsSave{2,side,sub} = eval(['total_' R.condname{2} '(' num2str(pind) ')']);
            [maxcoh1 ind] = max(mean(ON(frq>R.bandef(band,1) & frq<R.bandef(band,2),:),1));
            highbetacoh{1,side,sub} = [maxcoh1 ind maxcoh1>0.05];
            [maxcoh2 ind] = max(mean(OFF(frq>R.bandef(band,1) & frq<R.bandef(band,2),:),1));
            highbetacoh{2,side,sub} = [maxcoh2 ind maxcoh2>0.05];
            highbetacoh{3,side,sub} = maxcoh1>0.05 && maxcoh2>0.05;
            refidSave{side,sub} = refid;
            savefigure_v2([R.datapathr R.subname{sub} '\images\spectral\'],['STN_Source_Power_analysis_' R.subname{sub} '_' R.siden{side}  '_'  R.bandname{band}],[],[],'-r200');
            close all
        end
    end
    figure
    load('source_template.mat')
    A = min(source.pos); B = max(source.pos);
    X = A(1):.75:B(1);Y = A(2):.75:B(2);Z = A(3):.75:B(3);
    [Xz Yz Zz] = meshgrid(Y,X,Z);
    mask =  reshape(source.avg.pow,source.dim)>0;
    fv = isosurface(double(mask),0); %Yz,Xz,Zz,
    patch(fv,'FaceColor',[.1 .1 .1],'EdgeColor',[0 0 0],'FaceAlpha',0.05);
    hold on
    xyz = squeeze(locbank(:,1,1,:));
    clear a
    a(1) = scatter3(xyz(2,:),xyz(1,:),xyz(3,:),100,'b','x','LineWidth',2); hold on
    xyz = squeeze(locbank(:,1,2,:));
    a(2) = scatter3(xyz(2,:),xyz(1,:),xyz(3,:),100,'b','o','LineWidth',2);
    xyz = squeeze(locbank(:,2,1,:));
    a(3) = scatter3(xyz(2,:),xyz(1,:),xyz(3,:),100,'r','x','LineWidth',2); hold on
    xyz = squeeze(locbank(:,2,2,:));
    a(4) = scatter3(xyz(2,:),xyz(1,:),xyz(3,:),100,'r','o','LineWidth',2);
    legend(a,{'Left ON','Right ON','Left OFF','Right OFF'})
    
    % figure
    barplot_coh_groups(R.datapathr,highbetacoh)
    savefigure_v2([R.datapathr 'results\spectral\'],['STN_Source_Power_analysis_GroupAverage_boxplots_'  R.bandname{band}],[],[],'-r100'); close all
    
    spectralplots_groups(R.datapathr,powsubgrand,cohsubgrand,frq,R.titular,R.bandname{band})
    
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