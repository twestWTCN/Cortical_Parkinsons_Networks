function plot_phase_amp_analysis_PLIs_Cagnan(R)
if nargin<1
    R = makeHeader_SubCort_Cort_Networks();
end
close all
analynames = {'Segment Length','CTX High Beta Amp','STN High Beta Amp','STN Low Beta Amp','STN/CTX High Beta Amp Correlation','Causal Density'};
% WIP
% Percentage Change over what? Baseline or baseline of segments?
QX = 8; % Bin Size
% sub 3 side 2 - very peaked OFF
% sub 1 side 1 - very flat but amplified
for breg = 2:length(R.bregname)
    cmapint = linspecer(3);
    cmap = linspecer(5);
    cmapint(4,:) = cmap(5,:);
    
    for sub = 1:length(R.subname)
        cmapint = cmapint*0.95;
        for side = 1:2
            for cond = 1:length(R.condname)
                load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                if vc_clean.specanaly.flag ~= 1 % Check significant coherences                        load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\' idd '_ROIvoxel_phaseamp_PLI_' R.ipsicon '_' R.siden{side} '_' R.bandname{band} '.mat'])
                    phi = vc_clean.PA.pA_pli_dist_save';
                    relativePhiCol = phi'; %wrapToPi(phi-circ_mean(phi(~isnan(phi'))) )'; %circ_mean(pA_pli_dist_save{2,nrOFF}(amp_pli_dist_save{2,nrOFF}(3,:)>95)')); %wrapToPi(pA_pli_dist_save{cond,nr}); %-circ_mean([pA_pli_dist_save{2,nrOFF}]'));
                    segLCol =  vc_clean.PA.segL_pli_dist_save; %((segL_pli_dist_save{cond,nr} - mean([segL_pli_dist_save{cond,nr}],2))./mean([segL_pli_dist_save{cond,nr}],2)  )*100;
                    ampSegCol = vc_clean.PA.amp_pli_dist_save; %((amp_pli_dist_save{cond,nr}  - mean([amp_pli_dist_save{cond,nr}],2))./mean([amp_pli_dist_save{cond,nr}],2) )*100;
                    HdistSegCol = vc_clean.PA.H_dist_save(1,:);
                    tendtot = vc_clean.PA.timevec{1}(end)-vc_clean.PA.timevec{1}(1);
                    
                    figure(30)
                    panlist = [1 3 5 ; 2 4 6]; ylimlist = {{[-25 350];[-25 250];[-25 250]},{[-100 200];[-100 300];[-25 300]}};
                    obs = {[R.bregname{breg} ' ' R.bandinits{R.bregband{breg}}],['STN ' R.bandinits{R.bregband{breg}}],['STN ' R.bandinits{2}],[R.bregname{breg} '/STN ' R.bandinits{breg} ' Seg. Length']};
                    for i = 1:3
                        subplot(3,2,panlist(cond,i))
                        [A stat] = linplot_PD(log10(segLCol)',ampSegCol(i,:)','Seg Length (s)','Amplitude',cmapint(i,:)); xlim([-1.5 1])
                        ylim(ylimlist{breg}{i})
                        Rcoeff(:,:,i) = stat.modcoef.*(stat.p<0.05);
                        title(['STN-' R.bregname{breg} ' ' R.bandinits{R.bregband{breg}} ' Phase vs ' obs{i} ' Power'])
                        ylabel(['% Change in ' obs{i}]); xlabel('Segment Length (s)')
                        
                    end
                    ampsegColGroup{cond,side,sub} = ampSegCol;
                    segLColGroup{cond,side,sub} = segLCol;
                    RcoeffGroup{cond,side,sub} = Rcoeff;
                    phiBin = linspace(-pi,pi,QX);
                    phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);
                    [sub side cond]
                    for i = 1:3
                        [shiftPhiCol phipeak(i)]= findAmpPhi(R,ampSegCol(i,:),relativePhiCol,phiBin);
                        [ampBinMu(i,:) ampBinSEM(i,:)] = binstats(shiftPhiCol,ampSegCol(i,:),phiBin);
                    end
                    
                    [shiftPhiCol phipeak(4)]= findAmpPhi(R,segLCol,relativePhiCol,phiBin);
                    [segBinMu segBinSEM] = binstats(shiftPhiCol,segLCol,phiBin);
                    segBinSEM(isnan(segBinSEM)) = 0;
                    
                    
                    ampBinGroup{cond,side,sub} = ampBinMu;
                    segBinGroup{cond,side,sub} = segBinMu;
                    phipeakGroup{cond,side,sub} = phipeak;
                    
                    figure(1)
                    cmap = linspecer(3);
                    ylimlist = {{[-25 350];[-25 250];[-25 250]},{[-25 75];[-25 200];[-25 250]}};
                    for i = 1:3
                        subplot(4,2,panlist(cond,i))
                        hl = plot(phiBinMid', ampBinMu(i,:)','--','color',cmap(i,:)); hold on
                        % %                         [hl, hp] = boundedline(phiBinMid', ampBinMu(i,:)',ampBinSEM(i,:)','cmap',cmap(i,:)); hold on
                        % %                         if cond == 1; hl.LineStyle = '--'; end
                        % %                         hp.FaceAlpha = 0.4;
                        % %                         [xq yq R2 exitflag] = VMfit(phiBinMid',ampBinMu(i,:)',20,[0,-pi,-100],[100,pi,1e3],0);
                        % %                         if exitflag == 1
                        % %                             hold on; plot(xq,yq,'color',cmap(i,:));
                        % %                         end
                        % %                         [xq yq R2 exitflag] = sinfit(phiBinMid', ampBinMu(i,:)',20,[25; 2*pi; 0; -500],[300;2*pi; 2*pi; 500],0);
                        % %                         plot(xq,yq,'color',cmap(i,:));
                        ylim(ylimlist{breg}{i})
                        title(['STN-' R.bregname{breg} ' ' R.bandinits{R.bregband{breg}} ' Phase vs ' obs{i} ' Power'])
                        ylabel(['% Change in ' obs{i}]); xlabel('Relative Phase (rads)')
                        grid on
                    end
                    
                    cmap = linspecer(5);
                    panlist(:,4) = [7; 8];
                    subplot(4,2,panlist(cond,4))
                    hl = plot(phiBinMid', segBinMu(1,:)','--','color',cmap(5,:)); hold on
                    % %                     [hl, hp] = boundedline(phiBinMid', segBinMu(1,:)',segBinSEM(1,:)','cmap',cmap(5,:));
                    % %                     if cond == 1; hl.LineStyle = '--'; end
                    % %                     hp.FaceAlpha = 0.4;
                    % %                     [xq yq R2 exitflag] = VMfit(phiBinMid',segBinMu(1,:)',20,[0,-10,0.01],[1,10,1e3],0);
                    % %                     [xq yq R2 exitflag] = sinfit(phiBinMid', segBinMu(1,:)',20,[0.01; 2*pi; 0; -1],[0.5; 2*pi; 2*pi; 1],0);
                    % %                     if exitflag == 1
                    % %                     hold on; plot(xq,yq,'color',cmap(5,:));
                    % %                     end
                    ylimlist = {[0 1],[0 1]};
                    ylim(ylimlist{breg}); grid on
                    title(['STN-' R.bregname{breg} ' ' R.bandinits{R.bregband{breg}} ' Phase vs ' R.bandinits{R.bregband{breg}} ' Frame Length'])
                    
                    ylabel(['LB Segment Length (s)']); xlabel('Relative Phase (rads)')
                end
            end
        end
    end
    figure(30)
    set(gcf,'Position',[1190         330         627         750])
    x1 = vertcat(RcoeffGroup{1,:,:});
    x2 = vertcat(RcoeffGroup{2,:,:});
    xl = -1:.1:0.5; xl= [ones(length(xl),1) xl'];
    for cond = 1:2
        for i=1:3
            subplot(3,2,panlist(cond,i))
            if cond == 1; coefs = x1(:,:,i); else;   coefs = x2(:,:,i)'; end
            if size(coefs,2)>size(coefs,1); coefs = coefs'; end
            W = sum(coefs(:,1)~=0)./size(coefs,1);
            coefs(coefs(:,1)==0,:) = [];
            if ~isempty(coefs); coefs = (mean(coefs,1).*W)'; else; coefs = [0 0]'; end
            plot(xl(:,2),xl*coefs,'color',cmap(i,:),'LineWidth',2)
            a = gca; uistack(findobj(a.Children,'type','Line'),'top');
        end
    end
    a = gca;
    %     mean(x1)
    %     mean(x2)
    [h,p, chi2stat,df] = prop_test([sum(x2(:,1:2)>0)], [size(x2,1),size(x2,1)],'false');
    [h,p, chi2stat,df] = prop_test([sum(x1(:,1:2)>0)], [size(x1,1),size(x1,1)],'false');
    [h p]=ttest2(x1(:,2),x2(:,2));
    [h p]=ttest2(x1(:,1),x2(:,1));
    
    figure(1)
    set(gcf,'Position',[539    81   626   999])
    cmap = linspecer(3);
    for cond=1:2
        y =  horzcat(ampBinGroup{cond,:,:});
        for i = 1:3
            subplot(4,2,panlist(cond,i)); hold on
            y1 = reshape(y(i,:),size(phiBinMid,2),[])';
            ysave{cond,i} = y1;
            [hl hp] = boundedline(phiBinMid', nanmedian(y1)',abs([prctile(y1,16)' prctile(y1,84)']-nanmedian(y1)'),'cmap',cmap(i,:)); %
            hl.LineWidth = 2;
            L = get(gca,'children'); L(2) = L(end); L(end) = hp;set(gca,'children',L)
            xlim([-pi pi])
            if cond ==2
                x1 = sqres(ysave{1,i},nanmedian(ysave{1,i}))/1000; x2 = sqres(ysave{2,i},nanmedian(ysave{2,i}))/1000;
                p(i) = ranksum(x1,x2);
                a = gca;[figx figy] = dsxy2figxy(gca, -3,a.YLim(1)*0.95);
                annotation(gcf,'textbox',...
                    [figx figy 0.3 0.03],...
                    'String',sprintf('RS test = %.3f',p(i)),...
                    'LineStyle','none',...
                    'FontSize',8,...
                    'FontWeight','bold',...
                    'FontAngle','italic',...
                    'FitBoxToText','off');
            end
        end
    end
    cmapA = linspecer(3);
    cmap = linspecer(5);
    cmapA(4,:) = cmap(5,:);
    for cond=1:2
        y =  vertcat(segBinGroup{cond,:,:});
        ysave{cond,4} = y;
        subplot(4,2,panlist(cond,4)); hold on
        [hl hp] = boundedline(phiBinMid', nanmedian(y)',abs([prctile(y,16)' prctile(y,84)']-nanmedian(y)'),'cmap',cmap(5,:));
        hl.LineWidth = 2;
        L = get(gca,'children'); L(2) = L(end); L(end) = hp; set(gca,'children',L)
        xlim([-pi pi])
        if cond ==2
            x1 = sqres(ysave{1,i},nanmedian(ysave{1,i}))/1000; x2 = sqres(ysave{2,i},nanmedian(ysave{2,i}))/1000;
            p(4) = ranksum(x1,x2);
            a = gca;[figx figy] = dsxy2figxy(gca, -3,a.YLim(1)*0.95);
            annotation(gcf,'textbox',...
                [figx figy 0.3 0.03],...
                'String',sprintf('RS test = %.3f',p(4)),...
                'LineStyle','none',...
                'FontSize',8,...
                'FontWeight','bold',...
                'FontAngle','italic',...
                'FitBoxToText','off');
        end
    end
    clear p
    [ptt,pvart] = cagnan_stats_test(phiBinMid,ysave,1,panlist,cmapA);
    annotation(gcf,'textbox',...
        [ 0.1820    0.9570    0.2550    0.0350],...
        'String',{'OFF L-Dopa'},...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',18,...
        'FitBoxToText','off');
    
    annotation(gcf,'textbox',...
        [0.583 0.957 0.304 0.0349],...
        'String','ON L-Dopa',...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',18,...
        'FitBoxToText','off');
    savefigure_v2([R.resultspathr '\Group\seganalysis\'],...
        [R.bregname{breg} '_Group_CagnanAnalysis'],[],[],'-r100');
    close all
    clear ampBinGroup segBinGroup
    
    figure
    set(gcf,'Position',[823   421   785   677])
    x1 = vertcat(phipeakGroup{1,:,:});
    x2 = vertcat(phipeakGroup{2,:,:});
    
    legloc = legloc();
    for i = 1:4
        subplot(2,2,i)
        x11 = x1(:,i);
        x22 = x2(:,i);
        polarhistogram(x11,6,'Normalization','Probability','FaceColor',cmapA(i,:).*0.75,'LineStyle','none','FaceAlpha',.5); hold on
        a(i,1) =  polarhistogram(x11,6,'Normalization','Probability','DisplayStyle','stairs','EdgeColor',cmapA(i,:).*0.75,'LineWidth',1.2);
        polarhistogram(x22,6,'Normalization','Probability','FaceColor',cmapA(i,:),'LineStyle','none','FaceAlpha',.3); hold on
        a(i,2) = polarhistogram(x22,6,'Normalization','Probability','DisplayStyle','stairs','EdgeColor',cmapA(i,:),'LineWidth',1.5);
        title(['RP Dist. for Peak ' obs{i}]); rlim([0 0.5])
        legend(a(i,:),R.condname,'Position',legloc(i,:),'Orientation','Horizontal','Color','none','EdgeColor','none')
        %         [p,m,pv] = circ_cmtest([x11;x22]',[repmat(1,size(x11));repmat(2,size(x22))]');
        %circ_wwtest([x11;x22],[repmat(1,size(x11));repmat(2,size(x22))]); % ttest
        p(i,1) = circ_rtest(x11);rayleg(p(i,1),i,1);
        p(i,2) = circ_rtest(x22);rayleg(p(i,2),i,2);
        %         [sig(i)] = disp_anova(x11,x22)
    end
    savefigure_v2([R.resultspathr '\Group\seganalysis\'],...
        [R.bregname{breg} '_Group_CagnanAnalysis_RP_rose'],[],[],'-r100'); close all
    
end

function legloc = legloc();
legloc = [0.1710    0.5210    0.2500    0.0300;
    0.6110    0.5210    0.2500    0.0300;
    0.1710    0.0530    0.2500    0.0300;
    0.6110    0.0530    0.2500    0.0300];
function rayleg(P,i,cond)
posbank = {[0.1500    0.5000    0.1560    0.0250; 0.3400    0.5000    0.1560    0.0250];
    [0.5800    0.5000    0.1560    0.0250; 0.7700    0.5000    0.1560    0.0250];
    [0.1700    0.0210    0.1560    0.0250; 0.3400    0.0210    0.1560    0.0250];
    [0.5800    0.0210    0.1560    0.0250; 0.7700    0.0210    0.1560    0.0250]};

annotation(gcf,'textbox',...
    posbank{i}(cond,:),...
    'String',sprintf('Ray. P = %.3f',P),...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontAngle','italic',...
    'FitBoxToText','off');

%%[sprintf('OFF Res. = %.1f ',nanmean(x1)) sprintf('ON Res. = %.1f ',nanmean(x2))];
