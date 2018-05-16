function plot_phase_amp_analysis_PLIs_Cagnan(R)
if nargin<1
    R = makeHeader_SubCort_Cort_Networks();
end
close all
analynames = {'Segment Length','CTX High Beta Amp','STN High Beta Amp','STN Low Beta Amp','STN/CTX High Beta Amp Correlation','Causal Density'};

QX = 8 ; % Bin Size
for breg = 1:length(R.bregname)
    for sub = 1:length(R.subname)
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
                    
                    
                    phiBin = linspace(-pi,pi,QX);
                    phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);
                    [sub side cond]
                    for i = 1:3
                        shiftPhiCol = findAmpPhi(R,ampSegCol(i,:),relativePhiCol,phiBin);
                        [ampBinMu(i,:) ampBinSEM(i,:)] = binstats(shiftPhiCol,ampSegCol(i,:),phiBin);
                    end
                    
                    shiftPhiCol = findAmpPhi(R,segLCol,relativePhiCol,phiBin);
                    [segBinMu segBinSEM] = binstats(shiftPhiCol,segLCol,phiBin);
                    segBinSEM(isnan(segBinSEM)) = 0;
                    
                    
                    ampBinGroup{cond,sub} = ampBinMu;
                    segBinGroup{cond,sub} = segBinMu;
                    figure(1)
                    cmap = linspecer(3);
                    panlist = [1 3 5 ; 2 4 6]; ylimlist = {{[-50 100];[-50 100];[-50 100]},{[-50 100];[-50 150];[-50 150]}};
                    obs = {[R.bregname{breg} ' ' R.bandinits{R.bregband{breg}}],['STN ' R.bandinits{R.bregband{breg}}],['STN ' R.bandinits{2}]};
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
                        ylabel(['% Change in ' obs{i}]); xlabel('Relative Phase')
                        grid on
                    end
                    
                    cmap = linspecer(5);
                    panlistseg = [7; 8];
                    subplot(4,2,panlistseg(cond))
                    hl = plot(phiBinMid', segBinMu(1,:)','--','color',cmap(5,:)); hold on
                    % %                     [hl, hp] = boundedline(phiBinMid', segBinMu(1,:)',segBinSEM(1,:)','cmap',cmap(5,:));
                    % %                     if cond == 1; hl.LineStyle = '--'; end
                    % %                     hp.FaceAlpha = 0.4;
                    % %                     [xq yq R2 exitflag] = VMfit(phiBinMid',segBinMu(1,:)',20,[0,-10,0.01],[1,10,1e3],0);
                    % %                     [xq yq R2 exitflag] = sinfit(phiBinMid', segBinMu(1,:)',20,[0.01; 2*pi; 0; -1],[0.5; 2*pi; 2*pi; 1],0);
                    % %                     if exitflag == 1
                    % %                     hold on; plot(xq,yq,'color',cmap(5,:));
                    % %                     end
                    ylimlist = {[0 1],[0 0.3]};
                    ylim(ylimlist{breg}); grid on
                    title(['STN-' R.bregname{breg} ' ' R.bandinits{R.bregband{breg}} ' Phase vs ' R.bandinits{R.bregband{breg}} ' Frame Length'])
                    ylabel(['LB Segment Length (s)']); xlabel('Relative Phase')
                end
            end
        end
    end
    
    figure(1)
    set(gcf,'Position',[680    -3   864   999])
    for cond=1:2
        y =  horzcat(ampBinGroup{cond,:});
        for i = 1:3
            subplot(4,2,panlist(cond,i)); hold on
            y1 = reshape(y(i,:),size(phiBinMid,2),[])';
            [hl hp] = boundedline(phiBinMid', nanmean(y1)',nanstd(y1)','cmap',cmap(i,:)); %
            hl.LineWidth = 2;
            L = get(gca,'children'); L(2) = L(end); L(end) = hp;set(gca,'children',L)
            %             a = plot(phiBinMid,nanmean(reshape(y(i,:),size(phiBinMid,2),[])'),'color',cmap(i,:),'LineWidth',2);
            % %             if cond == 1; a.LineStyle = '--'; end
        end
    end
    
    cmap = linspecer(5);
    for cond=1:2
        y =  vertcat(segBinGroup{cond,:});
        subplot(4,2,panlistseg(cond)); hold on
        [hl hp] = boundedline(phiBinMid', nanmean(y)',nanstd(y),'cmap',cmap(5,:));
        hl.LineWidth = 2;
        L = get(gca,'children'); L(2) = L(end); L(end) = hp; set(gca,'children',L)
        %         a = plot(phiBinMid,nanmean(y),'color',cmap(5,:),'LineWidth',2);
        % %         if cond == 1; a.LineStyle = '--'; end
    end
    
    annotation(gcf,'textbox',...
        [0.187 0.956 0.205 0.035],...
        'String',{'OFF L-Dopa'},...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',18,...
        'FitBoxToText','off');
    
    annotation(gcf,'textbox',...
        [0.624 0.959 0.205 0.0349],...
        'String','ON L-Dopa',...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',18,...
        'FitBoxToText','off');
    savefigure_v2([R.resultspathr '\Group\seganalysis\'],...
        [R.bregname{breg} '_Group_CagnanAnalysis'],[],[],'-r100'); close all
    close all
    clear ampBinGroup segBinGroup
end




