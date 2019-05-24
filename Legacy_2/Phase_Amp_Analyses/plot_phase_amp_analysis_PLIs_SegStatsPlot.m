function plot_phase_amp_analysis_PLIs_SegStatsPlot(R)
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
                    tendtot = vc_clean.PA.timevec{1}(end)-vc_clean.PA.timevec{1}(1);
                    
                    % Length/Amp Correlations
                    for i = 1:3
                        x = segLCol'; y = ampSegCol(i,:)';
                        [x y] = remnan(x,y);
                        % %                         nleg = 50; nseg  = fix(length(x)/nleg);
                        % %                         x = reshape(x(1:nleg*nseg),nseg,nleg);
                        % %                         y = reshape(y(1:nleg*nseg),nseg,nleg);
                        % %                         for ip =1:nseg
                        % %                             [r(ip) p(ip)] = corr(log10(x(ip,:))',y(ip,:)','Type','Spearman');
                        % %                         end
                        % %                         rGroup(i,sub,side,cond) = mean(r);
                        % %                         pGroup(i,sub,side,cond) = mean(p);
                        
                        [r p] = corr(log10(x),y,'Type','Spearman');
                        rGroup(i,sub,side,cond) = r;
                        pGroup(i,sub,side,cond) = p;
                    end
                    
                    % Relative Phase
                    dp = pi/12;
                    phibins = -pi:dp:pi;
                    phibinmid = phibins(1:end-1)+(dp/2);
                    c = histcounts(relativePhiCol,phibins);
                    phiGroup{sub,side,cond} = c./tendtot;
                    % Segment Lengths
                    db = 0.1;
                    segbins = -1.5:db:1;
                    segbinmid = segbins(1:end-1)+(db/2);
                    c = histcounts(log10(segLCol),segbins); %,'Normalization','probability');
                    segLGroup{sub,side,cond} = c./tendtot;
                    %
                end
            end
        end
    end
    
    figure(1+(10*breg))
    PhiPlot(phiGroup,phibinmid,cmap,R)
    figure(2+(10*breg))
    segLPlot(segLGroup,segbinmid,cmap,R)
    figure(3+(10*breg))
    obs = {['Seg. L. Correlation ' R.bregname{breg} ' ' R.bandinits{R.bregband{breg}}],['Seg. L. Correlation ' 'STN ' R.bandinits{R.bregband{breg}}],['Seg. L. Correlation ' 'STN ' R.bandinits{2}],[R.bregname{breg} '/STN ' R.bandinits{breg} ' Seg. Length']};
    for i = 1:3
        r1 = squeeze(rGroup(i,:,:,1));
        r2 = squeeze(rGroup(i,:,:,2));
        p1 = squeeze(pGroup(i,:,:,1)); p1 = p1(:);
        p2 = squeeze(pGroup(i,:,:,2)); p2 = p2(:);
        s1 = sum(p1<0.05 & p1>0); t1 = sum(p1>0); 
        r1s =r1(p1<0.05 & p1>0) 
        mean(r1s)
        s2 = sum(p2<0.05 & p2>0); t2 = sum(p2>0); r2s =r2(p1<0.05 & p1>0) 
        
        [h,p, chi2stat,df] = prop_test([s1 s2], [t1 t2],'false')

        
        subplot(1,3,i)
        boxploter(r1,r2,cmap(i,:),R,obs{i})
    end
    clear segLGroup phiGroup rGroup
    close all
end
%% FUNCTION ENDS HERE

