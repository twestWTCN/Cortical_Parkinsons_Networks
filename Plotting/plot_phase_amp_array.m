function plot_phase_amp_array(R)
if nargin<2
    R = makeHeader_SubCort_Cort_Networks();
end
close all
analynames = {'Segment Length','CTX High Beta Amp','STN High Beta Amp','STN Low Beta Amp','STN/CTX High Beta Amp Correlation','Causal Density'};
QXN = 8 ; QYN = 8; %
QX = linspace(-pi,pi,QXN);
climits = [0 0.025];

ampylim = {{[-25 350];[-25 250];[-25 250]},{[-25 100];[-25 250];[-25 250]}};
segylim = {[0 1],[0 1]};

%%%
for breg = 1:length(R.bregname)
    for sub = 1:length(R.subname)
        for side = 1:2
            for cond = 1:length(R.condname)
                load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                if vc_clean.specanaly.flag ~= 1 % Check significant coherences                        load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\' idd '_ROIvoxel_phaseamp_PLI_' R.ipsicon '_' R.siden{side} '_' R.bandname{band} '.mat'])
                    relativePhiCol = vc_clean.PA.pA_pli_dist_save;
                    segLCol =  vc_clean.PA.segL_pli_dist_save;
                    ampSegCol = vc_clean.PA.amp_pli_dist_save;
                    HdistSegCol = vc_clean.PA.H_dist_save{cond}(1,:);
                    tendtot= vc_clean.PA.timevec{1}(end)-vc_clean.PA.timevec{1}(1);
                end
                % %                     for i = 1:3
                % %                         scatter(segLCol,ampSegCol(i,:)); hold on
                % %                         [yCalc b Rsq] = linregress(segLCol(~isnan(ampSegCol(i,:)))',ampSegCol(i,~isnan(ampSegCol(i,:)))');
                % %                         plot(segLCol(~isnan(ampSegCol(i,:))),yCalc)
                % %                     end
                % Plots of rel phase vs seg length
                phishift = findAmpPhi(R,segLCol,relativePhiCol,QX);
                f(1) = figure(1*10 + 4); % PA vs SegL
                f(2) = figure((2)*10 +1); set(gcf,'Position',[263 517 1403 352])
                subplot(1,3,cond);
                [L4(cond) segLbin(:,cond) binmid pdense(1,cond)] = plot_PA_InterpArray(cond,phishift,segLCol,...
                    'SRP Length',QX,linspace(segylim{breg}(1),segylim{breg}(2),QYN),f,tendtot,R.PA.interpolgrid); % ,logspace(log10(0.5),log10(1.2),6)
                %                     xlim([-3.5 3.5]); %ylim([0 1.75]);
                figure(f(2)); set(gcf,'Position',[263 517 1403 352]);  if cond>0; colorbar; end; title(R.condname{cond})
                caxis(climits);
                
                % Plots of rel phase vs
                phishift = findAmpPhi(R,ampSegCol(1,:),relativePhiCol,QX);
                f(1) = figure(1*10 + 5); % PA vs M1 Amp
                f(2) = figure((2)*10 +2);
                subplot(1,3,cond);
                [L5(cond) ampbin(:,cond,1) binmid pdense(2,cond)] = plot_PA_InterpArray(cond,phishift,ampSegCol(1,:),...
                    [R.bregname{breg} ' ' R.bandname{R.bregband{breg}} ' Power'],QX,linspace(ampylim{breg}{1}(1),ampylim{breg}{1}(2),QYN),f,tendtot,R.PA.interpolgrid); %linspace(0,7,8) ,linspace(-100,200,QY)
                figure(f(2)); set(gcf,'Position',[263 517 1403 352]); if cond>0; colorbar; end; title(R.condname{cond})
                caxis(climits);
                
                % Plots of rel phase vs STN HB Amp
                phishift = findAmpPhi(R,ampSegCol(2,:),relativePhiCol,QX);
                f(1) = figure(1*10 + 6); % PA vs STN HB
                f(2) = figure((2)*10 +3);
                subplot(1,3,cond)
                [L6(cond) ampbin(:,cond,2) binmid pdense(3,cond)] = plot_PA_InterpArray(cond,phishift,ampSegCol(2,:),...
                    ['STN ' R.bandname{R.bregband{breg}} ' Power'],QX,linspace(ampylim{breg}{2}(1),ampylim{breg}{2}(2),QYN),f,tendtot,R.PA.interpolgrid); %linspace(0,7,8) linspace(-100,200,QY)
                figure(f(2)); set(gcf,'Position',[263 517 1403 352]);  if cond>0; colorbar; end; title(R.condname{cond})
                caxis(climits);
                
                % Plots of rel phase vs STN LB Amp
                phishift = findAmpPhi(R,ampSegCol(3,:),relativePhiCol,QX);
                f(1) = figure(1*10 + 7); % PA vs STN LB
                f(2) = figure((2)*10 +4); set(gcf,'Position',[263 517 1403 352])
                subplot(1,3,cond)
                [L7(cond) ampbin(:,cond,3) binmid pdense(4,cond)] = plot_PA_InterpArray(cond,phishift,ampSegCol(3,:),...
                    'STN Low Beta Amp',QX,linspace(ampylim{breg}{3}(1),ampylim{breg}{3}(2),QYN),f,tendtot,R.PA.interpolgrid); % linspace(0,4,8) linspace(-100,200,QY)
                figure(f(2)); set(gcf,'Position',[263 517 1403 352]);  if cond>0; colorbar; end; title(R.condname{cond})
                caxis(climits);
                %                 % Plots of rel phase vs Segment Length DDT
                %                 figure(100)
                %                 [H1(cond) hdist(:,cond,1)] = plot_segL_histogram(relativePhiCol,segLCol,'Segment Length',cond);
                
                %Plots of rel phase vs Segment Length DDT
                f(1) = figure(1*10 + 8); % PA vs SegL
                f(2) = figure((2)*10 +5);
                subplot(1,3,cond)
                set(gcf,'Position',[263 517 1403 352])
                %                 [H2(cond) hdist(:,cond,2)] = plot_segL_histogram(pA_pli_dist_save{cond},segL_pli_dist_save{cond},'PLI Segment Length',cond);
                [L6(cond) hbin(:,cond,1) binmid pdense(5,cond)] = plot_PA_InterpArray(cond,relativePhiCol,HdistSegCol,...
                    ['STN/' R.bregname{breg} ' ' R.bandname{R.bregband{breg}} ' Amp Correlation'],QX,linspace(-1.1,1.1,QYN),f,tendtot,R.PA.interpolgrid); % linspace(0,4,8)
                caxis(climits);
                %% Phase Angle Rose DDT
                figure(200)
                pA_dist = relativePhiCol;
                pA_dist(isnan(pA_dist)) = [];
                R_1(cond) = polarhistogram(pA_dist,18,'FaceAlpha',0.75); hold on
                % %                 phaseAng_dist(cond).circ_mean_std = [circ_mean(pA_dist') circ_var(pA_dist')];
                % %                 [h mu] =circ_mtest(pA_dist',circ_mean(pA_dist'));
                % %                 phaseAng_dist(cond).circ_meantest = [h mu];
                % %                 %                     [pval, z] = circ_rtest(pA_dist');
                % %                 phaseAng_dist(cond).circ_RayTest = [1 1]; %[pval, z];
                
                title('Phase Angle Distribution Phase Lock','FontSize',18)
                figure(252)
                y =  HdistSegCol; x = segLCol;
                scatter(x,y(1,:));hold on;[r1 p1] = corrcoef(x,y(1,:)); length_amp_corr(:,cond,sub) = [r1(2) p1(2)];
                xlabel('Segment Length'); ylabel('Amplitude Correlation')
            end
            %             legend(H1,R.condname)
            legend(R_1,R.condname)
            %         legend(H2,R.condname)
            %         legend(R_2,R.condname)
            legend(L5,R.condname)
            legend(L6,R.condname)
            legend(L7,R.condname)
            for i = 1:5
                figure(20+i)
                subplot(1,3,3)
                N = pdense(i,2).N - pdense(i,1).N;
                x = pdense(i,2).Xedges;
                y = pdense(i,2).Yedges;
                pcolor(x,y,N); xlim([-pi pi]); ylim([y(1) y(end)])
                %                      contourf(x,y,N); xlim([-pi pi]); ylim([y(1) y(end)])
                %                     pcolor(x,y,N');
                % imagesc(x,y,av1');
                caxis([-0.01 0.01]); title(R.condname{1});
                xlabel('Phi_1 - Phi_2'); ylabel(analynames{i}); set(gca,'YDir','normal');
                h = colorbar; a = ylabel(h, 'Occurence (s^{-1})');
                set(a,'rotation',270); set(a,'Position',get(a,'Position') + [0.8 0 0]);
                title('OFF - ON')
            end
            mkdir([R.resultspathr '\' R.subname{sub} '\seganalysis\' R.bregname{breg} '_' R.siden{side} '\'])
%             savefigure_v2([R.resultspathr '\' R.subname{sub} '\seganalysis\' R.bregname{breg} '_' R.siden{side} '\'],['STN_Source_PA_Amp_Array_' R.bregname{breg} '_' R.siden{side}],[],[],'-r100'); close all
            close all
            segLsave{sub,side} = segLbin;
            densesave{sub,side} = pdense;
            clear segLbin ampbin pdense hdist phaseAng_dist gc_dist_c gc_dist_cd
        end
    end
    mkdir([R.resultspathr '\Group\seganalysis'])
    save([R.resultspathr '\Group\seganalysis\pdense_' R.bregname{breg}],'densesave')
    clear densesave segLsave
end
% save([R.datapathr '\results\seganalysis\groupseganaly_' R.bandname{band}],'segLsave','amp_pli_dist_save','densesave','hdistsave','phaseAng_dist_save'); %,'gc_dist_sub_save','GC_stat_table')
% ON = sum(squeeze(length_amp_corr(:,1,:,:)),3);
% OFF = sum(squeeze(length_amp_corr(:,2,:,:)),3);
% [h p] = ttest2(ON(1,:),OFF(1,:))
% mean(ON(1,:))
% mean(OFF(1,:))
% 




