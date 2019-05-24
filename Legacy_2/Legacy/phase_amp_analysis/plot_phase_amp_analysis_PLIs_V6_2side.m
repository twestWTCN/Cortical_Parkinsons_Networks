function plot_phase_amp_analysis_PLIs_V6_2side(R)
if nargin<1
    R = makeHeader_SubCort_Cort_Networks();
end
idd = '';
% Collate the points and then do paired ttest between bins? A lot of
% information is lost by cutting off in histograms. Either that or somehow
% increase the number of samples - perhaps include all of the STN channels?
% Maybe combine left and right into one histogram? Collapse all of one
% subjects sync segs into one histogram?? All should be sampling from the
% same distribution so build up better resolved image of underlying 2D PDF.
% Or look at doing a hotelling test between the collated data??

close all
analynames = {'Segment Length','CTX High Beta Amp','STN High Beta Amp','STN Low Beta Amp','STN/CTX High Beta Amp Correlation','Causal Density'};

QXN = 8 ; QYN = 8; % 6 6
% QX = [-pi -pi/8 pi/8 pi]';
QX = linspace(-pi,pi,QXN);
% load([R.datapathr 'subject_hbWPLI075'])
% subscreen = squeeze(sum(subject_hbcohscreen>R.PA.WPLIscreen)==2);
a = logspace(0.5,log10(150),4); logscalez = [-100 -50 -15 15 50 100]; %-40:10:40; %linspace(-50,50,QY); %logscalez = [-3.^(4:-1:2) 3.^(2:4)];
for band = 3%; %1:numel(R.bandname)
    for sub = 1:numel(R.subname)
        load([R.datapathr 'results\spectral\screen.mat'])
        for side =1:2
            if screen(side,sub) == 1
                for cond = 1:2
                    [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
                    [datafileN,pp_mark,nrOFF,senscheck] = data_fileguide(R.subname{sub},1);
                    clear relativePhi segL ampSeg HdistSeg tends
                    for nr = 1:nrep
                        if nr == 2;                            a = 1; end
                        load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\' idd '_ROIvoxel_phaseamp_PLI_' R.ipsicon '_' R.siden{side} '_' R.bandname{band} '.mat'])
                        relativePhi{nr} = pA_pli_dist_save{cond,nr}; %circ_mean(pA_pli_dist_save{2,nrOFF}(amp_pli_dist_save{2,nrOFF}(3,:)>95)')); %wrapToPi(pA_pli_dist_save{cond,nr}); %-circ_mean([pA_pli_dist_save{2,nrOFF}]'));
                        segL{nr} =  segL_pli_dist_save{cond,nr}; %((segL_pli_dist_save{cond,nr} - mean([segL_pli_dist_save{cond,nr}],2))./mean([segL_pli_dist_save{cond,nr}],2)  )*100;
                        ampSeg{nr} = amp_pli_dist_save{cond,nr}; %((amp_pli_dist_save{cond,nr}  - mean([amp_pli_dist_save{cond,nr}],2))./mean([amp_pli_dist_save{cond,nr}],2) )*100;
                        HdistSeg{nr} = H_dist_save{cond,nr}(1,:);
                        tends(nr) = timevec{cond,nr}(end)-timevec{cond,nr}(1);
                    end
                    tendtot = sum(tends(:));
                    relativePhiCol = [relativePhi{:,:}];
                    segLCol = [segL{:,:}]; %segLCol = log10(segLCol);
                    ampSegCol = [ampSeg{:,:}];
                    HdistSegCol = [HdistSeg{:,:}];
                    
%                     for i = 1:3
%                     scatter(segLCol,ampSegCol(i,:)); hold on
%                     [yCalc b Rsq] = linregress(segLCol(~isnan(ampSegCol(i,:)))',ampSegCol(i,~isnan(ampSegCol(i,:)))');
%                     plot(segLCol(~isnan(ampSegCol(i,:))),yCalc)
%                     end
                    % Plots of rel phase vs seg length
                    phishift = findAmpPhi(segLCol,relativePhiCol,QX);
                    f(1) = figure(1*10 + 4); % PA vs SegL
                    f(2) = figure((2)*10 +1); set(gcf,'Position',[263 517 1403 352])
                    subplot(1,3,cond);
                    [L4(cond) segLbin(:,nr,cond) binmid pdense(1,cond)] = plot_PA_Dep_relation_v6(cond,phishift,segLCol,...
                        'Segment Length',QX,linspace(0.1,1.2,QYN),f,tendtot,R.PA.interpolgrid); % ,logspace(log10(0.5),log10(1.2),6)
%                     xlim([-3.5 3.5]); %ylim([0 1.75]);
                    figure(f(2)); set(gcf,'Position',[263 517 1403 352]);  if cond>0; colorbar; end; title(R.condname{cond})
                    caxis([0 0.06]);
                    
                    % Plots of rel phase vs M1 HB Amp
                    phishift = findAmpPhi(ampSegCol(1,:),relativePhiCol,QX);
                    f(1) = figure(1*10 + 5); % PA vs M1 Amp
                    f(2) = figure((2)*10 +2);
                    subplot(1,3,cond);
                    [L5(cond) ampbin(:,nr,cond,1) binmid pdense(2,cond)] = plot_PA_Dep_relation_v6(cond,phishift,ampSegCol(1,:),...
                        'M1 High Beta Amp',QX,linspace(-100,300,QYN),f,tendtot,R.PA.interpolgrid); %linspace(0,7,8) ,linspace(-100,200,QY)
                    figure(f(2)); set(gcf,'Position',[263 517 1403 352]); if cond>0; colorbar; end; title(R.condname{cond})
                    caxis([0 0.06]);
                    
                    % Plots of rel phase vs STN HB Amp
                    phishift = findAmpPhi(ampSegCol(2,:),relativePhiCol,QX);
                    f(1) = figure(1*10 + 6); % PA vs STN HB
                    f(2) = figure((2)*10 +3);
                    subplot(1,3,cond)
                    [L6(cond) ampbin(:,nr,cond,2) binmid pdense(3,cond)] = plot_PA_Dep_relation_v6(cond,phishift,ampSegCol(2,:),...
                        'STN High Beta Amp',QX,linspace(-100,300,QYN),f,tendtot,R.PA.interpolgrid); %linspace(0,7,8) linspace(-100,200,QY)
                    figure(f(2)); set(gcf,'Position',[263 517 1403 352]);  if cond>0; colorbar; end; title(R.condname{cond})
                    caxis([0 0.06]);
                    
                    % Plots of rel phase vs STN LB Amp
                    phishift = findAmpPhi(ampSegCol(3,:),relativePhiCol,QX);
                    f(1) = figure(1*10 + 7); % PA vs STN LB
                    f(2) = figure((2)*10 +4); set(gcf,'Position',[263 517 1403 352])
                    subplot(1,3,cond)
                    [L7(cond) ampbin(:,nr,cond,3) binmid pdense(4,cond)] = plot_PA_Dep_relation_v6(cond,phishift,ampSegCol(3,:),...
                        'STN Low Beta Amp',QX,linspace(-100,300,QYN),f,tendtot,R.PA.interpolgrid); % linspace(0,4,8) linspace(-100,200,QY)
                    figure(f(2)); set(gcf,'Position',[263 517 1403 352]);  if cond>0; colorbar; end; title(R.condname{cond})
                    caxis([0 0.06]);
                    % Plots of rel phase vs Segment Length DDT
                    figure(100)
                    [H1(cond) hdist(:,cond,1)] = plot_segL_histogram(relativePhiCol,segLCol,'Segment Length',cond);
                    
                    %Plots of rel phase vs Segment Length DDT
                    f(1) = figure(1*10 + 8); % PA vs SegL
                    f(2) = figure((2)*10 +5);
                    subplot(1,3,cond)
                    set(gcf,'Position',[263 517 1403 352])
                    %                 [H2(cond) hdist(:,cond,2)] = plot_segL_histogram(pA_pli_dist_save{cond,nr},segL_pli_dist_save{cond,nr},'PLI Segment Length',cond);
                    [L6(cond) hbin(:,nr,cond,1) binmid pdense(5,cond)] = plot_PA_Dep_relation_v6(cond,relativePhiCol,HdistSegCol,...
                        'STN CTX H Beta Amp Corr',QX,linspace(-1.1,1.1,QYN),f,tendtot,R.PA.interpolgrid); % linspace(0,4,8)
                    caxis([0 0.2]);
                    %% Phase Angle Rose DDT
                    figure(200)
                    pA_dist = [pA_pli_dist_save{cond,nr}];
                    pA_dist(isnan(pA_dist)) = [];
                    R_1(cond) = polarhistogram(pA_dist,18,'FaceAlpha',0.75); hold on
                    phaseAng_dist(cond).circ_mean_std = [circ_mean(pA_dist') circ_var(pA_dist')];
                    [h mu] =circ_mtest(pA_dist',circ_mean(pA_dist'));
                    phaseAng_dist(cond).circ_meantest = [h mu];
                    %                     [pval, z] = circ_rtest(pA_dist');
                    
                    phaseAng_dist(cond).circ_RayTest = [1 1]; %[pval, z];
                    
                    ON = pA_pli_dist_save{1};  OFF = pA_pli_dist_save{2}; OFF(isnan(OFF)) = []; ON(isnan(ON)) = [];
                    sampsz = min([length(ON) length(OFF)]);
                    ON = ON(round(randbetween(1,length(ON),1,sampsz)));
                    OFF = OFF(round(randbetween(1,length(OFF),1,sampsz)));
                    [pval, table] =circ_wwtest(ON,OFF);
                    phaseAng_dist(3).Ftest = table{2,5};
                    phaseAng_dist(3).pval = pval; disp(pval);
                    
                    title('Phase Angle Distribution Phase Lock','FontSize',18)
                    %%% CUT OUT HERE
                    %%% PASTE HERE
                    %% correlation
                    %                     figure;
                    %                     x =  amp_dist_save{1,nr}; y = segL_dist_save{1,nr};
                    %                     scatter(x(3,:),y);hold on
                    %                     x =  amp_dist_save{2,nr}; y = segL_dist_save{2,nr};
                    %                     scatter(x(3,:),y);hold on
                    figure(252)
                    y =  HdistSegCol; x = segLCol;
                    scatter(x,y(1,:));hold on;[r1 p1] = corrcoef(x,y(1,:)); length_amp_corr(:,cond,sub) = [r1(2) p1(2)];
                    xlabel('Segment Length'); ylabel('Amplitude Correlation')
                    
                end
                legend(H1,R.condname)
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
                    
%                     Nbed = zeros(size(N)+1);
%                     Nbed(1:end-1,1:end-1) = N;
%                     Nbed(Nbed==0) = NaN;
                     pcolor(x,y,N); xlim([-pi pi]); ylim([y(1) y(end)])
%                      contourf(x,y,N); xlim([-pi pi]); ylim([y(1) y(end)])
%                     pcolor(x,y,N');
                    % imagesc(x,y,av1');
                    caxis([-0.05 0.05]); title(R.condname{1});
                    xlabel('Phi_1 - Phi_2'); ylabel(analynames{i}); set(gca,'YDir','normal');
                    h = colorbar; a = ylabel(h, 'Occurence (s^{-1})');
                    set(a,'rotation',270); set(a,'Position',get(a,'Position') + [0.8 0 0]);
                    title('OFF - ON')
                end
                
                
%                         legend(GH,R.condname);
                %                     pt = [R.datapathr R.subname{sub} '\images\seganalysis\'];
                %                     eval(['!del /F /Q /S C:' pt(3:end) '*'])
%                         savefigure_v2([R.datapathr 'results\images\seganalysis\' R.subname{sub} '\' R.bandname{band} '\'],['PLIseg_analysis_ONvsOFF_MaxCohM1_' R.bandname{band}],[],[],[]);
                close all
                % Now shift relative to OFF condition
                segLsave{sub,side} = segLbin;
                
                for i = 1:5 % 6 for granger
                    if i<6
                        [pdense(i,2).shiftN, pdense(i,1).shiftN] = circshift2centre_array(pdense(i,2).N,pdense(i,1).N);
                    else
                        pdense(i,2).shiftN = pdense(i,2).N;
                        pdense(i,1).shiftN = pdense(i,1).N;
                    end
                end
                densesave{sub,side} = pdense;
                hdistsave{sub,side} = hdist;
                phaseAng_dist_save{sub,side} = phaseAng_dist;
                %         gc_dist_sub_save{sub,side} = gc_dist_c;
                %         gc_cd_dist_sub_save{sub,side} =gc_dist_cd;
                clear segLbin ampbin pdense hdist phaseAng_dist gc_dist_c gc_dist_cd
            end
        end
    end
    
    save([R.datapathr '\results\seganalysis\groupseganaly_' R.bandname{band}],'segLsave','amp_pli_dist_save','densesave','hdistsave','phaseAng_dist_save'); %,'gc_dist_sub_save','GC_stat_table')
    
    ON = sum(squeeze(length_amp_corr(:,1,:,:)),3);
    OFF = sum(squeeze(length_amp_corr(:,2,:,:)),3);
    [h p] = ttest2(ON(1,:),OFF(1,:))
    mean(ON(1,:))
    mean(OFF(1,:))
end






%%%%%%%%%%%%%%%
%         if cond == 2
%             pval = circ_cmtest(pA_dist_save{1}, pA_dist_save{2});
%             annotation(gcf,'textbox',...
%                 [0.71 0.0690 .250 .16],...
%                 'String',{'Circular CM test',sprintf('P: %.3f',pval)},...
%                 'FitBoxToText','off');
%         end

%% Phase Angle Rose - PLI
%                 figure(400)
%                 pA_pli_dist = [pA_pli_dist_save{cond,nr}];
%                 pA_pli_dist(isnan(pA_pli_dist)) = [];
%                 R_2(cond) = polarhistogram(pA_pli_dist,18); hold on
%                 %             x = get(h,'Xdata');
%                 %             y = get(h,'Ydata');
%                 %             g=patch(x,y,'r');
%                 title('Phase Angle Distribution PLI')
%         if cond == 2
%
%             pval = circ_cmtest(pA_pli_dist_save{1}, pA_pli_dist_save{2});
%             annotation(gcf,'textbox',...
%                 [0.71 0.0690 .250 .16],...
%                 'String',{'Circular CM test',sprintf('P: %.3f',pval)},...
%                 'FitBoxToText','off');
%         end
%% Granger Histograms
%                 gc_dist = [gc_dist_save{cond,nr}];
%                 gc_dist = gc_dist(7:9,:);
%                 gc_dist_c_fs{cond,side,sub} = gc_dist;
%                 gc_dist_c{cond} = gc_dist(1,gc_dist(2,:)<0.05);
%                 gc_dist_cd{cond} = gc_dist(3,gc_dist(2,:)<0.05);
%                 figure(238);
%                 GH(cond) = histogram(gc_dist_c{cond},-.5:0.025:.5,'Normalization','probability'); hold on
%                 xlabel('G-Causality'); ylabel('P(x)');
%                 figure(239);
%                 GH(cond) = histogram(gc_dist_cd{cond},0:0.025:.5,'Normalization','probability'); hold on
%                 xlabel('Seths Causal Density'); ylabel('P(x)');
% %                 [f xi] = ksdensity(gc_dist_c{cond},@normpdf);
% %                 plot(xi,f)
%                 if cond == 2
%                      figure(238);
%                     [h,p(1),ci,stats1] = ttest2(gc_dist_c{1},gc_dist_c{2});
%                     [h,p(2),ci,stats2] = vartest2(gc_dist_c{1},gc_dist_c{2});
%                     GC_stat_table{side,sub} = [mean(gc_dist_c{1}),std(gc_dist_c{1}),mean(gc_dist_c{2}),std(gc_dist_c{2}),...
%                         stats1.df,stats1.tstat,p(1),stats2.fstat,stats2.df1,stats2.df2,p(2)];
%                     annotation(gcf,'textbox',...
%                         [0.1590    0.7870    0.2790    0.0960],...
%                         'String',{sprintf('2 Samp t-test P: %.3f, 2 Samp f-test P: %.3f',p(1),p(2))},...
%                         'FitBoxToText','off','LineStyle','none');
%                      figure(239);
%                     [h,p(1),ci,stats1] = ttest2(gc_dist_cd{1},gc_dist_cd{2});
%                     [h,p(2),ci,stats2] = vartest2(gc_dist_cd{1},gc_dist_cd{2});
%                     GC_cd_stat_table{side,sub} = [mean(gc_dist_cd{1}),std(gc_dist_cd{1}),mean(gc_dist_cd{2}),std(gc_dist_cd{2}),...
%                         stats1.df,stats1.tstat,p(1),stats2.fstat,stats2.df1,stats2.df2,p(2)];
%                     annotation(gcf,'textbox',...
%                         [0.6200    0.6020    0.2790    0.0960],...
%                         'String',{sprintf('2 Samp t-test P: %.3f, 2 Samp f-test P: %.3f',p(1),p(2))},...
%                         'FitBoxToText','off','LineStyle','none');
%                 end
%
%                 f(1) = figure(421); % PA vs M1 Amp
%                 f(2) = figure(422);
%                 subplot(1,3,cond);
%                 gcdir = gc_dist_save{cond,nr}(7,:); gcpval = gc_dist_save{cond,nr}(8,:); pA = pA_dist_save{cond,nr};
%                 gcdir(gcpval>0.05) = []; pA(gcpval>0.05) = [];
%                 [L5(cond) gcbin(:,cond,1) binmid pdense(6,cond)] = plot_PA_Dep_relation_v6(cond,pA,gcdir,'Causal Density',QX,[linspace(-.3,-0.05,3) linspace(0.05,.3,3)],f,tend); %linspace(0,7,8) %
%                 xlim([-3.5 3.5]); ylim([-0.5 0.5])
%                  set(f(2),'Position',[263 517 1403 352]); caxis([0 0.08]); if cond>0; colorbar; end; title(R.condname{cond})
