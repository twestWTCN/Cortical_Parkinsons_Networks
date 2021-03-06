function plot_phase_amp_analysis_PLIs(R)
idd = '';
% Collate the points and then do paired ttest between bins? A lot of
% information is lost by cutting off in histograms. Either that or somehow
% increase the number of samples - perhaps include all of the STN channels?
% Maybe combine left and right into one histogram? Collapse all of one
% subjects sync segs into one histogram?? All should be samplign from the
% same distribution so build up better resolved image of underlying 2D PDF.
% Or look at doing a hotelling test between the collated data??

close all
QX = 6 ; QY = 8;
% QX = [-pi -pi/8 pi/8 pi]';
QX = linspace(-pi,pi,QX);
% load([R.datapathr 'subject_hbWPLI075'])
% subscreen = squeeze(sum(subject_hbcohscreen>R.PA.WPLIscreen)==2);
a = logspace(0.5,log10(150),4); logscalez = [-100 -50 -15 15 50 100]; %-40:10:40; %linspace(-50,50,QY); %logscalez = [-3.^(4:-1:2) 3.^(2:4)];
for band = 2%; %1:numel(R.bandname)
    for sub = 1:numel(R.subname)
        for side =1:2
            if 1==1 %subscreen(side,sub)
                load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\' idd '_ROIvoxel_phaseamp_PLI_' R.ipsicon '_' R.siden{side} '_' R.bandname{band} '.mat'])
                for cond = 1:2
                    [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
                    for nr = 1:nrep
                        
                        amp_dist_fs{cond,side,sub} = amp_pli_dist_save{cond,nr};
                        segL_dist_fs{cond,side,sub} = segL_pli_dist_save{cond,nr};
                        relphase = wrapToPi(pA_pli_dist_save{cond,nr}-circ_mean([pA_pli_dist_save{:,:}]'));
                        f(1) = figure(1*10 + 4); % PA vs SegL
                        f(2) = figure((2)*10 +1); set(gcf,'Position',[263 517 1403 352])
                        subplot(1,3,cond);
                        tend = timevec{cond,nr}(end)-timevec{cond,nr}(1);
                        
                        [L4(cond) segLbin(:,nr,cond) binmid pdense(1,nr,cond)] = plot_PA_Dep_relation(cond,relphase,median(segL_pli_dist_save{cond,nr})-segL_pli_dist_save{cond,nr},...
                            'Segment Length',QX,linspace(-0.75,0.75,QY),f,tend); % ,logspace(log10(0.5),log10(1.2),6)
                        xlim([-3.5 3.5]); %ylim([0 1.75]);
                        figure(f(2)); set(gcf,'Position',[263 517 1403 352]);  if cond>0; colorbar; end; title(R.condname{cond})
                        caxis([0 0.04]);
                        
                        f(1) = figure(1*10 + 5); % PA vs M1 Amp
                        f(2) = figure((2)*10 +2);
                        subplot(1,3,cond);
                        [L5(cond) ampbin(:,nr,cond,1) binmid pdense(2,nr,cond)] = plot_PA_Dep_relation(cond,relphase,amp_pli_dist_save{cond,nr}(1,:),...
                            'M1 High Beta Amp',QX,linspace(-50,50,QY),f,tend); %linspace(0,7,8) ,linspace(-100,200,QY)
                        xlim([-3.5 3.5]); %ylim([0 8])
                        figure(f(2)); set(gcf,'Position',[263 517 1403 352]); if cond>0; colorbar; end; title(R.condname{cond})
                        caxis([0 0.04]);
                        
                        f(1) = figure(1*10 + 6); % PA vs SegL
                        f(2) = figure((2)*10 +3);
                        subplot(1,3,cond)
                        [L6(cond) ampbin(:,nr,cond,2) binmid pdense(3,nr,cond)] = plot_PA_Dep_relation(cond,relphase,amp_pli_dist_save{cond,nr}(2,:),...
                            'STN High Beta Amp',QX,linspace(-75,75,QY),f,tend); %linspace(0,7,8) linspace(-100,200,QY)
                        xlim([-3.5 3.5]); %ylim([0 8]);
                        figure(f(2)); set(gcf,'Position',[263 517 1403 352]);  if cond>0; colorbar; end; title(R.condname{cond})
                        caxis([0 0.04]);
                        
                        f(1) = figure(1*10 + 7); % PA vs SegL
                        f(2) = figure((2)*10 +4); set(gcf,'Position',[263 517 1403 352])
                        subplot(1,3,cond)
                        [L7(cond) ampbin(:,nr,cond,3) binmid pdense(4,nr,cond)] = plot_PA_Dep_relation(cond,relphase,amp_pli_dist_save{cond,nr}(3,:),...
                            'STN Low Beta Amp',QX,linspace(-20,20,QY),f,tend); % linspace(0,4,8) linspace(-100,200,QY)
                        xlim([-3.5 3.5]); %ylim([0 5]);
                        figure(f(2)); set(gcf,'Position',[263 517 1403 352]);  if cond>0; colorbar; end; title(R.condname{cond})
                        caxis([0 0.0]);
                        %% Segment Length DDT
                        figure(100)
                        [H1(cond) hdist(:,nr,cond,1)] = plot_segL_histogram(pA_pli_dist_save{cond,nr},segL_pli_dist_save{cond,nr},'Segment Length',cond);
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
                        %         if cond == 2
                        %             pval = circ_cmtest(pA_dist_save{1}, pA_dist_save{2});
                        %             annotation(gcf,'textbox',...
                        %                 [0.71 0.0690 .250 .16],...
                        %                 'String',{'Circular CM test',sprintf('P: %.3f',pval)},...
                        %                 'FitBoxToText','off');
                        %         end
                        %% CTX Cross Correlations
                        f(1) = figure(300); % PA vs SegL
                        f(2) = figure(301);
                        subplot(1,3,cond)
                        set(gcf,'Position',[263 517 1403 352])
                        %                 [H2(cond) hdist(:,cond,2)] = plot_segL_histogram(pA_pli_dist_save{cond,nr},segL_pli_dist_save{cond,nr},'PLI Segment Length',cond);
                        [L6(cond) hbin(:,nr,cond,1) binmid pdense(5,nr,cond)] = plot_PA_Dep_relation(cond,relphase,H_dist_save{cond,nr}(1,:),...
                            'STN CTX H Beta Amp Corr',QX,linspace(0.25,1,QY),f,tend); % linspace(0,4,8)
                        caxis([0 0.08]);
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
                        %                 [L5(cond) gcbin(:,cond,1) binmid pdense(6,cond)] = plot_PA_Dep_relation(cond,pA,gcdir,'Causal Density',QX,[linspace(-.3,-0.05,3) linspace(0.05,.3,3)],f,tend); %linspace(0,7,8) %
                        %                 xlim([-3.5 3.5]); ylim([-0.5 0.5])
                        %                  set(f(2),'Position',[263 517 1403 352]); caxis([0 0.08]); if cond>0; colorbar; end; title(R.condname{cond})
                        %% correlation
                        %                     figure;
                        %                     x =  amp_dist_save{1,nr}; y = segL_dist_save{1,nr};
                        %                     scatter(x(3,:),y);hold on
                        %                     x =  amp_dist_save{2,nr}; y = segL_dist_save{2,nr};
                        %                     scatter(x(3,:),y);hold on
                        figure(252)
                        y =  H_dist_save{cond,nr}; x = segL_pli_dist_save{cond,nr};
                        scatter(x,y(1,:));hold on;[r1 p1] = corrcoef(x,y(1,:)); length_amp_corr(:,cond,sub,nr) = [r1(2) p1(2)];
                        xlabel('Segment Length'); ylabel('Amplitude Correlation')
                    end
                end
                legend(H1,R.condname)
                legend(R_1,R.condname)
                %         legend(H2,R.condname)
                %         legend(R_2,R.condname)
                legend(L5,R.condname)
                legend(L6,R.condname)
                legend(L7,R.condname)
                %         legend(GH,R.condname);
                %                     pt = [R.datapathr R.subname{sub} '\images\seganalysis\'];
                %                     eval(['!del /F /Q /S C:' pt(3:end) '*'])
                %                     savefigure_v2([R.datapathr R.subname{sub} '\images\seganalysis\PLI\'],['PLIseg_analysis_ONvsOFF_MaxCohM1 '],[],[],[]);
                close all
                % Now shift relative to OFF condition
                segLsave{sub,side} = segLbin;
                
                %             for i = 1:3
                %                 x = squeeze(ampbin(:,:,:,i));
                %                 ampbin(:,:,:,i) = circshift2centre(x',2)';
                %             end
                %             ampsave{sub,side} = ampbin;
                for i = 1:5 % 6 for granger
                    for nr = 1:size(pdense,2)
                        if i>6
                            [pdense(i,nr,2).shiftN, pdense(i,nr,1).shiftN] = circshift2centre_array(pdense(i,nr,2).N,pdense(i,nr,1).N);
                        else
                            pdense(i,nr,2).shiftN = pdense(i,nr,2).N;
                            pdense(i,nr,1).shiftN = pdense(i,nr,1).N;
                        end
                    end
                end
                densesave{sub,side} = pdense;
                hdistsave{sub,side} = hdist;
                phaseAng_dist_save{sub,side} = phaseAng_dist;
                %         gc_dist_sub_save{sub,side} = gc_dist_c;
                %         gc_cd_dist_sub_save{sub,side} =gc_dist_cd;
            end
            clear segLbin ampbin pdense hdist phaseAng_dist gc_dist_c gc_dist_cd
        end
    end
    save([R.datapathr '\results\seganalysis\groupseganaly'],'segLsave','amp_pli_dist_save','densesave','hdistsave','phaseAng_dist_save'); %,'gc_dist_sub_save','GC_stat_table')
    
    ON = sum(squeeze(length_amp_corr(:,1,:,:)),3);
    OFF = sum(squeeze(length_amp_corr(:,2,:,:)),3);
    [h p] = ttest2(ON(1,:),OFF(1,:))
    mean(ON(1,:))
    mean(OFF(1,:))
end