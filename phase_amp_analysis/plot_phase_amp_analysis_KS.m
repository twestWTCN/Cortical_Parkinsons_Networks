%function plot_phase_amp_analysis
clear; close all

datapathr = 'C:\Users\Tim\Documents\Work\Cortical_Networks\Data\';
% subname = {'JN','MC','SW'}
subname = {'JN','MC','SW','DF','JB','MW','DP','DS','JA'};
condname = {'ON','OFF'};
sidenm = {'RSTN','LSTN'}; ipsicon = 'ipsi';
QY = 7; QX = 7;
BWY = [0.05 0.25]; BWX = [pi/7];
for sub = 1:numel(subname)
    for side =1:2
        load([datapathr subname{sub} '\ftdata\ROI_analy\ROIvoxel_phaseamp_' ipsicon '_' sidenm{side}])
        for cond = 1:2
            [datafileN,pp_mark,nrep,senscheck] = data_fileguide(subname{sub},cond-1);
            for nr =nrep
                
                amp_dist_fs{cond,side,sub} = amp_dist_save{cond,nr};
                segL_dist_fs{cond,side,sub} = segL_dist_save{cond,nr};
                f(1) = figure(1*10 + 4); % PA vs SegL
                f(2) = figure((2)*10 +1); set(gcf,'Position',[263 517 1403 352])
                subplot(1,3,cond);
                tend = timevec{cond,nr}(end)-timevec{cond,nr}(1);
                
                [L4(cond) segLbin(:,cond) binmid pdense(1,cond)] = plot_PA_Dep_relation_KS(cond,pA_dist_save{cond,nr},segL_dist_save{cond,nr},'Segment Length',QX,linspace(0.2,0.75,QY),f,tend,[BWX BWY(1)]); %linspace(0,1.5,8)
                xlim([-3.5 3.5]); %ylim([0 1.75]);
                figure(f(2)); set(gcf,'Position',[263 517 1403 352]);% caxis([0 0.005]); 
                if cond>0; colorbar; end; title(condname{cond})
                
                 f(1) = figure(1*10 + 5); % PA vs M1 Amp
                f(2) = figure((2)*10 +2);
                subplot(1,3,cond);
                [L5(cond) ampbin(:,cond,1) binmid pdense(2,cond)] = plot_PA_Dep_relation_KS(cond,pA_dist_save{cond,nr},amp_dist_save{cond,nr}(1,:),'M1 High Beta Amp',QX,linspace(0.25,5.5,QY),f,tend,[BWX BWY(2)]); %linspace(0,7,8)
                xlim([-3.5 3.5]); %ylim([0 8])
                figure(f(2)); set(gcf,'Position',[263 517 1403 352]);% caxis([0 0.005]); 
                if cond>0; colorbar; end; title(condname{cond})
                
                f(1) = figure(1*10 + 6); % PA vs SegL
                f(2) = figure((2)*10 +3);
                subplot(1,3,cond)
                [L6(cond) ampbin(:,cond,2) binmid pdense(3,cond)] = plot_PA_Dep_relation_KS(cond,pA_dist_save{cond,nr},amp_dist_save{cond,nr}(2,:),'STN High Beta Amp',QX,linspace(0.25,5.5,QY),f,tend,[BWX BWY(2)]); %linspace(0,7,8)
                xlim([-3.5 3.5]); %ylim([0 8]);
                figure(f(2)); set(gcf,'Position',[263 517 1403 352]);% caxis([0 0.005]); 
                if cond>0; colorbar; end; title(condname{cond})
                
                f(1) = figure(1*10 + 7); % PA vs SegL
                f(2) = figure((2)*10 +4); set(gcf,'Position',[263 517 1403 352])
                subplot(1,3,cond)
                [L7(cond) ampbin(:,cond,3) binmid pdense(4,cond)] = plot_PA_Dep_relation_KS(cond,pA_dist_save{cond,nr},amp_dist_save{cond,nr}(3,:),'STN Low Beta Amp',QX,linspace(1,3,QY),f,tend,[BWX BWY(2)]); % linspace(0,4,8)
                xlim([-3.5 3.5]); %ylim([0 5]);
                figure(f(2)); set(gcf,'Position',[263 517 1403 352]); %caxis([0 0.005]); 
                if cond>0; colorbar; end; title(condname{cond})
                %% Segment Length DDT
                figure(100)
                [H1(cond) hdist(:,cond,1)] = plot_segL_histogram(pA_dist_save{cond,nr},segL_dist_save{cond,nr},'Segment Length',cond);
                %% Phase Angle Rose DDT
                figure(200)
                pA_dist = [pA_dist_save{cond,nr}];
                pA_dist(isnan(pA_dist)) = [];
                R_1(cond) = polarhistogram(pA_dist,18,'FaceAlpha',0.75); hold on
                phaseAng_dist(cond).circ_mean_std = [circ_mean(pA_dist') circ_var(pA_dist')];
                [h mu] =circ_mtest(pA_dist',circ_mean(pA_dist'));
                phaseAng_dist(cond).circ_meantest = [h mu];
                [pval, z] = circ_rtest(pA_dist');
                phaseAng_dist(cond).circ_RayTest = [pval, z];
                
                ON = pA_dist_save{1};  OFF = pA_dist_save{2}; OFF(isnan(OFF)) = []; ON(isnan(ON)) = [];
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
                %% Segment Length PLI
%                 figure(300)
%                 [H2(cond) hdist(:,cond,2)] = plot_segL_histogram(pA_pli_dist_save{cond,nr},segL_pli_dist_save{cond,nr},'PLI Segment Length',cond);
                
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
%                 [L5(cond) gcbin(:,cond,1) binmid pdense(5,cond)] = plot_PA_Dep_relation(cond,pA_dist_save{cond,nr},gc_dist_save{cond,nr}(3,:),'Causal Density',8,logspace(log10(0.05),log10(0.3),8),f,tend); %linspace(0,7,8)
%                 xlim([-3.5 3.5]); ylim([0 0.3])
%                  set(f(2),'Position',[263 517 1403 352]); caxis([0 0.08]); if cond>0; colorbar; end; title(condname{cond})
                
            end
        end
        legend(H1,condname)
        legend(R_1,condname)
%         legend(H2,condname)
%         legend(R_2,condname)
        legend(L5,condname)
        legend(L6,condname)
        legend(L7,condname)
%         legend(GH,condname);
        %         pt = [datapathr subname{sub} '\images\seganalysis\'];
        %         eval(['!del /F /Q /S C:' pt(3:end) '*'])
%         savefigure_v2([datapathr subname{sub} '\images\seganalysis\'],['seg_analysis_ONvsOFF_MaxCohM1 '],[],[],[]);
        close all
        % Now shift relative to OFF condition
        segLbin = circshift2centre(squeeze(segLbin)',2)';
        for i = 1:3
            x = squeeze(ampbin(:,:,i));
            ampbin(:,:,i) = circshift2centre(x',2)';
        end
        segLsave{sub,side} = segLbin;
        ampsave{sub,side} = ampbin;
        for i = 1:4 % granger
            [pdense(i,2).shiftN, pdense(i,1).shiftN] = circshift2centre_array(pdense(i,2).N,pdense(i,1).N);
            %             pdense(i,2).shiftN = pdense(i,2).N;
            %             pdense(i,1).shiftN = pdense(i,1).N;
        end
        densesave{sub,side} = pdense;
        hdistsave{sub,side} = hdist;
        phaseAng_dist_save{sub,side} = phaseAng_dist;
%         gc_dist_sub_save{sub,side} = gc_dist_c;
%         gc_cd_dist_sub_save{sub,side} =gc_dist_cd;
        clear segLbin ampbin pdense hdist phaseAng_dist gc_dist_c gc_dist_cd
    end
end
save([datapathr '\results\seganalysis\groupseganaly'],'segLsave','ampsave','densesave','hdistsave','phaseAng_dist_save'); %,'gc_dist_sub_save','GC_stat_table')