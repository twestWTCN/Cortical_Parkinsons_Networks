close all; clear
datapathr = 'C:\Users\Tim\Documents\Work\Cortical_Networks\Data\';
subname = {'JN','MC','SW','DF','JB','MW','DP','DS','JA'};
condname = {'ON','OFF'};
sidenm = {'RSTN','LSTN'}; ipsicon = 'ipsi';

load([datapathr '\results\seganalysis\groupseganaly'])
% CORRELATION BETWEEN BETA POWER AND DURATION OF SEGMENTS
cmap = [0 0 1; 1 0 0];
analynames = {'Segment Length','M1 High Beta Amp','STN High Beta Amp','STN Low Beta Amp','Causal Density'};
%Phase Amp Analy
for k = 1:4 %5 for granger
i = 0;
for sub = [1:9]
    for side = 1:2
        for cond = 1:2
            denseav(:,:,cond,side,sub) = densesave{sub,side}(k,cond).shiftN;
        end
    end
end
dimz = size(densesave{sub,side}(k,cond).N);
x = densesave{sub,side}(k,cond).Xedges;
y = densesave{sub,side}(k,cond).Yedges;

A = reshape(squeeze(denseav(:,:,1,:,:)),dimz(1),dimz(2),2*numel(subname));
B = reshape(squeeze(denseav(:,:,2,:,:)),dimz(1),dimz(2),2*numel(subname));
% for i = 1:size(A,1)
%     for j = 1:size(B,2)
%         [dum pmat(i,j)] = ttest2(squeeze(A(i,j,:)),squeeze(B(i,j,:)));
%     end
% end
% Change plotting to run through clusters instead of simply masking!!
stat = clusterstat_matrix(B,A,1000,dimz);
statmat = squeeze(stat.stat).*squeeze(stat.mask);
av1 = squeeze(nanmean(nanmean(denseav(:,:,1,:,:),4),5));
av2 = squeeze(nanmean(nanmean(denseav(:,:,2,:,:),4),5));
avdiff = av1-av2;

figure(k); set(gcf,'Position',[263 517 1403 352])

subplot(1,3,1);

N = av1;
imagesc(x,y,N);
xlabel('Phi_1 - Phi_2')
set(gca,'YDir','normal')
hold on

% imagesc(x,y,av1');
%caxis([0 0.04]); title(condname{1});
xlabel('Phi_1 - Phi_2'); ylabel(analynames{k}); set(gca,'YDir','normal');
h = colorbar; a = ylabel(h, 'Occurence (s^{-1})');
set(a,'rotation',270); set(a,'Position',get(a,'Position') + [0.8 0 0]);

subplot(1,3,2);
N = av2;
imagesc(x,y,N);
xlabel('Phi_1 - Phi_2')
set(gca,'YDir','normal')
hold on

% imagesc(x,y,av2');
%caxis([0 0.04]); title(condname{2});
xlabel('Phi_1 - Phi_2'); ylabel(analynames{k}); set(gca,'YDir','normal');
h = colorbar; a = ylabel(h, 'Occurence (s^{-1})'); 
set(a,'rotation',270); set(a,'Position',get(a,'Position') + [0.8 0 0]);
subplot(1,3,3);
N = statmat;
imagesc(x,y,N);
xlabel('Phi_1 - Phi_2')
set(gca,'YDir','normal')
hold on 
%,'alphadata', ~(isnan(statmat')|isinf(statmat')|statmat'==0));
% imagesc(x,y,statmat', 'alphadata', ~(isnan(statmat')|isinf(statmat')|statmat'==0));
caxis([-3 3]); title('OFF-ON t-stat'); xlabel('Phi_1 - Phi_2'); ylabel(analynames{k}); h = colorbar; a = ylabel(h, 't-stat');
set(a,'rotation',270); set(a,'Position',get(a,'Position') + [0.8 0 0]);
set(gca,'YDir','normal')

% figure(k); subplot(1,3,1); imagesc(x,y,av1'); caxis([0 0.1]); subplot(1,3,2); imagesc(x,y,av2'); caxis([0 0.1]); subplot(1,3,3); imagesc(x,y,avdiff'); caxis([-0.2 0.2])
end
 savefigure_v2([datapathr '\results\seganalysis\'],['PhaseAng_group_seg_analysis_ONvsOFF_MaxCohM1 '],[],[],[]); close all

 figure(777);
 x = logspace(-1.7,0.3,50);
 for sub = [1:9]
    for side = 1:2
        for cond = 1:2
            distpar = hdistsave{sub,side}(:,cond,1);
            distparsave(:,side,sub,cond) = distpar;
        y = lognpdf(x,distpar(1),distpar(2));
        plot(x,y,'color',cmap(cond,:)); hold on
        end
    end
 end
 meanpar(:,1) = mean(reshape(squeeze(distparsave(:,:,:,1)),2,numel(subname)*2),2);
 meanpar(:,2) = mean(reshape(squeeze(distparsave(:,:,:,2)),2,numel(subname)*2),2);
        y = lognpdf(x,meanpar(1,1),meanpar(2,1));
        h(1) = plot(x,y,'color',cmap(1,:),'LineWidth',4); 
          y = lognpdf(x,meanpar(1,2),meanpar(2,2));      
 h(2) = plot(x,y,'color',cmap(2,:),'LineWidth',4);
 legend(h,condname)
  set(gca,'xscale','log'); grid on;
  xlabel('Segment Duration (s)','FontSize',14); ylabel('P(x)','FontSize',14); title('LogN Distributions of Segment Duration','FontSize',16)
  savefigure_v2([datapathr '\results\seganalysis\'],['PhaseAng_group_seg_analysis_ONvsOFF_segduration'],[],[],[]); close all

  figure(888);
 for sub = [1:9]
    for side = 1:2
        for cond = 1:2
            ppar = phaseAng_dist_save{sub,side}(cond).circ_mean_std;
            phapar(:,side,sub,cond) = ppar;
            pdiff(sub,side) = ppar(2)-ppar(1);
            if cond ==1; mid =0; else mid = 0; end
            y = circ_vmpdf(-pi:.1:pi,mid,1/ppar(2));
            plot(-pi:.1:pi,y,'color',cmap(cond,:)); hold on
        end
    end
 end
 
 meanpar(:,1) = mean(reshape(squeeze(phapar(:,:,:,1)),2,numel(subname)*2),2);
 y = circ_vmpdf(-pi:.1:pi, 0,1/ meanpar(2,1));
 h(1) = plot(-pi:.1:pi,y,'color',cmap(1,:),'LineWidth',4); hold on
 meanpar(:,2) = mean(reshape(squeeze(phapar(:,:,:,2)),2,numel(subname)*2),2);
 y = circ_vmpdf(-pi:.1:pi, 0,1/ meanpar(2,2));
 h(2) = plot(-pi:.1:pi,y,'color',cmap(2,:),'LineWidth',4); hold on
 legend(h,condname)
 grid on;
  xlabel('Relative Phase','FontSize',14); ylabel('P(x)','FontSize',14); title('von-Mises Distributions of Relative Phase','FontSize',16)

  [p h] = circ_vtest(pdiff(:)',mean(pdiff(:)));
  [p h] = circ_rtest(pdiff(:)');
  
  savefigure_v2([datapathr '\results\seganalysis\'],['PhaseAng_group_seg_analysis_ONvsOFF_phaseVMisses'],[],[],[]); close all

 
 return
 i = 0;
for sub = 1:4
    for side = 1:2
        i = i+1;
        for cond = 1:2
            figure(100)
            hold on
            plot(binmid,segLsave{sub,side}(:,cond),'color',cmap(cond,:))
            x1(:,i,cond) = segLsave{sub,side}(:,cond);
            figure(200)
            hold on
            plot(binmid,ampsave{sub,side}(:,cond,1),'color',cmap(cond,:))
            x2(:,i,cond) = ampsave{sub,side}(:,cond,1);
            figure(300)
            hold on
            plot(binmid,ampsave{sub,side}(:,cond,2),'color',cmap(cond,:))
            x3(:,i,cond) = ampsave{sub,side}(:,cond,2);
            figure(400)
            hold on
            plot(binmid,ampsave{sub,side}(:,cond,3),'color',cmap(cond,:))
            x4(:,i,cond) = ampsave{sub,side}(:,cond,3);
            figure(500)
            
        end
    end
end

for cond = 1:2
    figure(100); grid on
    plot(binmid,squeeze(median(x1(:,:,cond),2)),'color',cmap(cond,:),'LineWidth',3)
    xlabel('Relative Phase Angle'); ylabel('Segment Length')
    figure(200); grid on
    plot(binmid,squeeze(median(x2(:,:,cond),2)),'color',cmap(cond,:),'LineWidth',3);
    xlabel('Relative Phase Angle'); ylabel('M1 High Beta')
    figure(300); grid on
    plot(binmid,squeeze(median(x3(:,:,cond),2)),'color',cmap(cond,:),'LineWidth',3)
    xlabel('Relative Phase Angle'); ylabel('STN High Beta')
    figure(400); grid on
    plot(binmid,squeeze(median(x4(:,:,cond),2)),'color',cmap(cond,:),'LineWidth',3)
    xlabel('Relative Phase Angle'); ylabel('STN Low Beta')
    legend({'ON','OFF'})
end

savefigure_v2([datapathr '\results\seganalysis\'],['derPhaseAng_group_seg_analysis_ONvsOFF_MaxCohM1 '],[],[],[]);