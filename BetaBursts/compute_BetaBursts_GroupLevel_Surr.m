function compute_BetaBursts_GroupLevel_Surr(R)
%% Main function to compute subject level BB analysis + phase syncs
% TO do:
% 1) Implement the rate normalization step to divide by recording time.
% 2) Is there some relationship between the change in the coherence and the
% the change in power?
% 3) For group level- inspect individually and classify by eye.
% 4) Think of the amplitude by duration plots as a cone - how can you fit
% to data and examine the difference? You need to fit an arc sector
if nargin<2
    R = makeHeader_SubCort_Cort_Networks();
end
%%%
close all

% % Phase randomize surrogate
% surrtag = 'surrogates'; %
% R.condname{3} = 'Surr_{phi}';

% % Permuted data surrogate
surrtag = 'surrogates_shuff';
R.condname{3} = 'Surr_{perm}';
% 
R.condcmap(3,:) = [0.7 0.7 0.7];
for sub = 1:length(R.subname)
    for breg = 1:length(R.bregname)
        for side = 1:2
            try
                load([R.datapathr R.subname{sub} '\ftdata\BetaBursts\BetaBurstAnalysis_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg} '_org'],'BB')
                for cond = 1:3
                    if cond == 3
                        load([R.datapathr R.subname{sub} '\ftdata\BetaBursts\BetaBurstAnalysis_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg} '_' surrtag],'BB')
                        condr = 1; % Surrogate is encoded in space of the first condition
                    else
                        condr = cond;
                    end
                    segAGroup{cond,side,sub} = BB.segA_save{condr};
                    segAPrcGroup{cond,side,sub} = BB.segAPrc_save{condr};
                    segLGroup{cond,side,sub} = BB.segL_t_save{condr};
                    segPLVGroup{cond,side,sub} = squeeze(BB.segPLV{condr}(:,3,:))';
                    recLGroup(cond,side,sub) = diff(BB.DTvec{condr}([1 end]))/60;
                    rpAGroup(:,cond,side,sub) = squeeze(BB.Amp_binRP(:,3,condr,1));
                    rpAPrcGroup(:,cond,side,sub) = squeeze(BB.AmpPrc_binRP(:,3,condr,1));
                    rpPLVGroup(:,cond,side,sub) = squeeze(BB.PLV_binRP(:,3,condr,1));
                    plvAmpPrc(:,cond,side,sub) = squeeze(BB.AmpPrc_binPLV(:,3,condr,1));
                    plvAmp(:,cond,side,sub) = squeeze(BB.AmpPrc_binPLV(:,3,condr,1));
                    
                    % Binned Amp by Dur
                    seg_binAmpGroup(:,cond,side,sub) = squeeze(BB.Seg_binAmp(:,condr,1));
                    amp_binDurGroup(:,cond,side,sub) = squeeze(BB.Amp_binDur(:,condr,1));
                    
                    % bs,band,cond,1
                    stvec = statvec(BB.segA_save,BB.segL_t_save,2);
                    regstatsGroup(:,:,side,sub) =stvec; %[stind cond side sub]
                    
                    % correlations
                    ampdurSpearman(:,cond,side,sub) = BB.corrs.ampdur.spearman{condr};
                    ampdurRegress(:,cond,side,sub) = BB.corrs.ampdur.rest_regress{condr};
                    
                    % DFA
                    psdfaGroup(:,cond,side,sub) = BB.PSDFA(:,condr);
                    
                end
            catch
                disp('load error!')
                segAGroup{cond,side,sub} = [];
                segLGroup{cond,side,sub} = [];
                segPLVGroup{cond,side,sub} = [];
                recLGroup(cond,side,sub) = 0;
                rpAGroup(:,cond,side,sub) = NaN(size(squeeze(BB.Amp_binRP(:,1,condr,1))));
                rpAPrcGroup(:,cond,side,sub) = NaN(size(squeeze(BB.AmpPrc_binRP(:,1,condr,1))));
                rpPLVGroup(:,cond,side,sub) = NaN(size(squeeze(BB.PLV_binRP(:,1,condr,1))));
                seg_binAmpGroup(:,cond,side,sub) = NaN(size(squeeze(BB.Seg_binAmp(:,condr,1))));
            end
            
            
            disp([cond side breg sub])
        end
    end
end

% (for cell formatted variables)
segAGroup = reshape(segAGroup,3,[]);
segAPrcGroup = reshape(segAPrcGroup,3,[]);
segLGroup = reshape(segLGroup,3,[]);
segPLVGroup = reshape(segPLVGroup,3,[]);

for cond = [1 2 3]; %[3 1 2]; %[1 2 3]; %
    if cond<3; condr = cond; else; condr = 3; end
    totRecL = sum(reshape(recLGroup(condr,:,:),1,[]));
    
    %% Group Histogramas
    figure(1)
    % Group Amplitude Histogram
    subplot(1,3,1)
    A = segAGroup(cond,:);
    AsegA{cond} = [A{:}];
    BAmp(:,cond) = [mean(AsegA{cond}); std(AsegA{cond}); numel(AsegA{cond})/20];
    BB.lrange.Amp = linspace(log10(1e-2),log10(BB.range.Amp(end))+0.5,numel(BB.range.Amp));
    h = plotRateDists(BB.lrange.Amp,log10(AsegA{cond}),totRecL,R.condcmap(cond,:),1);
    hl1(cond) = h(1);
    
    % Group Segment Length Histogram
    subplot(1,3,2)
    A = segLGroup(cond,:);
    AsegLGroup{cond} = [A{:}];
    BB.lrange.segDur = linspace(log10(1e1),log10(BB.range.segDur(end))+0.5,numel(BB.range.segDur));
    h = plotRateDists(BB.lrange.segDur,log10(AsegLGroup{cond}),totRecL,R.condcmap(cond,:),1);
    hl2(cond) = h(1);
    
    % Group PLV Histogram
    subplot(1,3,3)
    A = segPLVGroup(cond,:);
    AsegPLVGroup{cond} = [A{:}];
    h = plotRateDists(BB.range.PLV,AsegPLVGroup{cond},totRecL,R.condcmap(cond,:),1);
    hl3(cond) = h(1);
    %% RP Plots
    % Group RP to Amp
    A = recentreDist(squeeze(rpAGroup(:,cond,:,:)),cond);
    [R2,s,stdeqv] = vmpars(binEdge2Mid(BB.range.RP),A);
    rpAVMfit(:,:,cond) = [R2; s(3,:)]; %std(A,[],1)];
    rpAGroupData(:,:,cond) = A;
    rpAGroupStat(:,cond,1) = nanmean(A,2); rpAGroupStat(:,cond,2) = nanstd(A,1,2)./sqrt(size(A,2));
    
    % Group RP to AmpPrc
    A = recentreDist(squeeze(rpAPrcGroup(:,cond,:,:)),cond);
    [R2,s,stdeqv] = vmpars(binEdge2Mid(BB.range.RP),A);
    figure(222)
    hlvm(cond) = VmPlot(BB.range.RP,A,s,R.condcmap(cond,:),'o');

    rpAPrcVMfit(:,:,cond) = [R2; s(3,:)]; %std(A,[],1)
    rpAPrcGroupData(:,:,cond) = A;
    rpAPrcGroupStat(:,cond,1) = nanmean(A,2); rpAPrcGroupStat(:,cond,2) = nanstd(A,1,2)./sqrt(size(A,2));
    
    % Group RP to PLV
    A = recentreDist(squeeze(rpPLVGroup(:,cond,:,:)),cond);
    rpPLVGroupData(:,:,cond) = A;
    rpPLVGroupStat(:,cond,1) = nanmean(A,2); rpPLVGroupStat(:,cond,2) = nanstd(A,1,2)./sqrt(size(A,2));
    
    %% Binned Amplitude/Duration
    A = squeeze(seg_binAmpGroup(:,cond,:,:));
    A = reshape(A,size(A,1),[]); A(isnan(A)) = 0;
    seg_binAmpGroupData(:,:,cond) = A;
    seg_binAmpGroupStat(:,cond,1) = nanmean(A,2); seg_binAmpGroupStat(:,cond,2) = nanstd(A,1,2)./sqrt(size(A,2));
    
    A = squeeze(amp_binDurGroup(:,cond,:,:));
    A = reshape(A,size(A,1),[]);
    amp_binDurGroupData(:,:,cond) = A; A(isnan(A)) = 0;
    amp_binDurGroupStat(:,cond,1) = nanmean(A,2); amp_binDurGroupStat(:,cond,2) = nanstd(A,1,2)./sqrt(size(A,2));
    
    
    
    %% PLV Plots
    A = squeeze(plvAmpPrc(:,cond,:,:));
    A = reshape(A,size(A,1),size(A,2)*size(A,3));
    plvAMPprcGroupData(:,:,cond) = A;
    plvAMPprcGroupStat(:,cond,1) = nanmean(A,2); plvAMPprcGroupStat(:,cond,2) =  nanstd(A,1,2)./sqrt(size(A,2));
    % ScatterGrams and Regression of Amplitude vs Duration
    
    %     figure(4)
    %     X = reshape(squeeze(regstatsGroup(1,:,:,:)),2,[]);
    %     Xbar(:,cond) = mean(X'); Xhat = std(X');
    %     HB = superbar(1:2,Xbar,'E',Xhat)
    %
    %
    %% Amp Seg Scatters
    figure(6)
    subplot(2,2,2)
    groupPars = reshape( ampdurRegress(:,cond,:,:),6,20);
    B(2) = mean(groupPars(2,:),2); 
    B(1) = mean(groupPars(2,:),2)-std(groupPars(2,:),[],2); 
    B(3) = mean(groupPars(2,:),2)+std(groupPars(2,:),[],2);
    X = 0:1:150; 
    Y = B.*X';
    plot(X',Y(:,2),'color',R.condcmap(cond,:),'LineWidth',2);
        hold on
    plot(X',Y(:,[1 3]),'color',R.condcmap(cond,:),'LineWidth',2,'LineStyle','--');
    
    B(1) = mean(groupPars(5,:),2);
    B(2) = mean(groupPars(6,:),2);
    
    Y = B(1) + B(2).*X';
%     plot(X',Y(:,1),'color',R.condcmap(cond,:),'LineWidth',2,'LineStyle',':');
%     scatter(X(1:10:end)',Y(1:10:end,1),75,R.condcmap(cond,:),'Marker','+','LineWidth',2.5)
    s2(cond) = scatter(AsegA{cond},AsegLGroup{cond},25,R.condcmap(cond,:),'filled')
    s2(cond).MarkerFaceAlpha = 0.5;

    
end

%% Group Distributions
figure(1)
% Group Amp
subplot(1,3,1);
legend(hl1,R.condname,'Location','West'); ylabel('Frequency (min^{-1})'); xlabel(['log10 Burst Amplitude']);
% xlim(BB.lrange.Amp([1 end])); %ylim([0 12]); box off
xlim([0.25 2.25]); ylim([0 20]);
title('Burst Amplitude ')
A = cellfun(@mean,segAGroup);
statv = statvec(log10(A(1,:)),log10(A(3,:)),1);
txtlab = [sprintf('OFF-ON: %.3f [%0.2f %0.2f] ',statv(end-3:end-1)) ' (' textstar(statv(end),0.05) ')'];
TS = text(1,18, txtlab); TS.FontWeight = 'bold'; TS.FontSize = 8; TS.FontWeight = 'bold'; box off

% Group Dur
subplot(1,3,2); legend(hl2,R.condname,'Location','East'); ylabel('Frequency (min^{-1})'); xlabel('log Burst Duration (ms)');
xlim(BB.range.segDur([1 end])); % box off
xlim([1.75 3.75]); ylim([0 20]);
title('Burst Duration ')
A = cellfun(@mean,segLGroup);
statv = statvec(log10(A(1,:)),log10(A(2,:)),1);
txtlab = [sprintf('OFF-ON: %.3f [%0.2f %0.2f] ',statv(end-3:end-1)) ' (' textstar(statv(end),0.05) ')'];
TS = text(2.5,18, txtlab); TS.FontWeight = 'bold'; TS.FontSize = 8; TS.FontWeight = 'bold'; box off

% Group PLV
subplot(1,3,3); legend(hl3,R.condname,'Location','East'); ylabel('Frequency (min^{-1})'); xlabel('Burst PPC');
xlim([-0.05 BB.range.PLV(end)]); ylim([0 20]);
title('Burst PPC ')
A = cellfun(@mean,segPLVGroup);
statv =statvec(A(1,:),A(2,:),1);
txtlab = [sprintf('OFF-ON: %.3f [%0.2f %0.2f] ',statv(end-3:end-1)) ' (' textstar(statv(end),0.05) ')'];
TS = text(0.1,18, txtlab); TS.FontWeight = 'bold'; TS.FontSize = 8; TS.FontWeight = 'bold'; box off

set(gcf,'Position',[500   680  1090  292])

figure(23)
subplot(2,1,1)
p = vecTTest(seg_binAmpGroupData);
barplot160818(R,BB.binAmp,seg_binAmpGroupStat,p(2,:),0,1,10);
title('Beta Burst Duration by Amplitude')
xlabel('Burst Amplitude'); ylabel(['Burst Duration']);

subplot(2,1,2)
amp_binDurGroupData(amp_binDurGroupData==0) = NaN;
p = vecTTest(amp_binDurGroupData);
barplot160818(R,BB.binDur,amp_binDurGroupStat,p(2,:),0,1,1);
title('Beta Burst Amplitude by Duration')
xlabel('Burst Duration (ms)'); ylabel(['Burst Amplitude']);
set(gcf,'Position',[680    310   688   572])

%% Now do RP plots
figure(2)
subplot(3,1,2)
p = vecTTest(rpAGroupData);
barplot160818(R,BB.range.RP,rpAGroupStat,p(2,:),1,0,0.5);
title('Beta Burst Amplitude by Relative Phase')
xlabel('Relative Phase'); ylabel(['Abs. Burst Amplitude']); ylim([0 7])

subplot(3,1,3)
p = vecTTest(rpAPrcGroupData);
barplot160818(R,BB.range.RP,rpAPrcGroupStat,p(2,:),1,0,1);
title('% \Delta Burst Amplitude by Relative Phase')
xlabel('Relative Phase'); ylabel(['% \Delta Burst Amplitude']); ylim([0 7])

subplot(3,1,1)
p = vecTTest(plvAMPprcGroupData);
barplot160818(R,BB.range.PLV,plvAMPprcGroupStat,p(2,:),0,0,1);
title('Beta Burst Amplitude by Phase Synchronization')
xlabel('% \Delta PPC'); ylabel(['% \Delta Burst Amplitude']); ylim([-2 10])
set(gcf,'Position',[680.0000   83.5000  907.5000  894.5000])

figure(3)
p = vecTTest(rpPLVGroupData);
barplot160818(R,BB.range.RP,rpPLVGroupStat,p(2,:),1,0,0.025);
title('Burst PPC by Relative Phase')
xlabel('Relative Phase'); ylabel(['\beta_{1/2} PPC']); ylim([-20 20])

set(gcf,'Position',[680.0000   83.5000  907.5000  894.5000])
%     BB.stats.DurCond  = statvec(BB.segL_t_save{1},BB.segL_t_save{2},1);

% binN x conds x stat(mu,std)

%% PS-DFA
figure(841)
X = reshape(psdfaGroup,2,3,20);
X = permute(X,[1 3 2]);
for i=1:2
    subplot(1,2,i)
    a3WayBarPlot(R,X(i,:,:),0.1)
end

%% VM Fits
figure(5)
for i = 1:2
    if i == 1; ofs = 0.1; else; ofs = 0.1; end
    subplot(2,1,i)
    a3WayBarPlot(R,rpAPrcVMfit(i,:,:),ofs)
    if i == 1
        ylim([0 1]);ylabel('Least Squares R^2');
        title('VM Fit')
    else
        ylim([0 1.3]);  ylabel('VM Distribution Width');
        title('VM Distribution Width')
    end
end
set(gcf,'Position',[402   148   443   700])

figure(222)
xlabel('Relative Phase')
ylabel('% \Delta Burst Amplitude')
legend(hlvm,R.condname)
ylim([0 8]); xlim([-pi pi]);
title('NL Regression of Von Mises to Burst Amplitude vs Relative Phase')
set(gcf,'Position',[848   145   638   700])

figure(6)
ylist(:,1) = [0 0.3]; ylist(:,2) = [0 35];
ylist(:,3) = [0 0.5]; ylist(:,4) = [0 1]; ylist(:,5) = [0 1];
ylablist = {'Regression R^2','B Coefficient','White''s Statistic','Spearman''s R','P-value'};
titlist = {'Amp vs Duration R^2','Amp = B*Duration','White Statistic','Amp vs Dur Spearman''s','White P-Value'};
idlist = [1 2 3 4]; itlist = [4 2 3 5];
for i = [1 3 4]
    if i>1
        X = reshape(ampdurRegress,6,3,20);
    elseif i<2
        X = reshape(ampdurSpearman,2,3,20);
    end
    X = permute(X,[1 3 2]);
    if i == 1; ofs = 0.1; else; ofs = 0.05; end
    subplot(2,2,i)
    a3WayBarPlot(R,X(idlist(i),:,:),ofs)
    ylim(ylist(:,itlist(i)))
    ylabel(ylablist{itlist(i)})
    title(titlist{itlist(i)})
end

subplot(2,2,2)
xlim([0 150]); ylim([0 5000])
xlabel('Burst Amplitude (a.u.)');
ylabel('Burst Duration (ms)')
legend(s2,R.condname,'Location','NorthWest')
set(gcf,'Position',[681   184 910 789])

a= 1;

%% Helper Functions:

function h = plotRateDists(edges,y,dur,cmap,histscat)
if histscat == 1
    h = rateHistogram(y,edges,dur); hold on
    h.FaceColor = cmap;
    h.FaceAlpha = 0.85;
    %     if cond == 3; h.FaceAlpha = 0.85; end
elseif histscat == 2
    h = rateScatter(y,edges,dur,cmap); hold on
end

function A = recentreDist(A,cond)
A = reshape(A,size(A,1),[]);

for i = 1:size(A,2)
    Acol = A(:,i);
    % %     if cond==3; Acol =Acol(randperm(size(Acol,1))); end
    [dum m] = max(Acol); % Max
    K = (median(1:size(A,1)))-m;
    A(:,i) = circshift(Acol,K,1);
end

function p = vecTTest(A)
for i = 1:size(A,1)
    try
        p(1,i) = ranksum(A(i,:,1),A(i,:,2));
        statv = statvec(squeeze(A(i,:,1)),squeeze(A(i,:,2)),1);
        [dum p(2,i)] = ttest(A(i,:,1),A(i,:,2));
    catch
        p(1,i) = NaN;
    end
    %     [dum p(2,i)] = ttest(A(i,:,1),A(i,:,2));
end

function [R2,s,stdeqv] = vmpars(X,Y)

for i = 1:size(Y,2)
%             figure(4)
    [dum dum R2(i) s(:,i) stdeqv(i)] = vmfit(X,Y(:,i)',50,0);
%     hold on
%             close all
end

function hl = VmPlot(xin,A,s,cmap,marktype)
    vm=  @(b,x) ((exp(b(1).*cos(x-b(2))))./(2*pi*besselj(0,b(1))) ) + b(3);
    x = linspace(-pi,pi,50);
    y = vm(mean(s,2),x); yhat = vm(std(s,[],2)./sqrt(size(s,2)),x);
    [hl hp] = boundedline(x,y,yhat);
    hl.Color = cmap;
    hl.LineWidth = 2;
    hp.FaceColor = cmap;
    hp.FaceAlpha = 0.2;
    hp.EdgeColor = cmap;
    hp.LineStyle = '--';
    hold on
    sc = scatter(repmat(binEdge2Mid(xin),1,size(A,2)),A(:),15,cmap,'filled','Marker',marktype);
    sc.MarkerFaceAlpha =1;
    
    
function a3WayBarPlot(R,X,ofs)
% assumes X has dimensions 1 x N x cond
p = [];
p(1) = ranksum(squeeze(X(1,:,1)),squeeze(X(1,:,2)));
p(2) = ranksum(squeeze(X(1,:,1)),squeeze(X(1,:,3)));
statv = statvec(squeeze(X(1,:,1)),squeeze(X(1,:,2)),1);

P = nan(3);
P(1,2) = p(1);
P(2,1) = p(1);
P(1,3) = p(2);
P(3,1) = p(2);

bptmp(:,:,1) = nanmean(squeeze(X(1,:,:)));
bptmp(:,:,2) = nanstd(squeeze(X(1,:,:)))/sqrt(sum(~isnan(squeeze(X))));

[HB,LEG] = barplot160818(R,1:4,bptmp,P,0,2,ofs);
a = gca;
a.XTickLabel = R.condname; a.XTickLabelRotation = 0;
delete(LEG);


%% SCRIPT GRAVE
%% 3 WAY ANOVA
% % X = plvAMPprcGroupData;
% % XG = [];GFlat = [];
% % for j = 1:length(size(X))
% %     XG{j} = nan(size(X));
% %     for i = 1:size(X,j)
% %         if j == 1
% %             XG{j}(i,:,:) = repmat(i,size(XG{j}(i,:,:)));
% %         elseif j == 2
% %             XG{j}(:,i,:) = repmat(i,size(XG{j}(:,i,:)));
% %         elseif j ==3
% %             XG{j}(:,:,i) = repmat(i,size(XG{j}(:,:,i)));
% %         end
% %     end
% %     GFlat(:,j) = reshape(XG{j},1,[])
% % end
% %
% % XFlat = flatten(X);
% %
% % figure
% %
% % [~,~,stats]  = anovan(XFlat,{GFlat(:,1),GFlat(:,2),GFlat(:,3)},'varnames',{'bin','subject','condition'});
% % figure
% % results = multcompare(stats,'Dimension',[1 3]);
% % a= 1;
% % % %
% % varnames = {'amp','cond'};
% % t = table(GFlat(:,3),XFlat',GFlat(:,1),GFlat(:,2),'VariableNames',{'cond','meas1','meas2','meas3'});
% % rm = fitrm(t,'meas1-meas3~1+cond','WithinDesign',table([1 2 3]','VariableNames',{'Measurements'}))
% % [~,~,stats]  = ranova(rm,'WithinModel','TestCond*Attention*TMS');
% % results = multcompare(stats,'Dimension',[1 3]);

%% SCATTER DENSITY
%     figure(42)
%     levels = 100;
%     [ map xax yax ]  = dataDensity(AsegA{cond},AsegLGroup{cond}, 200, 200,[0 100 0 6000]);
%     map = map - min(min(map));
%     map = floor(map ./ max(max(map)) * (levels-1));
%     
%     newpoints = 75;
%     [xaxq,yaxq] = meshgrid(...
%         linspace(min(min(xax,[],2)),max(max(xax,[],2)),newpoints ),...
%         linspace(min(min(yax,[],1)),max(max(yax,[],1)),newpoints )...
%         );
%     dmap = interp2(xax,yax,map,xaxq,yaxq,'cubic');
%     [h] = contourfcmap(xaxq,yaxq,dmap,100-[95 68 34]',linspecer(2),'method','calccontour');
%     h.p = h.p(end:-1:1)
%     cmap = repmat(R.condcmap(cond,:),3,1).*[1 0.8 0.4]';
%     k = 0;
%     for i = [5 4 2]
%         k = k+1;
%         h.p(i).FaceColor = cmap(k,:);
% %         h.p(i).FaceAlpha = 0.3;
%         h.LineWidth = 1;
%     end
%     
