function compute_BetaBursts_GroupLevel(R)
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
for sub = 1:length(R.subname)
    for breg = 1:length(R.bregname)
        for side = 1:2
            try
                load([R.datapathr R.subname{sub} '\ftdata\BetaBursts\BetaBurstAnalysis_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg} '_org'],'BB')
                for cond = 1:2
                    segAGroup{cond,side,sub} = BB.segA_save{cond};
                    segLGroup{cond,side,sub} = BB.segL_t_save{cond};
                    segPLVGroup{cond,side,sub} = squeeze(BB.segPLV{cond}(:,3,:))';
                    recLGroup(cond,side,sub) = diff(BB.DTvec{cond}([1 end]))/60;
                    rpAGroup(:,cond,side,sub) = squeeze(BB.Amp_binRP(:,3,cond,1));
                    rpPLVGroup(:,cond,side,sub) = squeeze(BB.PLV_binRP(:,3,cond,1));

                end
            catch
                disp('load error!')
                segAGroup{cond,side,sub} = [];
                segLGroup{cond,side,sub} = [];
                segPLVGroup{cond,side,sub} = [];
                recLGroup(cond,side,sub) = 0;
                rpAGroup(:,cond,side,sub) = NaN(size(squeeze(BB.Amp_binRP(:,3,cond,1))));
                rpPLVGroup(:,cond,side,sub) = NaN(size(squeeze(BB.PLV_binRP(:,3,cond,1))));
            end
            
            
            disp([cond side breg sub])
        end
    end
end
segAGroup = reshape(segAGroup,2,[]);
segLGroup = reshape(segLGroup,2,[]);
segPLVGroup = reshape(segPLVGroup,2,[]);

for cond = 1:2
    totRecL = sum(reshape(recLGroup(cond,:,:),1,[]));
    figure(1)
    % Group Amplitude Histogram
    subplot(3,1,1)
    A = segAGroup(cond,:);
    AsegA{cond} = [A{:}];
    h = rateHistogram(AsegA{cond},BB.range.Amp,totRecL); hold on
    h.FaceColor = R.condcmap(cond,:);
    h.FaceAlpha = 0.75;
    
    % Group Segment Length Histogram
    subplot(3,1,2)
    A = segLGroup(cond,:);
    AsegLGroup{cond} = [A{:}];
    h = rateHistogram(AsegLGroup{cond},BB.range.segDur,totRecL); hold on
    h.FaceColor = R.condcmap(cond,:);
    h.FaceAlpha = 0.75;
    
    % Group PLV Histogram
    subplot(3,1,3)
    A = segPLVGroup(cond,:);
    AsegPLVGroup{cond} = [A{:}];
    h = rateHistogram(AsegPLVGroup{cond},BB.range.PLV,totRecL); hold on
    h.FaceColor = R.condcmap(cond,:);
    h.FaceAlpha = 0.75;
    
    % Group RP to Amp
    A = squeeze(rpAGroup(:,cond,:,:));
    A = reshape(A,size(A,1),[]);
    for i = 1:size(A,2)
        Acol = A(:,i);
        [dum m] = max(Acol); % Max
        % or centre of mass
        %         [dum m] = min(abs(Acol-prctile(Acol,80)));
        %         m = round(sum(Acol.*1:12)./sum(Acol))
        K = (floor(size(Acol,1))/2)-m;
        A(:,i) = circshift(Acol,K,1);
    end
    rpAGroupData(:,:,cond) = A;
    rpAGroupStat(:,cond,1) = nanmean(A,2); rpAGroupStat(:,cond,2) = nanstd(A,1,2)./sqrt(size(A,2));

    % Group RP to PLV
    A = squeeze(rpPLVGroup(:,cond,:,:));
    A = reshape(A,size(A,1),[]);
    for i = 1:size(A,2)
        Acol = A(:,i);
        Acol(Acol==0) = NaN;
        [dum m] = max(Acol); % Max
        K = (floor(size(Acol,1))/2)-m;
        A(:,i) = circshift(Acol,K,1);
    end

    rpPLVGroupData(:,:,cond) = A;
    rpPLVGroupStat(:,cond,1) = nanmean(A,2); rpPLVGroupStat(:,cond,2) = nanstd(A,1,2)./sqrt(size(A,2));
end




subplot(3,1,1); legend(R.condname); ylabel('Frequency (min^{-1})'); xlabel('Amplitude');  xlim([0 BB.range.Amp(end)]); ylim([0 10]); box off
title('Burst Amplitude ')
statvec(AsegA{1},AsegA{2},1)

subplot(3,1,2); legend(R.condname); ylabel('Frequency (min^{-1})'); xlabel('Duration (ms)'); xlim([0 BB.range.segDur(end)]); ylim([0 10]); box off
title('Burst Duration ')
statvec(AsegLGroup{1},AsegLGroup{2},1)

subplot(3,1,3); legend(R.condname); ylabel('Frequency (min^{-1})'); xlabel('PLV'); xlim([0 BB.range.PLV(end)]); ylim([0 10]); box off
title('Burst PLV ')
statvec(AsegPLVGroup{1},AsegPLVGroup{2},1)

figure(2)
subplot(2,1,1)
for i = 1:size(rpAGroupData,1)
    p(1,i) = ranksum(rpAGroupData(i,:,1),rpAGroupData(i,:,2));
    [dum p(2,i)] = ttest(rpAGroupData(i,:,1),rpAGroupData(i,:,2));
end
barplot160818(R,BB.range.RP,rpAGroupStat,p(2,:),1,2)
title('Burst Amplitude by Relative Phase')
xlabel('Relative Phase')

subplot(2,1,2)
for i = 1:size(rpPLVGroupData,1)
    p(1,i) = ranksum(rpPLVGroupData(i,:,1),rpPLVGroupData(i,:,2));
    [dum p(2,i)] = ttest(rpPLVGroupData(i,:,1),rpPLVGroupData(i,:,2));
end
barplot160818(R,BB.range.RP,rpPLVGroupStat,p(2,:),1,2)
title('Burst PLV by Relative Phase')
xlabel('Relative Phase')

%     BB.stats.DurCond  = statvec(BB.segL_t_save{1},BB.segL_t_save{2},1);
