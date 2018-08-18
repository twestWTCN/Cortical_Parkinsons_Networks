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
for sub = 2:3; %length(R.subname)
    for breg = 2:length(R.bregname)
        for side = 1:2
            load([R.datapathr R.subname{sub} '\ftdata\BetaBursts\BetaBurstAnalysis_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'BB')
            for cond = 1:2
                segAGroup{cond,side,sub} = BB.segA_save{cond};
                segLGroup{cond,side,sub} = BB.segL_t_save{cond};
                segPLVGroup{cond,side,sub} = squeeze(BB.segPLV{cond}(:,3,:))';
                recLGroup(cond,side,sub) = diff(BB.DTvec{cond}([1 end]))/60;
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
a = 1;

%     BB.stats.DurCond  = statvec(BB.segL_t_save{1},BB.segL_t_save{2},1);
