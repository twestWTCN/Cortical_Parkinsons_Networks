function BB = computeBetaBurstPLVStats(R,BB,F,plotop)
cnt = 0;
for cond = 1:2
    for band = 2:size(R.bandef,1)
        % Segment PLV during Bursts
        if isequal([band cond],[2 1]); BB.segPLV = []; BB.segRP = []; end
        for ci = 1:size(BB.segA_save{cond},2)
            tind = BB.segTInds{cond}{ci};
            swind = find(BB.SWTvec{cond} >= tind(1) & BB.SWTvec{cond} < tind(2));
            BB.segPLV{cond}(:,band,ci) = mean(BB.PLV{cond}(band,swind));
            BB.segRP{cond}(band,ci) = circ_mean(BB.RP{cond}(band,swind),[],2);
        end
        
        % Bin data by PLV and find Amplitude
        BB.binPLVEd = BB.range.PLV;
        if isequal([band cond],[2 1]); BB.Amp_binPLV = []; end
        for bs = 1:numel(BB.binPLVEd)-1
            s = find(BB.segPLV{cond}(:,band,:)>=BB.binPLVEd(bs) & BB.segPLV{cond}(:,band,:)<BB.binPLVEd(bs+1));
            x = nanmean(BB.segA_save{cond}(s));
            xv = nanstd(BB.segA_save{cond}(s))/sqrt(numel(s));
            w = (numel(s)/numel(BB.segA_save{cond}));
            BB.Amp_binPLV(bs,band,cond,1) = x*w;
            BB.Amp_binPLV(bs,band,cond,2) = xv*w;
        end
        
        if plotop == 1
            figure(F(1))
            subplot(3,2,band-1)
            colormap(R.condcmap)
            h = rateHistogram(squeeze(BB.segPLV{cond}(:,band,:)),BB.range.PLV,diff(BB.DTvec{cond}([1 end]))/60); hold on
            h.FaceColor = R.condcmap(cond,:);
            h.FaceAlpha = 0.75;
            legend(R.condname); ylabel('Frequency (min^{-1})'); xlabel('PLV');  box off; ylim([0 14]); xlim([0 1])
            title([R.bandinits{band} ' SMA/STN Phase Sync'])
            
            subplot(3,2,(2)+(band-1))
            s = scatter(BB.segPLV{cond}(:,band,:),BB.segA_save{cond},0.075*BB.segL_t_save{cond},R.condcmap(cond,:),'filled'); hold on
            s.MarkerFaceAlpha = 0.7;
            ylabel('Amplitude'); xlabel('PLV'); ylim([0 BB.range.Amp(end)]); xlim([0 1]);
            title([R.bandinits{band} ' PLV vs Burst ' R.bandinits{2} ' Amp.'])
            
            szvec = [200 400 600];
            for sz = 1:3
                sp(sz) = plot(NaN,NaN,'o', 'MarkerSize', sqrt(0.075*szvec(sz)));
                sp(sz).MarkerFaceColor = R.condcmap(cond,:);
                sp(sz).MarkerEdgeColor	 = 'none';
                szleg{sz} = sprintf('%.0f ms',szvec(sz));
            end
            legend(sp,szleg)
        end
        % %         % Scatter PLV vs Duration (No Findings!)
        % %         figure(F(2))
        % %         cnt = cnt + 1;
        % %         subplot(2,2,cnt)
        % %         s = scatter(BB.segL_t_save{cond}, BB.segPLV{cond}(:,band,:),25,R.condcmap(cond,:),'filled'); hold on
        % %         s.MarkerFaceAlpha = 0.75;
        % %         xlabel('Burst Duration'); ylabel('Wavelet PLV');% xlim([0 20]); ylim([0 1]);
        % %         title([R.bandname{band} ' Amp vs PLV'])
    end
end

cnt = 4;
if plotop == 1
    for band = 2:3
        cnt = cnt + 1;
        figure(F(1))
        subplot(3,2,cnt)
        % Bar Chart
        barplot160818(R,BB.binPLVEd,squeeze(BB.Amp_binPLV(:,band,:,:)))
        title('Burst Amp by PLV')
        xlabel('PLV'); ylabel('Wghtd. Amplitude'); ylim([0 5]);
    end
end

%% script grave
% %  COMPUTE NON-SEGMENTED SCATTERGRAMS
% % figure; p = 0;
% % for cond = 1:2
% %     for band = 2:3
% %         p = p+1;
% %         tarC = 1:10:length(BB.SWTvec{cond});
% %         A = BB.SWTvec{cond}(tarC);           % Slow Time
% %         B = BB.DTvec{cond}(1:end); % Fast Time
% %
% %         % Resample the Fast time relative to Slow
% %         C = [];
% %         for ac = 1:length(A); [c index] = min(abs(A(ac)-B)); C(ac) = index; end
% %         AmpResamp = BB.A{cond}(2,C); % Low beta amp
% %         PLVResamp = BB.PLV{cond}(band,tarC); % variable b1/2 PLV
% %         subplot(2,2,p)
% %         s = scatter(AmpResamp,PLVResamp,25,R.condcmap(cond,:),'filled');
% %         s.MarkerFaceAlpha = 0.75;
% %     end
% % end