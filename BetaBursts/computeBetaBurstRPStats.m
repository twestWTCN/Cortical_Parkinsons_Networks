function BB = computeBetaBurstRPStats(R,BB,F,plotop)
cnt = 0;
for band = 2:3
    cnt = cnt +1;
    for cond = 1:2
        
        if isequal([2 1],[band cond]); BB.Seg_binRP = []; BB.PLV_binRP = []; end
        for bs = 1:numel(BB.range.RP)-1
            s = find(BB.segRP{cond}(band,:)>=BB.range.RP(bs) & BB.segRP{cond}(band,:)<BB.range.RP(bs+1));
            
            w = (numel(s)/numel(BB.segA_save{cond}));
            
             x  = nanmean(BB.segA_save{cond}(s));
             xv = nanstd(BB.segA_save{cond}(s));
             BB.Amp_binRP(bs,band,cond,1) = w*x;
             BB.Amp_binRP(bs,band,cond,2) = w*xv;
             
             x = nanmean(BB.segPLV{cond}(:,band,s));
             xv = nanstd(BB.segPLV{cond}(:,band,s));
             BB.PLV_binRP(bs,band,cond,1)  = w*x;
             BB.PLV_binRP(bs,band,cond,2) = w*xv;
        end
        
        
        if plotop == 1
            figure(F(1))
            subplot(3,2,cnt)
            p(cond) = polarhistogram(BB.segRP{cond}(band,:),BB.range.RP,'Normalization','Probability'); hold on
            p(cond).FaceColor = R.condcmap(cond,:);
            p(cond).FaceAlpha = 0.8;
            title([R.bandinits{band} ' Relative Phase']); 
            pax = gca; pax.ThetaLim = [-180 180];
            pax.ThetaZeroLocation = 'Top';
            if isequal([2 3],[cond band]);
                l = legend(p,R.condname);
                l.Position = [0.4702    0.8007    0.0826    0.0375];
            end
        end
        
    end
end


if plotop == 1
    figure(F(1))
    cnt = 0;
    for band = 2:3
        cnt = cnt +1;
        subplot(3,2,cnt+2)
        
        barplot160818(R,BB.range.RP,squeeze(BB.Amp_binRP(:,band,:,:)),1)
        title(['Burst Amplitude by ' R.bandinits{band} ' Relative Phase'])
        ylabel(['Wghtd. ' R.bandinits{2} ' Amplitude']);
        xlabel([R.bandinits{band} ' Relative Phase']); ylim([0 6])
        
        subplot(3,2,cnt+4)
        barplot160818(R,BB.range.RP,squeeze(BB.PLV_binRP(:,band,:,:)),1)
        title(['PLV by ' R.bandinits{band} ' Relative Phase'])
        ylabel(['Wghtd. ' R.bandinits{band} ' SMA/STN PLV']);
        xlabel([R.bandinits{band} ' Relative Phase']); ylim([0 0.2])
        
    end
end