function BB = computeBetaBurstRPStats(R,BB,F,plotop)
cnt = 0;
for band = 3
    cnt = cnt +1;
    for cond = 1:length(R.condname)
        
        if isequal([3 1],[band cond]); BB.Seg_binRP = []; BB.PLV_binRP = []; end
        for bs = 1:numel(BB.range.RP)-1
            if size(BB.segRP{cond}(:),1)>1
            s = find(BB.segRP{cond}(band,:)>=BB.range.RP(bs) & BB.segRP{cond}(band,:)<BB.range.RP(bs+1));
            else
                s = 1; band = 1;
            end
            % Unnormalized
            w = (numel(s)/numel(BB.segA_save{cond}));
            if numel(s)<2; w = 0; end
            x  = nanmean(BB.segA_save{cond}(s));
            xv = nanstd(BB.segA_save{cond}(s))/sqrt(numel(s));
            
            BB.Amp_binRP(bs,band,cond,1) = w*x;
            BB.Amp_binRP(bs,band,cond,2) = w*xv;
            BB.Amp_binRP_data{cond}{bs} = BB.segA_save{cond}(s);
            
            % Percentage Deviation
            w = numel(s)/numel(BB.segAPrc_save{cond});
            if numel(s)<2; w = 0; end
            x  = nanmean(BB.segAPrc_save{cond}(s));
            xv = nanstd(BB.segAPrc_save{cond}(s))/sqrt(numel(s));
            
            BB.AmpPrc_binRP(bs,band,cond,1) = w*x;
            BB.AmpPrc_binRP(bs,band,cond,2) = w*xv;
            BB.AmpPrc_binRP_data{cond}{bs} = BB.segAPrc_save{cond}(s);
            
            x = nanmean(BB.segPLV{cond}(:,band,s));
            xv = nanstd(BB.segPLV{cond}(:,band,s))/sqrt(numel(s));
            BB.PLV_binRP(bs,band,cond,1)  = w*x;
            BB.PLV_binRP(bs,band,cond,2) = w*xv;
            BB.PLV_binRP_data{cond}{bs} = BB.segPLV{cond}(s);
        end
        
        
        if plotop == 1
            figure(F(1))
            subplot(2,3,3)
            p(cond) = polarhistogram(BB.segRP{cond}(band,:),BB.range.RP); hold on
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
    for band = 3
        cnt = cnt +1;
        % test
        pvec = pvec_bin_TTest(BB.AmpPrc_binRP_data);
        %Plot
        subplot(2,3,2+3)
        barplot160818(R,BB.range.RP,squeeze(BB.AmpPrc_binRP(:,band,:,:)),pvec,1,1)
        title(['% Amplitude by ' R.bandinits{band} ' Relative Phase'])
        ylabel(['Wghtd. ' R.bandinits{2} ' % Amplitude']);
        xlabel([R.bandinits{band} ' Relative Phase']); ylim(BB.plot.lims.wAmpPrc)
        
        % Test
        pvec = pvec_bin_TTest(BB.PLV_binRP_data);
        %Plot
        subplot(2,3,3+3)
        barplot160818(R,BB.range.RP,squeeze(BB.PLV_binRP(:,band,:,:)),pvec,1,1)
        title(['PLV by ' R.bandinits{band} ' Relative Phase'])
        ylabel(['Wghtd. ' R.bandinits{band} ' SMA/STN PLV']);
        xlabel([R.bandinits{band} ' Relative Phase']); ylim(BB.plot.lims.PLV)

        
    end
end