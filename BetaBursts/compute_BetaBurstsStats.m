function [BB] = compute_BetaBurstsStats(R,BB,F,plotop)
BB.guide = {...
    'segL_t_save = segment lengths in ms'
    'segA_save = seg amplitudes'
    'AmpBin = frequencies (occurence) of binned amps'
    'binSgEd = set of bin edges'
    'segTInds = start end of time'
    };

fsamp = R.fsamp;
for cond = 1:2
    X = BB.A{cond};
    Xcd = X>BB.epsAmp;
    Xcd = double(Xcd);
    
    period = (2/R.bandef(2,1))*fsamp;
    consecSegs = SplitVec(find(Xcd(2,:)),'consecutive');
    
    % Work first with lengths
    segL = cellfun('length',consecSegs);
    segInds = find(segL>(period)); % segs exceeding min length
    segL_t = (segL/fsamp)*1000; % Segment lengths in ms
    segL_t(setdiff(1:length(segL),segInds)) = [];
    
    % Segment Time indexes
    segT_ind = [];
    for ci = 1:numel(segInds); segT_ind{ci} = (consecSegs{segInds(ci)}([1 end])/fsamp); end
    
    % Now do Amplitudes
    segA = [];
    for ci = 1:numel(segInds); segA(ci) = nanmean(X(consecSegs{segInds(ci)})); end
    
    BB.segL_t_save{cond}= segL_t;
    BB.segA_save{cond} = segA;
    BB.segTInds{cond} = segT_ind;
    
    % Bin data by Amp and find Length
    BB.binAmp = [BB.range.Amp inf];
    if cond == 1; BB.Seg_binAmp = []; end
    for bs = 1:numel(BB.binAmp)-1
        s = find(BB.segA_save{cond}>=BB.binAmp(bs) & BB.segA_save{cond}<BB.binAmp(bs+1));
        x = nanmean(BB.segL_t_save{cond}(s));
        xv = nanstd(BB.segL_t_save{cond}(s))/sqrt(numel(s));
        w = (numel(s)/numel(segA));
        
        BB.Seg_binAmp(bs,cond,1) = x*w;
        BB.Seg_binAmp(bs,cond,2) = xv*w;
        BB.Seg_binAmpData{bs,cond} = BB.segL_t_save{cond}(s);
    end
    
    % Bin data by length and find amplitude
    BB.binDur = [BB.range.segDur inf];
    if cond == 1; BB.Amp_binDur = []; end
    for bs = 1:numel(BB.binDur)-1
        s = find(BB.segL_t_save{cond}>=BB.binDur(bs) & BB.segL_t_save{cond}<BB.binDur(bs+1));
        x = nanmean(BB.segA_save{cond}(s));
        xv = nanstd(BB.segA_save{cond}(s))/sqrt(numel(s));
        w = (numel(s)/numel(segA));
        
        BB.Amp_binDur(bs,cond,1) = x*w;
        BB.Amp_binDur(bs,cond,2) = xv*w;
        BB.Amp_binDurData{bs,cond} = BB.segA_save{cond}(s);
    end
    
    
    
    % Save across conditions
    if plotop == 1
        figure(F(1))
        colormap(R.condcmap)
        subplot(3,3,1)
        h = rateHistogram(segL_t,BB.range.segDur,diff(BB.DTvec{cond}([1 end]))/60); hold on
        h.FaceColor = R.condcmap(cond,:);
        h.FaceAlpha = 0.75;
        
        subplot(3,3,2)
        h = rateHistogram(segA,BB.binAmp,diff(BB.DTvec{cond}([1 end]))/60); hold on
        h.FaceColor = R.condcmap(cond,:);
        h.FaceAlpha = 0.75;
        
        % Now do duration vs amplitude
        subplot(3,3,3)
        sp(cond) = scatter(segA,segL_t,25,R.condcmap(cond,:),'filled'); hold on
        sp(cond).MarkerFaceAlpha = 0.7;
        [xCalc yCalc b Rsq] = linregress(segA',segL_t',1)
        plot(xCalc,yCalc,'color',R.condcmap(cond,:),'LineWidth',2)
        
    end % Plotop
end
if plotop == 1
    subplot(3,3,1); legend(R.condname); ylabel('Frequency (min^{-1})'); xlabel('Duration (ms)'); xlim([0 BB.range.segDur(end)]); ylim([0 10]); box off
    title('Burst Duration ')
    BB.stats.DurCond  = statvec(BB.segL_t_save{1},BB.segL_t_save{2},1);
    
    subplot(3,3,2); legend(R.condname); ylabel('Frequency (min^{-1})'); xlabel('Amplitude');  xlim([0 BB.range.Amp(end)]); ylim([0 10]); box off
    title('Burst Amplitude ')
    BB.stats.AmpCond = statvec(BB.segA_save{1},BB.segA_save{2},1);
    
    subplot(3,3,3); legend(sp,R.condname); ylabel('Duration (ms)'); xlabel('Amplitude'); xlim([0 BB.range.Amp(end)]); ylim([0 BB.range.segDur(end)]); box off
    title('Burst Amplitude vs Duration ')
    BB.stats.AmpCond = statvec(BB.segA_save,BB.segL_t_save,2);
    
    subplot(3,1,2);
    %     BB.Seg_binAmpStats
    barplot160818(R,BB.binAmp,squeeze(BB.Seg_binAmp(:,:,:)),[],0)
    title('Burst Duration by Amplitude')
    xlabel('Mean Amplitude'); ylabel('Wghtd. Duration (ms)'); ylim([0 160]);
    
    subplot(3,1,3);
    %     BB.Seg_binAmpStats
    barplot160818(R,BB.binDur,squeeze(BB.Amp_binDur(:,:,:)),[],0)
    title('Burst Amplitude by Duration')
    xlabel('Duration (ms)'); ylabel('Wghtd. Amplitude'); ylim([0 8]);
end

% function binstats(XY)