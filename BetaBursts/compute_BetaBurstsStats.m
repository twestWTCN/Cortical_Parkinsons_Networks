function [BB] = compute_BetaBurstsStats(R,BB,F,plotop)
% This functio
BB.guide = {...
    'segL_t_save = segment lengths in ms'
    'segA_save = seg amplitudes'
    'AmpBin = frequencies (occurence) of binned amps'
    'binSgEd = set of bin edges'
    'segTInds = start end of time'
    };

fsamp = R.fsamp;
for cond = 1:length(R.condname)
    X = BB.A{cond}; % Copy amplitude data
    
    % DFA-AE
    A = remnan(X(2,:)')';
    DFAP = [];
    DFAP(1) = BB.fsamp; DFAP(2) = 8/14;  DFAP(3) = 8;  DFAP(4) = 50;  DFAP(5) = 0;
    [bmod win evi alpha] = peb_dfa_gen(A,DFAP,0);
    BB.AEDFA(:,cond) = [alpha evi(2)];

    Xcd = X>BB.epsAmp; % Threshold on eps
    Xcd = double(Xcd); % Convert from logical
    
    % Work first with lengths
    BB.period = (2/R.bandef(2,1))*fsamp;
    consecSegs = SplitVec(find(Xcd(2,:)),'consecutive');
    segL = cellfun('length',consecSegs);
    segInds = find(segL>(BB.period)); % segs exceeding min length
    Xnorm = X;
    % Define the Amplitudes
    % Normalised
    % % %     %  XPrc = BB.APrc{cond};
    % % %     %         Xnorm = (X-max(X(:,[consecSegs{:}]),[],2))./max(X(:,[consecSegs{:}]),[],2);
    % % %     %     Xnorm = (X-min(X(:,[consecSegs{:}]),[],2))./min(X(:,[consecSegs{:}]),[],2);
    % % %     %         Xnorm = (X-nanmedian(X(:,[consecSegs{:}]),2))./nanmedian(X(:,[consecSegs{:}]),2);
    % % %     XL = log10(X); %##
    %         XL = X;
    %     XL = X(:,[consecSegs{segInds}]);
    %     Xnorm = (X-nanmedian(XL,2))./nanmedian(XL,2);
    %     Xnorm = Xnorm.^2;
    
    % or prctile
    %     Xnorm = [];
    %     for i = 1:2
    %         %         XL = remnan(X(1,:));
    %         [~,~,rnk] = unique(X(i,:));
    %         Xnorm(i,:) = rnk./size(X,2);
    %     end
    %     Xnorm = (Xnorm-0.5).*100;
    
    % Work first with lengths
    BB.segL_t_save{cond} = (segL/fsamp)*1000; % Segment lengths in ms
    BB.segL_t_save{cond}(setdiff(1:length(segL),segInds)) = [];
    
    % Segment Time indexes
   BB.segInds{cond}{1} = NaN(1,2); BB.segTInds{cond}{1} = NaN(1,2);
    for ci = 1:numel(segInds);
        BB.segInds{cond}{ci} = consecSegs{segInds(ci)};
        BB.segTInds{cond}{ci} = (consecSegs{segInds(ci)}([1 end])/fsamp);
    end
    
    % Now do Amplitudes
   BB.segA_save{cond} = NaN;
    for ci = 1:numel(segInds)
        BB.segA_save{cond}(ci) = nanmean(X(2,consecSegs{segInds(ci)}));
%         BB.segAPrc_save{cond}(ci) = nanmean(Xnorm(2,consecSegs{segInds(ci)})); % Uses the Normed Amplitudes
    end
    %%%
    X = BB.segA_save{cond};
    Xnorm = 100*(X-median(X))/median(X);
%     Xnorm = Xnorm-median(Xnorm); % 2nd normalization brings conditions closer to zero mean
    BB.segAPrc_save{cond} = Xnorm; %nanmean(Xnorm(2,consecSegs{segInds(ci)})); % Uses the Normed Amplitudes
    %%%
    
    % Bin data by Amp and find Length
    BB.binAmp = [BB.range.Amp inf]; 
    if cond == 1; BB.Seg_binAmp = []; end
    for bs = 1:numel(BB.binAmp)-1
        s = find(BB.segA_save{cond}>=BB.binAmp(bs) & BB.segA_save{cond}<BB.binAmp(bs+1));
        x = nanmean(BB.segL_t_save{cond}(s));
        xv = nanstd(BB.segL_t_save{cond}(s))/sqrt(numel(s));
        w = (numel(s)/numel(BB.segA_save{cond}));
%         if numel(s)<2; w = 0; end
        
        BB.Seg_binAmp(bs,cond,1) = x*w;
        BB.Seg_binAmp(bs,cond,2) = xv*w;
        BB.Seg_binAmpData{bs,cond} = BB.segL_t_save{cond}(s);
        Seg_binAmpDataLOC{cond}{bs} = BB.segL_t_save{cond}(s); % LOCAL (different formatting)
        
    end
    
    % Bin data by length and find amplitude
    BB.binDur = [BB.range.segDur inf];
    if cond == 1; BB.Amp_binDur = []; end
    for bs = 1:numel(BB.binDur)-1
        s = find(BB.segL_t_save{cond}>=BB.binDur(bs) & BB.segL_t_save{cond}<BB.binDur(bs+1));
        x = nanmean(BB.segA_save{cond}(s));
        xv = nanstd(BB.segA_save{cond}(s))/sqrt(numel(s));
        w =  (numel(s)/numel(BB.segA_save{cond}));
%                     if numel(s)<2; w = 0; end
        BB.Amp_binDur(bs,cond,1) = x*w;
        BB.Amp_binDur(bs,cond,2) = xv*w;
        BB.Amp_binDurData{bs,cond} = BB.segA_save{cond}(s);
        Amp_binDurDataLOC{cond}{bs} = BB.segA_save{cond}(s);% LOCAL (different formatting)
    end
    
    % Do amplitude/length correlation stats
    if numel(BB.segA_save{cond})>2
        [Rho pval] = corr(BB.segA_save{cond}',BB.segL_t_save{cond}','type','Spearman');
        [xCalc yCalc b Rsq bHd RsqHd wPval] = linregress(BB.segA_save{cond}',BB.segL_t_save{cond}',1);
        BB.corrs.ampdur.spearman{cond} = [Rho pval];
        BB.corrs.ampdur.rest_regress{cond} = [Rsq b RsqHd wPval bHd'];
    else
        BB.corrs.ampdur.spearman{cond} = [NaN NaN];
        BB.corrs.ampdur.rest_regress{cond} = [NaN NaN];
    end
    % Save across conditions
    if plotop == 1
        figure(F(1))
        colormap(R.condcmap)
        subplot(3,3,1)
        h = rateHistogram(BB.segL_t_save{cond},BB.range.segDur,diff(BB.DTvec{cond}([1 end]))/60); hold on
        h.FaceColor = R.condcmap(cond,:);
        h.FaceAlpha = 0.75;
        %         a = gca;
        %         a.XTick = round( a.XTick(1:2:end) / 10) * 10;
        
        if cond == 2
            statv =statvec((BB.segL_t_save{1}),(BB.segL_t_save{2}),1);
            txtlab = [sprintf('OFF-ON: %.1f [%0.1f %0.1f] ',statv(end-3:end-1)) ' (' textstar(statv(end),0.05) ')'];
            TS = text(35,13, txtlab); TS.FontWeight = 'bold'; TS.FontSize = 7; TS.FontWeight = 'bold';
        end
        
        subplot(3,3,2)
        h = rateHistogram(BB.segA_save{cond},BB.binAmp,diff(BB.DTvec{cond}([1 end]))/60); hold on
        h.FaceColor = R.condcmap(cond,:); ylim(BB.plot.lims.wAmpPrc)
        h.FaceAlpha = 0.75;
        
        if cond == 2
            statv =statvec(BB.segA_save{1},BB.segA_save{2},1);
            txtlab = [sprintf('OFF-ON: %.1f [%0.1f %0.1f] ',statv(end-3:end-1)) ' (' textstar(statv(end),0.05) ')'];
            TS = text(5,9.5, txtlab); TS.FontWeight = 'bold'; TS.FontSize = 7; TS.FontWeight = 'bold';
        end
        
        % Now do duration vs amplitude
        subplot(3,3,3)
        sp(cond) = scatter(BB.segA_save{cond},BB.segL_t_save{cond},25,R.condcmap(cond,:),'filled'); hold on
        sp(cond).MarkerFaceAlpha = 0.7;
        plot(xCalc,yCalc,'color',R.condcmap(cond,:),'LineWidth',2); ylim(BB.plot.lims.Dur); xlim(BB.plot.lims.Amp)
        
    end % Plotop
end
if plotop == 1
    subplot(3,3,1); legend(R.condname,'Location','East'); ylabel('Frequency (min^{-1})'); xlabel('Duration (ms)'); xlim([0 BB.range.segDur(end)]); %ylim([0 15]); box off
    title('Burst Duration '); ylim([0 20])
    BB.stats.DurCond  = statvec(BB.segL_t_save{1},BB.segL_t_save{2},1);
    
    subplot(3,3,2); legend(R.condname,'Location','East'); ylabel('Frequency (min^{-1})'); xlabel('Amplitude');  xlim([BB.range.Amp(1) BB.range.Amp(end)]); %ylim([0 15]); box off
    title('Burst Amplitude '); ylim([0 25])
    BB.stats.AmpCond = statvec(BB.segA_save{1},BB.segA_save{2},1);
    
    subplot(3,3,3); legend(sp,R.condname,'Location','NorthEast'); ylabel('Duration (ms)'); xlabel('Amplitude'); xlim([BB.range.Amp(1) BB.range.Amp(end)]); %ylim([0 BB.range.segDur(end)]); box off
    title('Burst Amplitude vs Duration ');
    BB.stats.AmpCond = statvec(BB.segA_save,BB.segL_t_save,2);
    
    subplot(3,1,2);
    pvec = pvec_bin_TTest(Seg_binAmpDataLOC);
    barplot160818(R,BB.binAmp,squeeze(BB.Seg_binAmp(:,:,:)),pvec,0,1)
    title('Burst Duration by Amplitude')
    xlabel('Mean Amplitude'); ylabel('Wghtd. Duration (ms)'); ylim([0 140]); % ylim([BB.range.Amp(1) BB.range.Amp(end)]);
    
    subplot(3,1,3);
    %     BB.Seg_binAmpStats
    pvec = pvec_bin_TTest(Amp_binDurDataLOC);
    barplot160818(R,BB.binDur,squeeze(BB.Amp_binDur(:,:,:)),pvec,0,1)
    title('Burst Amplitude by Duration')
    xlabel('Duration (ms)'); ylabel('Wghtd. Amplitude'); ylim([0 7]); %ylim([BB.range.Amp(1) BB.range.Amp(end)]);
end

% function binstats(XY)