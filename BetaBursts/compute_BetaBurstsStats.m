function [segL_t_save segA_save AmpBin binSgEd] = compute_BetaBurstsStats(R,AmpTime,eps,plotop)
fsamp = R.fsamp;
for cond = 1:2
    X = AmpTime{cond};
    Xcd = X>eps;
    Xcd = double(Xcd);
    
    period = (2/R.bandef(2,1))*fsamp;
    consecSegs = SplitVec(find(Xcd(2,:)),'consecutive');
    
    % Work first with lengths
    segL = cellfun('length',consecSegs);
    segInds = find(segL>(period)); % segs exceeding min length
    segL_t = (segL/fsamp)*1000; % Segment lengths in ms
    segL_t(setdiff(1:length(segL),segInds)) = [];
    
        % Now do Amplitudes
        segA = [];
    for ci = 1:numel(segInds); segA(ci) = nanmean(X(consecSegs{segInds(ci)})); end

       % Bin data by length and find amplitude
    binSgEd = [0:50:1000 inf];
    for bs = 1:numel(binSgEd)-1
        s = find(segL_t>=binSgEd(bs) & segL_t<binSgEd(bs+1));
        AmpBin(bs,cond,1) = nanmean(segA(s));
        AmpBin(bs,cond,2) = nanstd(segA(s));
    end 
    
        % Save across conditions
    segL_t_save{cond}= segL_t;
    segA_save{cond} = segA;

    
    if plotop == 1
        figure(1)
        colormap(R.condcmap)
        subplot(2,3,1)
        h = histogram((segL_t),10:50:1000); hold on
        h.FaceColor = R.condcmap(cond,:);
        h.FaceAlpha = 0.75;
        
        
        subplot(2,3,2)
        h = histogram(segA,0:.5:12); hold on
        h.FaceColor = R.condcmap(cond,:);
        h.FaceAlpha = 0.75;
        
        % Now do duration vs amplitude
        subplot(2,3,3)
        s = scatter(segA,segL_t,25,R.condcmap(cond,:),'filled'); hold on
        s.MarkerFaceAlpha = 0.75;
    end % Plotop
end
if plotop == 1
    subplot(2,3,1); legend(R.condname); ylabel('Frequency'); xlabel('Burst Length (ms)'); xlim([0 1000]); ylim([0 30])
    subplot(2,3,2); legend(R.condname); ylabel('Frequency'); xlabel('Amplitude');  xlim([0 15]); ylim([0 30]);
    subplot(2,3,3); legend(R.condname); ylabel('Burst Length (ms)'); xlabel('Amplitude'); xlim([0 15]); ylim([0 1000])
    
    subplot(2,1,2);
    X = binSgEd-(binSgEd(2)/2);
    X = [X(2:end-1) X(end-1)+(binSgEd(2)/2)];
    for m = 1:2;for n = 1:length(X);  cmn(n,m,:) = R.condcmap(m,:); end;  end
    HB = superbar(1:length(X),AmpBin(:,:,1),'E',AmpBin(:,:,2),'BarFaceColor',cmn) %m-by-n-by-3
%     errorbar(1:length(X),AmpBin(:,1,1),AmpBin(:,1,2),'.')
    a  = gca;
    a.XTick = 1:1:length(X);
    a.XTickLabel = X(1:1:end);
    a.XTickLabelRotation = 45;
    legend(HB(1,:),R.condname); ylabel('Amplitude'); xlabel('Burst Length (ms)');ylim([0 15])
    set(gcf,'Position',[680 347 913 632])
end