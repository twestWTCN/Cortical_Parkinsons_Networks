function [AmpTime,data_tvec,PLVTime,sw_tvec,RPTime]  = compute_BetaBursts(R,BB,data,shuff)

if shuff == 1
    for i = 1:size(data.trial{1},1)
        x = data.trial{1}(i,:);
        x = randperm(length(x));
        data.trial{1}(i,:) = x;
    end
end

for i = 1:size(data.trial{1},1)
    x = data.trial{1}(i,:);
    x(abs(x)>(4*std(x))) = NaN;
    x= fillmissing(x,'spline');
    data.trial{1}(i,:) = x;
end

data.trial{1} = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[12 38], [], 'fir', 'twopass', 'reduce');

for i = 1:size(data.trial{1},1)
    x = data.trial{1}(i,:);
    x = (x-mean(x))./std(x);
    data.trial{1}(i,:) = x;
end

cfg = [];
cfg.method     = 'wavelet';
cfg.width      = 10; % spectral bandwidth at frequency F: (F/width)*2
cfg.gwidth     = 5;
cfg.output     = 'pow';
cfg.foi        = 8:0.5:35;
cfg.toi        = data.time{1};
cfg.pad = 'nextpow2';
TFRwave = ft_freqanalysis(cfg, data);

% cfg = []
% cfg.channel = 'STN_L01';
% cfg.zlim = [0 60]
% ft_singleplotTFR(cfg, TFRwave);

% Indices of foi
% % fInd = find(TFRwave.freq>=R.bandef(2,1) & TFRwave.freq<=R.bandef(2,2)); % Take all of the indices
% OR account for smoothing and dont stray outside band
% % [dum fInd] = min(abs(TFRwave.freq-median(R.bandef(2,:)))); % Take middle!
% OR use peak of power
[dum fInd] = min(abs(TFRwave.freq-BB.powfrq)); % Take middle!

fbwid = (median(R.bandef(2,:))/cfg.width)*2;
disp([R.bandname{2} ' amp estimate at mid freq ' num2str(TFRwave.freq(fInd)) ' Hz +/- ' num2str(fbwid)])

AmpTime = squeeze(mean(TFRwave.powspctrm(:,fInd,:),2));
data_tvec = data.time{1};
%% Now Recover Phase
cfg = [];
cfg.method     = 'wavelet';
cfg.width      = 10;
cfg.gwidth     = 5;
cfg.output     = 'fourier';
cfg.foi        = 8:0.5:35;
cfg.toi        = data.time{1};
cfg.pad = 'nextpow2';
TFRwave = ft_freqanalysis(cfg, data);

for band = 2:3
    % Take all Inds within band
    %     %     fInd = find(TFRwave.freq>=R.bandef(band,1) & TFRwave.freq<=R.bandef(band,2));
    %     %     P = squeeze(mean(TFRwave.fourierspctrm(1,:,fInd,:),3));
    %     % OR Take middle:
    %     [dum fInd] = min(abs(TFRwave.freq-median(R.bandef(band,:)))); % Take middle!
    
    if band == 2
        [dum fInd] = min(abs(TFRwave.freq-median(R.bandef(band,:)))); % Take middle!
    else
        [dum fInd] = min(abs(TFRwave.freq-BB.cohfrq)); % Take coh peak!
    end
    fbwid = (median(R.bandef(band,:))/cfg.width)*2;
    disp([R.bandname{2} ' phase estimate at mid freq ' num2str(TFRwave.freq(fInd)) ' Hz +/- ' num2str(fbwid)])
    P = squeeze(TFRwave.fourierspctrm(1,:,fInd,:));
    
    % Now compute phase
    PhiTime = angle(P);
    period = 1.*data.fsample; overlap = 0.99;
    for ci = 1:size(PhiTime,1); PhiTime(ci,:) = unwrap(PhiTime(ci,:)); end
    
    % Sliding Window
    [slide_dphi_12,sind] = slideWindow(diff(PhiTime,1,1)', floor(period), floor(period*overlap));
    %     %     slide_dphi_12 = wrapToPi(slide_dphi_12);
    
    % PLV
    RPTime(band,:) = circ_mean(slide_dphi_12, [], 1);
    PLVTime(band,:) = abs(mean(exp(-1i*slide_dphi_12),1));
    % Get indices of SW centres
    sw_tvec = (round(median(sind,1)));
    % Convert to data time
    sw_tvec = data.time{1}(sw_tvec);
end

%% Now Wavlet Coherence
% % cfg = [];
% % cfg.method     = 'wavelet';
% % cfg.width      = 6;
% % cfg.output     = 'powandcsd';
% % cfg.foi        = 4:0.5:35;
% % cfg.toi        = data.time{1};
% %  cfg.pad = 'nextpow2';
% % TFRwave = ft_freqanalysis(cfg, data);


% % P2 =squeeze(TFRwave.fourierspctrm(1,:,fInd,:));
% % FX = squeeze(P2(1,:,:));
% % FY = squeeze(P2(2,:,:));
% % FXX = abs(FX.*FX);
% % FYY = abs(FY.*FY);
% % FXY = FY.*conj(FX);
% % CXY = (abs(FXY.*FXY).^2)./sqrt(FXX.*FYY);
%
% %
% % P2 =squeeze(TFRwave.powspctrm(:,fInd,:));
% % FXX = squeeze(P2(1,:,:));
% % FYY = squeeze(P2(2,:,:));
% % FXY = squeeze(TFRwave.crsspctrm(1,fInd,:));
% %
% % CXY = (abs(FXY).^2)./(FXX.*FYY);
%
% % wcoherence(data.trial{1}(1,:),data.trial{1}(2,:),data.fsample)
