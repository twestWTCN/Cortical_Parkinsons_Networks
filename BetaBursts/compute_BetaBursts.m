function AmpTime = compute_BetaBursts(R,data)

for i = 1:size(data.trial{1},1)
    x = data.trial{1}(i,:);
    x = (x-mean(x))./std(x);
    data.trial{1}(i,:) = x;
end

cfg = [];
cfg.method     = 'wavelet';                
cfg.width      = 6; 
cfg.output     = 'pow';	
cfg.foi        = 4:0.5:35;	                
cfg.toi        = data.time{1};
 cfg.pad = 'nextpow2';
TFRwave = ft_freqanalysis(cfg, data);

% cfg = []
% cfg.channel = 'STN_L01';
% cfg.zlim = [0 60]
% ft_singleplotTFR(cfg, TFRwave);
% 
fInd = find(TFRwave.freq>=R.bandef(2,1) & TFRwave.freq<=R.bandef(2,2));

AmpTime = squeeze(mean(TFRwave.powspctrm(:,fInd,:),2));
% plot(TFRwave.time,AmpTime); ylim([0 60])