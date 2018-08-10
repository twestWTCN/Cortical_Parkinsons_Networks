function AmpTime = compute_BetaBursts_cohwip(R,data);

cfg = [];
cfg.method     = 'wavelet';                
cfg.width      = 6; 
cfg.output     = 'fourier';	
cfg.foi        = 4:0.5:35;	                
cfg.toi        = data.time{1};
 cfg.pad = 'nextpow2';
TFRwave = ft_freqanalysis(cfg, data);

x = squeeze(TFRwave.fourierspctrm(1,1,:,:));
xx = x.*conj(x);
y = squeeze(TFRwave.fourierspctrm(1,2,:,:));
yy = y.*conj(y);
xy = x.*conj(y);

coh = imag(xy./(xx.*yy)).^2;
TFRwave.powspctrm = reshape(coh,1,size(coh,1),size(coh,2));
%  powspctrm: [2×63×23873 double]
cfg = [];
cfg.method = 'coh';
TCR = ft_connectivityanalysis(cfg,TFRwave);

cfg = []
% cfg.parameter = 'cohspctrm';
% cfg.refchannel = 'STN_L01';
% cfg.channel = ;
cfg.zlim = [0 1]
ft_singleplotTFR(cfg, TFRwave);
% 
fInd = find(TFRwave.freq>=R.bandef(2,1) & TFRwave.freq<=R.bandef(2,2));

AmpTime = squeeze(mean(TFRwave.powspctrm(:,fInd,:),2));
% plot(TFRwave.time,AmpTime); ylim([0 60])