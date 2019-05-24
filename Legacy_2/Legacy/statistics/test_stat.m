load('C:\Users\Tim\Downloads\freqFIC_planar_cmb.mat')
load('C:\Users\Tim\Downloads\freqFC_planar_cmb.mat')
clear fakedata

D = reshape(A,7,1,7,12);
D = permute(D,[4 2 1 3]);
fakedata = freqFC_planar_cmb;
fakedata.label = 'A';
fakedata.freq  = 1:length(x)-1;
fakedata.time = 1:length(y)-1;
fakedata.powspctrm = D; %rand(12,1,7,7).*10;
fakedata = rmfield(fakedata,{'grad','cumtapcnt','cfg'});

D = reshape(B,7,1,7,12);
D = permute(D,[4 2 1 3]);
fakedata(2) =  fakedata(1);
fakedata(2).powspctrm = D; %rand(12,1,7,7); %D;

cfg = [];
cfg.latency          = 'all';
% cfg.frequency        = 20;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.1;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.minnbchan = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;

design = zeros(1,size(fakedata(1).powspctrm,1) + size(fakedata(2).powspctrm,1));
design(1,1:size(fakedata(1).powspctrm,1)) = 1;
design(1,(size(fakedata(1).powspctrm,1)+1):(size(fakedata(2).powspctrm,1)+...
    size(fakedata(2).powspctrm,1))) = 2;

cfg.design           = design;
cfg.ivar             = 1;

[stat] = ft_freqstatistics(cfg, fakedata(1), fakedata(2));

