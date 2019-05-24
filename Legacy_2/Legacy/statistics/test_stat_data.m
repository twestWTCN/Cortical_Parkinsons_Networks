load('C:\Users\Tim\Downloads\freqFIC_planar_cmb.mat')
load('C:\Users\Tim\Downloads\freqFC_planar_cmb.mat')
clear fakedata

fakedata = freqFC_planar_cmb;
fakedata.label = 'A';
fakedata.freq  = 1:6;
fakedata.time = 1:7;
fakedata.powspctrm = rand(6,1,6,7);
fakedata = rmfield(fakedata,{'grad','cumtapcnt','cfg'})

fakedata(2) =  fakedata(1);
fakedata(2).powspctrm = rand(6,1,6,7);

cfg = [];
cfg.latency          = 'all';
% cfg.frequency        = 20;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;

design = zeros(1,size(fakedata(1).powspctrm,1) + size(fakedata(2).powspctrm,1));
design(1,1:size(fakedata(1).powspctrm,1)) = 1;
design(1,(size(fakedata(1).powspctrm,1)+1):(size(fakedata(2).powspctrm,1)+...
    size(fakedata(2).powspctrm,1))) = 2;

cfg.design           = design;
cfg.ivar             = 1;

[stat] = ft_freqstatistics(cfg, fakedata(1), fakedata(2));