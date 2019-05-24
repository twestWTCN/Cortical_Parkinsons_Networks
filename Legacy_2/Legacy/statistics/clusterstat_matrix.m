function stat = clusterstat_matrix(A,B,N,dimz)
% 24     1     5     6
if ndims(A) == 4 && ndims(B) == 4
    D = A;
else
D = reshape(A,dimz(1),1,dimz(2),dimz(3)*2);
D = permute(D,[4 2 1 3]);
end
D(D==0) = NaN;

% fakedata = freqFC_planar_cmb;
ABdata.label = 'A';
% ABdata.dimord = 'rpt_chan_freq_time';
ABdata.freq  = 1:size(A,1);
ABdata.time = 1:size(A,2);
ABdata.powspctrm = D; %rand(12,1,7,7).*10;
% fakedata = rmfield(fakedata,{'grad','cumtapcnt','cfg'});

if ndims(A) == 4 && ndims(B) == 4
    D = B;
else
D = reshape(B,dimz(1),1,dimz(2),dimz(3)*2);
D = permute(D,[4 2 1 3]);
end
D(D==0) = NaN;

ABdata(2) =  ABdata(1);
ABdata(2).powspctrm = D; %rand(12,1,7,7); %D;

cfg = [];
cfg.latency          = 'all';
% cfg.frequency        = 20;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';  % ft_statfun_indepsamplesT
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.2;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 0;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = N;

% if ndims(dimz)<2
design = zeros(1,size(ABdata(1).powspctrm,1) + size(ABdata(2).powspctrm,1));
design(1,1:size(ABdata(1).powspctrm,1)) = 1;
design(1,(size(ABdata(1).powspctrm,1)+1):(size(ABdata(2).powspctrm,1)+...
    size(ABdata(2).powspctrm,1))) = 2;
design = [reshape(repmat(1:dimz(3),2,1),1,[]) reshape(repmat(1:dimz(3),2,1),1,[]); design];
% end

cfg.design           = design;
cfg.ivar             = 2;
cfg.uvar             = 1;
[stat] = ft_freqstatistics(cfg, ABdata(1), ABdata(2));

