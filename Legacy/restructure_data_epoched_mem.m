function restructure_data_epoched_mem(datapathr,subname,ref_chan)
% datapathr = 'C:\Users\Tim\Documents\Work\Cortical_Networks\Data\';
% subname = {'LM'};
% ref_chan = 'STN_R23';

i = 1; nr = 1
load([datapathr subname{i} '\ftdata\r' subname{i} '_LCMV_source'])
load([datapathr subname{i} '\ftdata\virtual_sources_' num2str(nr)])

ftdata_cell = cell(size(edgInd));
for x = 1:numel(edgInd)
    ftdata = [];
    ftdata.fsample = 200;
    ftdata.label = {sprintf('%.1f ',source.pos(x)); 'STN_23'};
    ftdata.trial = source_sens{x};
    ftdata.time = repmat({linspace(0,2.5,length(ftdata.trial{1}') )},1,size(ftdata.trial,2));
    ftdata_cell{x} = ftdata;
end

% cfg = [];
% cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
% ft_databrowser(cfg, ftdata_cell{6,4,7});
gASym = zeros(numel(ftdata_cell{1}.trial),numel(edgInd));
isum = sum(source.dim);
% delete(gcp('nocreate'))
% parpool()

% Compute connectivities
for x = 1:numel(edgInd)
    dataP = ftdata_cell{x};
%     ppm = ParforProgMon('Computing Source Space Metrics',size(cSurf,1));
    
    for t = 1:numel(dataP.trial)
        if size(dataP.trial{1},2)>1
            cfg = [];
            cfg.trials = t;
            data = ft_selectdata(cfg,dataP);
            
            cfg         = [];
            cfg.order   = 8;
            cfg.method = 'bsmart';
            mdata       = ft_mvaranalysis(cfg, data);
            
            cfg        = [];
            cfg.method = 'mvar';
            mfreq      = ft_freqanalysis(cfg, mdata);
            
            %             cfg           = [];
            %             cfg.method    = 'mtmfft';
            %             cfg.taper     = 'dpss';
            %             cfg.output    = 'fourier';
            %             cfg.pad       = 128;
            %             cfg.tapsmofrq = 1;
            %             freq          = ft_freqanalysis(cfg, ftdata_cell{x,y,z});
            
            cfg           = [];
            cfg.method    = 'coh';
            %             coh           = ft_connectivityanalysis(cfg, freq);
            cohm          = ft_connectivityanalysis(cfg, mfreq);
            %             cASym =
            %                                 figure
            %                                 cfg           = [];
            %                                 cfg.parameter = 'cohspctrm';
            %                                 cfg.zlim      = [0 1];
            %                                 ft_connectivityplot(cfg, cohm) %, coh);
            
            cfg           = [];
            cfg.method    = 'granger';
            granger       = ft_connectivityanalysis(cfg, mfreq);
            indz = find(granger.freq>=14 & granger.freq<30);
            gASym(t,x) = mean(squeeze(granger.grangerspctrm(2,1,indz) - granger.grangerspctrm(1,2,indz)));
            %             figure
            %             cfg           = [];
            %             cfg.parameter = 'grangerspctrm';
            %             cfg.zlim      = [0 0.1];
            %             ft_connectivityplot(cfg, granger);
        else
            gASym(t,x) = NaN;
        end
        disp(t); %ppm.increment();
    end
    disp(x)
end

