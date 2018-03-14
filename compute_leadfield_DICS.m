function compute_leadfield_DICS(datapathr,subname,steps,fresh,band)
if nargin<5
    band = [4 48 20 15]; % [lf hf cf ntap]
end

condname = {'ON','OFF'};
for i = 1:numel(subname)
mkdir([datapathr subname{i} '\images\sourcespace'])
    for cond = 1:2
    [datafileN,pp_mark,nrep,senscheck] = data_fileguide(subname{i},cond-1);
    for nr = 1:nrep
        datafilen = [pp_mark datafileN{nr}];
        %         if senscheck{nr} == 1
        %             nr = nr+1;
        %             if senscheck{nr} == 1
        %                 break
        %             end
        %         end
        if fresh == 1
            which_dir =  [datapathr subname{i} '\images'];
            dinfo = dir(which_dir);
            dinfo([dinfo.isdir]) = [];   %skip directories
            filenames = fullfile(which_dir, {dinfo.name});
            delete( filenames{:} )
        end
        load([datapathr subname{i} '\ftdata\meg_clean_data_' num2str(nr) '_' condname{cond}],'meg_clean_data')
        
        cfg = [];
        cfg.resamplefs = 128;
        megdata = ft_resampledata(cfg,meg_clean_data);
        
        cfg = [];
        cfg.lpfilter = 'yes';
        cfg.hpfilter = 'yes';
        cfg.lpfreq = band(2);
        cfg.hpfreq = band(1);
        megdata = ft_preprocessing(cfg,megdata);
        
        save([datapathr subname{i} '\ftdata\megdata_' num2str(nr) '_' condname{cond}],'megdata')
        
        % compute CSD
        load([datapathr subname{i} '\ftdata\megdata_' num2str(nr) '_' condname{cond}])
        
        cfg = [];
        cfg.method    = 'mtmfft';
        cfg.output    = 'powandcsd';
        cfg.tapsmofrq = band(4);
        cfg.foi   = band(3);
        cfg.keeptrials = 'yes';
        freqdata = ft_freqanalysis(cfg, megdata);
        
        
        cfg        = [];
        cfg.layout = 'CTF275_helmet.mat';
        cfg.ylim   = [14 30];
        cfg.parameter = 'powspctrm';
        ft_topoplotER(cfg, freqdata);
        
        save([datapathr subname{i} '\ftdata\freqdata_' num2str(nr) '_' condname{cond}],'freqdata')
%                 load([datapathr subname{i} '\ftdata\STN_ref_list'],'ref_list')
    ref_list = {'STN_R01','STN_R12','STN_R23','STN_L01','STN_L12','STN_L23'}

        for refN = 1:numel(ref_list)
            for meth = 1:numel(steps)
                switch steps{meth}
                    case 'leadfield'
                        ref_chan = ref_list{refN};
                        
                        %% Prepare Leadfield
                        load([datapathr subname{i} '\ftdata\coreg_headmodel'],'vol')
                        % compute leadfield
                        if senscheck{nr} == 1
                            sens = ft_read_sens([datapathr subname{i} '\' condname{cond} '\' datafilen],'senstype','meg');
                        else
                                [datafileNd] = data_fileguide(subname{i},0);
                                datafilend = [pp_mark datafileNd{nr}];

                            sens = ft_read_sens([datapathr subname{i} '\' condname{1} '\' datafilend],'senstype','meg');
%                             load([datapathr 'default_sens']);
                        end
                        
                        H = figure; shg
                        ft_plot_sens(sens,'style','*b'); %,'coil','true')
                        
                        hold on
                        ft_plot_vol(vol); shg
                        pause(3)
                        cfg                 = [];
%                         cfg.grad            = freqdata.grad;
                        cfg.grad            = sens; %% THIS ONE WORKED
%                         WELL
                        cfg.normalize       = 'yes';
                        cfg.headmodel       = vol;
                        % cfg.reducerank      = 2;
                        cfg.channel         = {'MEG',ref_chan};
                        cfg.grid.resolution = 1;   % use a 3-D grid with a 1 cm resolution
                        cfg.grid.unit       = 'cm';
                        grid = ft_prepare_leadfield(cfg);
                        
                        cfg              = [];
                        cfg.method       = 'dics';
                        cfg.grid         = grid;
                        cfg.headmodel    = vol;
                        cfg.refchan       = ref_chan;
                        cfg.frequency    = band(3);
                        cfg.dics.projectnoise = 'yes';
                        cfg.dics.lambda       = '5%';
                        cfg.dics.fixedori     = 'yes';
                        cfg.dics.keepfilter = 'yes';
                        cfg.dics.realfilter = 'yes';
                        source = ft_sourceanalysis(cfg,freqdata);
                        %                     Neural activity index ( avoid centre of head bias)
                        source.avg.pow = source.avg.pow ./ source.avg.noise;
                        
                        [cohmax, cohind] = max(source.avg.coh);
                        [d1 d2 d3] = ind2sub(size(source.avg.coh),cohmax);
                        DICS_peak = {[d1 d2 d3] cohind source.dim cohmax};
                        save([datapathr subname{i} '\ftdata\r' subname{i} '_DICS_peakinfo_source' condname{cond} 'nrep_' num2str(nr) '_' ref_chan],'DICS_peak')
                        save([datapathr subname{i} '\ftdata\r' subname{i} '_DICS_source' condname{cond} 'nrep_' num2str(nr) '_' ref_chan],'source','-v7.3')
                        load([datapathr subname{i} '\ftdata\r' subname{i} 'rs'],'mri')
                        %                 cfg = [];
                        %                 cfg.parameter = 'avg.pow';
                        %                 cfg.filename = [datapathr subname{i} '_tmp\MRI\source_data']
                        %                 ft_sourcewrite(cfg, source)
                        load([datapathr subname{i} '\ftdata\rs_transform_vox2ctf'])
                        
                        T = transform_vox2ctf; %/transform_vox2spm;%;
                        Tmri = ft_transform_geometry(T,mri);
                        Tmri = ft_convert_units(Tmri,'cm');
                        
                    case 'plotsourcepow'
                        
                        load([datapathr subname{i} '\ftdata\r' subname{i} '_DICS_source' condname{cond} 'nrep_' num2str(nr) '_' ref_chan])
                        load([datapathr subname{i} '\ftdata\rs_transform_vox2ctf'])
                        cfg            = [];
                        % cfg.downsample = 1;
                        cfg.parameter = 'avg.coh';
                        
                        sourceInt  = ft_sourceinterpolate(cfg, source , Tmri);
                        
                        %                 Z-score power
                        %                     sourceInt.coh = reshape(sourceInt.coh,size(sourceInt.anatomy));
                        %                     incoh = sourceInt.coh.*sourceInt.inside;
                        %                     incoh(incoh==0) = NaN;
                        %                     incoh = (incoh - nanmean(incoh(:)))/nanstd(incoh(:));
                        %                     %                                     inpow(isnan(inpow)) = 0;
                        %                     sourceInt.coh = incoh;
                        %                     sourceInt.coh = reshape(sourceInt.coh,[],1);
                        
                        cfg = [];
                        cfg.method        = 'slice';
                        cfg.funparameter  = 'avg.coh';
                        cfg.maskparameter = cfg.funparameter;
                        % cfg.funcolorlim   = [4.0 6.2];
                        % cfg.opacitylim    = [4.0 6.2];
                        cfg.opacitymap    = 'rampup';
                        figure
                        ft_sourceplot(cfg, sourceInt);
                        
                        sourceInt.coordsys = 'ctf';
                        cfg = [];
                        cfg.nonlinear     = 'no';
                        sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
%                         volumewrite_spm([datapathr subname{i} '\MRI\sourcespace\r' subname{i} '_' num2str(nr) '_DICS_sourcecoh_int_' ref_chan '.nii'], sourceIntNorm.coh, sourceIntNorm.transform, 'SPM12');
                        
                        cfg = [];
                        cfg.method         = 'surface';
                        cfg.funparameter   = 'avg.coh';
                        cfg.maskparameter  = cfg.funparameter;
                        % cfg.funcolorlim    = [0.0 1.2];
                        cfg.funcolormap    = 'jet';
%                         cfg.opacitylim     = [0. 0.05];
                        cfg.opacitymap     = 'auto';
                        cfg.projmethod     = 'nearest';
%                         cfg.downsample      = 2                 
                        cfg.surffile       = 'surface_white_both.mat';
%                         cfg.surfinflated   = 'surface_inflated_both_caret.mat';
                        ft_sourceplot(cfg, sourceIntNorm);
                        view ([0 90]); 
                        a = gca; a.Children(3).FaceAlpha = 0;
                         savefigure_v2([datapathr subname{i} '\images\sourcespace\'],['r' subname{i} '_rep_' num2str(nr) '_DICS_surf_avg_pow_zscore'  '_' condname{cond}],[],[],[]);
%                         savefig(gcf,[datapathr subname{i} '\images\sourcespace\r' subname{i} '_rep_' num2str(nr) '_surf_avg_DICS_coh_' ref_chan condname{cond}])
                close all
                end % reflist
            end
        end
        close all
    end
    end
end