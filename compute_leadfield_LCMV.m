function compute_leadfield_LCMV(datapathr,subname,steps,fresh,band)
% fresh = 0; band = [6 40];
if nargin<5
    band = [6 40];
end
condname = {'ON','OFF'};
for i = 1:numel(subname)
    mkdir([datapathr subname{i} '\images\sourcespace'])
    for cond = 1:2
        mriname = ['r' subname{i} '.nii'];
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(subname{i},cond-1);
        for nr = 1:nrep
            datafilen = [pp_mark datafileN{nr}];
            %         if senscheck{nr} == 1
            %             nr = nr+1;
            %             if se nscheck{nr} == 1
            %                 break
            %             end
            %         end
            if fresh == 1
                which_dir =  [datapathr subname{i} '\ftdata'];
                dinfo = dir(which_dir);
                dinfo([dinfo.isdir]) = [];   %skip directories
                filenames = fullfile(which_dir, {dinfo.name});
                delete( filenames{:} )
                
                which_dir =  [datapathr subname{i} '\MRI'];
                dinfo = dir(which_dir);
                dinfo([dinfo.isdir]) = [];   %skip directories
                filenames = fullfile(which_dir, {dinfo.name});
                delete( filenames{:} )
            end
            for meth = 1:numel(steps)
                switch steps{meth}
                    case 'leadfield'
                        load([datapathr subname{i} '\ftdata\meg_clean_data_' num2str(nr) '_' condname{cond}],'meg_clean_data')
                        %% Prepare Leadfield
                        cfg = [];
                        cfg.resamplefs = 128;
                        megdata = ft_resampledata(cfg,meg_clean_data);
                        
                        cfg = [];
                        cfg.lpfilter = 'yes';
                        cfg.hpfilter = 'yes';
                        cfg.channel         = {'MEG'};
                        cfg.lpfreq = band(2);
                        cfg.hpfreq = band(1);
                        megdata = ft_preprocessing(cfg,megdata);
                        
                        save([datapathr subname{i} '\ftdata\megdata_' num2str(nr) '_' condname{cond}],'megdata')
                        
                        % compute CSD
                        load([datapathr subname{i} '\ftdata\megdata_' num2str(nr) '_' condname{cond}])
                        
                        %                     cfg = [];
                        %                 cfg.method    = 'mtmfft';
                        %                 cfg.output    = 'powandcsd';
                        %                 cfg.channel         = {'MEG'};
                        %                 cfg.tapsmofrq = 2;
                        %                 cfg.foilim    = [14 18];
                        %                 freqdata = ft_freqanalysis(cfg, megdata);
                        
                        % Compute data covariance
                        cfg = [];
                        cfg.channel         = {'MEG'};
                        cfg.covariance='yes';
                        cfg.vartrllength       = 2;
                        %                 cfg.covariancewindow = [-.3 .3];
                        covardata = ft_timelockanalysis(cfg,megdata);
                        
                        % compute leadfield
                        if senscheck{nr} == 1
                            sens = ft_read_sens([datapathr subname{i} '\' condname{cond} '\' datafilen],'senstype','meg');
                            load([datapathr subname{i} '\ftdata\coreg_headmodel'],'vol')
                        else
                            load([datapathr subname{i} '\ftdata\coreg_headmodel'],'vol')
                            load([datapathr 'default_sens']);
%                             load([datapathr subname{i} '\ftdata\volsens' condname{cond}])
%                             sens = volsens.MEG.sens;
%                             vol  = volsens.MEG.vol;
% %                             load([datapathr subname{i} '\ftdata\coreg_headmodel'],'vol')
%                             sens = ft_convert_units(sens,'m');
%                             Q =pi/2;
%                             T = [cos(Q) sin(Q) 0 0;
%                                 -sin(Q) cos(Q) 0 0;
%                                 0       0  1 0;
%                                 0       0  0 1];
%                             vol = ft_transform_geometry(T,vol);
                        end
                        H = figure; shg
                        ft_plot_sens(sens,'style','*b'); %,'coil','true')
                        
                        hold on
                        ft_plot_vol(vol)
                        pause(5)
                        
                        cfg                 = [];
                        cfg.grad            = sens;
                        cfg.normalize       = 'yes';
                        cfg.headmodel       = vol;
                        % cfg.reducerank      = 2;
                        cfg.channel         = {'MEG'};
                        cfg.grid.resolution = 1;   % use a 3-D grid with a 1 cm resolution
                        cfg.grid.unit       = 'cm';
                        grid = ft_prepare_leadfield(cfg);
                        
                        cfg              = [];
                        cfg.method       = 'lcmv';
                        cfg.lcmv.reducerank  = 2;
                        cfg.grid         = grid;
                        cfg.headmodel    = vol;
                        %                 cfg.frequency    = 18;
                        %                 cfg.dics.projectnoise = 'yes';
                        %                 cfg.dics.lambda       = '5%';
                        cfg.lcmv.keepfilter = 'yes';
                        source = ft_sourceanalysis(cfg,covardata);
                        % Nueral activity index ( avoid centre of head bias)
%                         source.avg.pow = source.avg.pow ./ source.avg.noise;
                        save([datapathr subname{i} '\ftdata\r' subname{i} '_LCMV_source' '_' condname{cond} 'nrep_' num2str(nr)],'source','-v7.3')
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
                        %                 load([datapathr subname{i} '\ftdata\r' subname{i} '_source'])
                        load([datapathr subname{i} '\ftdata\rs_transform_vox2ctf'])
                        cfg            = [];
                        % cfg.downsample = 1;
                        cfg.parameter = 'avg.pow';
                        sourceInt  = ft_sourceinterpolate(cfg, source , Tmri);
                        
                        %                 Z-score power
                        sourceInt.pow = reshape(sourceInt.pow,size(sourceInt.anatomy));
                        inpow = sourceInt.pow.*sourceInt.inside;
                        inpow(inpow==0) = NaN;
                        %                         inpow = (inpow - nanmean(inpow(:)))/nanstd(inpow(:));
                        %                 inpow(isnan(inpow)) = 0;
                        sourceInt.pow = inpow;
                        sourceInt.pow = reshape(sourceInt.pow,[],1);
                        
                        cfg = [];
                        cfg.method        = 'slice';
                        cfg.funparameter  = 'avg.pow';
                        cfg.maskparameter = cfg.funparameter;
                        % cfg.funcolorlim   = [4.0 6.2];
                        % cfg.opacitylim    = [4.0 6.2];
                        cfg.opacitymap    = 'rampup';
                        figure
                        ft_sourceplot(cfg, sourceInt);
                        
                        sourceInt.coordsys = 'ctf';
                        cfg = [];
                        cfg.nonlinear     = 'no';
                        cfg.spmversion  = 'spm8';
                        sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
                        %                     volumewrite_spm([datapathr subname{i} '\MRI\sourcespace\r' subname{i} '_' num2str(nr) '_sourcepow_int.nii'], sourceIntNorm.pow, sourceIntNorm.transform, 'SPM12');
                        %
                        cfg = [];
                        cfg.method         = 'surface';
                        cfg.funparameter   = 'avg.pow';
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
                        view ([0 90])
                        %                         a = gca; a.Children(3).FaceAlpha = 0;
                         savefigure_v2([datapathr subname{i} '\images\sourcespace\'],['r' subname{i} '_rep_' num2str(nr) '_LCMV_surf_avg_pow_zscore'  '_' condname{cond}],[],[],[]);
%                         savefig(gcf,[datapathr subname{i} '\images\sourcespace\r' subname{i} '_rep_' num2str(nr) '_LCMV_surf_avg_pow_zscore'  '_' condname{cond}]); close all
                end
            end
            close all
        end
    end
end