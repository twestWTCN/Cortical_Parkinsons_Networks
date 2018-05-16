function compute_leadfield_LCMV_v2(R,steps,fresh)
for sub = 1:numel(R.subname)
    mkdir([R.datapathr R.subname{sub} '\images\sourcespace'])
    for cond = 1:2
        mriname = ['r' R.subname{sub} '.nii'];
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
        for nr = 1:nrep
            datafilen = [pp_mark datafileN{nr}];
            if fresh == 1
                which_dir =  [R.datapathr R.subname{sub} '\ftdata'];
                dinfo = dir(which_dir);
                dinfo([dinfo.isdir]) = [];   %skip directories
                filenames = fullfile(which_dir, {dinfo.name});
                delete( filenames{:} )
                
                which_dir =  [R.datapathr R.subname{sub} '\MRI'];
                dinfo = dir(which_dir);
                dinfo([dinfo.isdir]) = [];   %skip directories
                filenames = fullfile(which_dir, {dinfo.name});
                delete( filenames{:} )
            end
            for meth = 1:numel(steps)
                switch steps{meth}
                    case 'leadfield'
                        load([R.datapathr R.subname{sub} '\ftdata\meg_clean_data_' num2str(nr) '_' R.condname{cond}],'meg_clean_data')
                        %% Prepare Leadfield
                        cfg = [];
                        cfg.resamplefs = 128;
                        megdata = ft_resampledata(cfg,meg_clean_data);
                        
                        %                         cfg = [];
                        %                         cfg.lpfilter = 'yes';
                        %                         cfg.hpfilter = 'yes';
                        %                         cfg.channel         = {'MEG'};
                        %                         cfg.lpfreq = band(2);
                        %                         cfg.hpfreq = band(1);
                        %                         megdata = ft_preprocessing(cfg,megdata);
                        
                        save([R.datapathr R.subname{sub} '\ftdata\megdata_' num2str(nr) '_' R.condname{cond}],'megdata')
                        load([R.datapathr R.subname{sub} '\ftdata\megdata_' num2str(nr) '_' R.condname{cond}])
                        
                        % Compute data covariance
                        cfg = [];
                        cfg.channel         = {'MEG'};
                        cfg.covariance='yes';
                        cfg.vartrllength       = 2;
                        %                 cfg.covariancewindow = [-.3 .3];
                        covardata = ft_timelockanalysis(cfg,megdata);
                        
                        % compute leadfield
                        load([R.datapathr R.subname{sub} '\ftdata\volsens' R.condname{cond}])
                        sens = volsens.MEG.sens;
                        vol = volsens.MEG.vol;
                        sens = ft_convert_units(sens,'cm');
                        vol = ft_convert_units(vol,'cm');
                        %                             Q =-pi/2;
                        %                             T = [cos(Q) sin(Q) 0 0;
                        %                                 -sin(Q) cos(Q) 0 0;
                        %                                   0       0  1 0;
                        %                                 0       0  0 1];
                        %                             vol = ft_transform_geometry(T,sens);
                        %                         end
                        H = figure; shg
                        ft_plot_sens(sens,'style','*b');  %,'coil','true')
                        
                        hold on
                        ft_plot_vol(vol)
                        savefig(H,[R.datapathr R.subname{sub} '\images\r' R.subname{sub} '_' num2str(nr) '_VLsrcmod_coreg']); close all
                        %                         savefigure_v2([datapathr subname{1} '\images\sourcespace\'],[subname{1} '_VLsrcmod_coreg'],[],[],[]);
                        %                         pause
                        
                        cfg                 = [];
                        cfg.grad            = sens;
                        cfg.normalize       = 'yes';
                        cfg.headmodel       = vol;
                        % cfg.reducerank      = 2;
                        cfg.channel         = {'MEG'};
                        cfg.grid.resolution = 0.75;   % use a 3-D grid with a 1 cm resolution
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
                        save([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_LCMV_source' '_' R.condname{cond} 'nrep_' num2str(nr)],'source','-v7.3')
                        
                        
                    case 'plotsourcepow'
                        
                        %                         load([datapathr subname{sub} '\ftdata\r' subname{sub} '_LCMV_source' '_' condname{cond} 'nrep_' num2str(nr)])
                        %                         Tmri = ft_read_mri([datapathr subname{sub} '\MRI\orig\r' subname{sub} '.nii'],'dataformat','nifti_spm');
                        %                         Tmri = ft_convert_units(Tmri,'cm');
                        Tmri = ft_read_mri([R.datapathr 'template_MRI\single_subj_T1_1mm.nii'],'dataformat','nifti_spm');
                        Tmri = ft_convert_units(Tmri,'cm');
                        %                 cfg = [];
                        %                 cfg.parameter = 'avg.pow';
                        %                 cfg.filename = [datapathr subname{i} '_tmp\MRI\source_data']
                        %                 ft_sourcewrite(cfg, source)
                        %                         load([datapathr subname{i} '\ftdata\rs_transform_vox2ctf'])
                        %                         load([datapathr subname{i} '\ftdata\transform_vox2spm'])
                        %                         T = volsens.transforms.toMNI; %/transform_vox2spm;%;
                        %                         Tmri = ft_transform_geometry(T,Tmri);
                        %                                                 T = [1  0   0    0
                        %                                                     0   1   0    0
                        %                                                     0   0   1    2
                        %                                                     0   0   0    1];
                        %                                                     source = ft_transform_geometry(T,source);
                        
                        cfg            = [];
                        % cfg.downsample = 1;
                        cfg.parameter = 'avg.pow';
                        sourceInt  = ft_sourceinterpolate(cfg, source , Tmri);
                        sourceInt.coordsys = 'mni';
                        
                        %                 Z-score power
                        %                         sourceInt.pow = reshape(sourceInt.pow,size(sourceInt.anatomy));
                        %                         inpow = sourceInt.pow.*sourceInt.inside;
                        %                         inpow(inpow==0) = NaN;
                        %                         %                         inpow = (inpow - nanmean(inpow(:)))/nanstd(inpow(:));
                        %                         %                 inpow(isnan(inpow)) = 0;
                        %                         sourceInt.pow = inpow;
                        %                         sourceInt.pow = reshape(sourceInt.pow,[],1);
                        close all
                        cfg = [];
                        cfg.method        = 'slice';
                        cfg.funparameter  = 'avg.pow';
                        cfg.maskparameter = cfg.funparameter;
                        % cfg.funcolorlim   = [4.0 6.2];
                        % cfg.opacitylim    = [4.0 6.2];
                        cfg.opacitymap    = 'rampup';
                        figure
                        ft_sourceplot(cfg, sourceInt);
                        
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
                        cfg.funcolorlim    = [0.0 1.4];
                        cfg.funcolormap    = 'jet';
                        %                         cfg.opacitylim     = [0. 0.05];
                        cfg.opacitymap     = 'auto';
                        cfg.projmethod     = 'nearest';
                        %                         cfg.downsample      = 2
                        cfg.surffile       = [R.datapathr 'template_MRI\surface_white_both.mat'];
                        %                         cfg.surfinflated   = 'surface_inflated_both_caret.mat';
                        ft_sourceplot(cfg, sourceIntNorm);
                        view ([0 90])
%                         a = gca; a.Children(3).FaceAlpha = 0;
                        savefigure_v2([R.datapathr R.subname{sub} '\images\sourcespace\'],[R.subname{sub} '_LCMV_VLsrcmod_power'],[],[],[]);
                end
            end
            close all
        end
    end
end