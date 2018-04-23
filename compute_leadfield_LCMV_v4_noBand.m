function compute_leadfield_LCMV_v4_noBand(R,steps,fresh)
rmpath(genpath('C:\spm12\external\fieldtrip\external\spm8'))
for sub =1:numel(R.subname)
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
                        
                        %                             cfg = [];
                        %                             cfg.bpfilter = 'yes';
                        %                             cfg.channel  = {'MEG'};
                        %                             cfg.bpfreq = R.bandef(band,:);
                        %                             megdata = ft_preprocessing(cfg,megdata);
                        
                        %                         save([R.datapathr R.subname{sub} '\ftdata\megdata_' num2str(nr) '_' R.condname{cond}],'megdata')
                        %                         load([R.datapathr R.subname{sub} '\ftdata\megdata_' num2str(nr) '_' R.condname{cond}])
                        
                        % Compute data covariance
                        cfg = [];
                        cfg.channel         = {'MEG'};
                        cfg.covariance = 'yes';
                        cfg.vartrllength       = 2;
                        covardata = ft_timelockanalysis(cfg,megdata);
                        
                        % compute leadfield
                        load([R.datapathr R.subname{sub} '\ftdata\volsens' R.condname{cond}])
                        sens = volsens.MEG.sens;
                        vol = volsens.MEG.vol;
                        sens = ft_convert_units(sens,'mm');
                        vol = ft_convert_units(vol,'mm');
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
                        load([R.datapathr 'template_MRI/standard_sourcemodel3d10mm']);
                        template_grid = sourcemodel;
                        clear sourcemodel;
                        
                        Tmri = ft_read_mri([R.datapathr R.subname{sub} '\MRI\orig\r' R.subname{sub} '.nii'],'dataformat','nifti_spm');
                        Tmri = ft_convert_units(Tmri,'mm');
                        Tmri.coordsys = 'spm';
                        
                        cfg                = [];
                        cfg.grid.warpmni   = 'yes';
                        cfg.grid.template  = template_grid;
                        cfg.grid.nonlinear = 'yes';
                        cfg.mri            = Tmri;
                        cfg.grid.unit      ='mm';
                        grid               = ft_prepare_sourcemodel(cfg);
                        
                        cfg              = [];
                        cfg.method       = 'lcmv';
                        cfg.lcmv.reducerank  = 2;
                        cfg.lcmv.lambda      = '0.01%';
                        cfg.grid         = grid;
                        cfg.headmodel    = vol;
                        cfg.normalize = 'yes';
                        cfg.lcmv.keepfilter = 'yes';
                        source = ft_sourceanalysis(cfg,covardata);
                        source.pos = grid.pos;
                        % Nueral activity index ( avoid centre of head bias)
                        %                         source.avg.pow = source.avg.pow ./ source.avg.noise;
                        save([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_LCMV_source' '_' R.condname{cond} 'nrep_' num2str(nr)],'source','-v7.3')
                        
                        
                    case 'plotsourcepow'
                        
                        %                         load([datapathr subname{sub} '\ftdata\r' subname{sub} '_LCMV_source' '_' condname{cond} 'nrep_' num2str(nr)])
                        %                         Tmri = ft_read_mri([datapathr subname{sub} '\MRI\orig\r' subname{sub} '.nii'],'dataformat','nifti_spm');
                        %                         Tmri = ft_convert_units(Tmri,'cm');
                        Tmri = ft_read_mri([R.datapathr R.subname{sub} '\MRI\orig\r' R.subname{sub} '.nii'],'dataformat','nifti_spm');
                        Tmri = ft_convert_units(Tmri,'mm');
                       
                        cfg            = [];
                        % cfg.downsample = 1;
                        cfg.parameter = 'avg.pow';
                        sourceInt  = ft_sourceinterpolate(cfg, source , Tmri);
                        sourceInt.coordsys = 'spm';
                        
                        %                 Z-score power
                        sourceInt.pow = reshape(sourceInt.pow,size(sourceInt.anatomy));
                        inpow = sourceInt.pow.*sourceInt.inside;
                        inpow(inpow==0) = NaN;
                                                inpow = (inpow - nanmean(inpow(:)))/nanstd(inpow(:));
%                                         inpow(isnan(inpow)) = 0;
                        sourceInt.pow = inpow;
                        sourceInt.pow = reshape(sourceInt.pow,[],1);
                        close all
                        
                        cfg = [];
                        cfg.method        = 'slice';
                        cfg.funparameter  = 'pow';
                        cfg.funcolorlim    = [prctile(inpow(:),1) prctile(inpow(:),99)];
                        cfg.maskparameter = cfg.funparameter;
                        % cfg.funcolorlim   = [4.0 6.2];
                        % cfg.opacitylim    = [4.0 6.2];
                        cfg.opacitymap    = 'rampup';
                        ft_sourceplot(cfg, sourceInt);
                        
                        rmpath(genpath('C:\spm12\external\fieldtrip\external\spm8'))
                        cfg = [];
                        cfg.nonlinear     = 'no';
                        cfg.spmversion  = 'spm12';
                        cfg.template    = [R.datapathr 'template_MRI\single_subj_T1_1mm.nii'];
                        sourceIntNorm = ft_volumenormalise(cfg, sourceInt);                        
                        %                     volumewrite_spm([datapathr subname{i} '\MRI\sourcespace\r' subname{i} '_' num2str(nr) '_sourcepow_int.nii'], sourceIntNorm.pow, sourceIntNorm.transform, 'SPM12');
                        %
                        cfg = [];
                        cfg.method         = 'surface';
                        cfg.funparameter   = 'avg.pow';
                        cfg.maskparameter  = cfg.funparameter;
                        cfg.funcolorlim    = [prctile(inpow(:),1) prctile(inpow(:),99)];
                        cfg.funcolormap    = 'jet';
                        %                         cfg.opacitylim     = [0. 0.05];
                        cfg.opacitymap     = 'auto';
                        cfg.projmethod     = 'nearest';
                        %                         cfg.downsample      = 2
                        cfg.surffile       = [R.datapathr 'template_MRI\surface_white_both.mat'];
                        %                         cfg.surfinflated   = 'surface_inflated_both_caret.mat';
                        ft_sourceplot(cfg, sourceIntNorm);
                        view ([0 90])
                        a = gca; a.Children(3).FaceAlpha = 0.1;
                        savefigure_v2([R.datapathr R.subname{sub} '\images\sourcespace\'],[R.subname{sub} '_LCMV_VLsrcmod_power'],[],[],[]);
                end
            end
            close all
        end
    end
end
