function compute_leadfield_DICS_v4(R,steps,fresh)
rmpath(genpath(['C:\spm12\external\fieldtrip\external\spm8']))
for sub = 1:numel(R.subname)
    if fresh == 1 % Clear image folders
        for band = [1 3];
            which_dir =  [R.datapathr R.subname{sub} '\images\sourcespace\' R.bandname{band} '\'];
            dinfo = dir(which_dir);
            dinfo([dinfo.isdir]) = [];   %skip directories
            filenames = fullfile(which_dir, {dinfo.name});
            if ~isempty(filenames)
                delete( filenames{:} )
            end
        end
    end
    mkdir([R.datapathr R.subname{sub} '\images\sourcespace'])
    for cond = 1:2
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
        for nr = 1:nrep % Use the later trial
            for band = [1 3]
                % Load in preprocessed data
                load([R.datapathr R.subname{sub} '\ftdata\meg_clean_data_' num2str(nr) '_' R.condname{cond}],'meg_clean_data')
                %% Resample for efficiency
                cfg = [];
                cfg.resamplefs = 128;
                megdata = ft_resampledata(cfg,meg_clean_data);
                % Dont need to bandpass, this is done through
                % ft_freqanalysis
                %                 cfg = [];
                %                 cfg.lpfilter = 'yes';
                %                 cfg.hpfilter = 'yes';
                %                 cfg.lpfreq = R.dics.bp(band,2);
                %                 cfg.hpfreq = R.dics.bp(band,1);
                %                 megdata = ft_preprocessing(cfg,megdata);
                % save([R.datapathr R.subname{sub} '\ftdata\megdata_' num2str(nr) '_' R.condname{cond}],'megdata')
                
                %% compute CSD
                %  load([R.datapathr R.subname{sub} '\ftdata\megdata_' num2str(nr) '_' R.condname{cond}])
                cfg = [];
                cfg.method    = 'mtmfft';
                cfg.output    = 'powandcsd';
                cfg.tapsmofrq = R.dics.bp(band,4);
                cfg.foi   = R.dics.bp(band,3);
                cfg.keeptrials = 'yes';
                freqdata = ft_freqanalysis(cfg, megdata);
                
                % Optional topoplot for sensor space
                %             cfg        = [];
                %             cfg.layout = [R.datapathr 'template_MRI\CTF275_helmet.mat'];
                %             cfg.ylim   = R.dics.bp(1:2);
                %             cfg.parameter = 'powspctrm';
                %             ft_topoplotER(cfg, freqdata);
                
                save([R.datapathr R.subname{sub} '\ftdata\freqdata_' num2str(nr) '_' R.condname{cond} '_' R.bandname{band}],'freqdata')
                load([R.datapathr 'template_MRI/standard_sourcemodel3d10mm']);
                template_grid = sourcemodel;
                clear sourcemodel;
                % Get list of STN channels
                ref_list = {meg_clean_data.label{strmatch('STN',meg_clean_data.label)}};
                % Start looping through refs
                for refN =1:numel(ref_list)
                    for meth = 1:numel(steps)
                        switch steps{meth}
                            case 'leadfield'
                                ref_chan = ref_list{refN};
                                
                                %% Prepare Leadfield
                                %                             load([R.datapathr 'default_sens']);
                                load([R.datapathr R.subname{sub} '\ftdata\volsens' R.condname{cond}])
                                sens = volsens.MEG.sens;
                                vol = volsens.MEG.vol;
                                sens = ft_convert_units(sens,'mm');
                                vol = ft_convert_units(vol,'mm');
                                
                                % Plot sensors on model to check coreg
                                                                H = figure; shg
                                                                ft_plot_sens(sens,'style','*b'); %,'coil','true')
                                                                hold on
                                                                ft_plot_vol(vol); shg
                                %                                 %                                 savefig(H,[R.datapathr R.subname{sub} '\images\r' R.subname{sub} '_' num2str(nr) '_VLsrcmod_coreg']);
                                                                close all
                                
                                %% Prepare Leadfield
                                %                                 cfg                 = [];
                                %                                 cfg.grad            = sens;
                                %                                 cfg.normalize       = 'yes';
                                %                                 cfg.headmodel       = vol;
                                %                                 % cfg.reducerank      = 2;
                                %                                 cfg.channel         = {'MEG',ref_chan};
                                %                                 cfg.grid.resolution = 0.75;   % use a 3-D grid with a 1 cm resolution
                                %                                 cfg.grid.unit       = 'cm';
                                %                                 cfg.grid.pos = template_grid.pos;
                                %                                 cfg.grid.inside = template_grid.inside;
                                %                                 cfg.grid.dim = template_grid.dim;
                                %                                                                 cfg.grid.warpmni   = 'yes'; %%%
                                %                                                                 cfg.grid.template  = template_grid;
                                %                                                                 cfg.grid.nonlinear = 'yes';
                                %                                 grid = ft_prepare_leadfield(cfg);
                                
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
                                
                                %% Run DICs
                                cfg              = [];
                                cfg.method       = 'dics';
                                cfg.grid         = grid;
                                cfg.headmodel    = vol;
                                cfg.refchan       = ref_chan;
                                cfg.frequency    = freqdata.freq; %R.dics.bp(band,3);
                                cfg.dics.projectnoise = 'yes';
                                cfg.dics.lambda       = '0.01%';
                                cfg.dics.fixedori     = 'yes';
                                cfg.dics.keepfilter = 'yes';
                                cfg.dics.realfilter = 'yes';
                                source = ft_sourceanalysis(cfg,freqdata);
                                %                     Neural activity index ( avoid centre of head bias)
                                source.avg.pow = source.avg.pow ./ source.avg.noise;
                                source.pos = grid.pos;
                                [cohmax, cohind] = max(source.avg.coh);
                                [d1 d2 d3] = ind2sub(size(source.avg.coh),cohmax);
                                DICS_peak = {[d1 d2 d3] cohind source.dim cohmax};
                                save([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_DICSv2_peakinfo_source' R.condname{cond} 'nrep_' num2str(nr) '_' ref_chan '_' R.bandname{band}],'DICS_peak')
                                save([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_DICSv2_source' R.condname{cond} 'nrep_' num2str(nr) '_' ref_chan '_' R.bandname{band}],'source','-v7.3')
                            case 'plotsourcepow'
                                % Read in MRI and shift source in Z to
                                % match
                                Tmri = ft_read_mri([R.datapathr R.subname{sub} '\MRI\orig\r' R.subname{sub} '.nii'],'dataformat','nifti_spm');
                                Tmri = ft_convert_units(Tmri,'mm');
                                T =[0.9   0   0    0
                                     0   0.9  0    0
                                     0    0   1   2.2
                                     0    0   0    1];
                                source = ft_transform_geometry(T,source);
                                
%                                 Tmri = ft_read_mri([R.datapathr 'template_MRI\single_subj_T1_1mm.nii'],'dataformat','nifti_spm');
                                %% Interpolate functional to anatomical MRI
                                cfg            = [];
                                cfg.parameter = 'avg.coh';
                                sourceInt  = ft_sourceinterpolate(cfg, source , Tmri);
                                sourceInt.coordsys = 'spm';
                                
%                                 cfg            = [];
%                                 cfg.downsample = 5;
%                                 cfg.parameter = 'avg.coh';
%                                 sourceIntDS  = ft_sourceinterpolate(cfg, source , Tmri);
%                                 sourceIntDS.coordsys = 'mni';
                                
                                %% Plot source
                                cfg = [];
                                cfg.funparameter = 'avg.coh';
                                ft_sourceplot(cfg,sourceInt)
                                
                                %                 Z-score power
                                %                     sourceInt.coh = reshape(sourceInt.coh,size(sourceInt.anatomy));
                                %                     incoh = sourceInt.coh.*sourceInt.inside;
                                %                     incoh(incoh==0) = NaN;
                                %                     incoh = (incoh - nanmean(incoh(:)))/nanstd(incoh(:));
                                %                     %                                     inpow(isnan(inpow)) = 0;
                                %                     sourceInt.coh = incoh;
                                %                     sourceInt.coh = reshape(sourceInt.coh,[],1);
                                %                                 close all
                                %% Plot slice
                                cfg = [];
                                cfg.method        = 'slice';
                                cfg.funparameter  = 'avg.coh';
                                cfg.maskparameter = cfg.funparameter;
                                %                             cfg.funcolorlim    = [0 .3];
                                % cfg.opacitylim    = [4.0 6.2];
                                cfg.opacitymap    = 'rampup';
                                ft_sourceplot(cfg, sourceInt);
                                
                                rmpath(genpath('C:\spm12\external\fieldtrip\external\spm8'))
                                %% Now Normalise to Template MNI anatomical MRI
                                cfg = [];
                                cfg.nonlinear     = 'yes';
                                cfg.spmversion  = 'spm12';
                                cfg.template    = [R.datapathr 'template_MRI\single_subj_T1_1mm.nii'];
                                sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
                                % Transform interpolated source using
                                % Transformation matrix
%                                 sourceTrans = ft_transform_geometry(sourceIntNorm.transform,sourceIntDS);
                                
                                %                                 %% Downsample to make manageable
                                %                                 cfg = [];
                                %                                 cfg.downsample = 5;
                                %                                 cfg.smooth     = 2.5;
                                %                                 sourceIntNormDS = ft_volumedownsample(cfg,sourceTrans);
                                
                                % Plot 3D surface
                                cfg = [];
                                cfg.method         = 'surface';
                                cfg.funparameter   = 'avg.coh';
                                cfg.maskparameter  = cfg.funparameter;
                                %                             cfg.funcolorlim    = [0 .8];
                                cfg.funcolormap    = 'jet';
                                %                         cfg.opacitylim     = [0. 0.05];
                                cfg.opacitymap     = 'auto';
                                cfg.projmethod     = 'nearest';
                                %                         cfg.downsample      = 2
                                cfg.surffile       = [R.datapathr 'template_MRI\surface_white_both.mat'];
                                %                         cfg.surfinflated   = 'surface_inflated_both_caret.mat';
                                %                                 cfg.atlas =
                                % atlas = ft_read_atlas([R.datapathr 'template_MRI\ROI_MNI_V4.nii']);
                                % atlas.coordsys = 'spm';
                                % cfg.atlas =  atlas;
                                
                                ft_sourceplot(cfg, sourceIntNorm);
                                view ([0 90]);
                                a = gca; a.Children(3).FaceAlpha = 0.1;
                                
                                rmpath(genpath('C:\spm12\external\fieldtrip\external\spm8'))
                                savefigure_v2([R.datapathr R.subname{sub} '\images\sourcespace\' R.bandname{band} '\'],['r' R.subname{sub} '_rep_' num2str(nr) '_DICSv2_sourcecoh_' ref_chan '_' R.condname{cond} '_' R.bandname{band}],[],[],[]);
                                close all
                        end % reflist
                    end
                end
            end
        end
    end
end