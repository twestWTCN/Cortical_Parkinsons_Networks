% nii raw MRI processing to fieldtrip output.
clear; close all
% spm eeg
% close all
addpath(genpath('C:\Users\Tim\Documents\Work\Cortical_Networks'))
% addpath(genpath('C:\Users\Tim\Documents\Work\Cortical_Networks'))

datapathr = 'C:\Users\Tim\Documents\Work\Cortical_Networks\Data\';
subname = {'subA','subB','subC'};

steps = {'clearMRIfolder','loadmri','realign_coords',...
    'skullstrip','segment','reslice','headmodel',...
   'coregister','preproc_MEG','leadfield','plotsourcepow'};
for i = 2 %:numel(subname)
    mriname = ['r' subname{i}];
    datafilen = [subname{i} '_eData.ds'];
    for meth = 1:numel(steps)
        switch steps{meth}
            case 'clearMRIfolder'
                disp(['Step ' num2str(meth) ', cleaning MRI folder.'])
                % Start fresh MRI
                which_dir = [datapathr subname{i} '_tmp\MRI\'];
                dinfo = dir(which_dir);
                dinfo([dinfo.isdir]) = [];   %skip directories
                filenames = fullfile(which_dir, {dinfo.name});
                if ~isempty(filenames)
                    delete( filenames{:} )
                end
                % Start fresh FTDATA
                which_dir = [datapathr subname{i} '_tmp\ftdata\'];
                dinfo = dir(which_dir);
                dinfo([dinfo.isdir]) = [];   %skip directories
                filenames = fullfile(which_dir, {dinfo.name});
                if ~isempty(filenames)
                    delete( filenames{:} )
                end
                
                clear T vol vol2
            case 'loadmri'
                disp(['Step ' num2str(meth) ', loading raw MRI.'])
                cd([datapathr subname{i} '_tmp\MRI\'])
                mri = ft_read_mri([datapathr subname{i} '_tmp\MRI\orig\' mriname '.nii'],'dataformat','nifti_spm');
                save([datapathr subname{i} '_tmp\MRI\' mriname],'mri')
                volumewrite_spm([datapathr subname{i} '_tmp\MRI\' mriname '.nii'], mri.anatomy, mri.transform, 'SPM12');
            case 'realign_coords'
                disp(['Step ' num2str(meth) ', getting coreg transform.'])
                %% Realignment
                load([datapathr subname{i} '_tmp\MRI\' mriname],'mri')
                fidtab = fiducial_reader([datapathr subname{i} '_tmp\MRI\orig\' subname{i} '_smri_fid.txt']);
                
                mri = ft_determine_coordsys(mri,'interactive','no');
                %                 Compute Transforms
                cfg = [];
                cfg.method = 'fiducial';
                cfg.coordsys = 'acpc';
                cfg.fiducial.nas    = table2array(fidtab(1,2:4)); % position of nasion
                cfg.fiducial.lpa    = table2array(fidtab(2,2:4)); % position of LPA
                cfg.fiducial.rpa    = table2array(fidtab(3,2:4)); % position of RPA
                cfg.fiducial.zpoint = [45 45 75];
                mri_spm = ft_volumerealign(cfg,mri);
                transform_vox2spm = mri_spm.transform;
                %                 save(['r' subname{1} '_ra'],'mri_spm')
                cfg = [];
                cfg.method = 'interactive';
                cfg.coordsys = 'ctf';
                %                 cfg.fiducial.nas    = table2array(fidtab(1,2:4)); % position of nasion
                %                 cfg.fiducial.lpa    = table2array(fidtab(2,2:4)); % position of LPA
                %                 cfg.fiducial.rpa    = table2array(fidtab(3,2:4)); % position of RPA
                %                 cfg.fiducial.zpoint = [45 45 75];
                mri_ctf = ft_volumerealign(cfg,mri);
                transform_vox2ctf = mri_ctf.transform;
                volumewrite_spm([datapathr subname{i} '_tmp\MRI\' mriname '_ra.nii'], mri_ctf.anatomy, mri_ctf.transform, 'SPM12');
                save([datapathr subname{i} '_tmp\ftdata\mri_ctf'],'mri_ctf')
                save([datapathr subname{i} '_tmp\ftdata\transform_vox2ctf'],'transform_vox2ctf')
            case 'skullstrip'
                disp(['Step ' num2str(meth) ', skullstripping MRI.'])
                %% Skullstrip
                if exist([datapathr subname{i} '_tmp\MRI\' mriname '_ss.nii'])>0
                    delete([datapathr subname{i} '_tmp\MRI\' mriname '_ss.nii']);
                end
                system(['C:\Users\Tim\Documents\ROBEX\runROBEX.bat ' datapathr subname{i} '_tmp\MRI\orig\' mriname '.nii ' datapathr subname{i} '_tmp\MRI\' mriname '_ss.nii'])
            case 'segment'
                disp(['Step ' num2str(meth) ', segmenting MRI.'])
                %% Segmentations
                iterNum_outer=15;  % outer iteration
                iterCM=2;  % inner interation for C and M
                iter_b=1;  % inner iteration for bias
                q = 1.5;   % fuzzifier
                th_bg = 5;  %% threshold for removing background
                N_region = 3; %% number of tissue types, e.g. WM, GM, CSF
                tissueLabel=[1, 2, 3];
                
                MICO_3Dseq({[datapathr subname{i} '_tmp\MRI\' mriname '_ss.nii']}, N_region, q, th_bg, iterNum_outer, iter_b, iterCM, tissueLabel);
                %                 mriname = [mriname '_seg'];
                segmri = ft_read_mri([datapathr subname{i} '_tmp\MRI\' mriname '_ss_seg.nii'],'dataformat','nifti_spm');
                
                % get segmentation
                segmri.brain  = segmri.anatomy>0;
                segmri = rmfield(segmri,'anatomy');
                save([datapathr subname{i} '_tmp\ftdata\' mriname 'ss_seg'],'segmri')
            case 'reslice'
                disp(['Step ' num2str(meth) ', reslicing MRI.'])
                
                load([datapathr subname{i} '_tmp\ftdata\' mriname 'ss_seg'])
                %% Reslice
                cfg = [];
                cfg.resolution = 1;
                cfg.dim = [256 256 256];
                segmri = ft_volumereslice(cfg,segmri);
                
                save([datapathr subname{i} '_tmp\ftdata\' mriname 'ss_seg_rs'],'segmri')
                volumewrite_spm([datapathr subname{i} '_tmp\MRI\' mriname 'ss_seg_rs.nii'], segmri.brain, segmri.transform, 'SPM12');
            case 'headmodel'
                disp(['Step ' num2str(meth) ', preparing headmodel.'])
                %% Prepare Headmodel
                load([datapathr subname{i} '_tmp\ftdata\' mriname 'ss_seg_rs'])
                cfg = [];
                cfg.tissue = {'brain'};
                cfg.numvertices = [1200];
                mesh = ft_prepare_mesh(cfg,segmri);
                
                cfg = [];
                cfg.tissue = {'brain'};
                cfg.method = 'singleshell'; %'openmeeg';
                vol = ft_prepare_headmodel(cfg,segmri);
                save([datapathr subname{i} '_tmp\ftdata\r' subname{i} '_headmodel'],'vol')
                
            case 'coregister'
                disp(['Step ' num2str(meth) ', applying coregistration.'])
                %% Apply coregistration
                close all
                load([datapathr subname{i} '_tmp\ftdata\r' subname{i} '_headmodel'])
                load([datapathr subname{i} '_tmp\ftdata\transform_vox2ctf']);
                T = transform_vox2ctf;%;
                vol2 = ft_transform_geometry(T,vol);
                vol2 = ft_convert_units(vol2,'cm');
                vol2.coordsys = 'ctf';
%                 sens = ft_read_sens([datapathr subname{i} '_tmp\' datafilen],'senstype','meg');
                sens = ft_read_sens([datapathr subname{1} '_tmp\' subname{1} '_eData.ds'],'senstype','meg'); % use subject 1 sensor locations
                save([datapathr subname{i} '_tmp\ftdata\r' subname{i} '_coreg_sens_vol'],'sens','vol2')
                
                % PLOT
                % close all
                figure
                ft_plot_sens(sens,'style','*b'); %,'coil','true')
                
                hold on
                ft_plot_vol(vol2)
                drawnow
            case 'preproc_MEG'
                disp(['Step ' num2str(meth) ', preprocessing MEG data.'])
                %% Prepare Leadfield
                cfg = [];
                cfg.datafile = [datapathr subname{i} '_tmp\' datafilen];
                cfg.lpfilter = 'yes';
                cfg.hpfilter = 'yes';
                cfg.channel         = {'MEG'};
                cfg.lpfreq = 24;
                cfg.hpfreq = 5;
                megdata = ft_preprocessing(cfg);
                
                cfg = [];
                cfg.resamplefs = 400;
                megdata = ft_resampledata(cfg,megdata);
                save([datapathr subname{i} '_tmp\ftdata\megdata'],'megdata')
                
            case 'leadfield'
                disp(['Step ' num2str(meth) ', preparing leadfield.'])
                % compute CSD
                load([datapathr subname{i} '_tmp\ftdata\megdata'])
                
                cfg = [];
                cfg.method    = 'mtmfft';
                cfg.output    = 'powandcsd';
                cfg.channel         = {'MEG'};
                cfg.tapsmofrq = 2;
                cfg.foilim    = [4 30];
                freqdata = ft_freqanalysis(cfg, megdata);
                
                
                load([datapathr subname{i} '_tmp\ftdata\r' subname{i} '_coreg_sens_vol'])
                
                % compute leadfield
                cfg                 = [];
                cfg.grad            = freqdata.grad;
                cfg.normalize       = 'yes';
                cfg.headmodel       = vol2;
                % cfg.reducerank      = 2;
                cfg.channel         = {'MEG'};
                cfg.grid.resolution = 0.5;   % use a 3-D grid with a 1 cm resolution
                cfg.grid.unit       = 'cm';
                grid = ft_prepare_leadfield(cfg);
                
                cfg              = [];
                cfg.method       = 'dics';
                cfg.frequency    = 24;
                cfg.grid         = grid;
                cfg.headmodel    = vol2;
                cfg.dics.projectnoise = 'yes';
                cfg.dics.lambda       = '5%';
                cfg.keepfilter = 'yes';
                source = ft_sourceanalysis(cfg,freqdata);
                
                % Nueral activity index ( avoid centre of head bias)
                source.avg.pow = source.avg.pow ./ source.avg.noise;
                save([datapathr subname{i} '_tmp\ftdata\source_data'],'source')
                
                
            case 'plotsourcepow'
                disp(['Step ' num2str(meth) ', plotting source power.'])
                load([datapathr subname{i} '_tmp\ftdata\r' subname{i} '_coreg_sens_vol'])
                load([datapathr subname{i} '_tmp\ftdata\transform_vox2ctf']);
                load([datapathr subname{i} '_tmp\ftdata\mri_ctf'])
                load([datapathr subname{i} '_tmp\ftdata\source_data'])
                
                T = transform_vox2ctf; %/transform_vox2spm;%;
                Tmri = ft_convert_units(mri_ctf,'cm');
                
                cfg            = [];
                % cfg.downsample = 1;
                cfg.parameter = 'avg.pow';
                sourceInt  = ft_sourceinterpolate(cfg, source , Tmri);
                
                cfg = [];
                cfg.method        = 'slice';
                cfg.funparameter  = 'avg.pow';
                cfg.maskparameter = cfg.funparameter;sourc
                % cfg.funcolorlim   = [4.0 6.2];
                % cfg.opacitylim    = [4.0 6.2];
                cfg.opacitymap    = 'rampup';
                figure
                ft_sourceplot(cfg, sourceInt);
                
                sourceInt.coordsys = 'ctf';
                cfg = [];
                cfg.nonlinear     = 'no';
                sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
                volumewrite_spm([datapathr subname{i} '_tmp\MRI\r' subname{i} '_sourcepow.nii'], sourceIntNorm.pow, sourceIntNorm.transform, 'SPM12');
                
                cfg = [];
                cfg.method         = 'surface';
                cfg.funparameter   = 'avg.pow';
                cfg.maskparameter  = cfg.funparameter;
                % cfg.funcolorlim    = [0.0 1.2];
                cfg.funcolormap    = 'jet';
                % cfg.opacitylim     = [0.0 1.2];
                cfg.opacitymap     = 'auto';
                cfg.projmethod     = 'nearest';
                cfg.surffile       = 'surface_white_both.mat';
                cfg.surfdownsample = 2;
                ft_sourceplot(cfg, sourceIntNorm);
                view ([90 0])
        end
    end
end

