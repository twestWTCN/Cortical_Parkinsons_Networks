% nii raw MRI processing to fieldtrip output.
clear; close all

addpath(genpath('C:\Users\Tim\Documents\Work\Cortical_Networks'))

datapathr = 'C:\Users\Tim\Documents\Work\Cortical_Networks\Data\';
subname = {'subA','subB'};

steps = {'clearmri','loadmri','reslice','realign_coords',...
    'skullstrip','segment','headmodel','coregister','leadfield','plotsourcepow'};

for i = 2 %:numel(subname)
    mriname = ['r' subname{i} '.nii'];
    for meth = 1:numel(steps)
        switch steps(meth)
            case 'cleanmri'
                % Start fresh
                which_dir = [datapathr subname{i} '_tmp\MRI\'];
                dinfo = dir(which_dir);
                dinfo([dinfo.isdir]) = [];   %skip directories
                filenames = fullfile(which_dir, {dinfo.name});
                if ~isempty(filenames)
                    delete( filenames{:} )
                end
                clear T vol vol2
                
            case 'loadmri'
                cd([datapathr subname{i} '_tmp\MRI\'])
                mri = ft_read_mri([datapathr subname{i} '_tmp\MRI\orig\' mriname],'dataformat','nifti_spm');
                
            case 'reslice'
                %% Reslice
                cfg = [];
                cfg.resolution = 1;
                cfg.dim = [256 256 256];
                mri = ft_volumereslice(cfg,mri);
                
                mriname = [mriname(1:end-4) '_rs'];
                save([datapathr subname{i} '_tmp\MRI\' mriname],'mri')
                volumewrite_spm([datapathr subname{i} '_tmp\MRI\' mriname '.nii'], mri.anatomy, mri.transform, 'SPM12');
                
            case 'realign_coords'
                %% Realignment
                mri = ft_read_mri([datapathr subname{i} '_tmp\MRI\' mriname '.nii'],'dataformat','nifti_spm');
                fidtab = fiducial_reader([datapathr subname{i} '_tmp\MRI\orig\' subname{i} '_smri_fid.txt']);
                
                mri = ft_determine_coordsys(mri,'interactive','no');
                % Compute Transforms
%                 cfg = [];
                % cfg.method = 'interactive'
                % cfg.method = 'fiducial';
                % cfg.coordsys = 'spm';
                % cfg.fiducial.nas    = table2array(fidtab(1,2:4)); % position of nasion
                % cfg.fiducial.lpa    = table2array(fidtab(2,2:4)); % position of LPA
                % cfg.fiducial.rpa    = table2array(fidtab(3,2:4)); % position of RPA
                % cfg.fiducial.zpoint = [45 45 75];
                % mri_spm = ft_volumerealign(cfg,mri);
                % transform_vox2spm = mri_spm.transform;
                % save(['r' subname{1} '_ra'],'mri_spm')
                
                cfg = [];
                cfg.method = 'fiducial';
                cfg.coordsys = 'ctf';
                cfg.fiducial.nas    = table2array(fidtab(1,2:4)); % position of nasion
                cfg.fiducial.lpa    = table2array(fidtab(2,2:4)); % position of LPA
                cfg.fiducial.rpa    = table2array(fidtab(3,2:4)); % position of RPA
                cfg.fiducial.zpoint = [45 45 75];
                mri_ctf = ft_volumerealign(cfg,mri);
                transform_vox2ctf = mri_ctf.transform;
                
                mriname = [mriname(1:end-4) '_ra'];
                if exist([datapathr subname{i} '_tmp\MRI\' mriname '.nii'])>0
                    delete([datapathr subname{i} '_tmp\MRI\' mriname '.nii']);
                end
                save([[datapathr subname{i} '_tmp\MRI\' mriname],'mri_ctf')
                volumewrite_spm([datapathr subname{i} '_tmp\MRI\r' subname{i} '_ra.nii'], mri.anatomy, mri.transform, 'SPM12');
            case 'skullstrip'
                
                %% Skullstrip
                if exist([datapathr subname{i} '_tmp\MRI\r' subname{i} '_ss.nii'])>0
                    delete([datapathr subname{i} '_tmp\MRI\r' subname{i} '_ss.nii']);
                end
                
                system(['C:\Users\Tim\Documents\ROBEX\runROBEX.bat ' datapathr subname{i} '_tmp\MRI\r' subname{i} '_ra.nii ' datapathr subname{i} '_tmp\MRI\r' subname{i} '_ss.nii'])
                
            case 'segment'
                %% Segmentations
                iterNum_outer=15;  % outer iteration
                iterCM=2;  % inner interation for C and M
                iter_b=1;  % inner iteration for bias
                q = 1.5;   % fuzzifier
                th_bg = 5;  %% threshold for removing background
                N_region = 3; %% number of tissue types, e.g. WM, GM, CSF
                tissueLabel=[1, 2, 3];
                
                MICO_3Dseq({[datapathr subname{i} '_tmp\MRI\r' subname{i} '_ss.nii']}, N_region, q, th_bg, iterNum_outer, iter_b, iterCM, tissueLabel);
                
                segmri = ft_read_mri([datapathr subname{i} '_tmp\MRI\r' subname{i} '_ss.nii'],'dataformat','nifti_spm');
                
                % get segmentation
                segmri.brain  = segmri.anatomy>0;
                segmri = rmfield(segmri,'anatomy');
                save([datapathr subname{i} '_tmp\MRI\r' subname{i} '_seg'],'segmri')
            case 'headmodel'
                %% Prepare Headmodel
                load([datapathr subname{i} '_tmp\MRI\r' subname{i} '_seg'])
                cfg = [];
                cfg.tissue = {'brain'};
                cfg.numvertices = [1200];
                mesh = ft_prepare_mesh(cfg,segmri);
                
                cfg = [];
                cfg.tissue = {'brain'};
                cfg.method = 'singleshell'; %'openmeeg';
                vol = ft_prepare_headmodel(cfg,segmri);
            case 'coregister'
                
                %% Apply coregistration
                T = transform_vox2ctf; %/transform_vox2spm;%;
                vol2 = ft_transform_geometry(T,vol);
                vol2 = ft_convert_units(vol2,'cm');
                
                datafilen = [subname{i} '_eData.ds'];
                
                sens = ft_read_sens([datapathr subname{i} '_tmp\' datafilen],'senstype','meg');
                % PLOT
                % close all
                figure
                ft_plot_sens(sens,'style','*b'); %,'coil','true')
                
                hold on
                ft_plot_vol(vol2)
        end
    end
end

