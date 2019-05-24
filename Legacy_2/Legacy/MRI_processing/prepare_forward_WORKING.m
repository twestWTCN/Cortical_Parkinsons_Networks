function prepare_forward(datapathr,subname,steps,fresh)
% nii raw MRI processing to fieldtrip output.
for i = 1:numel(subname)
    mriname = ['r' subname{i} '.nii'];
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
            case 'loadmri'
                cd([datapathr subname{i} '\MRI\'])
                mri = ft_read_mri([datapathr subname{i} '\MRI\orig\' mriname],'dataformat','nifti_spm');
                
            case 'realign_coords'
                %% Realignment
                % Think that the 'r' prefixed MRIs are in SPM space,
                % convert to CTF, try the unprefixed (again), or manual
                % specification of the fiducials in SPM space to convert to
                % CTF.
                fidtab = fiducial_reader([datapathr subname{i} '\MRI\orig\' subname{i} '_smri_fid.txt']);
                mri = ft_determine_coordsys(mri,'interactive','no');
                %                 Compute Transforms
                cfg = [];
                %                 cfg.method = 'spm'
                cfg.method = 'fiducial';
                cfg.coordsys = 'spm';
                cfg.viewresult     = 'yes'
                cfg.fiducial.nas    = table2array(fidtab(1,2:4)); % position of nasion
                cfg.fiducial.lpa    = table2array(fidtab(2,2:4)); % position of LPA
                cfg.fiducial.rpa    = table2array(fidtab(3,2:4)); % position of RPA
                cfg.fiducial.zpoint = [45 45 75];
                mri_spm = ft_volumerealign(cfg,mri);
                transform_vox2spm = mri_spm.transform;
                save([datapathr subname{i} '\ftdata\transform_vox2spm'],'transform_vox2spm')
                %                 save(['r' subname{1} '_ra'],'mri_spm')
                
                cfg = [];
                cfg.method = 'fiducial';
                cfg.coordsys = 'ctf';
                cfg.fiducial.nas    = table2array(fidtab(1,2:4)); % position of nasion
                cfg.fiducial.lpa    = table2array(fidtab(2,2:4)); % position of LPA
                cfg.fiducial.rpa    = table2array(fidtab(3,2:4)); % position of RPA
                cfg.fiducial.zpoint = [45 45 75];
                mri_ctf = ft_volumerealign(cfg,mri);
                transform_vox2ctf = mri_ctf.transform;
                save([datapathr subname{i} '\ftdata\transform_vox2ctf'],'transform_vox2ctf')
                if exist([datapathr subname{i} '\MRI\' mriname '.nii'])>0
                    delete([datapathr subname{i} '\MRI\' mriname '.nii']);
                end
                save([datapathr subname{i} '\MRI\' mriname],'mri')
                volumewrite_spm([datapathr subname{i} '\MRI\r' subname{i} '_ra.nii'], mri.anatomy, mri.transform, 'SPM12');
                
                cfg = [];
                cfg.location = table2array(fidtab(2,2:4));          %%% fill in a fiducial point from the current subject by hand, to check if coregistration went correctly
                cfg.locationcoordinates = 'head';
                ft_sourceplot(cfg,mri); 
                
            case 'skullstrip'
                
                %% Skullstrip
                if exist([datapathr subname{i} '\MRI\r' subname{i} '_ss.nii'])>0
                    delete([datapathr subname{i} '\MRI\r' subname{i} '_ss.nii']);
                end
                system(['C:\Users\Tim\Documents\ROBEX\runROBEX.bat ' datapathr subname{i} '\MRI\orig\' mriname ' ' datapathr subname{i} '\MRI\r' subname{i} '_ss.nii'])
                
            case 'segment'
                %% Segmentations
                iterNum_outer=15;  % outer iteration
                iterCM=2;  % inner interation for C and M
                iter_b=1;  % inner iteration for bias
                q = 1.5;   % fuzzifier
                th_bg = 5;  %% threshold for removing background
                N_region = 3; %% number of tissue types, e.g. WM, GM, CSF
                tissueLabel=[1, 2, 3];
                
                MICO_3Dseq({[datapathr subname{i} '\MRI\r' subname{i} '_ss.nii']}, N_region, q, th_bg, iterNum_outer, iter_b, iterCM, tissueLabel);
                
                mri = ft_read_mri([datapathr subname{i} '\MRI\r' subname{i} '_ss.nii'],'dataformat','nifti_spm');
                
                % get segmentation
                mri.brain  = mri.anatomy>0;
%                 mri = rmfield(mri,'anatomy');
                save([datapathr subname{i} '\ftdata\r' subname{i} 'ss_seg'],'mri')
                
            case 'reslice'
                load([datapathr subname{i} '\ftdata\r' subname{i} 'ss_seg'])
                %% Reslice
                cfg = [];
                cfg.resolution = 1;
                cfg.dim = [256 256 256];
                mri = ft_volumereslice(cfg,mri);
                save([datapathr subname{i} '\MRI\' mriname '_rs'],'mri')
                volumewrite_spm([datapathr subname{i} '\MRI\' mriname '_rs.nii'], mri.anatomy, mri.transform, 'SPM12');
                cfg = [];
                %                 cfg.method = 'spm'
                cfg.method = 'fiducial';
                cfg.coordsys = 'ctf';
                cfg.viewresult     = 'yes'
                cfg.fiducial.nas    = table2array(fidtab(1,2:4)); % position of nasion
                cfg.fiducial.lpa    = table2array(fidtab(2,2:4)); % position of LPA
                cfg.fiducial.rpa    = table2array(fidtab(3,2:4)); % position of RPA
                cfg.fiducial.zpoint = [45 45 75];
                mri_spm = ft_volumerealign(cfg,mri);
                transform_vox2ctf = mri_spm.transform;
                save([datapathr subname{i} '\ftdata\transform_vox2spm'],'transform_vox2spm')

            case 'headmodel'
                %% Prepare Headmodel
                load([datapathr subname{i} '\MRI\' mriname '_rs'])
                cfg = [];
                cfg.tissue = {'brain'};
                cfg.numvertices = [1200];
                mesh = ft_prepare_mesh(cfg,mri);
                
                cfg = [];
                cfg.tissue = {'brain'};
                cfg.method = 'singleshell'; %'openmeeg';
                vol = ft_prepare_headmodel(cfg,mri);
                save([datapathr subname{i} '\ftdata\r' subname{i} '_headmodel'],'vol')
            case 'coregister'
                clear vol
                load([datapathr subname{i} '\ftdata\r' subname{i} '_headmodel'],'vol')
                %% Apply coregistration
%                 load([datapathr subname{i} '\ftdata\transform_vox2ctf'])
                load([datapathr subname{i} '\ftdata\transform_vox2spm'])
                
                T = transform_vox2ctf;%;
                vol = ft_transform_geometry(T,vol);
                vol = ft_convert_units(vol,'cm');
                
                [datafilen,pp_mark] = data_fileguide(subname{i},1)
                datafilen = ['OFF\' pp_mark datafilen{1}];
%                 save([datapathr subname{i} '\ftdata\coreg_headmodel'],'vol2')
                sens = ft_read_sens([datapathr subname{i} '\' datafilen],'senstype','meg');
                % PLOT
                % close all
                figure
                ft_plot_sens(sens,'style','*b'); %,'coil','true')
                
                hold on
                ft_plot_vol(vol)
        end
    end
end

