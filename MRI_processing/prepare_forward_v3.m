function prepare_forward_v3(datapathr,subname,steps,fresh)
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
                %                 Compute Transforms
                cfg = [];
                cfg.method = 'fiducial';
                cfg.coordsys = 'ctf';
                cfg.fiducial.nas    = table2array(fidtab(1,2:4)); % position of nasion
                cfg.fiducial.lpa    = table2array(fidtab(2,2:4)); % position of LPA
                cfg.fiducial.rpa    = table2array(fidtab(3,2:4)); % position of RPA
                cfg.fiducial.zpoint = [45 45 75];
                mri = ft_volumerealign(cfg,mri);
                transform_spm2ctf = mri.transform;
                save([datapathr subname{i} '\ftdata\transform_spm2ctf'],'transform_spm2ctf')
            case 'skullstrip'
                %% Skullstrip
                if exist([datapathr subname{i} '\MRI\r' subname{i} '_ss.nii'])>0
                    delete([datapathr subname{i} '\MRI\r' subname{i} '_ss.nii']);
                end
                %                 system(['C:\Users\Tim\Documents\ROBEX\runROBEX.bat ' datapathr subname{i} '\MRI\orig\' mriname ' ' datapathr subname{i} '\MRI\r' subname{i} '_ss.nii'])
                system(['C:\Users\Tim\Documents\ROBEX\runROBEX.bat ' datapathr subname{i} '\MRI\orig\r' subname{i} '.nii'  ' ' datapathr subname{i} '\MRI\r' subname{i} '_ss.nii'])
                
            case 'segment'
                %% Segmentations
                iterNum_outer=15;  % outer iteration
                iterCM=2;  % inner interation for C and M
                iter_b=1;  % inner iteration for bias
                q = 1.5;   % fuzzifier
                th_bg = 2;  %% threshold for removing background
                N_region = 3; %% number of tissue types, e.g. WM, GM, CSF
                tissueLabel=[1, 2, 3];
                
                MICO_3Dseq({[datapathr subname{i} '\MRI\r' subname{i} '_ss.nii']}, N_region, q, th_bg, iterNum_outer, iter_b, iterCM, tissueLabel);
                mri = ft_read_mri([datapathr subname{i} '\MRI\r' subname{i} '_ss_seg.nii'],'dataformat','nifti_spm');
                
                % get segmentation
                mri.brain  = mri.anatomy>1;
                %                 mri = rmfield(mri,'anatomy');
                save([datapathr subname{i} '\ftdata\r' subname{i} 'ss_seg'],'mri','-v7.3')
            case 'reslice'
                load([datapathr subname{i} '\ftdata\r' subname{i} 'ss_seg'])
                
                
%                 fidtab = fiducial_reader([datapathr subname{i} '\MRI\orig\' subname{i} '_smri_fid.txt']);
%                 cfg = [];
%                 cfg.method = 'fiducial';
%                 cfg.coordsys = 'ctf';
%                 cfg.fiducial.nas    = table2array(fidtab(1,2:4)); % position of nasion
%                 cfg.fiducial.lpa    = table2array(fidtab(2,2:4)); % position of LPA
%                 cfg.fiducial.rpa    = table2array(fidtab(3,2:4)); % position of RPA
%                 cfg.fiducial.zpoint = [45 45 75];
%                 mri = ft_volumerealign(cfg,mri);
%                 mri.unit = 'mm'; 
%                 %% Reslice
%                 cfg = [];
%                 cfg.resolution = 1;
%                 cfg.dim = [256 256 256];
%                 mri = ft_volumereslice(cfg,mri);
                save([datapathr subname{i} '\ftdata\r' subname{i} '_rs'],'mri','-v7.3')
                volumewrite_spm([datapathr subname{i} '\MRI\r' subname{i} '_rs.nii'], mri.anatomy, mri.transform, 'SPM12');
            case 'headmodel'
                %% Prepare Headmodel
                load([datapathr subname{i} '\ftdata\r' subname{i} '_rs'])
                T = transform_spm2ctf;%;
                mri = ft_transform_geometry(T,mri);
                mri.coordsys = 'ctf';
                
                cfg = [];
                ft_sourceplot(cfg,mri)
                cfg = [];
                cfg.tissue = {'brain'};
                cfg.numvertices = [1200];
                mesh = ft_prepare_mesh(cfg,mri);
                
                cfg = [];
                cfg.tissue = {'brain'};
                cfg.method = 'singleshell'; %'openmeeg';
                vol = ft_prepare_headmodel(cfg,mesh);
                save([datapathr subname{i} '\ftdata\r' subname{i} '_headmodel'],'vol')
            case 'coregister'
                clear vol
                load([datapathr subname{i} '\ftdata\r' subname{i} '_headmodel'],'vol')
                %                 %% Apply coregistration
                %                 %                 load([datapathr subname{i} '\ftdata\transform_vox2ctf'])
                load([datapathr subname{i} '\ftdata\transform_spm2ctf'])
                T = transform_spm2ctf;%;
                vol = ft_transform_geometry(T,vol);
                vol = ft_convert_units(vol,'mm');
                [datafileN,pp_mark,nrep,senscheck] = data_fileguide(subname{i},1);
                for nr = 1:nrep
                    if senscheck{nr} == 1
                        datafilen = ['OFF\' pp_mark datafileN{1}];
                        save([datapathr subname{i} '\ftdata\coreg_headmodel'],'vol')
                        sens = ft_read_sens([datapathr subname{i} '\' datafilen],'senstype','meg');
                    else
                        load([datapathr 'default_sens']);
                    end
                    % PLOT
                    % close all
                    sens = ft_convert_units(sens,'cm');

                    H = figure; shg
                    ft_plot_sens(sens,'style','*b'); %,'coil','true')
                    
                    hold on
                    ft_plot_vol(vol)
                    close all
                    %%
                    az = 0;
                    el = 90;
                    view(az, el);
                    drawnow
                    pause(1.5);
                    
                    az = 0;
                    el = 0;
                    view(az, el);
                    drawnow
                    pause(1.5);
                    
                    az = 90;
                    el = 0;
                    view(az, el);
                    drawnow
                    pause(1.5);
                    
                    az = -165;
                    el = 25;
                    view(az, el);
                    drawnow
                    pause(1.5);
                    
                    az = -45;
                    el = -45;
                    view(az, el);
                    drawnow
                    pause(1.5);
                    %
                    %%
                    mkdir([datapathr subname{i} '\images\'])
                    savefig(H,[datapathr subname{i} '\images\r' subname{i} '_' num2str(nr) '_coreg'])
                end
        end
    end
    close all
end

