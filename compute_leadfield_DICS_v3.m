function compute_leadfield_DICS_v3(R,steps,fresh)
for sub = 2:numel(R.subname)
    mkdir([R.datapathr R.subname{sub} '\images\sourcespace'])
    for cond = 1:2
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
        for nr = 1:nrep % Use the later trial
            for band = 1:2
                datafilen = [pp_mark datafileN{nr}];
                %         if senscheck{nr} == 1
                %             nr = nr+1;
                %             if senscheck{nr} == 1
                %                 break
                %             end
                %         end
                if fresh == 1
                    which_dir =  [R.datapathr R.subname{sub} '\images'];
                    dinfo = dir(which_dir);
                    dinfo([dinfo.isdir]) = [];   %skip directories
                    filenames = fullfile(which_dir, {dinfo.name});
                    delete( filenames{:} )
                end
                load([R.datapathr R.subname{sub} '\ftdata\meg_clean_data_' num2str(nr) '_' R.condname{cond}],'meg_clean_data')
                
                cfg = [];
                cfg.resamplefs = 128;
                megdata = ft_resampledata(cfg,meg_clean_data);
                
                cfg = [];
                cfg.lpfilter = 'yes';
                cfg.hpfilter = 'yes';
                cfg.lpfreq = R.dics.bp(band,2);
                cfg.hpfreq = R.dics.bp(band,1);
                megdata = ft_preprocessing(cfg,megdata);
                
                save([R.datapathr R.subname{sub} '\ftdata\megdata_' num2str(nr) '_' R.condname{cond}],'megdata')
                
                % compute CSD
                load([R.datapathr R.subname{sub} '\ftdata\megdata_' num2str(nr) '_' R.condname{cond}])
                
                cfg = [];
                cfg.method    = 'mtmfft';
                cfg.output    = 'powandcsd';
                cfg.tapsmofrq = R.dics.bp(band,4);
                cfg.foi   = R.dics.bp(band,3);
                cfg.keeptrials = 'yes';
                freqdata = ft_freqanalysis(cfg, megdata);
                
                
                %             cfg        = [];
                % %             cfg.layout = [R.datapathr 'template_MRI\CTF275_helmet.mat'];
                %             cfg.ylim   = R.dics.bp(1:2);
                %             cfg.parameter = 'powspctrm';
                %             ft_topoplotER(cfg, freqdata);
                
                save([R.datapathr R.subname{sub} '\ftdata\freqdata_' num2str(nr) '_' R.condname{cond} '_' R.bandname{band}],'freqdata')
                %                 load([R.datapathr R.subname{i} '\ftdata\STN_ref_list'],'ref_list')
                %             ref_list = {'STN_R01','STN_R12','STN_R23','STN_L01','STN_L12','STN_L23'}
                ref_list = {meg_clean_data.label{strmatch('STN',meg_clean_data.label)}};
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
                                sens = ft_convert_units(sens,'cm');
                                vol = ft_convert_units(vol,'cm');
                                
                                H = figure; shg
                                ft_plot_sens(sens,'style','*b'); %,'coil','true')
                                
                                hold on
                                ft_plot_vol(vol); shg
%                                 savefig(H,[R.datapathr R.subname{sub} '\images\r' R.subname{sub} '_' num2str(nr) '_VLsrcmod_coreg']); 
                            close all
                                
                                cfg                 = [];
                                %                         cfg.grad            = freqdata.grad;
                                cfg.grad            = sens; %% THIS ONE WORKED
                                %                         WELL
                                cfg.normalize       = 'yes';
                                cfg.headmodel       = vol;
                                % cfg.reducerank      = 2;
                                cfg.channel         = {'MEG',ref_chan};
                                cfg.grid.resolution = 0.75;   % use a 3-D grid with a 1 cm resolution
                                cfg.grid.unit       = 'cm';
                                grid = ft_prepare_leadfield(cfg);
                                
                                cfg              = [];
                                cfg.method       = 'dics';
                                cfg.grid         = grid;
                                cfg.headmodel    = vol;
                                cfg.refchan       = ref_chan;
                                cfg.frequency    = freqdata.freq; %R.dics.bp(band,3);
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
                                save([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_DICSv2_peakinfo_source' R.condname{cond} 'nrep_' num2str(nr) '_' ref_chan '_' R.bandname{band}],'DICS_peak')
                                save([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_DICSv2_source' R.condname{cond} 'nrep_' num2str(nr) '_' ref_chan '_' R.bandname{band}],'source','-v7.3')
                                
                            case 'plotsourcepow'
                                load([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_DICSv2_source' R.condname{cond} 'nrep_' num2str(nr) '_' ref_chan])
                                Tmri = ft_read_mri([R.datapathr R.subname{sub} '\MRI\orig\r' R.subname{sub} '.nii'],'dataformat','nifti_spm');
                                Tmri = ft_convert_units(Tmri,'cm');
                                
                                T = [1  0   0    0
                                    0   1   0    0
                                    0   0   1    2.2
                                    0   0   0    1];
                                source = ft_transform_geometry(T,source);
                                cfg            = [];
                                % cfg.downsample = 1;
                                cfg.parameter = 'avg.coh';
                                
                                sourceInt  = ft_sourceinterpolate(cfg, source , Tmri);
                                cfg = [];
                                cfg.funparameter = 'avg.coh';
                                ft_sourceplot(cfg,sourceInt)
                                
                                sourceInt.coordsys = 'mni';                            %                 Z-score power
                                %                     sourceInt.coh = reshape(sourceInt.coh,size(sourceInt.anatomy));
                                %                     incoh = sourceInt.coh.*sourceInt.inside;
                                %                     incoh(incoh==0) = NaN;
                                %                     incoh = (incoh - nanmean(incoh(:)))/nanstd(incoh(:));
                                %                     %                                     inpow(isnan(inpow)) = 0;
                                %                     sourceInt.coh = incoh;
                                %                     sourceInt.coh = reshape(sourceInt.coh,[],1);
                                close all
                                cfg = [];
                                cfg.method        = 'slice';
                                cfg.funparameter  = 'avg.coh';
                                cfg.maskparameter = cfg.funparameter;
                                %                             cfg.funcolorlim    = [0 .3];
                                % cfg.opacitylim    = [4.0 6.2];
                                cfg.opacitymap    = 'rampup';
                                figure
                                ft_sourceplot(cfg, sourceInt);
                                
                                sourceInt.coordsys = 'ctf';
                                %                             cfg = [];
                                %                             cfg.nonlinear     = 'no';
                                %                             cfg.spmversion  = 'spm12';
                                %                             sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
                                %                             %                         volumewrite_spm([R.datapathr R.subname{i} '\MRI\sourcespace\r' R.subname{i} '_' num2str(nr) '_DICS_sourcecoh_int_' ref_chan '.nii'], sourceIntNorm.coh, sourceIntNorm.transform, 'SPM12');
                                %
                                %                             cfg = [];
                                %                             cfg.method         = 'surface';
                                %                             cfg.funparameter   = 'avg.coh';
                                %                             cfg.maskparameter  = cfg.funparameter;
                                % %                             cfg.funcolorlim    = [0 .8];
                                %                             cfg.funcolormap    = 'jet';
                                %                             %                         cfg.opacitylim     = [0. 0.05];
                                %                             cfg.opacitymap     = 'auto';
                                %                             cfg.projmethod     = 'nearest';
                                %                             %                         cfg.downsample      = 2
                                %                         cfg.surffile       = [R.datapathr 'template_MRI\surface_white_both.mat'];
                                %                             %                         cfg.surfinflated   = 'surface_inflated_both_caret.mat';
                                %                             ft_sourceplot(cfg, sourceIntNorm);
                                %                             view ([0 90]);
                                %                             a = gca; a.Children(3).FaceAlpha = 0;
                                savefigure_v2([R.datapathr R.subname{sub} '\images\sourcespace\' R.bandname{band} '\'],['r' R.subname{sub} '_rep_' num2str(nr) '_DICSv2_sourcecoh_' ref_chan '_' R.condname{cond} '_' R.bandname{band}],[],[],[]);
                                %                         savefig(gcf,[R.R.datapathr R.R.subname{i} '\images\sourcespace\r' R.R.subname{i} '_rep_' num2str(nr) '_surf_avg_DICS_coh_' ref_chan R.R.R.condname{cond}])
                                close all
                        end % reflist
                    end
                end
            end
            close all
        end
    end
end