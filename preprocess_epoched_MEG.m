function preprocess_epoched_MEG(R)
for sub = 1:numel(R.subnames)
    for cond = 1:2
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subnames{sub},cond-1);
        for nr = 1:nrep
            %         datafilen = [pp_mark datafileN{nr}];
            %             % Epoch those fuckers
            %             cfg = [];
            %             cfg.datafile = ['C:\data\TimExtracts190516\' subname{sub} '\' subname{sub} '_R_' num2str(nr) '_' condnamelc{cond} '.mat'];
            %             megdata = ft_preprocessing(cfg);
            %             ref_list = find(strncmp(megdata.label,'STN_',4));
            %             save([datapathr subname{sub} '\ftdata\STN_ref_list' condname{cond}],'ref_list')
            %
            %             %        for p = 1:size(megdata.trial{1},1)
            %             %            t = megdata.trial{1}(p,:);
            %             %            if any(isnan(t))
            %             %                ai(p) = 1;
            %             %            end
            %             %        end
            %             cfg = [];
            %             cfg.resamplefs = 512;
            %             megdata = ft_resampledata(cfg,megdata);#
            load([R.datapathr R.subname{sub} '\ftdata\meg_clean_data_cont_' num2str(nr) '_' R.condname{cond}],'data_clean')
            %             load([datapathr subname{sub} '\ftdata\meg_clean_data_' num2str(nr) '_' condname{cond}])
            megdata = data_clean;
            
            cfg = [];
            cfg.length = 1.5;
            megdata = ft_redefinetrial(cfg,megdata);
            
            % Smooth dem badboiz
            %             cfg = [];
            %             cfg.lpfilter = 'yes';
            %             cfg.hpfilter = 'yes';
            %             %         cfg.channel         = {'MEG','eeg'};
            %             cfg.lpfreq = 98;
            %             cfg.hpfreq = 4;
            %             megdata = ft_preprocessing(cfg,megdata);
            
            %% JUMP
            % cut the crap
            cfg = [];
            cfg.artfctdef.zvalue.channel    = 'MEG';
            cfg.artfctdef.zvalue.cutoff     = 25;
            cfg.artfctdef.zvalue.trlpadding = 0;
            cfg.artfctdef.zvalue.artpadding = 0;
            cfg.artfctdef.zvalue.fltpadding = 0;
            
            % algorithmic parameters
            cfg.artfctdef.zvalue.cumulative    = 'yes';
            cfg.artfctdef.zvalue.medianfilter  = 'yes';
            cfg.artfctdef.zvalue.medianfiltord = 9;
            cfg.artfctdef.zvalue.absdiff       = 'yes';
            
            % make the process interactive
            cfg.artfctdef.zvalue.interactive = 'yes';
            [cfg, artifact_jump] = ft_artifact_zvalue(cfg,megdata);
            
            %% MUSCLE
            cfg            = [];
            % channel selection, cutoff and padding
            cfg.artfctdef.zvalue.channel = 'MRT*';
            cfg.artfctdef.zvalue.cutoff      = 4;
            cfg.artfctdef.zvalue.trlpadding  = 0;
            cfg.artfctdef.zvalue.fltpadding  = 0;
            cfg.artfctdef.zvalue.artpadding  = 0.1;
            
            % algorithmic parameters
            cfg.artfctdef.zvalue.bpfilter    = 'yes';
            cfg.artfctdef.zvalue.bpfreq      = [110 120];
            cfg.artfctdef.zvalue.bpfiltord   = 5;
            cfg.artfctdef.zvalue.bpfilttype  = 'but';
            cfg.artfctdef.zvalue.hilbert     = 'yes';
            cfg.artfctdef.zvalue.boxcar      = 0.2;
            
            % make the process interactive
            cfg.artfctdef.zvalue.interactive = 'yes';
            
            [cfg, artifact_muscle] = ft_artifact_zvalue(cfg,megdata);
            
            if any(strncmp(megdata.label,'HEOG',3))
                disp('Using EOG data to detect artefacts')
                % EOG
                cfg            = [];
                % channel selection, cutoff and padding
                cfg.artfctdef.zvalue.channel     = 'EOG';
                cfg.artfctdef.zvalue.cutoff      = 10;
                cfg.artfctdef.zvalue.trlpadding  = 0;
                cfg.artfctdef.zvalue.artpadding  = 0.1;
                cfg.artfctdef.zvalue.fltpadding  = 0;
                
                % algorithmic parameters
                cfg.artfctdef.zvalue.bpfilter   = 'yes';
                cfg.artfctdef.zvalue.bpfilttype = 'but';
                cfg.artfctdef.zvalue.bpfreq     = [1 15];
                cfg.artfctdef.zvalue.bpfiltord  = 4;
                cfg.artfctdef.zvalue.hilbert    = 'yes';
                
                % feedback
                cfg.artfctdef.zvalue.interactive = 'yes';
                
                [cfg, artifact_EOG] = ft_artifact_zvalue(cfg,megdata);
                
                cfg=[];
                cfg.artfctdef.reject = 'complete';
                cfg.artfctdef.eog.artifact = artifact_EOG; %
                cfg.artfctdef.jump.artifact = artifact_jump;
                cfg.artfctdef.muscle.artifact = artifact_muscle;
                meg_clean_data = ft_rejectartifact(cfg,megdata);
            else
                disp('No EOG data found, no artefact rejection for eye movement')
                cfg=[];
                cfg.artfctdef.reject = 'complete';
                cfg.artfctdef.jump.artifact = artifact_jump;
                cfg.artfctdef.muscle.artifact = artifact_muscle;
                meg_clean_data = ft_rejectartifact(cfg,megdata);
            end
            save([R.datapathr R.subname{sub} '\ftdata\meg_clean_data_' num2str(nr) '_' R.condname{cond}],'meg_clean_data')
        end
    end
end