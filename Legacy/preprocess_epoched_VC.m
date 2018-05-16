function preprocess_epoched_VC(R)
for sub = 1:numel(R.subname)
    for cond = 1:2
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
        for nr = 1:nrep
            for band = [1 3]
                for side = 1:2
                    vchansave = []
                    load([R.datapathr R.subname{sub} '\ftdata\virtual_sources_' num2str(nr) '_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bandname{band}])
                    %             load([datapathr subname{sub} '\ftdata\meg_clean_data_' num2str(nr) '_' condname{cond}])
                    clear concattrial
                    for  vc = 1:size(vchansave,2)
                        concattrial(:,vc) = vchansave(vc).trial{1}(1,:);
                    end
                    %                     [U,S,V] = svd(concattrial,'econ')
                    vchansave = vchansave(1);
                    vchansave.trial{1}(1,:) = mean(concattrial,2);
                    
                    for vc = 1:size(vchansave,2)
                        megdata = vchansave(vc);
                        megdata = rmfield(megdata,'time');
                        megdata.time{1} = vchansave(vc).time;
                        cfg = [];
                        cfg.length = 1;
                        megdata = ft_redefinetrial(cfg,megdata);
                        
                        %% JUMP
                        % cut the crap
                        cfg = [];
                        cfg.artfctdef.zvalue.channel    = 'all';
                        cfg.artfctdef.zvalue.cutoff     = 25;
                        cfg.artfctdef.zvalue.trlpadding = 0;
                        cfg.artfctdef.zvalue.artpadding = 0;
                        cfg.artfctdef.zvalue.fltpadding = 0;
                        
                        % algorithmic parameters
                        cfg.artfctdef.zvalue.cumulative    = 'yes';
                        cfg.artfctdef.zvalue.medianfilter  = 'yes';
                        cfg.artfctdef.zvalue.medianfiltord = 9;
                        %                     cfg.artfctdef.zvalue.absdiff       = 'yes';
                        
                        % make the process interactive
                        %                     cfg.artfctdef.zvalue.interactive = 'yes';
                        [cfg, artifact_jump] = ft_artifact_zvalue(cfg,megdata);
                        
                        %% MUSCLE
                        cfg            = [];
                        % channel selection, cutoff and padding
                        cfg.artfctdef.zvalue.channel = 'all';
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
                        %                     cfg.artfctdef.zvalue.interactive = 'yes';
                        
                        [cfg, artifact_muscle] = ft_artifact_zvalue(cfg,megdata);
                        
                        cfg=[];
                        cfg.artfctdef.reject = 'complete';
                        cfg.artfctdef.jump.artifact = artifact_jump;
                        vc_clean(vc) = ft_rejectartifact(cfg,megdata);
                    end
                    save([R.datapathr R.subname{sub} '\ftdata\virtual_sources_clean_' num2str(nr) '_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bandname{band}],'vc_clean')
                    clear vc_clean
                end
            end
        end
    end
end
end