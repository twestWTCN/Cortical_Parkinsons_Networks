function preprocess_cont_MEG_v6(R)
for sub = 1:numel(R.subname) %problems with 7(JB 8(KB
    for cond = 1:2 
            load(
            
            ref_list = find(strncmp(fulldata.label,'STN_',4));
            save([R.datapathr R.subname{sub} '\ftdata\STN_ref_list' R.condname{cond}],'ref_list')

            % decide where to truncate
            if exist([R.datapathr R.subname{sub} '\ftdata\meg_clean_DEDGE_' num2str(nr) '_' R.condname{cond} '.mat'])==0
                figure
                plot(repmat(fulldata.time{1},32,1)',fulldata.trial{1}(1:32,:)'); shg
                figure(100); close 100
                options = [];
                options.WindowStyle = 'normal';
                uinputer = inputdlg('Enter signal edges:',...
                    'Truncate',1,{''},options);
                dedge = str2num(uinputer{:});
                save([R.datapathr R.subname{sub} '\ftdata\meg_clean_DEDGE_' num2str(nr) '_' R.condname{cond}],'dedge'); close all
            end
            load([R.datapathr R.subname{sub} '\ftdata\meg_clean_DEDGE_' num2str(nr) '_' R.condname{cond}],'dedge')
            close all
            
            % Truncate
            cfg = [];
            cfg.latency = dedge;
            fulldata = ft_selectdata(cfg,fulldata);
            Z = fulldata.trial{1};
            % Fill in NaNs
            %             ppm = ParforProgMon('Example: ', size(fulldata.trial{1},1), 1, 300, 80);
            for j = 1:size(fulldata.trial{1},1)
                %                 close all
                X = Z(j,:);
                X = X-nanmean(X);
                %                 figure
                [X flag(j)] = removeNaNs(X,0); %title(fulldata.label{j});
                Z(j,:) = X;
                disp(sprintf(' %1.f percent complete',(j./size(fulldata.trial{1},1))*100 ))
                %                  ppm.increment()
            end
            %             ppm.delete()
            iremv = find(flag);
            Z(iremv,:) = [];
            for i = 1:length(iremv)
                fulldata.label{i} = [];
            end
            fulldata.trial{1} = Z;
            
            %%%%% Full data for granger analyses
            cfg = [];
            cfg.resamplefs = R.pp.cont.full.fs;
            fullsampdata = ft_resampledata(cfg,fulldata);
            
            % Preprocess
            cfg = [];
            cfg.dftfilter= 'yes';
            cfg.dftfreq = [50 100 150];
            cfg.lpfilter = 'yes';
            cfg.hpfilter = 'yes';
            cfg.lpfreq = R.pp.cont.full.bp(2);
            cfg.hpfreq = R.pp.cont.full.bp(1);
            cfg.demean = 'yes';
            %             cfg.polyremoval   = 'yes';
            %             cfg.polyorder = 4;
            fullsampdata = ft_preprocessing(cfg,fullsampdata);
            
            % Truncate to remove filter artefact
            cfg = [];
            cfg.latency = [5 fullsampdata.time{1}(end)-2.5];
            fullsampdata = ft_selectdata(cfg,fullsampdata);
            
            
            % Remove Jumps
            for j = 1:size(fullsampdata.trial{1},1)
                close all
                X = fullsampdata.trial{1}(j,:);
                figure
                X = removejumps(X,8,1); title(fullsampdata.label{j});
                fullsampdata.trial{1}(j,:) = X;
            end
            close all
            save([R.datapathr R.subname{sub} '\ftdata\meg_clean_fullsampdata_cont_' num2str(nr) '_' R.condname{cond}],'fullsampdata')
            
            %%%%
            % Now resample for standard data use
            % resample
            cfg = [];
            cfg.resamplefs = R.pp.cont.thin.fs;
            fulldata = ft_resampledata(cfg,fulldata);
            
            
            % Preprocess
            cfg = [];
            cfg.dftfilter= 'yes';
            cfg.dftfreq = [50];
            cfg.lpfilter = 'yes';
            cfg.hpfilter = 'yes';
            cfg.lpfreq = R.pp.cont.thin.bp(2);
            cfg.hpfreq = R.pp.cont.thin.bp(1);
            %             cfg.demean = 'yes';
            cfg.polyremoval   = 'yes';
            cfg.polyorder = 4;
            fulldata = ft_preprocessing(cfg,fulldata);
            
            % Truncate to remove filter artefact
            cfg = [];
            cfg.latency = [5 fulldata.time{1}(end)-2.5];
            fulldata = ft_selectdata(cfg,fulldata);
            
            
            % Remove Jumps
            for j = 1:size(fulldata.trial{1},1)
                close all
                X = fulldata.trial{1}(j,:);
                figure
                X = removejumps(X,8,1); title(fulldata.label{j});
                fulldata.trial{1}(j,:) = X;
            end
            
            
            
            %             cfg = [];
            %             cfg.channel    = {'meg'};
            %             megdata = ft_selectdata(cfg,fulldata);
            
            %             cfg            = [];
            %             cfg.method     = 'runica';
            %             comp           = ft_componentanalysis(cfg, megdata);
            %
            %             cfg           = [];
            %             cfg.component = [1:24];       % specify the component(s) that should be plotted
            %             cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
            %             cfg.comment   = 'no';
            %             ft_topoplotIC(cfg, comp)
            %
            %             figure
            %             cfg          = [];
            %             cfg.channel  = [1:24]; % components to be plotted
            %             cfg.viewmode = 'component';
            %             ft_databrowser(cfg, comp)
            %
            %             figure(100)
            %             options = [];
            %             options.WindowStyle = 'normal';
            %             uinputer = inputdlg('Enter components to remove:',...
            %                 'ICA Clean',1,{''},options);
            %             cmpntrm = str2num(uinputer{:});
            %
            %             % the original data can now be reconstructed, excluding those components
            %             cfg           = [];
            %             cfg.component = cmpntrm;
            %             data_clean    = ft_rejectcomponent(cfg, comp,fulldata);
            data_clean = fulldata;
            
            
            %             cfg          = [];
            %             cfg.channel  =  1:48;
            %             cfg.blocksize = 150;
            %             cfg.viewmode = 'vertical';
            %             ft_databrowser(cfg, megdata)
            %             close all
            save([R.datapathr R.subname{sub} '\ftdata\meg_clean_data_cont_' num2str(nr) '_' R.condname{cond}],'data_clean')
            disp([nr cond sub])
        end
    end
end