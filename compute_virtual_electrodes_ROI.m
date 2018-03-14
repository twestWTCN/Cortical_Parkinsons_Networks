function compute_virtual_electrodes_ROI(datapathr,subname,ref_chan)
% subname = {'LM'};
% datapathr = 'C:\Users\Tim\Documents\Work\Cortical_Networks\Data\';
% ref_chan = 'STN_R23';

normv = @(x) (x-mean(x))./std(x);

for i = 1:numel(subname)
    [datafileN,pp_mark,nrep,senscheck] = data_fileguide(subname{i},1);
    for nr = 1:nrep
%         cfg = [];
%         cfg.datafile = ['C:\data\TimExtracts190516\' subname{i} '\' subname{i} '_R_' num2str(nr) '_off.mat'];
%         cfg.lpfilter = 'yes';
%         cfg.hpfilter = 'yes';
%         cfg.channel         = {'MEG',ref_chan};
%         cfg.lpfreq = 40;
%         cfg.hpfreq = 4;
%         megdata = ft_preprocessing(cfg);
%         
%         cfg = [];
%         cfg.resamplefs = 200;
%         megdata = ft_resampledata(cfg,megdata);
%         
%         cfg = [];
%         cfg.length = 1.5;
%         cfg.overlap = 0.25;
%         megdata = ft_redefinetrial(cfg,megdata);
%         save([datapathr subname{i} '\ftdata\r' subname{i} '_meg_cont_pp'],'megdata')
        load([datapathr subname{i} '\ftdata\r' subname{i} '_meg_cont_pp'],'megdata')
        load([datapathr subname{i} '\ftdata\r' subname{i} '_LCMV_source'])
        source = source;
        chansel = ft_channelselection('MEG', megdata.label); % find MEG sensor names
        chansel = match_str(megdata.label, chansel);         % find MEG sensor indices
        
        refsel = match_str(megdata.label,ref_chan);         % find MEG sensor indices
        
        %         source_sens = [];
        %         source_sens.trial = {};
        %         source_sens.trial = [];
        %
        cnt = 0;
           ROI = [4 7 10];
           mask = makeSphereMask(source.dim,ROI,2);
           ROI_list = find(source.inside.*reshape(mask,[],1));
           
        for p = 1:size(ROI_list,1)
            iz = ROI_list(p); %sub2ind(source.dim,ROI_list(p,1),ROI_list(p,2),ROI_list(p,3));
%             MNI = source.pos(iz,:);
%             fprintf('You have chosen a point at %1.f %1.f %1.f in MNI',MNI)
            
            vchan = [];
            vchan.label = {'gam_pow_x', 'gam_pow_y', 'gam_pow_z'};
            vchan.time = megdata.time;
            sfilter = source.avg.filter{iz};
            mtrial = megdata.trial;
            trial = cell(1,numel(mtrial));
            for trln = 1:numel(mtrial)
                if ~isempty(sfilter)
                    trial{trln} = sfilter * mtrial{trln}(chansel,:);
                else
                    trial{trln} = NaN;
                end
                %                   ppm.increment();
            end
            vchan.trial = trial;
            x = cat(2,vchan.trial{:});
            x(isnan(x)) = [];
            [u,s,v] = svd(x,'econ');
            
            for trln = 1:numel(mtrial)
                if ~isempty(sfilter)
                    A = normv(u(:,1)' * sfilter * mtrial{trln}(chansel,:));
                    B = normv(mtrial{trln}(refsel,:));
                    vchan.trial{trln} = [A; B];
                else
                    vchan.trial{trln} = [NaN; NaN];
                end
            end
            vchan.label = {'src'};
            source_sens{p} = vchan.trial;
            disp(p)
%             ppm.increment();
            %         cfg = [];
            %         cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
            %         ft_databrowser(cfg, vchan);
            %         source_sens.label{p} = sprintf('%f ',source.pos(p,:));
            %         x = cat(2,vchan.trial{:});
            %
            %         source_sens.trial = cat(source_sens.trial,vchan.trial);
            %         source_sens.time  = cat(source_sens.time,vchan.time);
        end
        vtime = vchan.time;
        save([datapathr subname{i} '\ftdata\virtual_sources_' num2str(nr) '_ROI'],'source_sens','vtime','ROI_list')
    end
end