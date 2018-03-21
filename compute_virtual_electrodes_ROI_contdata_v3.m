function compute_virtual_electrodes_ROI_contdata_v3(R)
%%%
% This function computes virtual electrodes at the list of voxels given by
% decide_coh script. We compute all and pass to vchansave structure for
% selection in subsequent scripts.

%%%
normv = @(x) (x-mean(x)); %./std(x);
scatleg = {'bo','bx';'ro','rx'};
for band = 1 %:numel(R.bandname)
for sub = 1:numel(R.subname)
    for cond = 1:2
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
        [~,~,nrepOFF,~] = data_fileguide(R.subname{sub},1);
        for nr = 1:nrep
            if nrepOFF<nrep
                nroff = 1;
            else
                nroff = nr;
            end
%                         eval(['!del /q ' R.datapathr R.subname{sub} '\ftdata\ROI_analy\'])
%                         delete([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon]);
            load([R.datapathr R.subname{sub} '\ftdata\meg_clean_fullsampdata_cont_' num2str(nr) '_' R.condname{cond}])
            %             load([R.datapathr R.subname{i} '\ftdata\meg_clean_data_cont_' num2str(nr) '_' R.condname{cond}])
            load([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_DICSv2_source_channelselectioninfo_' R.bandname{band}])
            load([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_LCMV_source_' R.condname{cond} 'nrep_' num2str(nr)])
            for side = 1:2
                switch R.ipsicon
                    case 'ipsi'
                        seldet = STNselection{side,1,nrepOFF,2}; % selected for ipsi % Use OFF
                    case 'contra'
                        seldet = STNselection{side,2,nrepOFF,2}; % selected for contra % Use OFF
                end
                chansel = ft_channelselection('MEG', fullsampdata.label); % find MEG sensor names
                chansel = match_str(fullsampdata.label, chansel);         % find MEG sensor indices
                
                refsel = match_str(fullsampdata.label,seldet{5});         % find STN ref indice
                
                cnt = 0;
                %                 figure(200)
                %                 mask =  reshape(source.avg.pow,source.dim)>0;
                %                 hold on
                %                 fv = isosurface(mask,0);
                %                 patch(fv,'FaceColor',[.1 .1 .1],'EdgeColor',[0 0 0],'FaceAlpha',0.1); hold on
                
                ROI = seldet{7};
                ROIsave(:,cond,nr,side,sub) = ROI;
                %                 scatter3(ROI(1),ROI(2),ROI(3),125,scatleg{cond,side},'filled'); hold on
                mask = makeSphereMask(source.dim,ROI,R.ROI.maskrho_vc,0);
                mask = permute(mask,[2 1 3]);
                ROI_list = find(source.inside.*mask(:));
                POS_list = source.pos(source.inside.*mask(:)>0,:);
                clear vchansave
                
                for p = 1:size(ROI_list,1)
                    iz = ROI_list(p); size(ROI_list,1)
                    loc = zeros(1,3);
                    [loc(1),loc(2),loc(3)] = ind2sub(source.dim,ROI_list(p));
                    %                     scatter3(x,y,z,'rx'); hold on
                    
                    iz = ROI_list(p); %sub2ind(source.dim,ROI_list(p,1),ROI_list(p,2),ROI_list(p,3));
                    %             MNI = source.pos(iz,:);
                    %             fprintf('You have chosen a point at %1.f %1.f %1.f in MNI',MNI)
                    
                    vchan = [];
                    vchan.label = {'gam_pow_x', 'gam_pow_y', 'gam_pow_z'};
                    vchan.time = fullsampdata.time;
                    sfilter = source.avg.filter{iz};
                    mtrial = fullsampdata.trial;
                    winsz = fullsampdata.fsample*2;
                    mx = mtrial{1}; N = floor(size(mx,2)./winsz);
                    mx = mx(:,1:N*winsz);
                    mx = squeeze(num2cell(reshape(mx,size(mx,1),[],N),[1 2]));
                    trial = cell(1,numel(mx));
                    vchan.time = fullsampdata.time{1}(1:winsz*N);
                    %                     parpool('local',4);
                    N = numel(mx);
%                     progressStepSize = 1;
%                     ppm = ParforProgMon('Example: ', N, progressStepSize, 600, 160);
                    for trln = 1:numel(mx)
                        if ~isempty(sfilter)
                            trial{trln} = sfilter * mx{trln}(chansel,:);
                        else
                            trial{trln} = NaN;
                        end
%                         ppm.increment();
                    end
%                     ppm.delete()
                    vchan.trial = trial;
                    x = cat(2,vchan.trial{:});
                    x(isnan(x)) = [];
                    [u,s,v] = svd(x,'econ');
%                     ppm = ParforProgMon('Example: ', N, progressStepSize, 600, 160);
                    for trln = 1:numel(mx)
                        if ~isempty(sfilter)
                            A = normv(u(:,1)' * sfilter * mx{trln}(chansel,:));
                            B = normv(mx{trln}(refsel,:));
                            trial{trln} = [A; B];
                        else
                            trial{trln} = [NaN; NaN];
                        end
%                         ppm.increment();
                    end
%                     ppm.delete()
                    vchan.trial = {[trial{:}]};
                    vchan.label = {'src','ref'};
                    vchan.loc = loc;
                    vchan.pos = POS_list(p,:);
                    vchansave(p) = vchan;
                    %             source_sens{p} = vchan.trial;
                    disp(p)
                    %             ppm.increment();
                    %                             cfg = [];
                    %                             cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
                    %                             ft_databrowser(cfg, vchan);
                    %                             source_sens.label{p} = sprintf('%f ',source.pos(p,:));
                    %                             x = cat(2,vchan.trial{:});
                    %
                    %         source_sens.trial = cat(source_sens.trial,vchan.trial);
                    %         source_sens.time  = cat(source_sens.time,vchan.time);
                    vtime = vchan.time;
                end
                
                save([R.datapathr R.subname{sub} '\ftdata\virtual_sources_' num2str(nr) '_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bandname{band}],'vchansave','vtime','ROI_list')
                disp([side nr cond sub])
            end
        end
    end
end
end
a = 1;
