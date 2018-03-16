function decide_max_coh_v2(R,suborthoplot)
% suborthoplot = 0;
for sub = 8:numel(R.subname)
    clear source_avg_dics peak_loc peak_ind contra_peak_mag ipsi_peak_mag contra_peak_loc ipsi_peak_loc source_cohmag source_cohloc
    for cond = 1:2
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
        for nr = 1:nrep
            for band = 1:numel(R.bandname)
%             load([R.datapathr R.subname{sub} '\ftdata\megdata_' num2str(nr) '_' R.condname{cond}],'megdata')
%             R.ref_list = {megdata.label{strmatch('STN',megdata.label)}};
%             clear megdata
% you need to set indices for left and right to do LM which is missing two
% channels left i.e. 4 now becomes 2
            for refN = 1:numel(R.ref_list)
                ref_chan = R.ref_list{refN};
                premotor_R = [-20 -6 84]./10;
                premotor_L = [20 -6 84]./10;
                load([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_DICSv2_peakinfo_source' R.condname{cond} 'nrep_' num2str(nr) '_' ref_chan])
                load([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_DICSv2_source' R.condname{cond} 'nrep_' num2str(nr) '_' ref_chanz],'source')
                source_avg_dics(:,refN,nr,cond) = source.avg.coh;
                figure(400)
                mask =  reshape(source.avg.coh,source.dim)>0;
                hold on
                fv = isosurface(mask,0);
                patch(fv,'FaceColor',[.1 .1 .1],'EdgeColor',[0 0 0],'FaceAlpha',0.1);
                
                dm = pdist2(source.pos,repmat(premotor_L,size(source.pos,1),1));
                [x i] =min(dm(:,2));
                [max_xyz(1) max_xyz(2) max_xyz(3)] = ind2sub(source.dim,i);
                PM_mask_L = makeSphereMask(source.dim,max_xyz,R.ROI.maskrho_dic);
                PM_mask_L = permute(PM_mask_L,[2 1 3]);
                fv = isosurface(PM_mask_L,0);
                patch(fv,'FaceColor',[.1 .1 .1],'EdgeColor',[1 0 0],'FaceAlpha',0.1);
                
                dm = pdist2(source.pos,repmat(premotor_R,size(source.pos,1),1));
                [x i] =min(dm(:,2));
                [max_xyz(1) max_xyz(2) max_xyz(3)] = ind2sub(source.dim,i);
                PM_mask_R = makeSphereMask(source.dim,max_xyz,R.ROI.maskrho_dic);
                PM_mask_R = permute(PM_mask_R,[2 1 3]);
                fv = isosurface(PM_mask_R,0);
                patch(fv,'FaceColor',[.1 .1 .1],'EdgeColor',[0 0 1],'FaceAlpha',0.1);
                
                peak_loc{refN,nr,cond} = DICS_peak{1};
                peak_ind(:,refN,nr,cond) = DICS_peak{3};
                % This recovers the peak locations and magnitudes from the
                % DICS coh structure. For ipsi, set all contra vertices to
                % zero, vice versa for contra.
                if refN<4 % RHS Ref
                    LR = source.avg.coh;
                    LR(source.pos(:,1)<=0) = 0; % set RHS to zero
                    LR = LR.*PM_mask_L(:)
                    [contra_peak_mag(refN,nr,cond) contra_loc] = max(LR); % DICS_peak{4};
                    
                    LR = source.avg.coh;
                    LR(source.pos(:,1)>=0) = 0; % set LHS to zero
                    LR = LR.*PM_mask_R(:)
                    [ipsi_peak_mag(refN,nr,cond) ipsi_loc] = max(LR); % DICS_peak{4};
                else % LHS Ref
                    LR = source.avg.coh;
                    LR(source.pos(:,1)>=0) = 0; % set LHS to zero
                    LR = LR.*PM_mask_R(:);
                    [contra_peak_mag(refN,nr,cond) contra_loc] = max(LR); % DICS_peak{4};
                    
                    LR = source.avg.coh;
                    LR(source.pos(:,1)<=0) = 0; % set RHS to zero
                    LR = LR.*PM_mask_L(:);
                    [ipsi_peak_mag(refN,nr,cond) ipsi_loc] = max(LR); % DICS_peak{4};
                end
                [subx suby subz] = ind2sub(source.dim,contra_loc);
                contra_peak_loc(:,refN,nr,cond) = source.pos(contra_loc,:);
                
                [subx suby subz] = ind2sub(source.dim,ipsi_loc);
                ipsi_peak_loc(:,refN,nr,cond) = source.pos(ipsi_loc,:);
                
                disp([cond nr refN])
            end
        end
    end
    
    % %     load([R.datapathr R.subname{i} '\ftdata\r' R.subname{i} 'rs'],'mri')
    % %     load([R.datapathr R.subname{i} '\ftdata\rs_transform_vox2ctf'])
    % %
    % %     T = transform_vox2ctf; %/transform_vox2spm;%;
    % %     Tmri = ft_transform_geometry(T,mri);
    % %     Tmri = ft_convert_units(Tmri,'cm');
    
    %     load([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_DICSv2_source' R.condname{cond} 'nrep_' num2str(nr) '_' ref_chan])
    
    Tmri = ft_read_mri([R.datapathr R.subname{sub} '\MRI\orig\r' R.subname{sub} '.nii'],'dataformat','nifti_spm');
    Tmri = ft_convert_units(Tmri,'cm');
    
    T = [1  0   0    0
        0   1   0    0
        0   0   1    2.2
        0   0   0    1];
    source = ft_transform_geometry(T,source);
    %                             cfg = [];
    %                             cfg.method = 'glassbrain';
    %                             cfg.funparameter = 'avg.coh';
    %                             ft_sourceplot(cfg,Tmri)
    %     close all
    GA_sourceDICS_L = source;
    LR = source_avg_dics;
    LR(source.pos(:,1)>=0,:,:,:) = 0; % set all left to zero
    %     LR = LR.*PM_mask_R(:);
    GA_sourceDICS_L.avg.coh = reshape(nanmean(nanmean(nanmean(LR(:,1:3,:,:),4),2),3),source.dim);
    if suborthoplot == 1; generic_source_plot(R.datapathr,R.subname,GA_sourceDICS_L,Tmri,'avg.coh','left_DICSv2','mni'); end
    close all
    GA_sourceDICS_R = source;
    LR = source_avg_dics;
    LR(source.pos(:,1)<=0,:,:,:) = 0; % set all right to zero
    %     LR = LR.*PM_mask_L(:);
    GA_sourceDICS_R.avg.coh = reshape(nanmean(nanmean(nanmean(LR(:,4:6,:,:),4),2),3),source.dim);
    if suborthoplot == 1; generic_source_plot(R.datapathr,R.subname,GA_sourceDICS_R,Tmri,'avg.coh','right_DICSv2','mni'); end
    
    figure
    load([R.datapathr R.subname{sub} '\ftdata\volsens' R.condname{cond}])
    vol = volsens.MEG.vol;
    vol = ft_convert_units(vol,'cm');
    ft_plot_mesh(vol.bnd);
    hold on
    a = gca;
    a.Children(1).FaceAlpha = 0.2;
    switch R.ipsicon
        case 'ipsi'
            xpeak_loc = ipsi_peak_loc;
            xpeak_mag = ipsi_peak_mag;
        case 'contra'
            xpeak_loc = contra_peak_loc;
            xpeak_mag = contra_peak_mag;
    end
    % These are the ON peaks
    locs = reshape(xpeak_loc(:,:,:,1),3,[]); locs(:,sum(locs==0)==3) = [];
    mags = reshape(xpeak_mag(:,:,1),1,[]); mags(mags==0) = [];
    scatter3(locs(1,:),locs(2,:),locs(3,:),mags*5000,'r','filled','d')
    
    hold on
    % These are the OFF peaks
    locs = reshape(xpeak_loc(:,:,:,2),3,[]); locs(:,sum(locs,1)==0) = [];
    mags = reshape(xpeak_mag(:,:,2),1,[]); mags(mags==0) = [];
    scatter3(locs(1,:),locs(2,:),locs(3,:),mags*5000,'b','filled')
    
    %% Run again to define global ROIs
    % This just takes the peak from the mean DICs image left and right
    % Define Left ROI
    A = GA_sourceDICS_L.avg.coh.*PM_mask_R;
    [a1 xyz] = max(A(:));
    max_pos = GA_sourceDICS_L.pos(xyz,:);
    [max_xyz(1) max_xyz(2) max_xyz(3)] = ind2sub(GA_sourceDICS_L.dim,xyz);
    mask = makeSphereMask(GA_sourceDICS_L.dim,max_xyz,R.ROI.maskrho_dic);
    maskL = permute(mask,[2 1 3]);
    %  ROI_list_L = find(GA_sourceDICS_L.inside.*reshape(mask,[],1));
    
    % Define Right ROI
    A = GA_sourceDICS_R.avg.coh.*PM_mask_L;
    [a1 xyz] = max(A(:));
    max_pos = GA_sourceDICS_R.pos(xyz,:);
    [max_xyz(1) max_xyz(2) max_xyz(3)] = ind2sub(GA_sourceDICS_R.dim,xyz);
    mask = makeSphereMask(GA_sourceDICS_R.dim,max_xyz,R.ROI.maskrho_dic);
    maskR = permute(mask,[2 1 3]);
    ROI_list_R = find(GA_sourceDICS_R.inside.*reshape(mask,[],1));
    
    %% Channel Level ROIs
    % Now search within repeats and conditions for maximums that fall
    % within the    bounded sphere defined by the global maximum
    for cond = 1:2
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
        for nr = 1:nrep
            for refN = 1:numel(R.ref_list)
                ref_chan = R.ref_list{refN};
                load([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_DICSv2_source' R.condname{cond} 'nrep_' num2str(nr) '_' ref_chan],'source')
                ref_chan = R.ref_list{refN};
                %% Plot brain mesh
                figure(100)
                mask =  reshape(source.avg.coh,source.dim)>0;
                hold on
                fv = isosurface(mask,0);
                patch(fv,'FaceColor',[.1 .1 .1],'EdgeColor',[0 0 0],'FaceAlpha',0.1);
                
                %% source_cohmag(LM1/RM1,refn,nr,cond)
                cohvol = reshape(source.avg.coh,source.dim); cohvol = cohvol.*maskL;
                cohvol = cohvol(:); cohvol(source.pos(:,1)>=0) = NaN; %reshape(cohvol,source.dim)
                [a1 xyz] = max(cohvol(:));
                [max_xyz(1) max_xyz(2) max_xyz(3)] = ind2sub(source.dim,xyz);
                source_cohmag(1,refN,nr,cond) = a1;
                source_cohloc(1,refN,nr,cond) = xyz;
                % Plot Left mask
                mask =  reshape(cohvol,source.dim);
                hold on
                fv = isosurface(mask,0);
                patch(fv,'FaceColor',[0 0 .7],'EdgeColor',[0 0 1],'FaceAlpha',0.4);
                
                cohvol = reshape(source.avg.coh,source.dim); cohvol = cohvol.*maskR;
                cohvol = cohvol(:); cohvol(source.pos(:,1)<=0) = NaN;% reshape(cohvol,source.dim)
                [a1 xyz] = max(cohvol(:));
                [max_xyz(1) max_xyz(2) max_xyz(3)] = ind2sub(source.dim,xyz);
                source_cohmag(2,refN,nr,cond) = a1;
                source_cohloc(2,refN,nr,cond) = xyz;
                % Plot Right mask
                mask =  reshape(cohvol,source.dim);
                hold on
                fv = isosurface(mask,0);
                patch(fv,'FaceColor',[0.7 0 0],'EdgeColor',[1 0 0],'FaceAlpha',0.4);
                view(45,45);
                axis equal;
                drawnow; shg
            end
        end
    end
    for cond = 1:2
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
        for nr = 1:nrep
            %% STNselection{Right/Left,Ipsi/Contra,nr,cond}
            %% = {Cond, set, cohmax, ref_ind,refname,sourceind,sourcesub,sourcepos
            % Contra -Right STN, Left M1
            maxcohs = squeeze(source_cohmag(1,1:3,nr,cond));
            [a1 xyz] = max(maxcohs);
            maxlocs = squeeze(source_cohloc(1,1:3,nr,cond));
            [loc(1) loc(2) loc(3)] = ind2sub(source.dim,maxlocs(xyz));
            STNselection{1,2,nr,cond} = {R.condname{cond} 'contra_RSTN_LM1' a1 xyz R.ref_list{xyz} maxlocs(xyz)  loc source.pos(maxlocs(xyz),:)};
            
            % Ipsi -Left STN, Left M1
            maxcohs = squeeze(source_cohmag(1,4:6,nr,cond));
            [a1 xyz] = max(maxcohs);
            cohmax(1,nr,sub,cond) = a1;
            locmax(:,1,nr,sub,cond) = maxlocs(xyz);
            maxlocs = squeeze(source_cohloc(1,4:6,nr,cond));
            [loc(1) loc(2) loc(3)] = ind2sub(source.dim,maxlocs(xyz));
            STNselection{2,1,nr,cond} = {R.condname{cond} 'ipsi_LSTN_LM1' a1 xyz R.ref_list{xyz} maxlocs(xyz)  loc source.pos(maxlocs(xyz),:)};
            
            % Contra -Left STN, Right M1
            maxcohs = squeeze(source_cohmag(2,4:6,nr,cond));
            [a1 xyz] = max(maxcohs);
            maxlocs = squeeze(source_cohloc(2,4:6,nr,cond));
            [loc(1) loc(2) loc(3)] = ind2sub(source.dim,maxlocs(xyz));
            STNselection{2,2,nr,cond} = {R.condname{cond} 'contra_LSTN_RM1' a1 xyz R.ref_list{xyz} maxlocs(xyz)  loc source.pos(maxlocs(xyz),:)};
            
            % Ipsi -Right STN, Right M1
            maxcohs = squeeze(source_cohmag(2,1:3,nr,cond));
            [a1 xyz] = max(maxcohs);
            cohmax(2,nr,sub,cond) = a1;
            locmax(:,1,nr,sub,cond) = maxlocs(xyz);
            maxlocs = squeeze(source_cohloc(2,1:3,nr,cond));
            [loc(1) loc(2) loc(3)] = ind2sub(source.dim,maxlocs(xyz));
            STNselection{1,1,nr,cond} = {R.condname{cond} 'ipsi_RSTN_RM1' a1 xyz R.ref_list{xyz} maxlocs(xyz)  loc source.pos(maxlocs(xyz),:)};
        end
    end
    save([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_DICSv2_source_channelselectioninfo'],'STNselection')
    
    
    % Subject Effects
    % source_cohmag(left/rightM1,refn,nr,cond)
    % IPSIlateral sources
    OFF = [squeeze(source_cohmag(1,4:6,1,2)) squeeze(source_cohmag(2,1:3,1,2))]; OFF = OFF(:); OFF(OFF==0) = [];
    ON = [squeeze(source_cohmag(1,4:6,:,1)) squeeze(source_cohmag(2,1:3,:,1))]; ON = ON(:); ON(ON==0) = [];
    
    figure
    [centre height] = aboxplotTW([ON,OFF],'Labels',{'ON';'OFF'},'Colorgrad','blue_down');
    ylabel('DICs Max. Coherence','FontSize',14);
    title('STN/DICs Ipsilateral High Beta Source Coherence','FontSize',14)
    [H,P,CI,STATS] = ttest(ON,OFF);
    dum = [ON, OFF]; ylim([0 1.2*max(dum(:))])
    annotation(gcf,'textbox',...
        [0.442619519094767 0.786158401184308 0.149892857142857 0.0880952380952401],...
        'String',{sprintf('P = %.3f',P)},...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',12,...
        'FitBoxToText','off');
    dics_cohstattable.ispi.cohstat = [mean(ON) std(ON)/sqrt(numel(ON)) mean(OFF) std(OFF)/sqrt(numel(OFF)) mean(ON)-mean(OFF) STATS.tstat P];
    dics_cohstattable.ispi.cohmax = [max(ON) max(OFF)];
    % CONTRAlateral sources
    OFF = [squeeze(source_cohmag(1,1:3,:,2)) squeeze(source_cohmag(2,4:6,:,2))]; OFF = OFF(:); OFF(OFF==0) = [];
    ON = [squeeze(source_cohmag(1,1:3,:,1)) squeeze(source_cohmag(2,4:6,:,1))]; ON = ON(:); ON(ON==0) = [];
    figure
    [centre height] = aboxplotTW([ON,OFF],'Labels',{'ON';'OFF'},'Colorgrad','blue_down');
    ylabel('DICs Max. Coherence','FontSize',14);
    title('STN/DICs Contralateral High Beta Source Coherence','FontSize',14)
    [H P] = ttest(ON,OFF);
    dum = [ON, OFF]; ylim([0 1.2*max(dum(:))])
    annotation(gcf,'textbox',...
        [0.442619519094767 0.786158401184308 0.149892857142857 0.0880952380952401],...
        'String',{sprintf('P = %.3f',P)},...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',12,...
        'FitBoxToText','off');
    dics_cohstattable.contra.cohstat = [mean(ON) std(ON)/sqrt(numel(ON)) mean(OFF) std(OFF)/sqrt(numel(OFF)) mean(ON)-mean(OFF) STATS.tstat P];
    dics_cohstattable.contra.cohmax = [max(ON) max(OFF)];
    mkdir([R.datapathr R.subname{sub} '\stats\'])
    save([R.datapathr R.subname{sub} '\stats\DICS_ONvOFF'],'dics_cohstattable')
    savefigure_v2([R.datapathr R.subname{sub} '\images\sourcespace\'],['r' R.subname{sub} '_DICsv2_grandaverage_ON_OFF'],[],[],[]); close all
end
barplot_coh_groups_nocell(R.datapathr,cohmax)
savefigure_v2([R.datapathr 'results\images\'],['groupDICs_barplots'],[],[],[]); %close all


