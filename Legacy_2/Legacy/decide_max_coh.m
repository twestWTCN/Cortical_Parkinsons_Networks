function decide_max_coh(datapathr,subname)
ref_list = {'STN_R01','STN_R12','STN_R23','STN_L01','STN_L12','STN_L23'};
suborthoplot = 1;
condname = {'ON','OFF'}; ipsicon = 'ipsi';
for i = 1:numel(subname)
    for cond = 1:2
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(subname{i},cond-1);
        for nr = 1:nrep
            for refN = 1:numel(ref_list)
                ref_chan = ref_list{refN};
                load([datapathr subname{i} '\ftdata\r' subname{i} '_DICS_peakinfo_source' condname{cond} 'nrep_' num2str(nr) '_' ref_chan])
                load([datapathr subname{i} '\ftdata\r' subname{i} '_DICS_source' condname{cond} 'nrep_' num2str(nr) '_' ref_chan],'source')
                source_avg_dics(:,refN,nr,cond) = source.avg.coh;
                
                peak_loc{refN,nr,cond} = DICS_peak{1};
                peak_ind(:,refN,nr,cond) = DICS_peak{3};
                % This recovers the peak locations and magnitudes from the
                % DICS coh structure. For ipsi, set all contra vertices to
                % zero, vice versa for contra. 
                if refN<4 % RHS Ref
                    LR = source.avg.coh;
                    LR(source.pos(:,2)>=0) = 0; % set RHS to zero
                    [contra_peak_mag(refN,nr,cond) contra_loc] = max(LR); % DICS_peak{4};
                    
                    LR = source.avg.coh;
                    LR(source.pos(:,2)<=0) = 0; % set LHS to zero
                    [ipsi_peak_mag(refN,nr,cond) ipsi_loc] = max(LR); % DICS_peak{4};
                else % LHS Ref
                    LR = source.avg.coh;
                    LR(source.pos(:,2)<=0) = 0; % set LHS to zero
                    [contra_peak_mag(refN,nr,cond) contra_loc] = max(LR); % DICS_peak{4};
                    
                    LR = source.avg.coh;
                    LR(source.pos(:,2)>=0) = 0; % set RHS to zero
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
    
    load([datapathr subname{i} '\ftdata\r' subname{i} 'rs'],'mri')
    load([datapathr subname{i} '\ftdata\rs_transform_vox2ctf'])
    
    T = transform_vox2ctf; %/transform_vox2spm;%;
    Tmri = ft_transform_geometry(T,mri);
    Tmri = ft_convert_units(Tmri,'cm');
    
    close all
    GA_sourceDICS_L = source;
    LR = source_avg_dics;
    LR(source.pos(:,2)>=0,:,:,:) = 0; % set all left to zero
    GA_sourceDICS_L.avg.coh = reshape(nanmean(nanmean(nanmean(LR(:,1:3,:,:),4),2),3),source.dim);
    if suborthoplot == 1; generic_source_plot(datapathr,subname,GA_sourceDICS_L,Tmri,'avg.coh','left'); end
    close all
    GA_sourceDICS_R = source;
    LR = source_avg_dics;
    LR(source.pos(:,2)<=0,:,:,:) = 0; % set all right to zero
    GA_sourceDICS_R.avg.coh = reshape(nanmean(nanmean(nanmean(LR(:,4:6,:,:),4),2),3),source.dim);
    if suborthoplot == 1; generic_source_plot(datapathr,subname,GA_sourceDICS_R,Tmri,'avg.coh','right'); end
    
    figure
    load([datapathr subname{i} '\ftdata\coreg_headmodel'],'vol')
    cfg = [];
    ft_plot_mesh(vol.bnd);
    hold on
    a = gca;
    a.Children(1).FaceAlpha = 0.2
    switch ipsicon
        case 'ipsi'
            xpeak_loc = ipsi_peak_loc;
            xpeak_mag = ipsi_peak_mag;
        case 'contra'
            xpeak_loc = contra_peak_loc;
            xpeak_mag = contra_peak_mag;
    end
    % These are the ON peaks         
    locs = reshape(xpeak_loc(:,:,:,1),3,[]); locs(:,sum(locs==0)==3) = [];
    mags = reshape(xpeak_mag(:,:,1),1,[]); mags(sum(mags==0)==3) = [];
    scatter3(locs(1,:),locs(2,:),locs(3,:),mags*5000,'r','filled','d')
    
    hold on
    % These are the OFF peaks 
    locs = reshape(xpeak_loc(:,:,:,2),3,[]); locs(:,sum(locs==0)==3) = [];
    mags = reshape(xpeak_mag(:,:,2),1,[]); mags(sum(mags==0)==3) = [];
    scatter3(locs(1,:),locs(2,:),locs(3,:),mags*5000,'b','filled')
    
    %% Run again to define global ROIs
    % This just takes the peak from the mean DICs image left and right
    % Define Left ROI
    [a1 xyz] = max( GA_sourceDICS_L.avg.coh(:));
    max_pos = GA_sourceDICS_L.pos(xyz,:);
    [max_xyz(1) max_xyz(2) max_xyz(3)] = ind2sub(GA_sourceDICS_L.dim,xyz);
    mask = makeSphereMask(GA_sourceDICS_L.dim,max_xyz,3);
    maskL = permute(mask,[2 1 3]);
    %  ROI_list_L = find(GA_sourceDICS_L.inside.*reshape(mask,[],1));
    
    % Define Right ROI
    [a1 xyz] = max( GA_sourceDICS_R.avg.coh(:));
    max_pos = GA_sourceDICS_R.pos(xyz,:);
    [max_xyz(1) max_xyz(2) max_xyz(3)] = ind2sub(GA_sourceDICS_R.dim,xyz);
    mask = makeSphereMask(GA_sourceDICS_R.dim,max_xyz,3);
    maskR = permute(mask,[2 1 3]);
    ROI_list_R = find(GA_sourceDICS_R.inside.*reshape(mask,[],1));
    
    %% Channel Level ROIs
    % Now search within repeats and conditions for maximums that fall
    % within the bounded sphere defined by the global maximum
    for cond = 1:2
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(subname{i},cond-1);
        for nr = 1:nrep
            for refN = 1:numel(ref_list)
                ref_chan = ref_list{refN};
                load([datapathr subname{i} '\ftdata\r' subname{i} '_DICS_source' condname{cond} 'nrep_' num2str(nr) '_' ref_chan],'source')
                ref_chan = ref_list{refN};
                %% source_cohmag(LM1/RM1,refn,nr,cond)
                cohvol = reshape(source.avg.coh,source.dim); cohvol = cohvol.*maskL;
                cohvol = cohvol(:); cohvol(source.pos(:,2)>=0) = NaN; %reshape(cohvol,source.dim)
                [a1 xyz] = max(cohvol(:));
                [max_xyz(1) max_xyz(2) max_xyz(3)] = ind2sub(source.dim,xyz);
                source_cohmag(1,refN,nr,cond) = a1;
                source_cohloc(1,refN,nr,cond) = xyz;
                
                cohvol = reshape(source.avg.coh,source.dim); cohvol = cohvol.*maskR;
                cohvol = cohvol(:); cohvol(source.pos(:,2)<=0) = NaN;% reshape(cohvol,source.dim)
                [a1 xyz] = max(cohvol(:));
                [max_xyz(1) max_xyz(2) max_xyz(3)] = ind2sub(source.dim,xyz);
                source_cohmag(2,refN,nr,cond) = a1;
                source_cohloc(2,refN,nr,cond) = xyz;
            end
        end
    end
    for cond = 1:2
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(subname{i},cond-1);
        for nr = 1:nrep
            %% STNselection{Right/Left,Ipsi/Contra,nr,cond}
            %% = {Cond, set, cohmax, ref_ind,refname,sourceind,sourcesub,sourcepos
            % Contra -Right STN, Left M1
            maxcohs = squeeze(source_cohmag(1,1:3,nr,cond));
            [a1 xyz] = max(maxcohs);
            maxlocs = squeeze(source_cohloc(1,1:3,nr,cond));
            [loc(1) loc(2) loc(3)] = ind2sub(source.dim,maxlocs(xyz));
            STNselection{1,2,nr,cond} = {condname{cond} 'contra_RSTN_LM1' a1 xyz ref_list{xyz} maxlocs(xyz)  loc source.pos(maxlocs(xyz),:)}
            
            % Ipsi -Left STN, Left M1
            maxcohs = squeeze(source_cohmag(1,4:6,nr,cond));
            [a1 xyz] = max(maxcohs);
            maxlocs = squeeze(source_cohloc(1,4:6,nr,cond));
            [loc(1) loc(2) loc(3)] = ind2sub(source.dim,maxlocs(xyz));
            STNselection{2,1,nr,cond} = {condname{cond} 'ipsi_LSTN_LM1' a1 xyz ref_list{xyz} maxlocs(xyz)  loc source.pos(maxlocs(xyz),:)}
            
            % Contra -Left STN, Right M1
            maxcohs = squeeze(source_cohmag(2,4:6,nr,cond));
            [a1 xyz] = max(maxcohs);
            maxlocs = squeeze(source_cohloc(2,4:6,nr,cond));
            [loc(1) loc(2) loc(3)] = ind2sub(source.dim,maxlocs(xyz));
            STNselection{2,2,nr,cond} = {condname{cond} 'contra_LSTN_RM1' a1 xyz ref_list{xyz} maxlocs(xyz)  loc source.pos(maxlocs(xyz),:)}
            
            % Ipsi -Right STN, Right M1
            maxcohs = squeeze(source_cohmag(2,1:3,nr,cond));
            [a1 xyz] = max(maxcohs);
            maxlocs = squeeze(source_cohloc(2,1:3,nr,cond));
            [loc(1) loc(2) loc(3)] = ind2sub(source.dim,maxlocs(xyz));
            STNselection{1,1,nr,cond} = {condname{cond} 'ipsi_RSTN_RM1' a1 xyz ref_list{xyz} maxlocs(xyz)  loc source.pos(maxlocs(xyz),:)}
        end
    end
    save([datapathr subname{i} '\ftdata\r' subname{i} '_DICS_source_channelselectioninfo'],'STNselection')
    
    
    % Subject Effects
    % source_cohmag(left/rightM1,refn,nr,cond)
    % IPSIlateral sources
    OFF = [squeeze(source_cohmag(1,4:6,:,2)) squeeze(source_cohmag(2,1:3,:,2))]; OFF = OFF(:); OFF(OFF==0) = NaN;
    ON = [squeeze(source_cohmag(1,4:6,:,1)) squeeze(source_cohmag(2,1:3,:,1))]; ON = ON(:); ON(ON==0) = NaN;
    
    figure
    [centre height] = aboxplotTW([ON,OFF],'Labels',{'ON';'OFF'},'Colorgrad','blue_down');
    ylabel('DICs Max. Coherence','FontSize',14);
    title('STN/DICs Ipsilateral High Beta Source Coherence','FontSize',14)
    [H P] = ttest(ON,OFF);
    dum = [ON, OFF]; ylim([0 1.2*max(dum(:))])
    annotation(gcf,'textbox',...
        [0.442619519094767 0.786158401184308 0.149892857142857 0.0880952380952401],...
        'String',{sprintf('P = %.3f',P)},...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',12,...
        'FitBoxToText','off');
    
    % CONTRAlateral sources
    OFF = [squeeze(source_cohmag(1,1:3,:,2)) squeeze(source_cohmag(2,4:6,:,2))]; OFF = OFF(:); OFF(OFF==0) = NaN;
    ON = [squeeze(source_cohmag(1,1:3,:,1)) squeeze(source_cohmag(2,4:6,:,1))]; ON = ON(:); ON(ON==0) = NaN;
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
    savefigure_v2([datapathr subname{1} '\images\sourcespace\'],['r' subname{1} '_grandaverage_ON_OFF'],[],[],[]); 
end
