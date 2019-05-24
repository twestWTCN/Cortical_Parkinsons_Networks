%function compute_connectivity_metric_ROI(datapathr,subname)
clear
I = i;
% clear;
close all;
datapathr = 'C:\Users\Tim\Documents\Work\Cortical_Networks\Data\';
subname = {'JB'}; %{'MW'};
i = 1;
% ref_list = {'STN_L12','STN_R23'};
% nr = 2;
side =1;sidenm = {'RSTN','LSTN'};
ia = 1;
condname = {'ON','OFF'};
load('ref_ci.mat')
for J = 1; %:numel(ref_list)
    for  cond = 1:2 % USE just OFF to identify ROI for now
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(subname{i},cond-1);
        for nr = 1; %:nrep
            load([datapathr subname{i} '\ftdata\virtual_sources_' num2str(nr) '_ROI_' condname{cond} '_' sidenm{side} '_contra'])
            clear XD
            %         [U,S,V] = svd(XD,'econ'); V= V';
            fsamp = 256;
            for x = 1:numel(vchansave)
                Xdata.fsample = fsamp;
                Xdata.label = vchansave(x).label;
                Xdata.trial{1} = vchansave(x).trial{1};
                Xdata.time{1} = vchansave(x).time{1};
                Xdata.trial{1} = Xdata.trial{1}(:,2*fsamp:(end-2*fsamp));
                Xdata.time{1} = Xdata.time{1}(:,2*fsamp:(end-2*fsamp));
                
                cfg = [];
                cfg.length = 1;
                Xseg = ft_redefinetrial(cfg,Xdata);
                
                % cfg = [];
                % cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
                % ft_databrowser(cfg, Xdata);
                
                cfg           = [];
                cfg.method    = 'mtmfft';
                cfg.taper     = 'dpss';
                cfg.output    = 'fourier';
                cfg.foilim = [4 48];
                % cfg.pad       = 128;
                cfg.tapsmofrq  = 1.5;
                freq         = ft_freqanalysis(cfg, Xseg);
                
                figure(1+10*cond)
                plot(repmat(freq.freq,2,1)',normaliseV(squeeze(mean(abs(freq.fourierspctrm(:,:,:)),1)))')
                hold on
                freqstore(:,:,x,cond) = normaliseV(squeeze(mean(abs(freq.fourierspctrm(:,:,:)),1)))';
                figure(2+10*cond)
                cfg           = [];
                cfg.method    = 'wpli_debiased';
                %         cfg.complex = 'imag';
                coh           = ft_connectivityanalysis(cfg, freq);
                
                icoh = abs(squeeze(coh.wpli_debiasedspctrm(2,1,:)));
                plot(coh.freq,icoh)
                hold on
                cohstore(:,x,cond) = icoh;
                
                [maxcoh(x) fi] = max(icoh(coh.freq>24 & coh.freq<40));
                frq(x) = 24+fi;
                
                stn_pow = squeeze(mean(abs(freq.fourierspctrm(:,2,:)),1));
                [dum fi] = max(icoh(coh.freq>13 & coh.freq<24));
                stn_lb_frq(x) = 13+fi;
                %                 close all
            end
            maxcoh(maxcoh>0.8) = NaN;
            [dum id] = max(maxcoh);
            idbank(nr,cond) = id;
            frqbank(nr,cond) = frq(id);
            stn_lb_frqbank(nr,cond) = stn_lb_frq(id);
            clear maxcoh frq
            
        end
    end
        figure(101)
        ON_STN = squeeze(freqstore(:,2,:,1)); ON_M1 = squeeze(freqstore(:,1,:,1));
        ON_STN(:,sum(ON_STN,1)==0) = [];     ON_M1(:,sum(ON_M1,1)==0) = [];
    
        OFF_STN = squeeze(freqstore(:,2,:,2)); OFF_M1 = squeeze(freqstore(:,1,:,2));
        OFF_STN(:,sum(OFF_STN,1)==0) = [];     OFF_M1(:,sum(OFF_M1,1)==0) = [];
    
        hold on
        plot(repmat(freq.freq,size(ON_STN,2),1)',ON_STN,'b--'); plot(freq.freq',mean(ON_STN,2),'b');
        plot(repmat(freq.freq,size(OFF_STN,2),1)',OFF_STN,'r--'); plot(freq.freq',mean(OFF_STN,2),'r');
    
    for cond = 1:2
        id = idbank(nr,2);
        frq = frqbank(nr,cond);
        stn_lb_frq = frqbank(nr,cond);
        load([datapathr subname{i} '\ftdata\virtual_sources_' num2str(nr) '_ROI_' condname{cond}  '_' sidenm{side} '_contra'])
        load([datapathr subname{i} '\ftdata\r' subname{i} '_LCMV_source_' condname{cond} 'nrep_' num2str(nr)])
        % cfg           = [];
        % cfg.parameter = 'wplispctrm';
        % cfg.zlim      = [0 1];
        % ft_connectivityplot(cfg, coh) %, coh);
        %
        % for x=1
        %% Dynamic
        %                 Xdata = Xdata;
        Xdata.fsample = fsamp;
        %
        Xdata.fsample = fsamp;
        Xdata.trial{1} = vchansave(id).trial{1};
        Xdata.time{1} = vchansave(id).time{1};
        Xdata.trial{1} = Xdata.trial{1}(:,2*fsamp:(end-2*fsamp));
        Xdata.time{1} = Xdata.time{1}(:,2*fsamp:(end-2*fsamp));
        %
        bwid = 1.5;
        % bandpass
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [frq-bwid frq+bwid];
        Xdata = ft_preprocessing(cfg,Xdata);
        fsamp = Xdata.fsample;
        
        amp1 = abs(hilbert(Xdata.trial{1}(1,:)));
        amp2 = abs(hilbert(Xdata.trial{1}(2,:)));
        
        phi1 = angle(hilbert(Xdata.trial{1}(1,:)));
        phi2 = angle(hilbert(Xdata.trial{1}(2,:)));
        dphi_12 = unwrap(phi1-phi2);
        ampw = (amp1.*amp2)./(max(amp1)*max(amp2));  %%
        dphi_12 = angle(ampw.*exp(I.*(phi1-phi2)));%%
        dphi_12 = (dphi_12-(1/sqrt(length(phi1))))./(1-(1/sqrt(length(phi1)))); %%
        dphi_12_dt = diff(dphi_12);
        
        % lb_stn_bandpass
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [stn_lb_frq-bwid stn_lb_frq+bwid];
        stn_lb_data = ft_preprocessing(cfg,Xdata);
        stn_lb_amp1 = abs(hilbert(stn_lb_data.trial{1}(1,:)));

        
        
        mwid = 4; % Minimum number of cycles to consider sync
        
        period = (mwid/20)*fsamp;
        [slide_dphi_12,sind] = slideWindow(phi1-phi2, floor(period), floor(period*0.75));
        %         PLV = abs(mean(exp(1i*slide_dphi_12),1));
        %         PLV_tvec = Xdata.time{1}(round(median(sind,1)));
        
        PLI = abs(mean(sign(slide_dphi_12),1));  % PLI
        %         PLI= abs(mean(abs(slide_dphi_12).*sign(slide_dphi_12),1))./ mean(abs(slide_dphi_12),1);  %wPLI
        PLI_tvec = Xdata.time{1}(round(median(sind,1)));
        
        
        %%% Find segments
        dphi_12_dt_sm = smooth(dphi_12_dt,period/mwid)';
        qstable = find(abs(dphi_12_dt_sm)<0.035); % The points that are below threshold
        consecSegs = SplitVec(qstable,'consecutive');
        % lengths
        segL_ddt = cellfun('length',consecSegs);
        seg_ddt = [consecSegs{segL_ddt>(period/mwid)}];
        segL_ddt = segL_ddt(segL_ddt>(period/mwid))/fsamp;
        segL_ddt_save{cond} = segL_ddt;
        % phase angles of segments
        dphi_12_uw = phi1-phi2;
        clear pA_dist amp_dist
        for j = 1:numel(consecSegs)
            if length(consecSegs{j}) > period/mwid
                pA_dist(j) = circ_mean(wrapTo2Pi(dphi_12_uw(consecSegs{j}))');
                amp_dist(j) = max(stn_lb_amp1(consecSegs{j}))/mean(stn_lb_amp1);
            else
                pA_dist(j) = NaN;
                amp_dist(j) = NaN;
            end
        end
        pA_dist(isnan(pA_dist)) = [];
        pA_dist_save{cond} = pA_dist;
        amp_dist(isnan(amp_dist)) = [];
        %PLI segments
        qstable = find(abs(PLI)>wpli_ci); % The points that are below threshold
        consecSegs = SplitVec(qstable,'consecutive');
        % lengths
        segL_pli = cellfun('length',consecSegs);
        segL_pli = segL_pli(segL_pli>mwid);
        segL_pli = (segL_pli/mwid)*(1/20);
        segL_pli_save{cond} = segL_pli;
        % phase angles of segments
        dphi_12_uw = phi1-phi2;
        clear pA_pli_dist
        for j = 1:numel(consecSegs)
            if length(consecSegs{j}) > mwid
                pA_pli_dist(j) = circ_mean(wrapTo2Pi(dphi_12_uw(consecSegs{j}))');
            else
                pA_pli_dist(j) = NaN;
            end
        end
        pA_pli_dist(isnan(pA_pli_dist)) = [];
        
        pA_pli_dist_save{cond} = pA_pli_dist;
        %% Plots
        figure(cond*10 + 1)
        ax(1) = subplot(4,1,1);
        plot(Xdata.time{1},normaliseV(Xdata.trial{1}(1,:)),Xdata.time{1},normaliseV(Xdata.trial{1}(2,:)));%xlim([60 70])
        ylabel('\beta activity'); title(['MaxCohM1 ' condname{cond}])
        
        
        ax(2) = subplot(4,1,2);
        plot(Xdata.time{1},phi1,Xdata.time{1},phi2);%xlim([60 70])
        ylabel('\phi_{1/2}');
        
        
        % ddt Time Series
        ax(3) = subplot(4,1,3);
        plot(Xdata.time{1}(2:end),dphi_12_dt); %xlim([60 70]);
        
        tvec = nan(size(Xdata.time{1}(2:end))); tvec(seg_ddt) = Xdata.time{1}(seg_ddt);
        yvec = nan(size(dphi_12_dt)); yvec(seg_ddt) = dphi_12_dt(seg_ddt);
        hold on; plot(tvec,yvec,'LineWidth',2)
        hold on; plot([0 Xdata.time{1}(end)],[ddphi_ci ddphi_ci],'k--');
        plot([0 Xdata.time{1}(end)],[-ddphi_ci -ddphi_ci],'k--')
        
        ylabel('$$ d \frac{(\phi_1 - \phi_2)}{dt} $$','Interpreter','latex')
        
        %         % PA Time Series
        %         ax(3) = subplot(4,1,3);
        %         plot(Xdata.time{1}(:),wrapToPi(unwrap(phi1-phi2))); %xlim([60 70]);
        %         ylabel('\phi_1 - \phi_2')
        
        % PLI Time Series
        ax(4) = subplot(4,1,4);
        %             plot(PLV_tvec,PLV);hold on
        plot(PLI_tvec,PLI); ylim([0 1]); % xlim([60 70]);
        
        tvec = nan(size(PLI_tvec)); tvec([consecSegs{:}]) = PLI_tvec([consecSegs{:}]);
        yvec = nan(size(PLI)); yvec([consecSegs{:}]) = PLI([consecSegs{:}]);
        hold on; plot(tvec,yvec,'LineWidth',2)
        ylabel('PLV/PLI'); xlabel('Time (s)'); legend({'PLV','PLI'})
        hold on; plot([0 Xdata.time{1}(end)],[wpli_ci wpli_ci],'k--')
        %             hold on; plot([0 Xdata.time{1}(end)],[plv_ci plv_ci],'k--')
        linkaxes(ax,'x');
        set(gcf,'Position',[287 72 1446 932]);
        xlim([80 85])
        shg
        
        figure(cond*10 + 4 )
        cmap = [0 0 0; 1 0 0];
        [N,Xedges,Yedges] = histcounts2([pA_dist],segL_ddt,-pi:(pi/3):pi,0.1:.05:.8)
        N = N./sum(N(:));
%         imagesc(Xedges,Yedges,N');
        hold on
        scatter(pA_dist,segL_ddt,[],cmap(cond,:),'filled'); set(gca,'YDir','normal');
        xlabel('phi_1 - phi_2'); ylabel('segment length (s)')
        ylim([0.1 0.8]); xlim([-3.5 3.5]); colorbar; caxis([0 0.05])
        title(['MaxCohM1 ' condname{cond}])
        bedge = linspace(-pi,pi,7);  clear segLbinned
        for be = 1:numel(bedge)-1
            segLbinned(:,be) = [mean(segL_ddt(pA_dist>=bedge(be) & pA_dist<bedge(be+1))) prctile(segL_ddt(pA_dist>=bedge(be) & pA_dist<bedge(be+1)),90) prctile(segL_ddt(pA_dist>=bedge(be) & pA_dist<bedge(be+1)),9)];
        end
        binmid = bedge(1:end-1) + diff(bedge(1:2));
        plot(binmid,segLbinned(1,:),'Color',cmap(cond,:),'LineWidth',4); plot(binmid,segLbinned(2,:),'Color',cmap(cond,:),'LineWidth',4,'linestyle','--'); plot(binmid,segLbinned(3,:),'Color',cmap(cond,:),'LineWidth',4,'linestyle','--')
        
        figure(cond*10 + 5)
        [N,Xedges,Yedges] = histcounts2([pA_dist],amp_dist,-pi:(pi/3):pi,0.1:.25:3)
        N = N./sum(N(:));
%         imagesc(Xedges,Yedges,N');
        hold on
        scatter(pA_dist,amp_dist,[],cmap(cond,:),'filled'); set(gca,'YDir','normal');
        xlabel('phi_1 - phi_2'); ylabel('cortical amp')
        ylim([0.1 3]); xlim([-3.5 3.5]); colorbar; caxis([0 0.05])
        title(['MaxCohM1  ' condname{cond}])
        bedge = linspace(-pi,pi,7);  clear segLbinned
        for be = 1:numel(bedge)-1
            segLbinned(:,be) = [mean(amp_dist(pA_dist>=bedge(be) & pA_dist<bedge(be+1))) prctile(amp_dist(pA_dist>=bedge(be) & pA_dist<bedge(be+1)),90) prctile(amp_dist(pA_dist>=bedge(be) & pA_dist<bedge(be+1)),10)];
        end
        binmid = bedge(1:end-1) + diff(bedge(1:2));
        plot(binmid,segLbinned(1,:),'Color',cmap(cond,:),'LineWidth',4);  plot(binmid,segLbinned(2,:),'Color',cmap(cond,:),'LineWidth',4,'linestyle','--'); plot(binmid,segLbinned(3,:),'Color',cmap(cond,:),'LineWidth',4,'linestyle','--')
        % end
        %                 savefigure_v2([datapathr subname{ia} '\images\seganalysis\'],['seg_analysis_' condname{cond} '_MaxCohM1 '],[],[],[]); close all
        
        % Segment Length DDT
        figure(1)
        X = segL_ddt_save{cond};
        pA = pA_dist_save{cond};
        intvl = [circ_mean(wrapTo2Pi(pA)')-(pi/4) circ_mean(wrapTo2Pi(pA)')+(pi/4)]
        X = X(intvl(1) < pA & pA<intvl(2))
        H1(cond) = histogram(X,logspace(-1.7,0.3,20),'Normalization','probability'); hold on
        set(gca,'xscale','log')
        X = (H1(cond).BinEdges(2:end)); Y = (H1(cond).Values); %-(H1(cond).BinWidth/2))
        f = fit(X',Y','exp1');
        plot(f,X,Y);
        R2 = 1 - sum((Y'-f(X)).^2) ./ sum((Y'-mean(Y)).^2);
        model = sprintf('f(x) = %.2f*exp(%.2f*x)',f.a,f.b);
        annotation(gcf,'textbox',...
            [0.520 .26-((cond-1)*0.1) 0.37 0.24],...
            'String',{sprintf('R2 = %.2f',R2); model},...
            'LineStyle','none',...
            'HorizontalAlignment','right',...
            'FontSize',12,...
            'FitBoxToText','off');
        xlabel('Segment Duration (s)'); ylabel('P(X)'); title('Phase Lock')
        
        %             histogram(segL_pli,50,'Normalization','probability','BinMethod','sturges');
        xlabel('Seg Length'); ylabel('P(x)');
        %         xlim([0 0.8]); ylim([0 0.3])
        
        % Phase Angle Rose DDT
        figure(2)
        pA_dist = [pA_dist_save{cond}];
        pA_dist(isnan(pA_dist)) = [];
        R_1(cond) = polarhistogram(pA_dist,36,'FaceAlpha',0.75); hold on
        %             x = get(h,'Xdata');
        %             y = get(h,'Ydata');
        %             g=patch(x,y,'r');
        title('Phase Angle Distribution Phase Lock','FontSize',18)
        if cond == 2
            pval = circ_cmtest(pA_dist_save{1}, pA_dist_save{2});
            annotation(gcf,'textbox',...
                [0.71 0.0690 .250 .16],...
                'String',{'Circular CM test',sprintf('P: %.3f',pval)},...
                'FitBoxToText','off');
        end
        % Segment Length PLI
        figure(3)
        H2(cond) = histogram(segL_pli_save{cond},logspace(-1.5,0.3,20)); hold on; %,'BinMethod','sturges'
        set(gca,'xscale','log')
        hold on
        X = (H2(cond).BinEdges(2:end)); Y = (H2(cond).Values); %-(H2(cond).BinWidth/2))
        f = fit(X',Y','exp1');
        plot(f,X,Y);
        R2 = 1 - sum((Y'-f(X)).^2) ./ sum((Y'-mean(Y)).^2);
        model = sprintf('f(x) = %.2f*exp(%.2f*x)',f.a,f.b);
        annotation(gcf,'textbox',...
            [0.520 .26-((cond-1)*0.1) 0.37 0.24],...
            'String',{sprintf('R2 = %.2f',R2); model},...
            'LineStyle','none',...
            'HorizontalAlignment','right',...
            'FontSize',12,...
            'FitBoxToText','off');
        xlabel('Segment Duration (s)'); ylabel('P(X)'); title('PLI')
        
        %             histogram(segL_pli,50,'Normalization','probability','BinMethod','sturges');
        xlabel('Seg Length'); ylabel('P(x)');
        %         xlim([0 0.8]); ylim([0 0.3])
        
        % Phase Angle Rose - PLI
        figure(4)
        pA_pli_dist = [pA_pli_dist_save{cond}];
        pA_pli_dist(isnan(pA_pli_dist)) = [];
        R_2(cond) = polarhistogram(pA_pli_dist,36); hold on
        %             x = get(h,'Xdata');
        %             y = get(h,'Ydata');
        %             g=patch(x,y,'r');
        title('Phase Angle Distribution PLI')
        if cond == 2
            
            pval = circ_cmtest(pA_pli_dist_save{1}, pA_pli_dist_save{2});
            annotation(gcf,'textbox',...
                [0.71 0.0690 .250 .16],...
                'String',{'Circular CM test',sprintf('P: %.3f',pval)},...
                'FitBoxToText','off');
        end
        
    end
    figure(1)
    legend(H1,condname)
    figure(2)
    legend(R_1,condname)
    figure(3)
    legend(H2,condname)
    figure(4)
    legend(R_2,condname)
    savefigure_v2([datapathr subname{ia} '\images\seganalysis\'],['seg_analysis_ONvsOFF_MaxCohM1 '],[],[],[]); close all
end

% save([datapathr subname{ia} '\ftdata\ROI_source_granger_' ref_chan '_' num2str(nr)],'gASym')
% save([datapathr subname{ia} '\ftdata\ROI_source_coh_' ref_chan '_' num2str(nr)],'cASym')

% load('C:\Users\Tim\Documents\Work\Cortical_Networks\tmp_files\source_granger.mat')
% load([datapathr subname{i} '\ftdata\r' subname{i} 'rs'],'mri')
% %                 cfg = [];
% %                 cfg.parameter = 'avg.pow';
% %                 cfg.filename = [datapathr subname{i} '_tmp\MRI\source_data']
% %                 ft_sourcewrite(cfg, source)
% load([datapathr subname{i} '\ftdata\rs_transform_vox2ctf'])
%
% T = transform_vox2ctf; %/transform_vox2spm;%;
% Tmri = ft_transform_geometry(T,mri);
% Tmri = ft_convert_units(Tmri,'cm');
%
% load('C:\Users\Tim\Documents\Work\Cortical_Networks\tmp_files\source_LM.mat')
% source.avg.granger = reshape(gASym,[],1);
% load([datapathr subname{i} '\ftdata\rs_transform_vox2ctf'])
% cfg            = [];
% % cfg.downsample = 1;
% cfg.parameter = 'granger';
% sourceInt  = ft_sourceinterpolate(cfg, source , Tmri);
%
% %                 Z-score power
% %                     sourceInt.coh = reshape(sourceInt.coh,size(sourceInt.anatomy));
% %                     incoh = sourceInt.coh.*sourceInt.inside;
% %                     incoh(incoh==0) = NaN;
% %                     incoh = (incoh - nanmean(incoh(:)))/nanstd(incoh(:));
% %                     %                                     inpow(isnan(inpow)) = 0;
% %                     sourceInt.coh = incoh;
% %                     sourceInt.coh = reshape(sourceInt.coh,[],1);
%
% cfg = [];
% cfg.method        = 'slice';
% cfg.funparameter  = 'granger';
% cfg.maskparameter = cfg.funparameter;
% % cfg.funcolorlim   = [4.0 6.2];
% % cfg.opacitylim    = [4.0 6.2];
% cfg.opacitymap    = 'rampup';
% figure
% ft_sourceplot(cfg, sourceInt);
%
% sourceInt.coordsys = 'ctf';
% cfg = [];
% cfg.nonlinear     = 'no';
% sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
% volumewrite_spm([datapathr subname{i} '\MRI\sourcespace\r' subname{i} '_' num2str(nr) '_DICS_sourcegranger_int_' ref_chan '.nii'], sourceIntNorm.granger, sourceIntNorm.transform, 'SPM12');
%
% cfg = [];
% cfg.method         = 'surface';
% cfg.funparameter   = 'granger';
% cfg.maskparameter  = cfg.funparameter;
% % cfg.funcolorlim    = [0.0 1.2];
% cfg.funcolormap    = 'jet';
% % cfg.opacitylim     = [0.0 1.2];
% cfg.opacitymap     = 'auto';
% cfg.projmethod     = 'nearest';
% cfg.surffile       = 'surface_white_both.mat';
% cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% H = figure;
% ft_sourceplot(cfg, sourceIntNorm);
% view ([90 0])
%
% savefig(gcf,[datapathr subname{i} '\images\sourcespace\r' subname{i} '_rep_' num2str(nr) '_surf_avg_DICS_granger_' ref_chan])
%


%             cfg = [];
%             cfg.length = 1.5;
% %             cfg.overlap = 0.25;
%             data = ft_redefinetrial(cfg,Xdata);
%
%             cfg = [];
%             cfg.trials = t;
%             data = ft_selectdata(cfg,data);
%
%             cfg         = [];
%             cfg.order   = 8;
%             cfg.method = 'bsmart';
%             mdata       = ft_mvaranalysis(cfg, data);
%
%             cfg        = [];
%             cfg.method = 'mvar';
%             mfreq      = ft_freqanalysis(cfg, mdata);
%
%             cfg           = [];
%             cfg.method    = 'mtmfft';
%             cfg.taper     = 'dpss';
%             cfg.output    = 'fourier';
%             cfg.pad       = 128;
%             cfg.tapsmofrq = 1;
%             freq          = ft_freqanalysis(cfg, data);
%
%             cfg           = [];
%             cfg.method    = 'coh';
%             cfg.complex     = 'imag';
%             %             coh           = ft_connectivityanalysis(cfg, freq);
%             cohm          = ft_connectivityanalysis(cfg, freq);
%             indz = find(cohm.freq>=20 & cohm.freq<30);
%             cA(t) = mean(squeeze(cohm.cohspctrm(2,1,indz)));
%
%
%             figure
%             cfg           = [];
%             cfg.parameter = 'cohspctrm';
%             cfg.zlim      = [0 1];
%             ft_connectivityplot(cfg, cohm) %, coh);

%             cfg           = [];
%             cfg.method    = 'granger';
%             granger       = ft_connectivityanalysis(cfg, mfreq);
%             indz = find(granger.freq>=20 & granger.freq<30);
%             gA(t) = mean(squeeze(granger.grangerspctrm(2,1,indz) - granger.grangerspctrm(1,2,indz)));
%             figure
%             cfg           = [];
%             cfg.parameter = 'grangerspctrm';
%             cfg.zlim      = [0 0.1];
%             ft_connectivityplot(cfg, granger);
