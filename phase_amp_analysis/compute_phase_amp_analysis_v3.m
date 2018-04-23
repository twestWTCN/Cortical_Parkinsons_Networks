function compute_phase_amp_analysis_v3(R,idd)
if nargin<2
    R = makeHeader_SubCort_Cort_Networks();
end
if nargin<2
    idd = '';
end
%%%
% This script computes the instantaneous phase sync analyses. We partition
% the recording into stable segments using threshold, can compute Granger
% etc within frames. DFA is also computed PS/AE versions. Barplots for DFA
% at end of script.
%%%
for band = 3; %[1 3] %1:numel(R.bandname)
    for sub = 1:numel(R.subname)
        if exist([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon '_' R.bandname{band} '.mat']) ==0
            [idbank frqbank stn_lb_frqbank] = find_voxel_pow_coh_v3(R,sub,band)
        else
            load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon '_' R.bandname{band}])
        end
        for side = 1:2
            for cond = 1:2
                [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
                [~,~,nrepOFF,~] = data_fileguide(R.subname{sub},cond-1);
                for nr = 1:nrepOFF
                    id = idbank(nr,side,cond);
                    frq = frqbank(band,nr,side,cond); %% THIS USES THE 3rd Band (HB)
                    stn_lb_frq = stn_lb_frqbank(nr,side,cond);
%                     load([R.datapathr R.subname{sub} '\ftdata\virtual_sources_' num2str(nr) '_ROI_' R.condname{cond}  '_' R.siden{side} '_' R.ipsicon])
                    load([R.datapathr R.subname{sub} '\ftdata\virtual_sources_' num2str(nr) '_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bandname{band}])
                    
                    %                 load([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_LCMV_source_' R.condname{cond} 'nrep_' num2str(nr)])
                    
                    % Create data structure
                    Xdata.fsample = R.pp.cont.full.fs;
                    Xdata.label = vchansave(id).label;
                    Xdata.trial{1} = vchansave(id).trial{1};
                    Xdata.time{1} = vchansave(id).time; %{1}
                    
                    %                     for ch = 1:2
                    %                         X = Xdata.trial{1}(ch,:);
                    %                         Xdata.trial{1}(ch,:) = (X-mean(X))./std(X);
                    %                     end
                    signalEnvAmp = median(abs(hilbert(Xdata.trial{1})),2);
                    if 1 %exist([R.datapathr R.subname{sub} '\ftdata\ROI_analy\' R.bandname{band} '_' R.condname{cond} '_nr_' num2str(nr) '_' R.siden{side} '_optfreq.mat']) ==0
                        [maxfrq maxPLV] = PLV_compute_optimalFrq(Xdata,R,band); % THIS USES THE HIGH BETA BAND FOR OPTIMIZATION
                        save([R.datapathr R.subname{sub} '\ftdata\ROI_analy\' R.bandname{band} '_' R.condname{cond} '_nr_' num2str(nr) '_' R.siden{side} '_optfreq'],'maxfrq','maxPLV')
                    else
                        load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\' R.bandname{band} '_' R.condname{cond} '_nr_' num2str(nr) '_' R.siden{side} '_optfreq'],'maxfrq','maxPLV')
                    end
                    % Compute data transforms (Hilbert)
                    [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Xdata,maxfrq,stn_lb_frq,R.PA.bwid(band),Xdata.fsample,R.PA.LowAmpFix);
                    
                    % Sliding Window PLV
                    if R.PA.SType == 1
                        WinSize = R.PA.slidingwindow*R.pp.cont.full.fs;
                        [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,R.PA.WinOver);
                        PLV_tvec = Xdata.time{1}(PLV_tvec);
                        amp_sw = cont2slidingwindow(amp,WinSize,floor(R.PA.WinOver*WinSize));
                        clear snr_sw
                        snr_sw(:,1) = log10(amp_sw(:,1)./signalEnvAmp(1));
                        snr_sw(:,2) =  log10(amp_sw(:,2)./signalEnvAmp(2));
                        snr_sw(:,3) =  log10(amp_sw(:,3)./signalEnvAmp(2));
                        dphi_12_sw = cont2slidingwindow(dphi_12,WinSize,round(R.PA.WinOver*WinSize));
                        
                        SW_sampr = max(diff(PLV_tvec));
                        tseries = dphi_12_sw; qstable = find(PLV>R.PA.PLVeps);
                        mwid = R.PA.mwid;
                        period = (mwid/frq)/SW_sampr;
                        % Minimum number of cycles to consider sync
                        [phi_dist amp_dist seg_ddt segL_ddt consecSegs H] = analysestablesegs(qstable,tseries,snr_sw,period,mwid,1/SW_sampr,R.PA.SNR(band));
                        
                    elseif R.PA.SType == 2 % Sliding Window PhaseAng. Stability
                        clear snr_sw
                        snr_sw(:,1) = log10(amp(:,1)./signalEnvAmp(1));
                        snr_sw(:,2) =  log10(amp(:,2)./signalEnvAmp(2));
                        snr_sw(:,3) =  log10(amp(:,3)./signalEnvAmp(2));
                        
                        mwid = R.PA.mwid; % Minimum number of cycles to consider sync
                        period = (mwid/frq)*R.pp.cont.full.fs;
                        cycle = (1/frq)*R.pp.cont.full.fs;
                        %%% Find segments
                        %                 dphi_12_dt_sm = smooth(dphi_12_dt,period/6)';
                        qstable = find(abs(dphi_12_dt')<R.PA.SRPeps); % The points that are below threshold
                        tseries = wrapToPi(dphi_12(2:end));
                        %                 tseries(:,2) = phi(:,2);
                        %                 tseries(:,1) = phi(:,1);
                        [phi_dist amp_dist seg_ddt1 segL_ddt consecSegs H] = analysestablesegs(qstable,tseries,amp,period,mwid,R.pp.cont.full.fs,R.PA.SNR(band),cycle);
%                         plot_example_phaseanalysis_trace(betaS,amp,phi,dphi_12_dt,seg_ddt1,R.PA.SRPeps);
                    end
                    
                    %                 [gc_dist] = analysestablesegs_granger(qstable,tseries,amp,period,mwid,R.pp.cont.full.fs);
                    gc_dist_save{cond,nr} = NaN; %gc_dist;
                    H_dist_save{cond,nr} = H;
                    pA_pli_dist_save{cond,nr} = phi_dist;
                    amp_pli_dist_save{cond,nr} = amp_dist;
                    segL_pli_dist_save{cond,nr} = segL_ddt;
                    timevec{cond,nr} = vchansave(id).time; %{1}
                    
%                     threshnetwork.PLV = PLV;
%                     threshnetwork.PLV_tvec = PLV_tvec;
                    threshnetwork.stn_lb_frq = stn_lb_frq;
                    threshnetwork.maxfrq = maxfrq;
                    threshnetwork.betaS = betaS;
                    threshnetwork.amp = amp;
                    threshnetwork.Xdata = Xdata;
                    threshnetwork.consecSegs = consecSegs;
                    mkdir([R.datapathr R.subname{sub} '\ftdata\networks\'])
                    save([R.datapathr R.subname{sub} '\ftdata\networks\threshnetwork_' R.ipsicon '_' R.siden{side} '_' R.bandname{band} '_' R.condname{cond}],'threshnetwork')
                    %                     figure
                    %                     PLV_sw_plot(Xdata,betaS,amp,phi,snr_sw,seg_ddt,PLV,PLV_tvec,consecSegs,R)
                    %                     %                 savefigure_v2([R.datapathr 'results\images\seganalysis\'],['example_seg_subject1_ON'],[],[],[]);
                    %                                         close all
                end
            end
            save([R.datapathr R.subname{sub} '\ftdata\ROI_analy\' idd '_ROIvoxel_phaseamp_PLI_' R.ipsicon '_' R.siden{side} '_' R.bandname{band}],...
                'pA_pli_dist_save','amp_pli_dist_save','segL_pli_dist_save','H_dist_save','timevec','gc_dist_save') %
        end
    end
end
% ! shutdown /h