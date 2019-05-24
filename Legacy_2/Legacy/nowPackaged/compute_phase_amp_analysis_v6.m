function compute_phase_amp_analysis_v6(R)
if nargin<2
    R = makeHeader_SubCort_Cort_Networks();
end
%%%
close all
for sub = 1:length(R.subname)
    for breg = 1:length(R.bregname)
        for side = 1:2
            for cond = 1:length(R.condname)
                load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                % Clear previous computations
                if isfield(vc_clean,'PA')
                    vc_clean = rmfield(vc_clean,'PA');
                    vc_clean.history = {vc_clean.history{~strncmp(vc_clean.history,mfilename,numel(mfilename))}};
                end
                
                if vc_clean.specanaly.flag ~= 1 % Check significant coherences
                    band = R.bregband{breg};
                    frq = vc_clean.specanaly.cohstats.frqcoh(R.bregband{breg}); % Alpha and high beta frequencies
                    stn_lb_frq = vc_clean.specanaly.cohstats.frqcoh(strncmp(R.bandname,'Low Beta',6)); % STN frequencies
                    
                    % Create data structure
                    Xdata = vc_clean;
                    Xdata = rmfield(Xdata,'specanaly');
                    % Get the overal signal amplitudes
                    sig =  Xdata.trial{1};
                    signalEnvAmp = median(abs(hilbert(sig)),2);
                    
                    % Compute optimal frequencies based on PLV
                    [maxfrq,maxPLV] = compute_optimal_PhaseLockFrq(R,Xdata,band,stn_lb_frq); %
                    [SRPeps Ampeps SNR_eps PLVeps] = phase_amp_surrComp(R,Xdata,band,maxfrq,stn_lb_frq,R.PA.AmpSurrN);
                    SNR_eps_z(1) = 10*log10(SNR_eps(1)./signalEnvAmp(1));
                    SNR_eps_z(2) =  10*log10(SNR_eps(2)./signalEnvAmp(2));
                    SNR_eps_z(3) =  10*log10(SNR_eps(3)./signalEnvAmp(2));
                    % Compute data transforms (Hilbert)
                    [amp phi dphi_12 dphi_12_dt Xdata] = comp_instant_angle_phase(Xdata,maxfrq,stn_lb_frq,R.PA.bwid,band);
                    
                    % Sliding Window PLV
                    if R.PA.SType == 1
                        fsamp = Xdata.fsample;
%                         swsize = R.PA.slidingwindow;
                        swsize = floor((10/maxfrq)*fsamp);
                        WinSize = swsize;
                        [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,R.PA.WinOver,R.PA.optimalPLFrqMeth);
                        PLV_tvec = Xdata.time{1}(PLV_tvec);
                        amp_sw = cont2slidingwindow(amp,WinSize,floor(R.PA.WinOver*WinSize));
                        clear SNR_sw
                        SNR_sw(:,1) = 10.*log10(amp_sw(:,1)./signalEnvAmp(1));
                        SNR_sw(:,2) =  10.*log10(amp_sw(:,2)./signalEnvAmp(2));
                        SNR_sw(:,3) =  10.*log10(amp_sw(:,3)./signalEnvAmp(2));
%                                             figure(1)
%                                             SNR_Inspector(R,PLV_tvec,SNR_sw,SNR_eps_z,band)
                        dphi_12_sw = cont2slidingwindow(dphi_12,WinSize,round(R.PA.WinOver*WinSize));
                        SW_sampr = max(diff(PLV_tvec));
                        tseries = dphi_12_sw; qstable = find(PLV>PLVeps);
                        mwid = R.PA.mwid;
                        period = (mwid/frq)/SW_sampr;
                        % Minimum number of cycles to consider sync
                        [phi_dist amp_dist seg_ddt1 segL_ddt consecSegs H] = analysestablesegs(qstable,tseries,amp,period,mwid,fsamp,SNR_eps_z(1:2),[],[],Ampeps,SNR_sw);
                
%                         any(amp_dist(:)>300)
%                         plot_example_phaseanalysis_SW(Xdata,amp,phi,PLV,seg_ddt1,PLVeps,PLV_tvec);
                    elseif R.PA.SType == 2 % Sliding Window PhaseAng. Stability
                        clear SNR_sw
                        SNR_sw(:,1) = 10.*log10(amp(:,1)./signalEnvAmp(1));
                        SNR_sw(:,2) =  10.*log10(amp(:,2)./signalEnvAmp(2));
                        SNR_sw(:,3) =  10.*log10(amp(:,3)./signalEnvAmp(2));
                        % Plot SNR
%                         figure(1)
%                         SNR_Inspector(R,Xdata.time{1},SNR_sw,SNR_eps_z,breg)
                        
                        %%% Set Length Constraints
                        fsamp = Xdata.fsample;
                        mwid = R.PA.mwid; % Minimum number of cycles to consider sync
                        period = (mwid/maxfrq)*fsamp;
                        
                        %%% Find segments
                        %                 dphi_12_dt_sm = smooth(dphi_12_dt,period/6)';
                        qstable = find(abs(dphi_12_dt')<SRPeps); % The points that are below threshold
                        tseries = wrapToPi(dphi_12(2:end));
                        [phi_dist amp_dist seg_ddt1 segL_ddt consecSegs H] = analysestablesegs(qstable,tseries,amp,period,mwid,fsamp,SNR_eps_z(1:2),[],[],Ampeps,SNR_sw);
                        figure(2)
                        plot_example_phaseanalysis_trace(Xdata,amp,phi,dphi_12_dt,seg_ddt1,SRPeps,fsamp);
                        clf(1); clf(2)
                    end
                else
                    phi_dist = NaN;
                    amp_dist =NaN;
                    segL_ddt = NaN;
                end
                % Save to data
                vc_clean.PA.pA_pli_dist_save = phi_dist;
                vc_clean.PA.amp_pli_dist_save = amp_dist;
                vc_clean.PA.segL_pli_dist_save = segL_ddt;
                vc_clean.PA.timevec = Xdata.time;
                vc_clean.PA.H_dist_save{cond} = H;
                vc_clean.history = [vc_clean.history{:} {[mfilename '_' date]}];
                save([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                
                % % Possible Granger stuff and Network
                % %                 [gc_dist] = analysestablesegs_granger(qstable,tseries,amp,period,mwid,R.pp.cont.full.fs);
                % %                 gc_dist_save{cond} = NaN; %gc_dist;
                % %                 H_dist_save{cond,nr} = H;
                % %                 threshnetwork.PLV = PLV;
                % %                 threshnetwork.PLV_tvec = PLV_tvec;
                % %                 threshnetwork.stn_lb_frq = R.PA.stn_lb_frq;
                % %                 threshnetwork.maxfrq = maxfrq;
                % %                 threshnetwork.betaS = betaS;
                % %                 threshnetwork.amp = amp;
                % %                 threshnetwork.Xdata = Xdata;
                % %                 threshnetwork.consecSegs = consecSegs;
                % %                 mkdir([R.datapathr R.subname{sub} '\ftdata\networks\'])
                % %                 save([R.datapathr R.subname{sub} '\ftdata\networks\threshnetwork_' R.ipsicon '_' R.siden{side} '_' R.bandname{band} '_' R.condname{cond}],'threshnetwork')
                disp([cond side breg sub])
            end
        end
    end
end
% ! shutdown /h