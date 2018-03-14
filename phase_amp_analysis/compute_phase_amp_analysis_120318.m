function compute_phase_amp_analysis_120318(R,idd)
if nargin<2
     idd = '';
end
%%%
% This script computes the instantaneous phase sync analyses. We partition
% the recording into stable segments using threshold, can compute Granger
% etc within frames. DFA is also computed PS/AE versions. Barplots for DFA
% at end of script.
%%%

for sub = 1:numel(R.subname)
    if 1; %exist([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon '.mat']) ==0
        [idbank frqbank stn_lb_frqbank] = find_voxel_pow_coh_050118(R.datapathr,R.subname{sub},R.condname,R.siden,R.pp.cont.full.fs,R.ipsicon)
    else
        load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon])
    end
    for side = 1:2
        for cond = 1:2
            [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
            [~,~,nrepOFF,~] = data_fileguide(R.subname{sub},cond-1);
            for nr = 1:nrepOFF
                id = idbank(1,side,cond);
                frq = frqbank(1,side,cond);
                stn_lb_frq = frqbank(1,side,cond);
                load([R.datapathr R.subname{sub} '\ftdata\virtual_sources_' num2str(nr) '_ROI_' R.condname{cond}  '_' R.siden{side} '_' R.ipsicon])
                %                 load([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_LCMV_source_' R.condname{cond} 'nrep_' num2str(nr)])
                
                % Create data structure
                Xdata.fsample = R.pp.cont.full.fs;
                Xdata.label = vchansave(id).label;
                Xdata.trial{1} = vchansave(id).trial{1};
                Xdata.time{1} = vchansave(id).time; %{1}
                
                for ch = 1:2
                    X = Xdata.trial{1}(ch,:);
                    Xdata.trial{1}(ch,:) = (X-mean(X)); %std(X);
                end
                
                % Compute data transforms (Hilbert)
                [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Xdata,frq,stn_lb_frq,R.PA.bwid); % bandwidth 1 is good
                WinSize = 1.5*R.pp.cont.full.fs;
                [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,0.95);
                amp = cont2slidingwindow(amp,WinSize,round(0.95*WinSize));
                dphi_12 = cont2slidingwindow(dphi_12,WinSize,round(0.95*WinSize));
                
                SW_sampr = max(diff(PLV_tvec./R.pp.cont.full.fs));
                tseries = dphi_12; qstable = find(PLV>R.PA.PLVeps);
                mwid = R.PA.mwid;
                period = (mwid/frq)/SW_sampr;
                % Minimum number of cycles to consider sync
                [phi_dist amp_dist seg_ddt segL_ddt consecSegs H] = analysestablesegs(qstable,tseries,amp,period,mwid,1/SW_sampr);
                %                 [gc_dist] = analysestablesegs_granger(qstable,tseries,amp,period,mwid,R.pp.cont.full.fs);
                gc_dist_save{cond,nr} = NaN; %gc_dist;
                H_dist_save{cond,nr} = H;
                pA_pli_dist_save{cond,nr} = phi_dist;
                amp_pli_dist_save{cond,nr} = amp_dist;
                segL_pli_dist_save{cond,nr} = segL_ddt;
                timevec{cond,nr} = vchansave(id).time; %{1}
                %                 plot_example_phaseanalysis_trace(betaS,amp,phi,dphi_12_dt,seg_ddt1,0.005,PLI_tvec,PLI,consecSegs);
                %                 savefigure_v2([R.datapathr 'results\images\seganalysis\'],['example_seg_subject1_ON'],[],[],[]);
                close all
            end
        end
        save([R.datapathr R.subname{sub} '\ftdata\ROI_analy\' idd '_ROIvoxel_phaseamp_PLI_' R.ipsicon '_' R.siden{side}],...
            'pA_pli_dist_save','amp_pli_dist_save','segL_pli_dist_save','H_dist_save','timevec','gc_dist_save') %
    end
end
% ! shutdown /h