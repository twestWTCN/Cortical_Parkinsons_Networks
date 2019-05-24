function compute_phase_amp_analysis_050118(R)
%%%
% This script computes the instantaneous phase sync analyses. We partition
% the recording into stable segments using threshold, can compute Granger
% etc within frames. DFA is also computed PS/AE versions. Barplots for DFA
% at end of script.
%%%

load('ref_ci.mat')
load([R.datapathr 'subject_hbWPLI075'])
subscreen = squeeze(sum(subject_hbcohscreen>R.PA.WPLIscreen)==2);
for sub = 1:numel(R.subname)
    if exist([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon '.mat']) ==0
        [idbank frqbank stn_lb_frqbank] = find_voxel_pow_coh_050118(R.datapathr,R.subname{sub},R.condname,R.sidenm,R.pp.cont.full.fs,R.ipsicon)
    else
        load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon])
    end
    for side = 1:2
        for cond = 1:2
            [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
            for nr = nrep
                id = idbank(1,side,2);
                frq = frqbank(1,side,2);
                stn_lb_frq = frqbank(1,side,2);
                load([R.datapathr R.subname{sub} '\ftdata\virtual_sources_' num2str(nr) '_ROI_' R.condname{cond}  '_' R.siden{side} '_' R.ipsicon])
%                 load([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_LCMV_source_' R.condname{cond} 'nrep_' num2str(nr)])
                
                % Create data structure
                Xdata.fsample = R.pp.cont.full.fs;
                Xdata.label = vchansave(id).label;
                Xdata.trial{1} = vchansave(id).trial{1};
                Xdata.time{1} = vchansave(id).time; %{1}
                % Truncate
                Xdata.trial{1} = Xdata.trial{1}(:,2*R.pp.cont.full.fs:(end-2*R.pp.cont.full.fs));
                Xdata.time{1} = Xdata.time{1}(:,2*R.pp.cont.full.fs:(end-2*R.pp.cont.full.fs)); %
                
                for ch = 1:2
                    X = Xdata.trial{1}(ch,:);
                    Xdata.trial{1}(ch,:) = (X-mean(X))/std(X);
                end
               
%                 [cxx fxx] = mscohere(Xdata.trial{1}(1,:),Xdata.trial{1}(2,:),hanning(R.pp.cont.full.fs),[],R.pp.cont.full.fs,R.pp.cont.full.fs);
                if 1==1 %subscreen(side,sub) % look for significant coherences
                    
                % Compute data transforms (Hilbert)
                [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Xdata,frq,stn_lb_frq,0.5); % bandwidth 1 is good
                
                % Variance
                var_amp_save(:,cond,nr,side,sub) = std(amp);
                var_dphidt_save(:,cond,nr,side,sub) = std(dphi_12_dt);
                
                % DFA business
                DFAP = [R.pp.cont.full.fs,(1/14)*12,(length(dphi_12)/8)/R.pp.cont.full.fs,50,0];
                [dum1,dum2,evi,alpha] = peb_dfa_TW_2018(dphi_12_dt',DFAP,-4,0);
                dfa_amp_save(:,1,cond,nr,side,sub) = [alpha evi];
                [dum1,dum2,evi,alpha] = peb_dfa_TW_2018(amp(:,1)',DFAP,-4,0);
                dfa_amp_save(:,2,cond,nr,side,sub) = [alpha evi];
                [dum1,dum2,evi,alpha] = peb_dfa_TW_2018(amp(:,2)',DFAP,-4,0);
                dfa_amp_save(:,3,cond,nr,side,sub) = [alpha evi];
                [dum1,dum2,evi,alpha] = peb_dfa_TW_2018(amp(:,3)',DFAP,-4,0);
                dfa_amp_save(:,4,cond,nr,side,sub) = [alpha evi];                
                
                mwid = 6; % Minimum number of cycles to consider sync
                period = (mwid/frq)*R.pp.cont.full.fs;
                [PLI PLI_tvec] = slidingwindowPLI(period,phi,Xdata.time{1});
                
                %%% Find segments
%                 dphi_12_dt_sm = smooth(dphi_12_dt,period/6)';
                qstable = find(abs(dphi_12_dt')<0.005); % The points that are below threshold
                tseries = wrapToPi(dphi_12(2:end));
%                 tseries(:,2) = phi(:,2);
%                 tseries(:,1) = phi(:,1);
                [phi_dist amp_dist seg_ddt1 segL_ddt consecSegs H] = analysestablesegs(qstable,tseries,amp,period,mwid,R.pp.cont.full.fs);
                pA_dist_save{cond,nr} = phi_dist;
                amp_dist_save{cond,nr} = amp_dist;
                segL_dist_save{cond,nr} = segL_ddt;
                H_dist_save{cond,nr} = H;
                % granger subanalysis
                tseries = Xdata.trial{1}';
%                 [gc_dist] = analysestablesegs_granger(qstable,tseries,amp,period,mwid,R.pp.cont.full.fs);
                gc_dist_save{cond,nr} = NaN; %gc_dist;
%                 
                %PLI segments
                qstable = find(abs(PLI)>wpli_ci); % The points that are below threshold
                tseries = PLI';
                [phi_dist amp_dist seg_ddt segL_ddt consecSegs] = analysestablesegs(qstable,tseries,amp,period/mwid,mwid,R.pp.cont.full.fs);
                pA_pli_dist_save{cond,nr} = phi_dist;
                amp_pli_dist_save{cond,nr} = amp_dist;
                segL_pli_dist_save{cond,nr} = segL_ddt;
                timevec{cond,nr} = vchansave(id).time; %{1}
%                 plot_example_phaseanalysis_trace(betaS,amp,phi,dphi_12_dt,seg_ddt1,0.005,PLI_tvec,PLI,consecSegs); 
%                 savefigure_v2([R.datapathr 'results\images\seganalysis\'],['example_seg_subject1_ON'],[],[],[]);
                close all
                else
                  fprintf('Coherence is too weak at %.3f maximum in high beta',max(cxx(fxx>=20 & fxx<36)))
                    pA_dist_save{cond,nr} = NaN;
                    amp_dist_save{cond,nr} = NaN;
                    segL_dist_save{cond,nr} = NaN;
                    H_dist_save{cond,nr} = NaN;
                    gc_dist_save{cond,nr} = NaN;
                    pA_pli_dist_save{cond,nr} = NaN;
                    amp_pli_dist_save{cond,nr} = NaN;
                    segL_pli_dist_save{cond,nr} = NaN;
                    timevec{cond,nr} = vchansave(id).time; %{1};
                    
                end
                    
            end
        end
        save([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_phaseamp_' R.ipsicon '_' R.siden{side}],'pA_dist_save','amp_dist_save','segL_dist_save',...
            'pA_pli_dist_save','amp_pli_dist_save','segL_pli_dist_save','H_dist_save','timevec','var_amp_save','var_dphidt_save','dfa_amp_save','gc_dist_save') %
    end
end
save([R.datapathr '\results\dfa_amp_save'],'dfa_amp_save')

load([R.datapathr '\results\dfa_amp_save'],'dfa_amp_save')
barplots_DFA
 savefigure_v2([R.datapathr 'results\images\DFA\'],['DFA_boxplots'],[],[],[]);
% Result on the DFA for amplitude envelopes in cortex
% ! shutdown /h