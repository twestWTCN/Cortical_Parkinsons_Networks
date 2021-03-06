function compute_phase_amp_analysis_v8(R)
if nargin<2
    R = makeHeader_SubCort_Cort_Networks();
end
%%%
close all
for sub = 1:length(R.subname)
    for breg = 2:length(R.bregname)
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
                    cohfrq = vc_clean.specanaly.plvstats.frqplv; % Alpha and high beta frequencies
                    powfrq = vc_clean.specanaly.powstats.frqpow; % STN frequencies
                    
                    % Create data structure
                    Xdata = vc_clean;
                    Xdata = rmfield(Xdata,'specanaly');
                    [phi_dist amp_dist segL_ddt ampT maxPLV PLVtime Ampeps surr OVL] = compute_dynPhaseLocking(R,Xdata,band,cohfrq,powfrq);
                else
                    phi_dist = NaN;
                    amp_dist =NaN;
                    segL_ddt = NaN;
                end
                % Save to data
                vc_clean.PA.maxplv = maxPLV;
                vc_clean.PA.pureAmp = ampT;
                vc_clean.PA.pA_pli_dist_save = phi_dist;
                vc_clean.PA.amp_pli_dist_save = amp_dist;
                vc_clean.PA.segL_pli_dist_save = segL_ddt;
                vc_clean.PA.timevec = PLVtime;
                vc_clean.PA.surrAmp = Ampeps;
                vc_clean.PA.OVL = OVL;
%                 vc_clean.PA.surr = surr;
                vc_clean.history = [vc_clean.history{:} {[mfilename '_' date]}];
                save([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                
                disp([cond side breg sub])
            end
        end
    end
end
% ! shutdown /h