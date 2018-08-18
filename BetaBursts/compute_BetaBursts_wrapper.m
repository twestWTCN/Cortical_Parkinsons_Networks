function compute_BetaBursts_wrapper(R)
if nargin<2
    R = makeHeader_SubCort_Cort_Networks();
end
%%%
close all
for sub = 1:length(R.subname)
    for breg = 2:length(R.bregname)
        for side = 1:2
            AmpTime = [];
            for cond = 1:length(R.condname)
                load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                AmpTime{cond} =compute_BetaBursts(R,vc_clean);
            end
            
            fsamp = vc_clean.fsample;
            R.fsamp = fsamp;
            T = linspace(0,length([AmpTime{:}])/fsamp,length([AmpTime{:}]));
            eps = prctile([AmpTime{:}],75,2);
            
            % Plot Bursting
            figure(4)
            a = plot(T,[AmpTime{1} AmpTime{2}]); hold on
            a(1).Color = [0 0 0]; a(2).Color = [0.2 0.2 0.6]
            plot(T([1 end]),[eps eps],'--');
            ylim([0 35]);
            xlabel('Time'); ylabel('Beta Amplitude')
            set(gcf,'Position',[824.0000  669.0000  911.5000  178.5000])
            [segL_t_save segA_save AmpBin binSgEd] = compute_BetaBurstsStats(R,AmpTime,eps,1)
            close all
            % Save to data
            vc_clean.BB.segL_t = segL_t_save;
            vc_clean.BB.segA = segA_save;
            vc_clean.BB.AmpBin = AmpBin;
            vc_clean.BB.binSgEd = binSgEd;
            
            vc_clean.history = [vc_clean.history{:} {[mfilename '_' date]}];
            save([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
            
            disp([cond side breg sub])
        end
    end
end
end
% ! shutdown /h