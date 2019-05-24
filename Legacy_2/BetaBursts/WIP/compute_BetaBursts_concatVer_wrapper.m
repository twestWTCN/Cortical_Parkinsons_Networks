function compute_BetaBursts_wrapper(R)
%% Main function to compute subject level BB analysis + phase syncs
% TO do:
% 1) Implement the rate normalization step to divide by recording time.
% 2) Is there some relationship between the change in the coherence and the
% the change in power?
% 3) For group level- inspect individually and classify by eye.
% 4) Think of the amplitude by duration plots as a cone - how can you fit
% to data and examine the difference? You need to fit an arc sector
if nargin<2
    R = makeHeader_SubCort_Cort_Networks();
end
%%%
close all
for sub = 1:length(R.subname)
    for breg = 2:length(R.bregname)
        BB.AmpTime = [];
        for cond = 1:length(R.condname)
            for side = 1:2
                load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                data{side} = vc_clean;
            end
            % Concatanate Hemispheres
            cfg = [];
            catdat = ft_appenddata(cfg, data{1}, data{2});
            catdat.trial{1} = [catdat.trial{1}(1:2,:) catdat.trial{1}(3:4,:)];
            catdat.label = catdat.label(1:2);
            catdat.fsample = data{1}.fsample;
            catdat.time = {linspace(0,size(catdat.trial{1},2)/catdat.fsample,size(catdat.trial{1},2))};
            [BB.A{cond},BB.DTvec{cond},BB.PLV{cond},BB.SWTvec{cond},BB.RP{cond}] = compute_BetaBursts(R,catdat);
        end
        BB.fsamp = vc_clean.fsample;
        BB.fsamp_sw = 1/(BB.SWTvec{1}(5)-BB.SWTvec{1}(4));
        R.fsamp = BB.fsamp;
        BB.T = linspace(0,length([BB.A{:}])/BB.fsamp,length([BB.A{:}]));
        BB.TSw = linspace(0,length([BB.PLV{:}])/BB.fsamp_sw,length([BB.PLV{:}]));
        
        % ThresholdAmplitude
        BB.epsAmp = prctile([BB.A{:}],65,2);
        BB.epsPLV = prctile([BB.PLV{:}],65,2);
        
        % Plot Bursting
        figure(1)
        plotExampleBurstPLV(R,BB)
        set(gcf,'Position',[679.5000  621.0000  913.5000  358.5000])
        
        %% Now Do Stats and Plot Single Example (if desired)
        % Amplitude Only
        BB.range.Amp = 0:2:35;
        BB.range.segDur = 10:100:2000;
        F = figure(2);
        BB = compute_BetaBurstsStats(R,BB,F,1);
        set(gcf,'Position',[680 347 913 632])
        
        % Amplitude and Sync
        BB.range.PLV = 0:0.1:1;
        F(1) = figure(3);
        BB = computeBetaBurstPLVStats(R,BB,F,1);
        set(F(1),'Position',[680.0000   83.5000  907.5000  894.5000])
        
        % Relative Phase
        F(1) = figure(4);
        set(F(1),'Position',[680.0000   9  907.5000  969.0000])
        BB.range.RP = -pi:pi/6:pi;
        BB = computeBetaBurstRPStats(R,BB,F,1);
        
        % Save to data structure
        vc_clean.BB = BB;
        vc_clean.history = [vc_clean.history{:} {[mfilename '_' date]}];
        save([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
        close all
        disp([cond side breg sub])
    end
end
% ! shutdown /h