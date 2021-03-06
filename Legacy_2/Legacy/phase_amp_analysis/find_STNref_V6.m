function find_STNref_V6(R)
% clear; close all
% R = makeHeader_SubCort_Cort_Networks();
for sub = 1:numel(R.subname)
    for cond = 1:2
        for band = [1 3]
            for side = 1:2
                load([R.datapathr R.subname{sub} '\ftdata\virtualV6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bandname{band}])
                
                cfg = [];
                cfg.length = 0.5;
                Xseg = ft_redefinetrial(cfg,vc_clean);
                
                % cfg = [];
                % cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
                % ft_databrowser(cfg, Xdata);
                
                cfg           = [];
                cfg.method    = 'mtmfft';
                cfg.taper     = 'hanning';
                cfg.output    = 'fourier';
                %                 cfg.foilim = [4 48];
                % cfg.pad       = 128;
                cfg.pad         =  'nextpow2';
                %                 cfg.tapsmofrq  = 1.5;
                freq         = ft_freqanalysis(cfg, Xseg);
                
                tpx = squeeze(mean(abs(freq.fourierspctrm(:,:,:)),1));
                % Normalisation
                for i = 1:numel(Xseg.label)
                    %                     tpx(i,:) = tpx(i,:)./max(tpx(i,freq.freq>4 & freq.freq<46));
                    tpx(i,:) = tpx(i,:)./sum(tpx(i,freq.freq>4 & freq.freq<45));
                end
                
                figure(1)
%                 subplot(1,2,cond)
                title(R.condname{cond})
                xlabel('Freq'); ylabel('Norm. Power')
                plot(freq.freq',tpx(1,:)','b')
                hold on
                plot(repmat(freq.freq,3,1)',tpx(2:end,:)','r')
                legend(vc_clean.label); grid on
                set(gcf,'Position',[1310         677         560         420])
                
                cfg           = [];
                cfg.method    = 'mtmfft';
                cfg.output    = 'fourier';
                cfg.tapsmofrq  = 3;
                cfg.pad         =  'nextpow2';
                freq         = ft_freqanalysis(cfg, Xseg);
                cfg           = [];
%                 cfg.method    = 'wpli';
                cfg.method    = 'coh';
%                 cfg.complex = 'absimag';
                coh           = ft_connectivityanalysis(cfg, freq);
                                icoh = abs(squeeze(coh.cohspctrm(2:end,1,:)));
%                 icoh = abs(squeeze(coh.wplispctrm(2:end,1,:)));
                
                figure(2)
%                 subplot(1,2,cond)
                title(R.condname{cond})
                xlabel('Freq'); ylabel('iCoh')
                plot(repmat(coh.freq,3,1)',icoh'); hold on
                ylim([0 0.8]); xlim([2 65]); legend({'STN01','STN02','STN03'}); grid on
                set(gcf,'Position',[1310         677         560         420])
                figure(2); shg
                gx = inputdlg(['Choose STN for ' R.bandname{band} ' band']);
                gx = str2num(gx{1});
                
                
                if gx<4
                    vc_clean.trial{1} = vc_clean.trial{1}([1 1+gx],:);
                    vc_clean.label = {vc_clean.label{[1 1+gx]}};
                % MaxCohs
                % Alpha
                [maxcoha fi] = max(icoh(gx,coh.freq>=R.bandef(1,1) & coh.freq<=R.bandef(1,2)));
                frqa = R.bandef(1,1)+(fi.*min(diff(freq.freq)));
                % Low Beta
                [maxcohb fi] = max(icoh(gx,coh.freq>=R.bandef(2,1) & coh.freq<=R.bandef(2,2)));
                frqb = R.bandef(2,1)+(fi.*min(diff(freq.freq)));
                % High Beta
                [maxcohc fi] = max(icoh(gx,coh.freq>=R.bandef(3,1) & coh.freq<=R.bandef(3,2)));
                frqc = R.bandef(3,1)+(fi.*min(diff(freq.freq)));
                maxcoh = [maxcoha maxcohb maxcohc];
                frqcoh = [frqa frqb frqc];
                
                %MaxPows
                [maxpowa fi] = max(tpx(1+gx,freq.freq> R.bandef(1,1) & freq.freq< R.bandef(1,2)));
                frqa = R.bandef(1,1)+(fi.*min(diff(freq.freq)));
                [maxpowb fi] = max(tpx(1+gx,freq.freq> R.bandef(2,1) & freq.freq< R.bandef(2,2)));
                frqb = R.bandef(2,1)+(fi.*min(diff(freq.freq)));
                [maxpowc fi] = max(tpx(1+gx,freq.freq> R.bandef(3,1) & freq.freq< R.bandef(3,2)));
                frqc = R.bandef(3,1)+(fi.*min(diff(freq.freq)));
                maxpow = [maxpowa maxpowb maxcohc];
                frqpow = [frqa frqb frqc];

                vc_clean.specanaly.icoh = icoh(gx,:);
                vc_clean.specanaly.normpow = tpx([1 1+gx],:);
                vc_clean.specanaly.frq = freq.freq;
                vc_clean.specanaly.cohstats.frqcoh = frqcoh;
                vc_clean.specanaly.cohstats.maxcoh = maxcoh;
                vc_clean.specanaly.powstats.frqpow = frqpow;
                vc_clean.specanaly.powstats.maxpow = maxpow;
                vc_clean.specanaly.flag = 0;
                else
                    vc_clean.specanaly.flag = 1;
                end  
                vc_clean.history = [vc_clean.history{:} {[mfilename '_' date]}];
                save([R.datapathr R.subname{sub} '\ftdata\virtualV6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bandname{band}],'vc_clean')
                close all
            end
        end
    end
end
