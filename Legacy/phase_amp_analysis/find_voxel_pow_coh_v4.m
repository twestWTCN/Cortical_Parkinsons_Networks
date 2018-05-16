function [] = find_STNref(R)
for sub = 1:numel(R.subname)
    for cond = 1:2
        for band = [1 3]
            for side = 1:2
                load([R.datapathr R.subname{sub} '\ftdata\virtualV6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bandname{band}],'vc_clean')
                %                 cfg = [];
                %                 cfg.length = 0.5;
                %                 Xseg = ft_redefinetrial(cfg,Xdata(x));
                Xseg = Xdata(x);
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
                subplot(1,2,cond)
                title(R.condname{cond})
                xlabel('Freq'); ylabel('Norm. Power')
                plot(freq.freq',tpx(1,:)','b')
                hold on
                plot(repmat(freq.freq,3,1)',tpx(2:end,:)','r')
                legend({'CTX','STN'}); grid on
                
                %                 cfg = [];
                %                 cfg.length = 1;
                %                 Xseg = ft_redefinetrial(cfg,Xdata(x));
                
                cfg           = [];
                cfg.method    = 'mtmfft';
                %                 cfg.taper     = 'hanning';
                cfg.output    = 'fourier';
                %                 cfg.foilim = [4 48];
                % cfg.pad       = 128;
                cfg.tapsmofrq  = 3;
                cfg.pad         =  'nextpow2';
                freq         = ft_freqanalysis(cfg, Xseg);
                cfg           = [];
                cfg.method    = 'wpli';
%                 cfg.method    = 'coh';
%                                         cfg.complex = 'absimag';
                coh           = ft_connectivityanalysis(cfg, freq);
%                 icoh = abs(squeeze(coh.cohspctrm(2:end,1,:)));
                icoh = abs(squeeze(coh.wplispctrm(2:end,1,:)));
                
                figure(2)
                subplot(1,2,cond)
                title(R.condname{cond})
                xlabel('Freq'); ylabel('iCoh')
                plot(repmat(coh.freq,3,1)',icoh'); hold on
                ylim([0 0.8]); xlim([2 65]); legend({'STN01','STN02','STN03'}); grid on
                
                figure(2); shg
                %                 gx = inputdlg(['Choose STN for ' R.bandname{band} ' band']);
                %                 gx = str2num(gx{1});
                gx = 2;
                powsave_vox(:,:,x) = tpx([1 gx],:);
                cohsave_vox(:,x) = icoh(gx,:);
                
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
                    maxcoh(:,x) = [maxcoha maxcohb maxcohc];
                    
                    frq(:,x) = [frqa frqb frqc];
                    %MaxPows
                    stn_pow = squeeze(mean(abs(freq.fourierspctrm(:,1+gx,:)),1));
                    [dum fi] = max(stn_pow(freq.freq> R.bandef(2,1) & freq.freq< R.bandef(2,2)));
                    stn_lb_frq(x) =  R.bandef(2,1)+(fi.*min(diff(freq.freq)));
                    frqsave_vox(:,x) = coh.freq;
                    %                 ppm.increment();
                close all
            end
            
            %             ppm.delete()
            refsave{nr,side,cond} = gx;
            powsave{nr,side,cond} = powsave_vox;
            cohsave{nr,side,cond} = cohsave_vox;
            frqsave{nr,side,cond} = frqsave_vox;
            maxcoh(maxcoh>0.9) = NaN;
            [dum id] = max(maxcoh(band,:));
            idbank(nr,side,cond) = id;
            frqbank(:,nr,side,cond) = frq(:,id);
            stn_lb_frqbank(nr,side,cond) = stn_lb_frq(id);
            clear maxcoh frq
        end
    end
end
mkdir([R.datapathr R.subname{sub} '\ftdata\ROI_analy\'])
save([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_power_' R.ipsicon '_' R.bandname{band}],'powsave','frqsave')
save([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_coh_' R.ipsicon '_' R.bandname{band}],'cohsave','frqsave')
save([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon '_' R.bandname{band}],'frqbank','idbank','stn_lb_frqbank','refsave')