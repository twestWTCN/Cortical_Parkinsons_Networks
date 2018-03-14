function [idbank frqbank stn_lb_frqbank] = find_voxel_pow_coh_050118(datapathr,subname,condname,sidenm,fsamp,ipsicon)
for  cond = 1:2 % USE just OFF to identify ROI for now
    [datafileN,pp_mark,nrep,senscheck] = data_fileguide(subname,cond-1);
    for side = 1:2
        for nr = 1:nrep
            vchansave = [];
            load([datapathr subname '\ftdata\virtual_sources_' num2str(nr) '_ROI_' condname{cond} '_' sidenm{side} '_' ipsicon])
            for x = 1:numel(vchansave)
                Xdata(x).fsample = fsamp;
                Xdata(x).label = vchansave(x).label;
                Xdata(x).trial = vchansave(x).trial; %{1};
                Xdata(x).time = {vchansave(x).time};
                Xdata(x).trial{1} = Xdata(x).trial{1}(:,2*fsamp:(end-2*fsamp));
                Xdata(x).time{1} = linspace(0,size( Xdata(x).trial{1},2)/fsamp, size( Xdata(x).trial{1},2));  %Xdata(x).time{1}(:,2*fsamp:(end-2*fsamp));
                Xdata(x).cfg = [];
                cfg = [];
                cfg.resamplefs = 256;
                X = ft_resampledata(cfg,Xdata(x));
                Xdata(x) = X;
            end
            N = numel(vchansave);
            progressStepSize = 1;
            %             if ~is_in_parallel; parpool; end
            ppm = ParforProgMon('Power Estimation: ', N, progressStepSize, 800, 300);
            parfor x = 1:numel(vchansave)
                cfg = [];
                cfg.length = 1.5;
                Xseg = ft_redefinetrial(cfg,Xdata(x));
                
                % cfg = [];
                % cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
                % ft_databrowser(cfg, Xdata);
                
                cfg           = [];
                cfg.method    = 'mtmfft';
                cfg.taper     = 'dpss';
                cfg.output    = 'fourier';
                %                 cfg.foilim = [4 48];
                % cfg.pad       = 128;
                cfg.pad         =  'nextpow2';
                cfg.tapsmofrq  = 2;
                freq         = ft_freqanalysis(cfg, Xseg);
                
                tpx = squeeze(mean(abs(freq.fourierspctrm(:,:,:)),1));
                % Normalisation
                for i = 1:2
                    %                     tpx(i,:) = tpx(i,:)./max(tpx(i,freq.freq>4 & freq.freq<46));
                    tpx(i,:) = tpx(i,:)./sum(tpx(i,freq.freq>6 & freq.freq<45));
                end
                powsave_vox(:,:,x) = tpx;
                %                                 figure(1+10*cond)
                %                                 plot(repmat(freq.freq,2,1)',tpx')
                %                                 hold on
                %                                 freqstore(:,:,x,cond) = normaliseV(squeeze(mean(abs(freq.fourierspctrm(:,:,:)),1)))';
                %                                 figure(2+10*cond)
                cfg           = [];
                cfg.method    = 'mtmfft';
                cfg.taper     = 'dpss';
                cfg.output    = 'fourier';
                %                 cfg.foilim = [4 48];
                % cfg.pad       = 128;
                cfg.tapsmofrq  = 2.5;
                cfg.pad         =  'nextpow2';
                freq         = ft_freqanalysis(cfg, Xseg);
                cfg           = [];
                cfg.method    = 'coh';
                %                         cfg.complex = 'imag';
                coh           = ft_connectivityanalysis(cfg, freq);
                
                icoh = abs(squeeze(coh.cohspctrm(2,1,:)));
                cohsave_vox(:,x) = icoh;
                
                [maxcoh(x) fi] = max(icoh(coh.freq>24 & coh.freq<34));
                frq(x) = 24+fi;
                
                %                 [maxcoh(x) fi] = max(icoh(coh.freq>13 & coh.freq<24));
                %                 frq(x) = 13+fi;
                
                stn_pow = squeeze(mean(abs(freq.fourierspctrm(:,2,:)),1));
                [dum fi] = max(icoh(coh.freq>13 & coh.freq<24));
                stn_lb_frq(x) = 13+fi;
                frqsave_vox(:,x) = coh.freq;
                ppm.increment();
                %                 close all
            end
            ppm.delete()
            powsave{nr,side,cond} = powsave_vox;
            cohsave{nr,side,cond} = cohsave_vox;
            frqsave{nr,side,cond} = frqsave_vox;
            maxcoh(maxcoh>0.8) = NaN;
            [dum id] = max(maxcoh);
            idbank(nr,side,cond) = id;
            frqbank(nr,side,cond) = frq(id);
            stn_lb_frqbank(nr,side,cond) = stn_lb_frq(id);
            clear maxcoh frq
            
        end
    end
end
mkdir([datapathr subname '\ftdata\ROI_analy\'])
save([datapathr subname '\ftdata\ROI_analy\ROIvoxel_power_' ipsicon],'powsave','frqsave')
save([datapathr subname '\ftdata\ROI_analy\ROIvoxel_coh_' ipsicon],'cohsave','frqsave')
save([datapathr subname '\ftdata\ROI_analy\ROIvoxel_bank_' ipsicon],'frqbank','idbank','stn_lb_frqbank')