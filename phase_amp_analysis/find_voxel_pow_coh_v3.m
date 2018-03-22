function [idbank frqbank stn_lb_frqbank] = find_voxel_pow_coh_v3(R,subname,band,nr)
for  cond = 1:2 % USE just OFF to identify ROI for now
    [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{subname},cond-1);
    for side = 1:2
        for nr = 1:nrep
            vchansave = [];
            load([R.datapathr R.subname{subname} '\ftdata\virtual_sources_' num2str(nr) '_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon]); % '_' R.bandnames{band}])
            for x = 1:numel(vchansave)
                Xdata(x).fsample = R.pp.cont.full.fs;
                Xdata(x).label = vchansave(x).label;
                Xdata(x).trial = vchansave(x).trial; %{1};
                Xdata(x).time = {vchansave(x).time};
                Xdata(x).trial{1} = Xdata(x).trial{1}(:,2*R.pp.cont.full.fs:(end-2*R.pp.cont.full.fs));
                Xdata(x).time{1} = linspace(0,size( Xdata(x).trial{1},2)/R.pp.cont.full.fs, size( Xdata(x).trial{1},2));  %Xdata(x).time{1}(:,2*fsamp:(end-2*fsamp));
                Xdata(x).cfg = [];
                cfg = [];
                cfg.resamplefs = R.pp.cont.thin.fs;
                X = ft_resampledata(cfg,Xdata(x));
                Xdata(x) = X;
            end
            N = numel(vchansave);
%             progressStepSize = 1;
            %             if ~is_in_parallel; parpool; end
%             ppm = ParforProgMon('Power Estimation: ', N, progressStepSize, 800, 300);
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
                    tpx(i,:) = tpx(i,:)./sum(tpx(i,freq.freq>4 & freq.freq<45));
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
                cfg.tapsmofrq  = 2;
                cfg.pad         =  'nextpow2';
                freq         = ft_freqanalysis(cfg, Xseg);
                cfg           = [];
                cfg.method    = 'coh';
                %                         cfg.complex = 'imag';
                coh           = ft_connectivityanalysis(cfg, freq);
                
                icoh = abs(squeeze(coh.cohspctrm(2,1,:)));
                cohsave_vox(:,x) = icoh;
                % MaxCohs
                % Alpha
                
                [maxcoha fi] = max(icoh(coh.freq>=R.bandef(1,1) & coh.freq<=R.bandef(1,2)));
                frqa = R.bandef(1,1)+(fi.*min(diff(freq.freq)));
                % Low Beta
                [maxcohb fi] = max(icoh(coh.freq>=R.bandef(2,1) & coh.freq<=R.bandef(2,2)));
                frqb = R.bandef(2,1)+(fi.*min(diff(freq.freq)));
                % High Beta
                [maxcohc fi] = max(icoh(coh.freq>=R.bandef(3,1) & coh.freq<=R.bandef(3,2)));
                frqc = R.bandef(3,1)+(fi.*min(diff(freq.freq)));
                maxcoh(:,x) = [maxcoha maxcohb maxcohc];
                
                frq(:,x) = [frqa frqb frqc];
                %MaxPows
                stn_pow = squeeze(mean(abs(freq.fourierspctrm(:,2,:)),1));
                [dum fi] = max(stn_pow(freq.freq> R.bandef(2,1) & freq.freq< R.bandef(2,2)));
                stn_lb_frq(x) =  R.bandef(2,1)+(fi.*min(diff(freq.freq)));
                frqsave_vox(:,x) = coh.freq;
%                 ppm.increment();
                %                 close all
            end
%             ppm.delete()
            powsave{nr,side,cond} = powsave_vox;
            cohsave{nr,side,cond} = cohsave_vox;
            frqsave{nr,side,cond} = frqsave_vox;
            maxcoh(maxcoh>0.8) = NaN;
            [dum id] = max(maxcoh(band,:)); % THIS NEEDS TO BE CHANGED FOR THE BAND OF INTEREST!
            idbank(nr,side,cond) = id;
            frqbank(:,nr,side,cond) = frq(:,id);
            stn_lb_frqbank(nr,side,cond) = stn_lb_frq(id);
            clear maxcoh frq
            
        end
    end
end
mkdir([R.datapathr R.subname{subname} '\ftdata\ROI_analy\'])
save([R.datapathr R.subname{subname} '\ftdata\ROI_analy\ROIvoxel_power_' R.ipsicon '_' R.bandname{band}],'powsave','frqsave')
save([R.datapathr R.subname{subname} '\ftdata\ROI_analy\ROIvoxel_coh_' R.ipsicon '_' R.bandname{band}],'cohsave','frqsave')
save([R.datapathr R.subname{subname} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon '_' R.bandname{band}],'frqbank','idbank','stn_lb_frqbank')