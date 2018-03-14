function [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Xdata,frq,stn_lb_frq,bwid)
        bwid = 1.5;
        % bandpass
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfilttype    = 'fir';
        cfg.bpfreq = [frq-bwid frq+bwid];
        Xdata = ft_preprocessing(cfg,Xdata);
        fsamp = Xdata.fsample;
        
        amp(:,1) = abs(hilbert(Xdata.trial{1}(1,:)));
        amp(:,2) = abs(hilbert(Xdata.trial{1}(2,:)));
        
        phi(:,1) = angle(hilbert(Xdata.trial{1}(1,:)));
        phi(:,2) = angle(hilbert(Xdata.trial{1}(2,:)));
        
        dphi_12 = unwrap(phi(:,1)-phi(:,2));
        % optional amp weighting
%         ampw = unwrap((amp(:,1).*amp(:,2))./(max(amp(:,1))*max(amp(:,2))));  %%
%         dphi_12 = angle(ampw.*exp(i.*(phi(:,1)-phi(:,2))));%%
%         dphi_12 = (dphi_12-(1/sqrt(length(phi(:,1)))))./(1-(1/sqrt(length(phi(:,1))))); %%
        dphi_12_dt = diff(dphi_12);
        
        % lb_stn_bandpass
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfilttype    = 'fir';
        cfg.bpfreq = [stn_lb_frq-bwid stn_lb_frq+bwid];
        stn_lb_data = ft_preprocessing(cfg,Xdata);
        amp(:,3) = abs(hilbert(stn_lb_data.trial{1}(1,:)));
    
        betaS = Xdata;