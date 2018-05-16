function [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Odata,frq,stn_lb_frq,bwid,fsamp,LowAmpFix)
if nargin<6
    LowAmpFix = 0;
end
%         % bandpass
%         cfg = [];
%         cfg.bpfilter = 'yes';
%         cfg.bpfilttype    = 'fir';
%         cfg.bpfreq = [frq-bwid frq+bwid];
%         Xdata = ft_preprocessing(cfg,Xdata);
%         fsamp = Xdata.fsample;
Xdata.trial{1} = ft_preproc_bandpassfilter(Odata.trial{1}, fsamp, [frq-bwid frq+bwid], [], 'but', 'twopass', 'reduce');
% Xdata.trial{1} = Xdata.trial{1}(:,fsamp:end-fsamp);
% if size(Odata.trial{1},1)>size(Odata.trial{1},2)
%     Xdata.trial{1} = filterEEG(Odata.trial{1},fsamp,frq-bwid,frq+bwid,floor(8*(frq+bwid)))';
% else
%     Xdata.trial{1} = filterEEG(Odata.trial{1},fsamp,frq-bwid,frq+bwid,floor(8*(frq+bwid)))';
%     Xdata.trial{1} = Xdata.trial{1}';
% end
betaS.trial{1}= Xdata.trial{1} ;

amp(:,1) = abs(hilbert(Xdata.trial{1}(1,:)));
amp(:,2) = abs(hilbert(Xdata.trial{1}(2,:)));
amp = amp(fsamp:end-fsamp,:);
phi(:,1) = angle(hilbert(Xdata.trial{1}(1,:)));
phi(:,2) = angle(hilbert(Xdata.trial{1}(2,:)));
phi = phi(fsamp:end-fsamp,:);

dphi_12 = unwrap(phi(:,1)-phi(:,2));
% optional amp weighting
% ampw = unwrap((amp(:,1).*amp(:,2))./(max(amp(:,1))*max(amp(:,2))));  %%
% dphi_12 = angle(ampw.*exp(i.*(phi(:,1)-phi(:,2))));%%
% dphi_12 = (dphi_12-(1/sqrt(length(phi(:,1)))))./(1-(1/sqrt(length(phi(:,1))))); %%
dphi_12_dt = diff(dphi_12);

% 
% if LowAmpFix == 1
%     TAmp =  amp(:,1).* amp(:,2);
%     TAmpNorm = (TAmp./mean(TAmp));
%     TAmpNormNeg = TAmpNorm;
%     TAmpNormNeg(TAmpNorm>1) = 1;
%     dphi_12_dt = dphi_12_dt.*TAmpNormNeg(2:end);
% end
% lb_stn_bandpass
%         cfg = [];
%         cfg.bpfilter = 'yes';
%         cfg.bpfilttype    = 'fir';
%         cfg.bpfreq = [stn_lb_frq-bwid stn_lb_frq+bwid];
%         stn_lb_data = ft_preprocessing(cfg,Xdata);

Xdata.trial{1} = ft_preproc_bandpassfilter(Odata.trial{1}, fsamp, [stn_lb_frq-bwid stn_lb_frq+bwid], [], 'but', 'twopass', 'no');
% Xdata.trial{1} = Xdata.trial{1}(:,fsamp:end-fsamp);
% if size(Odata.trial{1},1)>size(Odata.trial{1},2)
%     Xdata.trial{1} = filterEEG(Odata.trial{1},fsamp,stn_lb_frq-bwid,stn_lb_frq+bwidfloor(6*(frq-bwid)))';
% else
%     Xdata.trial{1} = filterEEG(Odata.trial{1},fsamp,stn_lb_frq-bwid,stn_lb_frq+bwid,floor(6*(frq-bwid)))';
%     Xdata.trial{1} =  Xdata.trial{1}';
% end
betaS.time{1} = Odata.time{1} ;
betaS.trial{1}(3,:) = Xdata.trial{1}(2,:) ;
az = abs(hilbert(Xdata.trial{1}(2,:)));
amp(:,3) = az(fsamp:end-fsamp);
%         betaS = Xdata;