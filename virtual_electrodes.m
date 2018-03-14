cfg = [];
cfg.datafile = [datapathr subname{i} '_tmp\' datafilen];
cfg.lpfilter = 'yes';
cfg.hpfilter = 'yes';
cfg.channel         = {'MEG'};
cfg.lpfreq = 48;
cfg.hpfreq = 2;
bbdata = ft_preprocessing(cfg);

[x roi_ind] = min(abs(sum((source.pos - [1.2 2.1 -3.2]),2)));
roi_ind = 4432;
bf_leadfield = source.avg.filter{roi_ind}
chansel = ft_channelselection('MEG',bbdata.label)
chansel = match_str(bbdata.label,chansel)

vc_data = [];
vc_data.label = {'vc_x','vc_y','vc_z'};
vc_data.time = [bbdata.time{:}];
for i = 1:length(bbdata.trial)
    vc_data.trial{i} = bf_leadfield * bbdata.trial{i}(chansel,:);
end
vc_data.trial =  {[vc_data.trial{:}]};
vc_data.time = linspace(0,max(bbdata.time{1})*numel(vc_data.time),numel(bbdata.time)*numel(bbdata.time{1}))
figure
for i = 1; %:3
    plot(vc_data.time,abs(vc_data.trial{1}(i,:))+((i-1)*0.01))
    hold on
end
% ROI 84 102 186