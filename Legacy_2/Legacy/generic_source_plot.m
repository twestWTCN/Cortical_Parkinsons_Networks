function generic_source_plot(datapathr, subname,sourcestruc,mri,funpar,appender,coordsys)
if nargin<7
    coordsys = 'ctf';
end
cfg            = [];
% cfg.downsample = 1;
cfg.parameter = funpar;
sourceInt  = ft_sourceinterpolate(cfg, sourcestruc , mri);

sourceInt.coordsys = coordsys;
cfg = [];
cfg.nonlinear     = 'no';
sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
%                         volumewrite_spm([datapathr subname{i} '\MRI\sourcespace\r' subname{i} '_' num2str(nr) '_DICS_sourcecoh_int_' ref_chan '.nii'], sourceIntNorm.coh, sourceIntNorm.transform, 'SPM12');
close all
scfg = [];
scfg.method         = 'surface';
scfg.funparameter   = 'avg.coh';
scfg.maskparameter  = scfg.funparameter;
% cfg.funcolorlim    = [0.0 1.2];
scfg.funcolormap    = 'jet';
%                         cfg.opacitylim     = [0. 0.05];
scfg.opacitymap     = 'auto';
scfg.projmethod     = 'nearest';
%                         cfg.downsample      = 2
scfg.surffile       = 'surface_white_both.mat';
%                         cfg.surfinflated   = 'surface_inflated_both_caret.mat';
viewang = {[0 0],[0 90],[90 0],[45 45]};
for ip = 1:4
    ft_sourceplot(scfg, sourceIntNorm);
    view (viewang{ip});
    a = gca; a.Children(3).FaceAlpha = 0;
end


surface = ft_read_headshape('surface_white_both.mat');
cfg = [];
cfg.parameter = 'coh';
cfg.intermethod = 'nearest';
cfg.projvec = 1;
cfg.projcomb = 'mean';
cfg.projweight = 1;
surfcoh = ft_sourceinterpolate(cfg, sourceIntNorm, surface)
[a1 xyz] = max(surfcoh.coh);
xyz = surfcoh.pos(xyz,:);


[X Y Z] = sphere(20); sprad= 25;
hold on
XX = X * sprad + xyz(1);
YY = Y * sprad + xyz(2);
ZZ = Z * sprad + xyz(3);
for ip = 1:4
    figure(ip); hold on
    surf(XX,YY,ZZ,'FaceAlpha', 0.2, 'FaceColor',[1 0 0],'edgecolor',[0.5 0.5 0.5])
%     savefig(gcf,[datapathr subname{1} '\images\sourcespace\r' subname{1} '_grandaverage_' appender '_DICs_OFF_view_' num2str(ip)])
    savefigure_v2([datapathr subname{1} '\images\sourcespace\'],['r' subname{1} '_grandaverage_' appender '_OFF_view_' num2str(ip)],ip,[],[]); 
end

end