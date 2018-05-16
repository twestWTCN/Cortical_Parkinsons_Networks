%function plot_timevarying_source(datapathr,subname)

% Plot time varying granger
nr = 1; subname = {'LM'}; i = 1; ref_chan = 'STN_R23';
load([datapathr subname{i} '\ftdata\source_granger_' ref_chan '_' num2str(nr)])
load([datapathr subname{i} '\ftdata\source_coh_' ref_chan '_' num2str(nr)])
load([datapathr subname{i} '\ftdata\r' subname{i} 'rs'],'mri')
load([datapathr subname{i} '\ftdata\source_loc_time_' num2str(nr)])
load([datapathr subname{i} '\ftdata\r' subname{i} '_LCMV_source'])
%                 cfg = [];
%                 cfg.parameter = 'avg.pow';
%                 cfg.filename = [datapathr subname{i} '_tmp\MRI\source_data']
%                 ft_sourcewrite(cfg, source)
load([datapathr subname{i} '\ftdata\rs_transform_vox2ctf'])

T = transform_vox2ctf; %/transform_vox2spm;%;
Tmri = ft_transform_geometry(T,mri);
Tmri = ft_convert_units(Tmri,'cm');


gSource = nan([size(gASym,2) size(source.inside,1)]);
for x = 1:numel(edgInd)
    gSource(:,edgInd(x)) = gASym(x,:)';
end
gSource = reshape(gSource,[size(gSource,1) source.dim]);
cSource = nan([size(cASym,2) size(source.inside,1)]);
for x = 1:numel(edgInd)
    cSource(:,edgInd(x)) = cASym(x,:)';
end
cSource = reshape(cSource,[size(cSource,1) source.dim]);

save([datapathr subname{i} '\ftdata\cSource4D'],'cSource')
save([datapathr subname{i} '\ftdata\gSource4D'],'gSource')

%%%%%%%%%%%%%%%
% gSource = cSource;
load([datapathr subname{i} '\ftdata\rs_transform_vox2ctf'])
close all
for t = 1:size(gSource,1)
    gtSource = squeeze(gSource(t,:));
    source.avg.granger = gtSource;
    cfg            = [];
    % cfg.downsample = 1;
    cfg.parameter = 'granger';
    sourceInt  = ft_sourceinterpolate(cfg, source , Tmri);
    
    %                 Z-score power
    %                     sourceInt.coh = reshape(sourceInt.coh,size(sourceInt.anatomy));
    %                     incoh = sourceInt.coh.*sourceInt.inside;
    %                     incoh(incoh==0) = NaN;
    %                     incoh = (incoh - nanmean(incoh(:)))/nanstd(incoh(:));
    %                     %                                     inpow(isnan(inpow)) = 0;
    %                     sourceInt.coh = incoh;
    %                     sourceInt.coh = reshape(sourceInt.coh,[],1);
    
%     cfg = [];
%     cfg.method        = 'slice';
%     cfg.funparameter  = 'granger';
%     cfg.maskparameter = cfg.funparameter;
%     cfg.funcolorlim   = [-0.2 0.2];
%     % cfg.opacitylim    = [4.0 6.2];
%     cfg.opacitymap    = 'rampup';
%     figure
%     ft_sourceplot(cfg, sourceInt);
    
    sourceInt.coordsys = 'ctf';
    cfg = [];
    cfg.nonlinear     = 'no';
    cfg.spmversion = 'spm12'
    sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
    volumewrite_spm([datapathr subname{i} '\MRI\sourcespace\anim\r' subname{i} '_' num2str(nr) '_DICS_sourcegranger_int_' ref_chan '_' sprintf('%03d',t) '.nii'], sourceIntNorm.granger, sourceIntNorm.transform, 'SPM12');
    
    cfg = [];
    cfg.method         = 'surface';
    cfg.funparameter   = 'granger';
    cfg.maskparameter  = cfg.funparameter;
    cfg.funcolorlim    = [0 0.1];
    cfg.funcolormap    = 'jet';
    cfg.opacitylim     = [0 0.1];
    cfg.opacitymap     = 'auto';
    cfg.projmethod     = 'nearest';
    cfg.surffile       = 'surface_white_both.mat';
    cfg.surfinflated   = 'surface_inflated_both_caret.mat';
    H = figure;
    ft_sourceplot(cfg, sourceIntNorm);
    view ([0 90])
    drawnow;
    M(t) = getframe;
    savefig(gcf,[datapathr subname{i} '\images\sourcespace\anim\r' subname{i} '_rep_' num2str(nr) '_surf_avg_DICS_granger_' ref_chan '_' sprintf('%03d',t)])
    close all
end

figure
axes('Position',[0 0 1 1])
movie(M,1,10)

% Duration, location, assymetry
% Threshold Granger - and compute time to resolve
% Threshold. MPI coding
