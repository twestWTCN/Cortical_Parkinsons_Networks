function plot_timevarying_source(

% Plot time varying granger
nr = 1; subname = {'LM'}; i = 1;
datapathr = 'C:\Users\Tim\Documents\Work\Cortical_Networks\Data\'; ref_chan = 'STN_R12';


load('C:\Users\Tim\Documents\Work\Cortical_Networks\tmp_files\gASym.mat')
load([datapathr subname{i} '\ftdata\r' subname{i} 'rs'],'mri')
%                 cfg = [];
%                 cfg.parameter = 'avg.pow';
%                 cfg.filename = [datapathr subname{i} '_tmp\MRI\source_data']
%                 ft_sourcewrite(cfg, source)
load([datapathr subname{i} '\ftdata\rs_transform_vox2ctf'])

T = transform_vox2ctf; %/transform_vox2spm;%;
Tmri = ft_transform_geometry(T,mri);
Tmri = ft_convert_units(Tmri,'cm');

load('C:\Users\Tim\Documents\Work\Cortical_Networks\tmp_files\source_LM.mat')
gSource = NaN([size(gASym,1) size(source.inside,1)])
for x = 1:numel(edgInd)
    gSource(:,edgInd(x)) = gASym(:,x);
end
gSource = reshape(gSource,[size(gSource,1) source.dim]);
load([datapathr subname{i} '\ftdata\rs_transform_vox2ctf'])
close all
for t = 75:size(gSource,1)
    gtSource = squeeze(gSource(t,:,:,:));
    source.avg.granger = reshape(gtSource,[],1);
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
    sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
    volumewrite_spm([datapathr subname{i} '\MRI\sourcespace\anim\r' subname{i} '_' num2str(nr) '_DICS_sourcegranger_int_' ref_chan '_' sprintf('%03d',t) '.nii'], sourceIntNorm.granger, sourceIntNorm.transform, 'SPM12');
    
    cfg = [];
    cfg.method         = 'surface';
    cfg.funparameter   = 'granger';
    cfg.maskparameter  = cfg.funparameter;
    cfg.funcolorlim    = [-0.2 0.2];
    cfg.funcolormap    = 'jet';
    cfg.opacitylim     = [-0.2 0.2];
    cfg.opacitymap     = 'auto';
    cfg.projmethod     = 'nearest';
    cfg.surffile       = 'surface_white_both.mat';
    cfg.surfinflated   = 'surface_inflated_both_caret.mat';
    H = figure;
    ft_sourceplot(cfg, sourceIntNorm);
    view ([0 90])
    drawnow;
    M(t) = getframe;
    savefig(gcf,[datapathr subname{i} '\images\sourcespace\r' subname{i} '_rep_' num2str(nr) '_surf_avg_DICS_granger_' ref_chan '_' sprintf('%03d',t)])
    close all
end

figure
axes('Position',[0 0 1 1])
movie(M,1,10)

% Duration, location, assymetry
% Threshold Granger - and compute time to resolve
% Threshold. MPI coding
