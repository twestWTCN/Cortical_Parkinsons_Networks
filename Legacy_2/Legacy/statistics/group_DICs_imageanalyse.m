function [] = group_DICs_imageanalyse(R)
load([datapathr 'results\images\groupDICsresults'],'source_avg_dics','source_dims')

a = subname{1}; b = condname{1};

ref_chan = R.reflist{1};
dataload = load([R.datapathr a '\ftdata\r' a '_DICSv2_source' b 'nrep_' num2str(1) '_' ref_chan],'source')
source = dataload.source;
T = [1  0   0    0
    0   1   0    0
    0   0   1    2
    0   0   0    1];
source = ft_transform_geometry(T,source);
Tmri = ft_read_mri([datapathr subname{1} '\MRI\orig\r' subname{1} '.nii'],'dataformat','nifti_spm');
Tmri = ft_convert_units(Tmri,'cm');

cfg            = [];
cfg.downsample = 2;
cfg.parameter = 'avg.coh';
sourceInt  = ft_sourceinterpolate(cfg, source , Tmri);

sider = [1:3;4:6];
for cond = 1:2
    for side =1:2
        x = nanmean(nanmean(source_avg_dics(:,sider(side,:),cond,:),2),4);
        X = reshape(x,squeeze(source_dims(:,1,1,1))');
        
        sourceInt.coh = X;
        a = 1;
        cfg = [];
        cfg.funparameter  = 'coh';
        cfg.method  = 'ortho';
         cfg.location      = 'max';
         cfg.axis = 'off';
         cfg.funcolorlim = [0 0.05];
        cfg.surffile       = 'surface_white_both.mat';
        cfg.funcolormap    = 'jet';
        %                         cfg.opacitylim     = [0. 0.05];
%         cfg.opacitymap     = 'auto';
        ft_sourceplot(cfg,sourceInt)
        set(gcf,'Position',[680 108 1070 870])
        set(gcf,'Name',['Grand Average ' condname{cond} ' ' siden{side} ' Hemi'])
        savefigure_v2([datapathr 'results\images\'],['groupDICsresults_' condname{cond} ' ' siden{side} ' Hemi'],[],[],[]); close all
    end
end
