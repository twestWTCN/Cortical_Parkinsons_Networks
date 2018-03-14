clear; close all
datapathr = 'C:\Users\Tim\Documents\Work\Cortical_Networks\Data\';
subname = {'LM'};
ref_chan = 'STN_R12';

i = 1
load([datapathr subname{i} '\ftdata\r' subname{i} '_LCMV_source'])
load('C:\Users\Tim\Documents\Work\Cortical_Networks\tmp_files\source_sens.mat')
source_sens = reshape(source_sens,source.dim);
ftdata_cell = cell(source.dim);
for x = 1:source.dim(1)
    for y = 1:source.dim(2)
        for z = 1:source.dim(3)
            ftdata = [];
            ftdata.fsample = 200;
            ftdata.label = {sprintf('%.1f ',[x y z]); 'STN_R12'};
            ftdata.trial = source_sens{x,y,z};
            ftdata.time = repmat({linspace(0,2.5,length(ftdata.trial{1}') )},1,size(ftdata.trial,2));
            ftdata_cell{x,y,z} = ftdata;
        end
    end
end

% cfg = [];
% cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
% ft_databrowser(cfg, ftdata_cell{6,4,7});
gASym = zeros(source.dim);
isum = sum(source.dim);
% Compute connectivities
for x = 1:source.dim(1)
    for y = 1:source.dim(2)
        for z = 1:source.dim(3)
            if size(ftdata_cell{x,y,z}.trial{1},2)>1
                cfg         = [];
                cfg.order   = 15;
                cfg.method = 'bsmart';
                mdata       = ft_mvaranalysis(cfg, ftdata_cell{x,y,z});
                
                cfg        = [];
                cfg.method = 'mvar';
                mfreq      = ft_freqanalysis(cfg, mdata);
                
                %             cfg           = [];
                %             cfg.method    = 'mtmfft';
                %             cfg.taper     = 'dpss';
                %             cfg.output    = 'fourier';
                %             cfg.pad       = 128;
                %             cfg.tapsmofrq = 1;
                %             freq          = ft_freqanalysis(cfg, ftdata_cell{x,y,z});
                
                cfg           = [];
                cfg.method    = 'coh';
                %             coh           = ft_connectivityanalysis(cfg, freq);
                cohm          = ft_connectivityanalysis(cfg, mfreq);
                
                %             figure
                %             cfg           = [];
                %             cfg.parameter = 'cohspctrm';
                %             cfg.zlim      = [0 1];
                %             ft_connectivityplot(cfg, cohm) %, coh);
                
                cfg           = [];
                cfg.method    = 'granger';
                granger       = ft_connectivityanalysis(cfg, mfreq);
                indz = find(granger.freq>=14 & granger.freq<30);
                gASym(x,y,z) = mean(squeeze(granger.grangerspctrm(2,1,indz) - granger.grangerspctrm(1,2,indz)));
%                 figure
%                 cfg           = [];
%                 cfg.parameter = 'grangerspctrm';
%                 cfg.zlim      = [0 0.1];
%                 ft_connectivityplot(cfg, granger);
            else
                gASym(x,y,z) = NaN;
            end
        end
        disp([x y z])
    end
end
load('C:\Users\Tim\Documents\Work\Cortical_Networks\tmp_files\source_granger.mat')
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
source.avg.granger = reshape(gASym,[],1);
load([datapathr subname{i} '\ftdata\rs_transform_vox2ctf'])
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

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'granger';
cfg.maskparameter = cfg.funparameter;
% cfg.funcolorlim   = [4.0 6.2];
% cfg.opacitylim    = [4.0 6.2];
cfg.opacitymap    = 'rampup';
figure
ft_sourceplot(cfg, sourceInt);

sourceInt.coordsys = 'ctf';
cfg = [];
cfg.nonlinear     = 'no';
sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
volumewrite_spm([datapathr subname{i} '\MRI\sourcespace\r' subname{i} '_' num2str(nr) '_DICS_sourcegranger_int_' ref_chan '.nii'], sourceIntNorm.granger, sourceIntNorm.transform, 'SPM12');

cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'granger';
cfg.maskparameter  = cfg.funparameter;
% cfg.funcolorlim    = [0.0 1.2];
cfg.funcolormap    = 'jet';
% cfg.opacitylim     = [0.0 1.2];
cfg.opacitymap     = 'auto';
cfg.projmethod     = 'nearest';
cfg.surffile       = 'surface_white_both.mat';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
H = figure;
ft_sourceplot(cfg, sourceIntNorm);
view ([90 0])

savefig(gcf,[datapathr subname{i} '\images\sourcespace\r' subname{i} '_rep_' num2str(nr) '_surf_avg_DICS_granger_' ref_chan])

