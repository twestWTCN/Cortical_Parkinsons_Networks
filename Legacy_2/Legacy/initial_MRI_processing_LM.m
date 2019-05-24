% nii raw MRI processing to fieldtrip output.
clear; close all

addpath(genpath('C:\Users\Tim\Documents\Work\Cortical_Networks'))

datapathr = 'C:\Users\Tim\Documents\Work\Cortical_Networks\Data\';
datafilen = 'ls050651_vladimirDBS_20071103_06.ds';
subname = {'LM'};
cd([datapathr subname{1} '_tmp\MRI\'])
mri = ft_read_mri([datapathr subname{1} '_tmp\MRI\r' subname{1} '.nii'],'dataformat','nifti_spm');


%% Reslice
cfg = [];
cfg.resolution = 1;
cfg.dim = [256 256 256];
mri = ft_volumereslice(cfg,mri);
save(['r' subname{1} 'rs'],'mri')
volumewrite_spm([datapathr subname{1} '_tmp\MRI\r' subname{1} '_rs.nii'], mri.anatomy, mri.transform, 'SPM12');

%% Realignment
load(['r' subname{1} 'rs'])
if exist([datapathr subname{1} '_tmp\MRI\r' subname{1} 'ra.nii'])>0
    delete([datapathr subname{1} '_tmp\MRI\r' subname{1} 'ra.nii']);
end
fidtab = fiducial_reader([datapathr subname{1} '_tmp\MRI\' subname{1} '_smri_fid.txt']);

mri = ft_determine_coordsys(mri,'interactive','no');
% Compute Transforms
cfg = [];
% cfg.method = 'interactive'
% cfg.method = 'fiducial';
% cfg.coordsys = 'spm';
% cfg.fiducial.nas    = table2array(fidtab(1,2:4)); % position of nasion
% cfg.fiducial.lpa    = table2array(fidtab(2,2:4)); % position of LPA
% cfg.fiducial.rpa    = table2array(fidtab(3,2:4)); % position of RPA
% cfg.fiducial.zpoint = [45 45 75];
% mri_spm = ft_volumerealign(cfg,mri);
% transform_vox2spm = mri_spm.transform;
% save(['r' subname{1} '_ra'],'mri_spm')

cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'ctf';
cfg.fiducial.nas    = table2array(fidtab(1,2:4)); % position of nasion
cfg.fiducial.lpa    = table2array(fidtab(2,2:4)); % position of LPA
cfg.fiducial.rpa    = table2array(fidtab(3,2:4)); % position of RPA
cfg.fiducial.zpoint = [45 45 75];
mri_ctf = ft_volumerealign(cfg,mri);
transform_vox2ctf = mri_ctf.transform;

save(['r' subname{1} '_ra'],'mri_ctf')
volumewrite_spm([datapathr subname{1} '_tmp\MRI\r' subname{1} '_ra.nii'], mri.anatomy, mri.transform, 'SPM12');
%% Skullstrip
if exist([datapathr subname{1} '_tmp\MRI\r' subname{1} '_ss.nii'])>0
    delete([datapathr subname{1} '_tmp\MRI\r' subname{1} '_ss.nii']);
end

system(['C:\Users\Tim\Documents\ROBEX\runROBEX.bat ' datapathr subname{1} '_tmp\MRI\r' subname{1} '_ra.nii ' datapathr subname{1} '_tmp\MRI\r' subname{1} '_ss.nii'])

%% Segmentations
iterNum_outer=15;  % outer iteration
iterCM=2;  % inner interation for C and M
iter_b=1;  % inner iteration for bias
q = 1.5;   % fuzzifier
th_bg = 5;  %% threshold for removing background
N_region = 3; %% number of tissue types, e.g. WM, GM, CSF
tissueLabel=[1, 2, 3];

MICO_3Dseq({[datapathr subname{1} '_tmp\MRI\r' subname{1} '_ss.nii']}, N_region, q, th_bg, iterNum_outer, iter_b, iterCM, tissueLabel);

segmri = ft_read_mri([datapathr subname{1} '_tmp\MRI\r' subname{1} '_ss.nii'],'dataformat','nifti_spm');

% get segmentation
segmri.brain  = segmri.anatomy>0;
segmri = rmfield(segmri,'anatomy');
save(['r' subname{1} '_seg'],'segmri')

%% Prepare Headmodel
cfg = [];
cfg.tissue = {'brain'};
cfg.numvertices = [1200];
mesh = ft_prepare_mesh(cfg,segmri);

cfg = [];
cfg.tissue = {'brain'};
cfg.method = 'singleshell'; %'openmeeg';
vol = ft_prepare_headmodel(cfg,segmri);

%% Apply coregistration
T = transform_vox2ctf; %/transform_vox2spm;%;
vol2 = ft_transform_geometry(T,vol);
vol2 = ft_convert_units(vol2,'cm');

sens = ft_read_sens([datapathr subname{1} '_tmp\' datafilen],'senstype','meg');
% PLOT
% close all
% figure
% ft_plot_sens(sens,'style','*b'); %,'coil','true')
% 
% hold on
% ft_plot_vol(vol2)

%% Prepare Leadfield
cfg = [];
cfg.datafile = [datapathr subname{1} '_tmp\' datafilen];
cfg.lpfilter = 'yes';
cfg.hpfilter = 'yes';
cfg.channel         = {'MEG'};
cfg.lpfreq = 18;
cfg.hpfreq = 14;
megdata = ft_preprocessing(cfg);

cfg = [];
cfg.resamplefs = 400;
megdata = ft_resampledata(cfg,megdata);
save([datapathr subname{1} '_tmp\ftdata\megdata'],'megdata')

% compute CSD
load([datapathr subname{1} '_tmp\ftdata\megdata'])

cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.channel         = {'MEG'};
cfg.tapsmofrq = 2;
cfg.foilim    = [18 18];
freqdata = ft_freqanalysis(cfg, megdata);


% compute leadfield
cfg                 = [];
cfg.grad            = freqdata.grad;
cfg.normalize       = 'yes';
cfg.headmodel       = vol2;
% cfg.reducerank      = 2;
cfg.channel         = {'MEG'};
cfg.grid.resolution = 0.5;   % use a 3-D grid with a 1 cm resolution
cfg.grid.unit       = 'cm';
grid = ft_prepare_leadfield(cfg);

cfg              = []; 
cfg.method       = 'dics';
cfg.frequency    = 18;  
cfg.grid         = grid; 
cfg.headmodel    = vol2;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = 0;
source = ft_sourceanalysis(cfg,freqdata);

% Nueral activity index ( avoid centre of head bias)
source.avg.pow = source.avg.pow ./ source.avg.noise;

T = transform_vox2ctf; %/transform_vox2spm;%;
Tmri = ft_transform_geometry(T,mri);
Tmri = ft_convert_units(Tmri,'cm');

cfg            = [];
% cfg.downsample = 1;
cfg.parameter = 'avg.pow';
sourceInt  = ft_sourceinterpolate(cfg, source , Tmri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'avg.pow';
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
volumewrite_spm([datapathr subname{1} '_tmp\MRI\r' subname{1} '_sourcepow.nii'], sourceIntNorm.pow, sourceIntNorm.transform, 'SPM12');

cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'avg.pow';
cfg.maskparameter  = cfg.funparameter;
% cfg.funcolorlim    = [0.0 1.2];
cfg.funcolormap    = 'jet';
% cfg.opacitylim     = [0.0 1.2]; 
cfg.opacitymap     = 'auto';  
cfg.projmethod     = 'nearest'; 
cfg.surffile       = 'surface_white_both.mat';
cfg.surfdownsample = 2; 
ft_sourceplot(cfg, sourceIntNorm);
view ([90 0])

