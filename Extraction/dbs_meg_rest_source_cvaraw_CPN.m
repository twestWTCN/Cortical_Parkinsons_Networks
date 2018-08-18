function Ds = dbs_meg_rest_source_cvaraw_CPN(initials, drug, prefix, ROI,R,breg,siden,CVAfrqrng)

druglbl = R.condnamelc;

try
    [files, seq, root, details] = dbs_subjects(initials, drug);
catch
    return;
end

if nargin<3
    prefix = '';
end
% From dbs_meg_rest_extract.. [prefix initials '_' seq{f}(1) '_' num2str(cnt) '_' druglbl{drug+1}]
try
    %     D = spm_eeg_load(fullfile(root, 'SPMrest', [prefix initials '_' druglbl{drug+1} '.mat']));
    D = spm_eeg_load([R.origpath initials '\' initials '_R_1_' druglbl{drug+1}]);
catch
    error('not found');
    %     D = dbs_meg_rest_prepare_spm12(initials, drug);
end

% Interpolate Missing Data
if siden == 'L'
    stnKO = find(strncmp(D.chanlabels,'STN_R',5));
elseif siden == 'R'
    stnKO = find(strncmp(D.chanlabels,'STN_L',5));
end
chlist = setdiff(1:numel(D.chanlabels),stnKO);

S = [];
S.D = D;
S.channels = D.chanlabels(chlist);
D = spm_eeg_crop(S);

% Missing Data
a = D(:,:,:);
for x = 1:size(a,1)
    if sum(isnan(a(x,:)))> 1
        warning([D.chanlabels{x} ' has ' num2str(sum(isnan(a(x,:)))) ' missing values, interpolating'])
        [F,TF] = fillmissing( a(x,:),'linear','SamplePoints',1:length(a(x,:)));
        a(x,:) = F;
    end
end
D(:,:,:) = a;
clear a

% Jumps
a = D(:,:,:);
for x = 1:size(a,1)
    dat = a(x,:);
    dat(abs(dat)>(6*std(dat))) = NaN;
    dat= fillmissing(dat,'spline');
    a(x,:) = dat;
end
D(:,:,:) = a;
clear a




% S = [];
% S.D = D; S.fsample_new = 128;
% D = spm_eeg_downsample(S);

% cd(fullfile(root, 'SPMrest'));
cd([R.datapathr initials]);

res = mkdir('BF');


roisum = 'max';%'svd';%'keep' 'max'
scramblestn = 'no';
keep = 0;
%%
source{1}.spm.tools.beamforming.data.dir = {[R.datapathr initials '\BF']};
source{1}.spm.tools.beamforming.data.D = {fullfile(D)};
source{1}.spm.tools.beamforming.data.val = 1;
source{1}.spm.tools.beamforming.data.gradsource = 'inv';
source{1}.spm.tools.beamforming.data.space = 'MNI-aligned';
source{1}.spm.tools.beamforming.data.overwrite = 1;
source{2}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
source{2}.spm.tools.beamforming.sources.reduce_rank = [2 3];
source{2}.spm.tools.beamforming.sources.keep3d = 1;

source{2}.spm.tools.beamforming.sources.plugin.voi.radius = R.voi.radius;
source{2}.spm.tools.beamforming.sources.plugin.voi.resolution = R.voi.resolution;
source{2}.spm.tools.beamforming.sources.visualise = 1;
source{3}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
source{3}.spm.tools.beamforming.features.whatconditions.all = 1;
source{3}.spm.tools.beamforming.features.woi = [-Inf Inf];
source{3}.spm.tools.beamforming.features.modality = {'MEG'};
source{3}.spm.tools.beamforming.features.fuse = 'no';
source{3}.spm.tools.beamforming.features.plugin.cov.foi = [0 98];
source{3}.spm.tools.beamforming.features.plugin.cov.taper = 'none';
if details.oxford
    source{3}.spm.tools.beamforming.features.regularisation.manual.lambda = 0;
else
    source{3}.spm.tools.beamforming.features.regularisation.manual.lambda = 0.01;
end
source{3}.spm.tools.beamforming.features.bootstrap = false;
source{4}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
source{4}.spm.tools.beamforming.inverse.plugin.lcmv.orient = true;
source{4}.spm.tools.beamforming.inverse.plugin.lcmv.keeplf = false;
source{5}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
source{5}.spm.tools.beamforming.output.plugin.montage.method = roisum;
source{5}.spm.tools.beamforming.output.plugin.montage.voidef = struct('label', {}, 'pos', {}, 'radius', {});
source{6}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
source{6}.spm.tools.beamforming.write.plugin.spmeeg.mode = 'write';
source{6}.spm.tools.beamforming.write.plugin.spmeeg.modality = 'MEG';
source{6}.spm.tools.beamforming.write.plugin.spmeeg.addchannels.channels{1}.type = 'LFP';
source{6}.spm.tools.beamforming.write.plugin.spmeeg.prefix = [R.bregname{breg} '_VC04052018_'];

for i = 1:size(ROI, 1)
    source{2}.spm.tools.beamforming.sources.plugin.voi.vois{i}.voidef.label = ROI{i, 1};
    source{2}.spm.tools.beamforming.sources.plugin.voi.vois{i}.voidef.pos = ROI{i, 4};
    %     source{2}.spm.tools.beamforming.sources.plugin.voi.vois{i}.maskdef.label = ROI{i, 1};
    %     source{2}.spm.tools.beamforming.sources.plugin.voi.vois{i}.maskdef.mask =  ROI(i, 2);
end

%************
% source = source(5:6);
% source{1}.spm.tools.beamforming.output.plugin.montage.method = roisum;
% source{1}.spm.tools.beamforming.output.BF = {fullfile(root, 'SPMrest', 'BF', 'BF.mat')};
% source{2}.spm.tools.beamforming.write.BF = {fullfile(root, 'SPMrest', 'BF', 'BF.mat')};
%*************
spm_jobman('run', source);
% spm_jobman('interactive', source)

Ds = spm_eeg_load(fullfile(D.path, [R.bregname{breg} '_VC04052018_' D.fname]));
Ds = chantype(Ds, ':', 'PHYS');
Ds = chantype(Ds, Ds.indchannel(details.chan), 'LFP');
mkdir([R.datapathr initials '\SPMdata\']);
Ds = Ds.path([R.datapathr initials '\SPMdata\']);
save(Ds);
% Ds = spm_eeg_load('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\Data\DF\SPMdata\STG_VC04052018_dDF_R_1_off.mat')
%%
% if isequal(scramblestn, 'yes')
% stnchans = find(strncmp(Ds.chanlabels,'STN',3));
% for i = 1:numel(stnchans)
%     x = Ds(stnchans(i),:);
%     x = x(randperm(length(x)));
%     Ds(stnchans(i),:) = x;
% end
% save(Ds);
% end
if isequal(roisum, 'keep')

    
    S = [];
    S.D = Ds;
    S.prefix = 'R';
    S.channels = details.chan;
    S.method = 'cva';
    
    for i = 1:size(ROI, 1)
        S.settings.chanset(i).cvachan.channels{1}.regexp = ['^' ROI{i, 1}];
        S.settings.chanset(i).refchan.channels{1}.type = 'LFP';
        S.settings.chanset(i).ncomp = 1;
        S.settings.chanset(i).outlabel = ROI{i, 1};
        S.settings.chanset(i).foi = ROI{i, 3}; %CVAfrqrng; %
        S.settings.chanset(i).tshiftwin = [-80 80];
        S.settings.chanset(i).tshiftres = 10;
    end
    
    S.keepothers = false;
    S.conditions.all = 1;
    S.timewin = [-Inf Inf];
    Ds = spm_eeg_reduce(S);
    save(Ds);
    if ~keep, delete(S.D); end;
end
%%


