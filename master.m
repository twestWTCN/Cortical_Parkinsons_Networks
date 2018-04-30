% MASTER
clear; close all
% Curr list subname = {'JN','MC','SW','DF','JB','MW','DP','DS','JA'};
%% TO DO:
% 1) Theshold segments as ratio of amplitude of envelope of bandpass vs total
% to yield an SNR
% 2) Look at Alpha networks - how do they coincide with beta frames? 
% 3) Adapt scripts to cope with LM - missing two STN channels left
add_corticalnetworks_paths()
R = makeHeader_SubCort_Cort_Networks();
% R = makeHeader_SubCort_Cort_Networks_SFLAP();
%% Compute Forward
% subname = {'DF'};
% % steps = {'loadmri','realign_coords','skullstrip','segment','reslice','headmodel','coregister'}; 
% steps = {'reslice','headmodel','coregister'}; %  'loadmri','realign_coords','skullstrip','segment'
% prepare_forward(datapathr,subname,steps,0)
% prepare_forward_v2(datapathr,subname,steps,0)
% prepare_forward_v3(datapathr,subname,steps,0)
%% Preprocess Data 
R.subname = {'LN01','LN02','LN03'};
preprocess_cont_MEG(R)
preprocess_epoched_MEG(R)

%% Compute Beamformer Weights
% subname = {'DF','DP','DS','JA','JB','JN','JP','KB','LM','MC','MW','SW','WB'}; 
% subname = {'LM'};
% compute_leadfield_LCMV(datapathr,subname,steps,0,[6 40 NaN NaN])
steps = {'leadfield','plotsourcepow'}; %,'leadfield',};
compute_leadfield_LCMV_v4(R,steps,0)

% compute_leadfield_DICS(datapathr,subname,steps,0,[24 30 27 3])
steps = {'leadfield','plotsourcepow'};
% compute_leadfield_DICS_v2(R,steps,0)
compute_leadfield_DICS_v4(R,steps,0)

% subname = {'DF'};
% decide_max_coh(datapathr,subname)
%%% subname = {'JN','MC','SW','DF','JB','MW','DP','DS','JA'}; %%%
R.subname = {'JN','MC','SW','DF','JB','MW','DP','DS','JA','LN01','LN02','LN03'};
% group_DICS_imagemean(R)
% group_DICs_imageanalyse % not sure what this does
decide_max_coh_v3(R,0)
plot_DICs_LocalROIdist(R)
compute_virtual_electrodes_ROI_contdata_v3(R)
preprocess_epoched_VC(R)
plot_subject_cohspectra(R)
grandAverageSpectra_Screened(R)

compute_phase_amp_analysis_v6(R)
getHists_phase_amp_analysis_PLIs(R)
plot_phase_amp_analysis_PLIs_V3_2side(R)
plot_group_level_phaseamp_v3(R)
network_plot_270317


% LEGACY:
% decide_max_coh_v2(R,0)
% compute_phase_amp_analysis_050118(R)
% compute_phase_amp_analysis_120318(R)
% plot_phase_amp_analysis(R)
% plot_phase_amp_analysis_PLIs(R)
% plot_group_level_phaseamp
