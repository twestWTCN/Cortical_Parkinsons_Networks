%%% MASTER PROGRAM FOR ANALYSIS OF PATIENT M1/STN PHASE COUPLING AND
%%% COHERENCE
clear; close all

% Initialise
add_corticalnetworks_paths()
R = makeHeader_SubCort_Cort_Networks();

% Extract Virtual Channels and STN Data
sourcextractloop_cortnet(R)
preprocess_cont_VC_STN_V6(R)
find_STNref_V6(R)
plot_subject_cohspectra_V6(R)
compute_phase_amp_analysis_v6(R)
getHists_phase_amp_analysis_PLIs(R)
plot_phase_amp_array(R)
plot_group_level_phaseamp_v6(R)

plot_phase_amp_analysis_PLIs_Cagnan(R)

%%
%%network_plot_270317


