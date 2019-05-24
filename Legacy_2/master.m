%%% MASTER PROGRAM FOR ANALYSIS OF PATIENT M1/STN PHASE COUPLING AND
%%% COHERENCE
clear; close all
%%
% TO DO:
% 1)Redo the phase-amp analysis but splitting the data up into smaller
% sections. Use the sample means to generate better statistics?
% Initialise
% 2)Consider doing the 2D plots using:
% https://uk.mathworks.com/matlabcentral/fileexchange/16458-bullseye-polar-data-plot
% 3) Consider plotting sigbar only on one side - i.e. the side that has
% shown the increase-
% 4) Consider simple model to test whether the signal power has any effect
% on the PA metrics.
% 5) Rethink the power metric - % over median scrambled? Does this mean
% anything...
add_corticalnetworks_paths()
R = makeHeader_SubCort_Cort_Networks();

% Extract Virtual Channels and STN Data
sourcextractloop_cortnet(R)
preprocess_cont_VC_STN_V6(R)
find_STNref_V6(R)
plot_subject_cohspectra_V6(R)
% compute_phase_amp_analysis_v6(R)
compute_phase_amp_analysis_v8(R)
getHists_phase_amp_analysis_PLIs(R)
plot_phase_amp_array(R)
plot_group_level_phaseamp_v6(R)

plot_phase_amp_analysis_PLIs_CagnanV8(R)
plot_phase_amp_analysis_PLIs_Cagnan_ExampleSub(R)
%%
%%network_plot_270317


