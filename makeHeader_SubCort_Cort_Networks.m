function R = makeHeader_SubCort_Cort_Networks()
if strcmp(getenv('COMPUTERNAME'),'SFLAP-2')
    R.origpath = 'C:\data\TimExtracts280716\';
    R.datapathr = 'C:\Users\Tim\Documents\Work\GIT\Cortical_Parkinsons_Networks\Data\';
    R.resultspathr = 'C:\Users\Tim\Documents\Work\GitHub\Cortical_Parkinsons_Networks\Results';
elseif strcmp(getenv('COMPUTERNAME'),'FREE')
    R.origpath = 'C:\home\data\TimExtracts280716\';
    R.datapathr = 'C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\Data\';
    R.resultspathr = 'C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\Results';
end
R.subname = {'DF','DP','DS','JA','JB','JN','LM','LN01','LN02','LN03','MC','MW','WB'};
% BAD SUBJS:
% 'JP' - Issue with 'spm_cfg_eeg_inv_headmodel' - 'spm_eeg_artefact' is
% marking all MEG channels bad. [~ 1 2 2] specifically
% 'KB' - Only has 1 Left STN channel
% 'SW' - Spectra look flat/weird

R.condname = {'OFF','ON'};
R.condnamelc = {'off','on'};

R.siden = {'Left','Right'};
R.titular = {'CTX','STN'}; % Possibly outdated
R.bregname = {'STG','SMA'};
R.bregband = {[1],[3]};
R.bregROI = {[-46 -30 -2;46 -30 -2];[-18 -6 58; 18 -6 58]};
% Virtual Channel Construction
R.voi.radius = 15;
R.voi.resolution = 5;

R.specanaly.epochL = 1; % epoch length (s)
R.specanaly.tapsmofrq = 3; % tapers
R.specanaly.npd.taper = 1;
R.specanaly.npd.taperwid = 'M3';

R.ipsicon = 'ipsi';
R.ref_list = {'STN_L01','STN_L12','STN_L23','STN_R01','STN_R12','STN_R23'};
R.bandname = {'Alpha','Low Beta','High Beta'};
R.bandinits = {'\alpha','\beta_1','\beta_2'};
R.bandef = [8 13; 14 20; 24 34];

R.condcmap = linspecer(2);

% Preprocessing
R.pp.cont.thin.fs = 128;
R.pp.cont.thin.bp = [4 48];
% R.pp.cont.full.bp = [4 400];
% R.pp.cont.full.fs = 1024;

% Phase Analy
R.PA.optimalPLFrqMeth = 'WPLV'; % Use PLI or PLV or WPLV to determine frequency for bandpass
R.PA.AmpSurrN = 100; % Number of draws to compute surrogate distributions for stats.
R.PA.SRPeps_prctile = 1; % Stable relative phase Percentile
R.PA.SNReps_prctile = 50; % Signal Noise Percentile
R.PA.PLVeps_prctile = 85; % PLV Percentile
R.PA.plotting.realignMeth = 'WghtedPrctleAmp75'; % Method to align phases
% 'MaxLBAmp' - Maximum Low Beta Amp
% 'WghtedMeanAmp' -Wghted Mean Low Beta Amp
% 'WghtedPrctleAmp##' -Wghted Prctile Low Beta at ==%
% 'PhiMean'  - Centre to Phi Mean
% 'PhiHist'  - Bin with highest frequency
% 'noshift'  - Dont do any shifting
R.PA.SType = 1; % 1 = sliding window PLI and 2 = SRP
R.PA.bwid = [1.5 1.5 1.5];
R.PA.mwid = 3; % minimum SRP length (cycles)
R.PA.LowAmpFix = 0; % 1 if SRP is adjusted to account for low amplitude

R.PA.interpolgrid = 1; % For interpolating the grids
R.PA.frqrange{1} = R.bandef(1,1):0.5: R.bandef(1,2);
R.PA.frqrange{2} =  R.bandef(2,1):0.5: R.bandef(2,2);
R.PA.frqrange{3} =  R.bandef(3,1):0.5: R.bandef(3,2);


% R.PA.PLVeps =  0.50;
% R.PA.WinOver = 0.98;
% R.PA.SRPeps = 0.0015; %0.006;
% R.PA.stn_lb_frq = 14;
% R.PA.SNR_eps = -1;
% R.PA.bwid = 0.75;
R.PA.slidingwindow = 1;
R.PA.WinOver = 0.99;
% R.PA.stn_lb_frq = 14;
% R.PA.frqrange{1} = 24:0.5:34;
% R.PA.frqrange{2} = 8:0.5:12;
% R.PA.SNR = -1;
% R.PA.SNR = [-8 -8 -8]; %-1.5;
