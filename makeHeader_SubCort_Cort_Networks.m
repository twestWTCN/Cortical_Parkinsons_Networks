function R = makeHeader_SubCort_Cort_Networks()
if strcmp(getenv('COMPUTERNAME'),'SFLAP-2')
    R.origpath = 'C:\data\TimExtracts190516\';
    R.datapathr = 'C:\Users\Tim\Documents\Work\GIT\Cortical_Parkinsons_Networks\Data\';
elseif strcmp(getenv('COMPUTERNAME'),'FREE')
    R.origpath = 'C:\home\data\TimExtracts190516\';
    R.datapathr = 'C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\Data\';
end
R.subname = {'DF','DP','DS','JA','JB','JN','LN01','LN02','LN03','MC','MW','SW'};
R.condname = {'ON','OFF'};
R.bandname = {'Alpha','Low Beta','High Beta'};
R.condnamelc = {'on','off'};
R.ref_list = {'STN_L01','STN_L12','STN_L23','STN_R01','STN_R12','STN_R23'};
R.ipsicon = 'ipsi';
R.siden = {'Left','Right'};
R.titular = {'CTX','STN'};
R.bandef = [8 12; 14 21; 22 34];
R.condcmap = linspecer(2);
% Preprocessing
% Cont
R.pp.cont.full.bp = [4 400];
R.pp.cont.full.fs = 1024;
R.pp.cont.thin.fs = 256;
R.pp.cont.thin.bp = [4 98];

% Beamforming
% DICs
R.dics.bp = [8 12  10 2;
             NaN NaN NaN NaN;
             24 34 29 5;
              ];

R.ROI.maskrho_global = 2;
R.ROI.maskrho_local  = 1.5;
R.ROI.maskrho_vc = 1;
R.ROI.bandROI(:,1,1) = [-46 -30 2]./10; % R_STG
R.ROI.bandROI(:,1,2) = [46 -30 2]./10; % L_STG
R.ROI.bandROI(:,2,1) = [NaN NaN NaN]; % L_STG
R.ROI.bandROI(:,2,2) = [NaN NaN NaN]; % L_STG
R.ROI.bandROI(:,3,1) = [-20 -6 84]./10; % premotor_R
R.ROI.bandROI(:,3,2) = [20 -6 84]./10; % premotor_L

R.VC.normalize = 1; % Normalize
% NPD
% % % R.NPD.multitaper = 'M0.5'; % NO MULTITAPER"!!
R.NPD.windowlength = 9;

% Phase Analy
% R.PA.bwid = 0.75;
% R.PA.slidingwindow = 2.5;
% R.PA.PLVeps =  0.35;
% R.PA.mwid = 12;
% R.PA.WinOver = 0.95;
% R.PA.stn_lb_frq = 14;
% R.PA.frqrange{1} = 24:0.5:34;
% R.PA.frqrange{2} = 8:0.5:12;
% R.PA.SNR = -1;
R.PA.SType = 2; % 1 = sliding window PLI and 2 = SRP
R.PA.bwid = [1 1 1];
R.PA.mwid = 6; % minimum SRP length (cycles)
R.PA.SRPeps = 0.0015; %0.006;
R.PA.LowAmpFix = 0; % 1 if SRP is adjusted to account for low amplitude

R.PA.slidingwindow = 1;
R.PA.PLVeps =  0.50;
R.PA.WinOver = 0.98;
% R.PA.stn_lb_frq = 14;
% R.PA.SNR_eps = -1;

R.PA.frqrange{1} = R.bandef(1,1):0.5: R.bandef(1,2);
R.PA.frqrange{2} =  R.bandef(2,1):0.5: R.bandef(2,2);
R.PA.frqrange{3} =  R.bandef(3,1):0.5: R.bandef(3,2);
R.PA.SNR = [-8 -8 -8]; %-1.5;
