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
R.bandname = {'High Beta','Alpha'}
R.condnamelc = {'on','off'};
R.ref_list = {'STN_L01','STN_L12','STN_L23','STN_R01','STN_R12','STN_R23'};
R.ipsicon = 'ipsi';
R.siden = {'Left','Right'};
R.titular = {'CTX','STN'};
R.condcmap = linspecer(2);
% Preprocessing
% Cont
R.pp.cont.full.bp = [4 400];
R.pp.cont.full.fs = 1024;
R.pp.cont.thin.fs = 256;
R.pp.cont.thin.bp = [4 98];

% Beamforming
% DICs
R.dics.bp = [24 34 29 5;
              8 12  10 2];

R.ROI.maskrho_dic = 2;
R.ROI.maskrho_vc = 1;

% NPD
R.NPD.multitaper = 'M2';
R.NPD.windowlength = 10;

% Phase Analy
R.PA.bwid = 0.75;
R.PA.slidingwindow = 2.5;
R.PA.PLVeps =  0.35;
R.PA.mwid = 12;
R.PA.WinOver = 0.95;
R.PA.stn_lb_frq = 14;
R.PA.frqrange{1} = 24:0.5:34;
R.PA.frqrange{2} = 8:0.5:12;
R.PA.SNR = -1;