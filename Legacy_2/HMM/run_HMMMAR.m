%% LOAD IN 
% load('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\HMM\6_state_save')
%%%%%%%%%%%%%%%%%%%%%%%


addpath(('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks'))
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\bplot')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linregress')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\neurospec21')

% add_corticalnetworks_paths
R = makeHeader_SubCort_Cort_Networks;
addpath(genpath('C:\Users\twest\Documents\Work\GitHub\HMM-MAR'))
% This is a template script that specifies the usual options for a handful
% of data modalities, runs the HMM-MAR, and obtains some basic information
% from the estimation (mainly, the spectral information of the states).
%
% This script is just a starting point, and does not guarantee to be the best 
% configuration for the user's particular data.
%
% Also, it assumes that the data has been loaded into the right format
% (there are multiple possible input formats - see Wiki Documentation),
% that the variables data_modality, no_states ,Hz and stochastic_inference 
% (see below) are correctly set, 
% and that the toolbox paths are in the right place. 
%
% Author: Diego Vidaurre, OHBA, University of Oxford (2016)
load([R.datapathr 'HMMdata\HMMdataconcat'],'data','T')
% data = data';
% T = T';
fsamp = 128;
list = 1:size(data,1);
data = data(setdiff(list,[3 12 19 16 33 34 38]));
T = T(setdiff(list,[3 12 19 16 33 34 38]));
if ~exist('data','var') || ~exist('T','var')
    error('You need to load the data (data and T - see Documentation)')
end

data_modality = 'LFP' ; % one of: 'fMRI', 'M/EEG', 'M/EEG power' or 'LFP' 
no_states = 8; % the number of states depends a lot on the question at hand
Hz = fsamp; % the frequency of the data
stochastic_inference = 0; % set to 1 if a normal run is too computationally expensive (memory or time)
N = length(T); % number of subjects

% getting the number of channels
if iscellstr(data) 
    dfilenames = data;
    if ~isempty(strfind(dfilenames{1},'.mat')), load(dfilenames{1},'X');
    else X = dlmread(dfilenames{1});
    end
elseif iscell(data)
    X = data{1};
end
ndim = size(X,2); 

% Setting the options
options = struct();
options.K = no_states;
options.standardise = 1;
options.verbose = 1;
options.Fs = Hz;

if iscell(T), sumT = 0; for j = 1:N, sumT = sumT + sum(T{j}); end
else, sumT = sum(T); 
end

if strcmp(data_modality,'fMRI') % Gaussian observation model
    options.order = 0;
    options.zeromean = 0;
    options.covtype = 'full';     

elseif strcmp(data_modality,'M/EEG power') % Gaussian observation model
    options.order = 0;
    options.zeromean = 0;
    options.covtype = 'full';     
    options.onpower = 1; 
    
elseif strcmp(data_modality,'M/EEG') && ndim > 10 % Embedded observation model
    options.order = 0;
    options.zeromean = 1;
    options.covtype = 'full';
    options.embeddedlags = -7:7;  
    options.pca = ndim*2;
    options.standardise_pc = options.standardise;

elseif strcmp(data_modality,'M/EEG') && ndim <= 10 % MAR observation model
    options.order = 5;
    options.zeromean = 1;
    options.covtype = 'diag';

elseif strcmp(data_modality,'LFP') && ndim > 5 % MAR observation model
    options.order = 5;
    options.zeromean = 1;
    options.covtype = 'diag';
    if ndim > 2, options.pca = 0.9; end
    options.DirichletDiag = round(sumT/10); 

elseif strcmp(data_modality,'LFP') && ndim <= 5  % MAR observation model
    options.order = 5; %11;
    options.zeromean = 1;
    options.covtype = 'diag';
    if ndim > 2, options.pca = 0.9; end
    options.DirichletDiag = round(sumT/10); 
    
else
    error('Option data_modality not recognised')

end

% stochastic options
if stochastic_inference
    options.BIGNbatch = max(round(N/30),5);
    options.BIGtol = 1e-7;
    options.BIGcyc = 500;
    options.BIGundertol_tostop = 5;
    options.BIGforgetrate = 0.7;
    options.BIGbase_weights = 0.9;
end

% HMM computation
[hmm, Gamma, Xi, vpath] = hmmmar(data,T,options);

data

% Spectral estimation
options_spectra = struct(); 
options_spectra.Fs = Hz; % Sampling rate 
options_spectra.fpass = [6 45];  % band of frequency you're interested in
options_spectra.p = 0; %0.01; % interval of confidence  
options_spectra.to_do = [1 0]; % turn off pdc
options_spectra.K = hmm.K;
if strcmp(data_modality,'LFP') % MAR spectra
    options_spectra.to_do = [1 1];
    options_spectra.order = 5; 
    options_spectra.Nf = 90; 
    spectra = hmmspectramar(data,T,[],Gamma,options_spectra);
else % Multi-taper spectra
    options_spectra.to_do = [1 0];
    options_spectra.tapers = [6 7]; % internal multitaper parameter
    options_spectra.win = 10 * Hz;
    spectra = hmmspectramt(data,T,Gamma,options_spectra);
end

% Some useful information about the dynamics
maxFO = getMaxFractionalOccupancy(Gamma,T,options); % useful to diagnose if the HMM 
            % is capturing dynamics or grand between-subject 
            % differences (see Wiki)
FO = getFractionalOccupancy (Gamma,T,options); % state fractional occupancies per session
LifeTimes = getStateLifeTimes (Gamma,T,options); % state life times
Intervals = getStateIntervalTimes (Gamma,T,options); % interval times between state visits
SwitchingRate =  getSwitchingRate(Gamma,T,options); % rate of switching between stats
onsets = getStateOnsets(vpath,T,Hz);
S = getStartFinish_StateVisits (vpath,T,options);

TP =getTransProbs(hmm);
TP(TP<0.3) = 0;


Scor = cell(size(T,1),size(S,2));
Scor(find(~cellfun(@isempty,T)),:) = S
Ts = get_state_tseries(data,Scor,idx(1:10)');

%% Subject Specific
T_adj = 25276;
Xi_subj = Xi(1:T_adj,:); 
Gamma_subj = Gamma(1:T_adj,:);
T_subj = T{1};
data_subj = data{1};
% 25276
hmm_subj = subject_hmm({data_subj},{T_subj},hmm,Gamma_subj); %,Xi_subj
 TP_subj =getTransProbs(hmm_subj);
TP_subj(TP_subj<0.3) = 0