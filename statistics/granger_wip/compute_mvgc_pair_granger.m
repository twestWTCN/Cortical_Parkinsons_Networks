function spgc = compute_mvgc_pair_granger(X,fsamp,blim,diagnos)
% clear; close all
% addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\mvgc_v1.0')
% run startup.m
% load testdata.mat
% blim = [24 34];
% fsamp = 1024;
[y, index, reject] = slideWindow(1:size(X,1)',128, 32);
for i = 1:size(index,2)-1
    XB = X(index(:,i),:)';
    XA(:,:,i) = XB-mean(XB)/std(XB);
end
X = XA;
ntrials   = size(X,3);     % number of trials
nobs      = size(X,2);   % number of observations per tria

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'BIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 32;     % maximum model order for model order estimation

acmaxlags = fsamp*5;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'chi2';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = fsamp; %200;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)

nvars = 2; % number of variables

% Residuals covariance matrix.

SIGT = eye(nvars);
% X is time domain (Ch x Time x Samples)


% Calculate information criteria up to specified maximum model order.
ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.
if diagnos
figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');
end
% amo = size(AT,3); % actual model order

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);
% fprintf('actual model order     = %d\n',amo);

% Select model order.

if     strcmpi(morder,'actual')
    morder = amo;
    fprintf('\nusing actual model order = %d\n',morder);
elseif strcmpi(morder,'AIC')
    morder = moAIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = moBIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end
% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

% NOTE: at this point we have a model and are finished with the data! - all
% subsequent calculations work from the estimated VAR parameters A and SIG.
% The autocovariance sequence drives many Granger causality calculations (see
% next section). Now we calculate the autocovariance sequence G according to the
% VAR model, to as many lags as it takes to decay to below the numerical
% tolerance level, or to acmaxlags lags if specified (i.e. non-empty).

ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

% The above routine does a LOT of error checking and issues useful diagnostics.
% If there are problems with your data (e.g. non-stationarity, colinearity,
% etc.) there's a good chance it'll show up at this point - and the diagnostics
% may supply useful information as to what went wrong. It is thus essential to
% report and check for errors here.

var_info(info,true); % report results (and bail out on error)
% Calculate time-domain pairwise-conditional causalities - this just requires
% the autocovariance sequence.

ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.
if diagnos
figure(2); clf;
subplot(1,3,1);
plot_pw(F);
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval);
title('p-values');
subplot(1,3,3);
plot_pw(sig);
title(['Significant at p = ' num2str(alpha)])
end
% For good measure we calculate Seth's causal density (cd) measure - the mean
% pairwise-conditional causality. We don't have a theoretical sampling
% distribution for this.

cd = mean(F(~isnan(F)));

fprintf('\ncausal density = %f\n',cd);

% Calculate spectral pairwise-conditional causalities at given frequency
% resolution - again, this only requires the autocovariance sequence.

ptic('\n*** autocov_to_spwcgc... ');
f = autocov_to_spwcgc(G,fres); % f is Ch x Ch x Hz
ptoc;

% Check for failed spectral GC calculation

assert(~isbad(f,false),'spectral GC calculation failed');

% Plot spectral causal graph.
if diagnos
figure(3); clf;
plot_spw(f,fs); %
end
fvec = linspace(0,fs/2,size(f,3));
spgc(:,1) = [mean(f(1,2,fvec>=blim(1) & fvec<blim(2))) pval(1,2) cd];
spgc(:,2) = [mean(f(2,1,fvec>=blim(1) & fvec<blim(2))) pval(2,1) cd];
gc_diff = (mean(f(2,1,fvec>=blim(1) & fvec<blim(2))).*(pval(1,2)<0.05)) - (mean(f(1,2,fvec>=blim(1) & fvec<blim(2))).*(pval(2,1)<0.05));
if pval(1,2)<0.05; a = pval(1,2); else a=1; end
if pval(2,1)<0.05; b = pval(2,1); else b=1; end
pvaldiff = a*b;
if gc_diff==0; gc_diff = NaN; end
spgc(:,3) = [gc_diff pvaldiff cd];

% Check that spectral causalities average (integrate) to time-domain
% causalities, as they should according to theory.

fprintf('\nchecking that frequency-domain GC integrates to time-domain GC... \n');
Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
mad = maxabs(F-Fint);
madthreshold = 1e-5;
if mad < madthreshold
    fprintf('maximum absolute difference OK: = %.2e (< %.2e)\n',mad,madthreshold);
else
    fprintf(2,'WARNING: high maximum absolute difference = %e.2 (> %.2e)\n',mad,madthreshold);
end