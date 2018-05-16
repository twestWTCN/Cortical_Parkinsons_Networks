clear; close all
load('Xdata.mat')
load('surdata.mat')
Odata= Xdata;
% [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Odata,37,16,0.5);
% qstable = find(abs(dphi_12_dt')<0.01); % The points that are below threshold
% consecSegs = SplitVec(qstable,'consecutive');
% % lengths
% segL_ddt = cellfun('length',consecSegs);
% histogram(log(segL_ddt/256),'FaceColor','b'); hold on



surfft1 = fft(surdata(1,:));
Amp1 = abs(surfft1);
Phi1 = angle(surfft1);

surfft2 = fft(surdata(2,:));
Amp2 = abs(surfft2);
Phi2 = angle(surfft2);


N = 6;
eps = logspace(-3,-2,N)
for i = 1:N
    figure
[amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Odata,37,16,0.5);
qstable = find(abs(dphi_12_dt')<eps(i)); % The points that are below threshold
consecSegs = SplitVec(qstable,'consecutive');
% lengths
segL_ddt = cellfun('length',consecSegs);
histogram(log(segL_ddt/256),'FaceColor','b'); hold on

Phishuff = Phi1(randperm(length(Phi1)));
phasescm1 = surfft1(randperm(length(surfft1))); % real(ifft(Amp1.*exp(sqrt(-1)*Phishuff)));
Phishuff = Phi2(randperm(length(Phi2)));
phasescm2 = surfft2(randperm(length(surfft2))); %real(ifft(Amp2.*exp(sqrt(-1)*Phishuff)));

Xdata.trial{1} = [phasescm1; phasescm2];

[amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Xdata,37,16,0.5);
qstable = find(abs(dphi_12_dt')<eps(i)); % The points that are below threshold
consecSegs = SplitVec(qstable,'consecutive');
% lengths
segL_ddt = cellfun('length',consecSegs);
% segL(:,i) = segL_ddt ;
histogram(log(segL_ddt/256),'FaceColor','r')

end
