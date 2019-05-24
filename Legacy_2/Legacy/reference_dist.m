% permute to get reference distribution
close all
N = 1000;
load('phi_permute.mat'); fsamp = 200;
phi_can1 = phi1; phi_can2 = phi2;
clear PLV PLI ddphi ddphi_seg
parfor n = 1:N
%     p1 = randperm(numel(phi_can1));
%     p2 = randperm(numel(phi_can2));
    
%     H1 = hilbert(Xdata.trial{1}(1,:));
%     H2 = hilbert(Xdata.trial{1}(2,:));
% 
%     phi_1 = phi_can1(p1); phi_2 = phi_can1(p2);
    phi_1 = randbetween(-pi,pi,1,numel(phi_can1));
    phi_2 = randbetween(-pi,pi,1,numel(phi_can1));
    dphi_12 = unwrap(phi1-phi2);
    dphi_12_dt = diff(dphi_12);
    
    qstable = find(abs(dphi_12_dt)<0.025); % The points that are below threshold
    consecSegs = SplitVec(qstable,'consecutive');
    ddphi_seg(n,:) = cellfun('length',consecSegs);
    
    period = (4/16)*fsamp;
    [slide_ddphi,sind] = slideWindow(dphi_12_dt, floor(period),0);
    ddphi(n,:) = mean(slide_ddphi,1); % Use this for detecting sync segments greater than chance.
    
    [slide_dphi_12,sind] = slideWindow(phi1-phi2, floor(period), 0);
    PLV(n,:) = abs(mean(exp(1i*slide_dphi_12),1));
%     PLI(n,:) = abs(mean(sign(slide_dphi_12),1));  % PLI
    PLI(n,:) = abs(mean(abs(1i*slide_dphi_12).*sign(slide_dphi_12),1))./ mean(abs(1i*slide_dphi_12));  %wPLI
%    PLI(n,:)  = abs(mean(    abs(imag(slide_dphi_12)) .* sign(imag(slide_dphi_12)) )) ./ mean(abs(imag(slide_dphi_12)));
    disp(n)
end

figure
subplot(1,3,1)
pli_ref = reshape(PLI,1,[]);
histogram(pli_ref,50,'Normalization','probability','BinMethod','sturges')
wpli_ci = prctile(pli_ref,95);
hold on; plot([wpli_ci wpli_ci],[0 1],'k--')
xlabel('wPLI'); ylabel('P(x)'); ylim([0 0.5])

% subplot(1,3,2)
% plv_ref = reshape(PLV,1,[]);
% histogram(plv_ref,50,'Normalization','probability','BinMethod','sturges');
% plv_ci = prctile(plv_ref,95);
% hold on; plot([plv_ci plv_ci],[0 1],'k--')
% xlabel('plv'); ylabel('P(x)'); ylim([0 0.3])

subplot(1,3,2)
ddphi_seg_ref = reshape(ddphi_seg,1,[]);
histogram(ddphi_seg_ref,50,'Normalization','probability','BinMethod','sturges');
ddphi_seg_ci = prctile(ddphi_seg_ref,95);
hold on; plot([ddphi_seg_ci ddphi_seg_ci],[0 1],'k--')
xlabel('Seg Length'); ylabel('P(x)'); ylim([0 0.3])

subplot(1,3,3)
ddphi_ref = reshape(ddphi,1,[]);
ddphi_ci = prctile(ddphi_ref,95);
histogram(ddphi_ref,50,'Normalization','probability','BinMethod','sturges');
hold on; plot([ddphi_ci ddphi_ci],[0 1],'k--');
hold on; plot([-ddphi_ci -ddphi_ci],[0 1],'k--')
xlabel('ddphi'); ylabel('P(x)'); ylim([0 0.2])

set(gcf,'Position',[265 469 1230 458])
