function [L Depbinned binmid pdense] = plot_PA_Dep_relation_v6(cond,pA_dist,Dep,Depname,Xedge,Yedge,f,tend,interpolgrid)
condname = {'ON','OFF'};
Dep(isnan(pA_dist)) = []; pA_dist(isnan(pA_dist)) = [];
pA_dist(isnan(Dep)) = []; Dep(isnan(Dep)) = [];
cmap = [0 0 1; 1 0 0];
bedge = Xedge; %linspace(-pi,pi,nsplit);  
clear segLbinned
figure(f(2))
pA_dist(isnan(Dep)) = [];  Dep(isnan(Dep)) = [];

[N,Xedges,Yedges] = histcounts2(pA_dist,Dep,Xedge,Yedge); %,'Normalization','probability');
N = N./tend;
% [N2] = histcountsn([pA_dist; Dep]',{Xedge(2:end),Yedge},[1 0])
N2 = [N;N;N]';
if interpolgrid == 1
    XMid = Xedges(1:end-1) + diff(Xedges(1:2))/2; XMid2 = [XMid XMid XMid]; XMid2 = unwrap(XMid2);
    YMid = Yedges(1:end-1) + diff(Yedges(1:2))/2;
    [XE YE] = meshgrid(XMid2,YMid);
    % N(N==0) = NaN;
    % h = fspecial('gaussian',[3 3]);
    % N = filter2(h, N);
    Xq = linspace(Xedges(1),Xedges(end),floor(length(Xedges)*2)); Xq2 = [Xq(1:end) Xq(2:end) Xq(2:end)]; Xq2 = unwrap(Xq2);
    Yq = linspace(YMid(1),YMid(end),floor(length(YMid)*2));
    [XQ YQ] = meshgrid(Xq2,Yq);
    Nint = interp2(XE,YE,N2,XQ,YQ);
    N = Nint(:,length(Xq):(length(Xq)+1)+length(Xq));
    Xedges = wrapToPi(Xq2(length(Xq):(length(Xq)+1)+length(Xq))); 
    Xedges(1) = -Xedges(1);
    Xedges(end-1:end) = [pi+diff(Xedges(1:2)) pi+2*diff(Xedges(1:2))];
    Xedges = unwrap(Xedges);
    Yedges = Yq;
end

pdense.N = N; pdense.Xedges = Xedges; pdense.Yedges = Yedges;
% Xmidges =Xedges(1:end-1)+(diff(Xedges)/2)
% Ymidges = Yedges(1:end-1)+(diff(Yedges)/2);
% Nbed = zeros(size(N)+1);
% Nbed(1:end-1,1:end-1) = N;
pcolor(Xedges,Yedges,N);xlim([-pi pi]); ylim([Yedges(1) Yedges(end)])
% contourf(Xedges,Yedges,N);xlim([-pi pi]); ylim([Yedges(1) Yedges(end)])
% imagesc(Xedges,Yedges,N');
 ylabel(Depname); xlabel('Phi_1 - Phi_2')
set(gca,'YDir','normal')
hold on
scatter(pA_dist,Dep,30,cmap(cond,:),'x','LineWidth',2); set(gca,'YDir','normal');

%          pA_dist = pA_dist(Dep>mean(Dep)); Dep = Dep(Dep>mean(Dep));%% ####
figure(f(1))
hold on
scatter(pA_dist,Dep,20,cmap(cond,:),'filled'); set(gca,'YDir','normal');
xlabel('phi_1 - phi_2'); ylabel(Depname)
%         colorbar; caxis([0 0.05])
title(['MaxCohM1 ' condname{cond}])
for be = 1:numel(bedge)-1
    Depbinned(:,be) = [mean(Dep(pA_dist>=bedge(be) & pA_dist<bedge(be+1))) prctile(Dep(pA_dist>=bedge(be) & pA_dist<bedge(be+1)),90) prctile(Dep(pA_dist>=bedge(be) & pA_dist<bedge(be+1)),10)];
end
binmid = bedge(1:end-1) + diff(bedge(1:2))/2;
L = plot(binmid,Depbinned(1,:),'Color',cmap(cond,:),'LineWidth',4); plot(binmid,Depbinned(2,:),'Color',cmap(cond,:),'LineWidth',4,'linestyle','--'); plot(binmid,Depbinned(3,:),'Color',cmap(cond,:),'LineWidth',4,'linestyle','--')

% Now shift to max
% x = Depbinned(2,:);
% p = find(x==max(x));
% pshift = p - (length(x)/2);
% Depshift = circshift(x,p);
% 
Depbinned = Depbinned(2,:);