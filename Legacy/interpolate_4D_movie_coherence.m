clear all;
datapathr = 'C:\Users\Tim\Documents\Work\Cortical_Networks\Data\';
nr = 1; subname = {'LM'}; i = 1; ref_chan = 'STN_R23';

load([datapathr subname{i} '\ftdata\cSource4D'])

gSource = cSource;
load([datapathr subname{i} '\ftdata\source_loc_time_' num2str(nr)])
% Compute timings
for xi = 1:numel(vtime)
    vt(xi) = median(vtime{xi});
end

x = 1:size(gSource,2); y = 1:size(gSource,3); z = 1:size(gSource,4); t = 1:size(gSource,1);
[t,x,y,z] = ndgrid(t,x,y,z);

N = 2;
xq = 1:1/N:size(gSource,2); %linspace(1,16,16*N);
yq = 1:1/N:size(gSource,3); %linspace(1,13,13*N);
zq = 1:1/N:size(gSource,4); %linspace(1,12,12*N);
tq = 1:1/N:size(gSource,1); %linspace(1,125,125*N);
[xr,yr,zr] = ndgrid(xq,yq,zq);
xr = reshape(xr,1,[]); yr = reshape(yr,1,[]); zr = reshape(zr,1,[]);
[tq,xq,yq,zq] = ndgrid(tq,xq,yq,zq);
% Interpolate time vector
vt = interp1(1:size(gSource,1),vt,1:0.2:size(gSource,1));

% gSource(isnan(gSource)) = 0;
Vq = interpn(t,x,y,z,gSource,tq,xq,yq,zq);

f = figure('renderer','zbuffer');
set(gcf,'Position',[680.0000  233.5000  867.0000  744.5000])
nframes = size(tq, 1);

% v = VideoWriter('granger_LM.avi','Motion JPEG 2000');
% v.FrameRate = 20;
% open(v);


for j = 1:nframes
    gtSource = squeeze(Vq(j,:,:,:));
    Vr = reshape(gtSource,1,[]);
    subplot(2,2,1)
    scatter3(xr,yr,zr,(Vr.^2).*150,Vr);
    view([0 0]); caxis([0 0.5]);
    
    subplot(2,2,2)
    scatter3(xr,yr,zr,(Vr.^2).*150,Vr);
    view([0 90]); caxis([0 0.5]);
    
    subplot(2,2,3)
    scatter3(xr,yr,zr,(Vr.^2).*150,Vr);
    view([90 0]); caxis([0 0.5]);
    
    subplot(2,2,4)
    scatter3(xr,yr,zr,(Vr.^2).*150,Vr);
    view([45 45]); caxis([0 0.5]);
    title(sprintf('Time %.1fs , frame %.0f',vt(j),j))
    
    colorbar('peer',gca,'Position',...
    [0.93 0.34 0.031 0.34]);

    M(j) = getframe(gcf);
%     writeVideo(v,M(j));
    %     contourslice(gtSource,1:2.5:16,1:2.5:13,1:2.5:12)
end
save([datapathr subname{i} '\ftdata\coherence_movie_' num2str(nr)],'M','-v7.3')

% close(v);
shg;movie(f,M); 