function [pdfxy] = KS_2D_DensityEst(x,y,Qx,Qy,KSwid)
% Generate some normally distributed data
% x = randn(50,1);
% y = randn(50,1);

% Estimate a continuous pdf from the discrete data
[pdfx xi]= ksdensity(x,Qx,'Bandwidth',KSwid(1));
[pdfy yi]= ksdensity(y,Qy,'Bandwidth',KSwid(2));

% Create 2-d grid of coordinates and function values, suitable for 3-d plotting
[xxi,yyi]     = meshgrid(xi,yi);
[pdfxx,pdfyy] = meshgrid(pdfx,pdfy);

% Calculate combined pdf, under assumption of independence
pdfxy = pdfxx.*pdfyy; 

% Plot the results
% mesh(xxi,yyi,pdfxy)
% set(gca,'XLim',[min(xi) max(xi)])
% set(gca,'YLim',[min(yi) max(yi)])
