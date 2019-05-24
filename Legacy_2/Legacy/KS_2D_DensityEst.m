function [pdfxy,xi,yi] = KS_2D_DensityEst(x,y,Qx,Qy,KSwid)
% Generate some normally distributed data
% x = randn(50,1);
% y = randn(50,1);
%% HOW TO TRANSLATE CONTINUOUS PDF TO EXPECTED RATES? 
% Estimate a continuous pdf from the discrete data
[pdfx xi]= ksdensity(x,Qx,'Bandwidth',0.25);
pdfx = (pdfx*diff(xi(1:2))); % Convert from density to probability
[pdfy yi]= ksdensity(y,Qy,'Bandwidth',0.25);
pdfy = (pdfy*diff(yi(1:2))); % Convert from density to probability
% Create 2-d grid of coordinates and function values, suitable for 3-d plotting
[xxi,yyi]     = meshgrid(xi,yi);
[pdfxx,pdfyy] = meshgrid(pdfx,pdfy);

% Calculate combined pdf, under assumption of independence
pdfxy = pdfxx.*pdfyy; 
pdfxy = pdfxy.*length(x); % Go from probability to expected count
% Plot the results
% mesh(xxi,yyi,pdfxy)
% set(gca,'XLim',[min(xi) max(xi)])
% set(gca,'YLim',[min(yi) max(yi)])
