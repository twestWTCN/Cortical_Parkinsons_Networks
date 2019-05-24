load('xy_test_Circ.mat')
xdata= x;
sin_fit = @(par,xdata) par(1).*cos(xdata + par(2)) + par(3);
% ypos = y(1,y(1,:)>0); xpos = xdata(1,y(1,:)>0);

for i = 1:3
ypos = y(i,:); xpos = xdata(1,:);

start = [5 0 0];
% [paramEsts,paramCIs] = mle(X, 'pdf',pdf_trunclognorm, 'start',start, 'lower',[-1000 0.1]);

[x,resnorm,~,exitflag,output] = lsqcurvefit(sin_fit,start,xpos,ypos);
figure
xfit = -pi:.1:pi;
scatter(xpos,ypos,'r'); hold on
plot(xfit,sin_fit(x,xfit),'b')
end