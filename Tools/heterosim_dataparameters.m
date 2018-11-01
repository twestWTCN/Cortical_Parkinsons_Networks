close all; clear
B = [28.6843   27.7277   26.4419];
% bR2 = [0.2423    0.1273    0.0181]#
bi =  [-79.9151   -9.6010   -8.1336
    15.3887   10.4056    6.9765];
bampstats=     [17.8697   12.8778    7.1261;
                3.2050   10.8054    1.8719;
                179.7500  115.5000   33.1500];
cmap = linspecer(2);
cmap(3,:) =[0.7 0.7 0.7];
% Series
N = 100;
for i = 1:3
x = randnbetween(log10(bampstats(1,i)),log10(bampstats(2,i)),floor(bampstats(3,i)))';
x = 10.^x;
x(x<5) = [];
x = sort(x);
n = 1:size(x,2);

% Error Model
eps = x.*bi(2,i).*randn(1,size(x,2)) + bi(1,i);

% Data Model
yhat = B(i).*x;

y = abs(yhat + eps);

[xCalc yCalc b(:,i) Rsq(i) bHd(:,i) RsqHd(i) wPval(i)] = linregress(x',y',1);

s = scatter(x,y,25,cmap(i,:),'filled'); 
s.MarkerFaceAlpha = 0.7
hold on
xlim([0 85])

end
RsqHd(end)
% subplot(2,1,1)
% scatter(rangei,RsqHd)

figure
scatter(x,y); 
hold on
plot(x,yhat,'b')
plot(xCalc,yCalc,'r')
xlim([0 100])
shg
% 
figure
scatter(x,abs(y-yhat))
[dum dum b Rsqr] = linregress(x',abs(y-yhat)')
