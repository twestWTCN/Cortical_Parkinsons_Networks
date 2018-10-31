close all; clear
% bhd = [0.2423    0.1273    0.0181]#
% R2 = 


% Series
N = 100;
rangei = [300 150 100];
bi = [1.2 1.15 1.05]; % linspace(0,2,N);
for i = 1:3
x = randbetween(10,100,round(rangei(i)))';
x = sort(x);
n = 1:size(x,2);

% Error Model
eps = (x.^bi(i)).*randn(1,size(x,2));

% Data Model
yhat = 25*x;

y = abs(yhat + eps);

[xCalc yCalc b(:,i) Rsq(i) bHd(:,i) RsqHd(i) wPval(i)] = linregress(x',y',1);

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
