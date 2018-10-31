close all; clear
% Series
N = 100;
rangei = linspace(100,1e5,N);
bi = linspace(0,1,N);
for i = 1:N
x = randbetween(0,1000,round(rangei(i)))';
x = sort(x);
n = 1:size(x,2);

% Error Model
eps = (n.^bi(end)).*randn(1,size(x,2));


% Data Model
yhat = 3*x;

y = 3*x + eps;

[xCalc yCalc b(:,i) Rsq(i) bHd(:,i) RsqHd(i) wPval(i)] = linregress(x',y',1);

end
subplot(2,1,1)
scatter(rangei,RsqHd)
subplot(2,1,2)

scatter(x,y); 
hold on
plot(x,yhat,'b')
plot(xCalc,yCalc,'r')
shg

figure
scatter(x,abs(y-yhat))
[dum dum b Rsqr] = linregress(x',abs(y-yhat)')
