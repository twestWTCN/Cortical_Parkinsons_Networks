function [H,paramEsts] = plot_segL_histogram(pA,X,Xtit,cond)
% intvl = [circ_mean(wrapTo2Pi(pA)')-(pi/4) circ_mean(wrapTo2Pi(pA)')+(pi/4)];
% X = X(intvl(1) < pA & pA<intvl(2));
H = histogram(X,logspace(-1.7,0.3,20),'Normalization','pdf'); hold on
set(gca,'xscale','log')
ed = (H.BinEdges(2:end)); den = (H.Values); %-(H1(cond).BinWidth/2))
trunc = ed(min(find(den))-1);
pdf_trunclognorm = @(x,mu,sigma) lognpdf(x,mu,sigma) ./ (1-logncdf(trunc,mu,sigma));
start = [mean(X) 10^std(X)];
% [paramEsts,paramCIs] = mle(X, 'pdf',pdf_trunclognorm, 'start',start, 'lower',[-1000 0.1]);
% Y = lognpdf(logspace(-1.7,0.3,50),paramEsts(1),paramEsts(2));
% plot(logspace(-1.7,0.3,50),Y);
paramEsts = [1 1];
% [0.0579]    [0.0193]
% f = fit(X',Y','exp1');
% plot(f,X,Y);
% R2 = 1 - sum((Y'-f(X)).^2) ./ sum((Y'-mean(Y)).^2);
% model = sprintf('f(x) = %.2f*exp(%.2f*x)',f.a,f.b);
% annotation(gcf,'textbox',...
%     [0.520 .26-((cond-1)*0.1) 0.37 0.24],...
%     'String',{sprintf('R2 = %.2f',R2); model},...
%     'LineStyle','none',...
%     'HorizontalAlignment','right',...
%     'FontSize',12,...
%     'FitBoxToText','off');
xlabel(Xtit); ylabel('P(X)'); title('Phase Lock')
