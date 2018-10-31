function UPDRS_Scatter(a,b,alab,xlab,ylab,anotoff,statloc)
a(isnan(b)) = []; alab(isnan(b)) = []; b(isnan(b)) = [];
b(isnan(a)) = []; alab(isnan(a)) = []; a(isnan(a)) = [];
[ar stat] = linplot_PD(a,b,xlab,ylab,[1 0 0]);
text(a+anotoff(1), b+anotoff(2), alab);
text(statloc(1),statloc(2),...
{sprintf('R = %.2f',stat.Rcoeff);
    sprintf('P = %.3f',stat.p);
    sprintf('R_{rob} = %.2f',stat.R_rob);
    sprintf('P_{rob} = %.3f',stat.p_rob)},...
'FontWeight','bold','Color','b')
