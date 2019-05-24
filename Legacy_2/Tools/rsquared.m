function Rsq = rsquared(Y,yCalc)
 
Rsq = 1 - sum((Y - yCalc).^2)/sum((Y - mean(Y)).^2);