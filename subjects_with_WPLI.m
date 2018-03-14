R.subname = {'DF','DP','DS','JA','JB','JN','LN01','LN02','LN03','MC','MW','SW'};
% JP KB LM WB?  both ON and OFF
% CP JS PB RP
%screen(side,sub) 
screen(:,1)  = [1 0]; %DF % right is borderline
screen(:,2)  = [0 1]; %DP
screen(:,3)  = [1 1]; %DS % good example
screen(:,4)  = [0 0]; %JA % flat spectra- no beta
screen(:,5)  = [1 1]; %JB % nice example
screen(:,6)  = [1 1]; %JN 
screen(:,7)  = [1 1]; %LN01
screen(:,8)  = [1 1]; %LN02 % crazy npd?
screen(:,9)  = [0 1]; %LN03
screen(:,10) = [0 1]; %MC
screen(:,11) = [1 1]; %MW
screen(:,12) = [0 0]; %SW % OFF recordings show no connectivity





