function barplot160818(R,binedge,bindata,binstats,degflag)
if nargin<5
    degflag = 0;
end
X = binedge(1:end-1) + ((binedge(2)-binedge(1))/2);
for m = 1:2;for n = 1:length(X);  cmn(n,m,:) = R.condcmap(m,:); end;  end
HB = superbar(1:length(X),bindata(:,:,1),'E',bindata(:,:,2),'BarFaceColor',cmn); %m-by-n-by-3
%     errorbar(1:length(X),AmpBin(:,1,1),AmpBin(:,1,2),'.')
a  = gca;
a.XTick = 1:1:length(X);
if degflag == 0
    a.XTickLabel = sprintfc('%.2f',X(1:1:end));
elseif degflag == 1
    a.XTickLabel =  sprintfc('%.f',rad2deg(X(1:1:end)));
end
a.XTickLabelRotation = 45;
legend(HB(1,:),R.condname);
xlim([0 length(X)+1]); ylim([0 inf])
box off

