function [BinMu,BinSEM,binDat,BinPr2] = binstats(PhiCol,SegCol,Bin)
for i = 1:length(Bin)-1
    binDat{i} = SegCol(:,PhiCol>=Bin(i) & PhiCol<=Bin(i+1))';
    binDatC = binDat{i}(binDat{i}>25);
    binW(i) =  nansum(binDatC);
    binNumel(i) = numel(binDatC);
    BinMu(:,i) = prctile(binDatC,50); %nanmean(binDat{i}); %*(numel(binDat{i})/numel(SegCol)));
    BinSEM(:,i) = nanstd(binDatC)/size(binDat{i},1);
    BinPr2(:,i) = [prctile(binDat{i},50)  prctile(binDat{i},85) prctile(binDat{i},15)];
end
% W = binNumel./sum(binNumel);
% BinMu = (W.*BinMu)./sum(W);
% BinPr2 =(W.*BinPr2)./sum(W);


% W = binNumel/max(binNumel);
% BinMu = BinMu.*W; %(BinMu.*W);
% BinPr2= (BinPr2.*W);
a= 1;