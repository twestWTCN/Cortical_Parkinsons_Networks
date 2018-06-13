function [BinMu,BinSEM,binDat,BinPr2] = binstats(PhiCol,SegCol,Bin)
for i = 1:length(Bin)-1
    binDat{i} = SegCol(:,PhiCol>=Bin(i) & PhiCol<=Bin(i+1))';
    binNumel(i) = numel(binDat{i});
    BinMu(:,i) = nanmean(binDat{i});%*(numel(binDat)/numel(SegCol));
    BinSEM(:,i) = nanstd(binDat{i})/size(binDat{i},1);
    BinPr2(:,i) = [prctile(binDat{i},50)  prctile(binDat{i},99) prctile(binDat{i},1)];
end
BinMu = BinMu.*(binNumel./max(binNumel));
BinPr2 = BinPr2.*(binNumel./max(binNumel));