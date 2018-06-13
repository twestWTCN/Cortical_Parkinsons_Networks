function [BinMu,BinSEM,binDat] = binstats(PhiCol,SegCol,Bin)
for i = 1:length(Bin)-1
    binDat{i} = SegCol(:,PhiCol>=Bin(i) & PhiCol<=Bin(i+1))';
    binNumel(i) = numel(binDat{i});
    BinMu(:,i) = nanmean(binDat{i});%*(numel(binDat)/numel(SegCol));
    BinSEM(:,i) = nanstd(binDat{i})/size(binDat{i},1);
end
BinMu = BinMu.*(binNumel./max(binNumel));