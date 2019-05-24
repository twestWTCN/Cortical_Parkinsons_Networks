function [Amp AmpBinDat Ampprc binNumel] = binAmp(PhiCol,SegCol,DurCol,Bin,thresh,tottime,mampi)
for i = 1:length(Bin)-1
    binDat = SegCol(:,PhiCol>=Bin(i) & PhiCol<=Bin(i+1))';
    durDat = DurCol(:,PhiCol>=Bin(i) & PhiCol<=Bin(i+1))';
    binNumel(i) = sum(durDat); %numel(durDat);
    AmpBinDat{i} = (((binDat-mampi)./mampi)*100)';
    Ampprc(i) = median(AmpBinDat{i}); % median(AmpBinDat{i}.*(durDat/max(durDat))'); 
    Amp(i,1) = (sum(durDat(binDat<=thresh(1)))./tottime)*100;
%     AmpBinDat{i,1} = binDat(binDat<=thresh(1));
    Amp(i,2) = (sum(durDat(binDat>=thresh(2)))./tottime)*100;
%     AmpBinDat{i,2} = binDat(binDat>=thresh(2));
end
W = binNumel./max(binNumel);
Ampprc = (W.*Ampprc); %./sum(W);
% BinPr2 =(W.*BinPr2)./sum(W);


% W = binNumel/max(binNumel);
% BinMu = BinMu.*W; %(BinMu.*W);
% BinPr2= (BinPr2.*W);
a= 1;