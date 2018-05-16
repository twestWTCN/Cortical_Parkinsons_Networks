function phishift = findAmpPhi(amp,phi,phiBin)
% phiBin = linspace(-pi,pi,N);
phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);
for i = 1:length(phiBin)-1
    binDat = amp(:,phi>=phiBin(i) & phi<=phiBin(i+1))';
%     ampBinMu(:,i) = nanmean(binDat);
    ampBinMu(:,i) = prctile(binDat,85)*(numel(binDat)/numel(amp)); %%     Weighted mean
    binN(:,i) = numel(binDat);
    ampBinSEM(:,i) = nanstd(binDat)/size(binDat,1);
end

phishift = phiBinMid(ampBinMu==max(ampBinMu));
% wrapToPi(phiBinMid - phishift)
phishift = wrapToPi(phi - phishift(1));